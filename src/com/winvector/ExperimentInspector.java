package com.winvector;

import java.net.URI;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.SortedMap;
import java.util.TreeMap;

import com.winvector.linalg.DenseVec;
import com.winvector.linalg.colt.ColtMatrix;
import com.winvector.linalg.jblas.JBlasMatrix;
import com.winvector.lp.LPEQProb;
import com.winvector.lp.LPSoln;
import com.winvector.lp.LPSolver;
import com.winvector.lp.impl.RevisedSimplexSolver;
import com.winvector.util.BurstMap;
import com.winvector.util.TrivialReader;

public class ExperimentInspector {
	private static final String CONDITIONCOL = "condition";
	private static final String OUTCOMECOL = "Decedents";
	private static final String TOTALCOL = "total";
	
	/**
	 * Right now expects there is at least one group of condition names with dots in them and they are 
	 * a mutually exclusive and complete set of levels.  All other conditions are thought of
	 * as stand-alone indicator variables.  The CONDITIONCOL and OUTCOMECOL and TOTALCOL are
	 * all expected to have the names we specify.
	 * 
	 * Right now working with lib/muscleData.csv ( file:lib/muscleData.csv )
	 *  which is related to http://www.win-vector.com/blog/2013/04/worry-about-correctness-and-repeatability-not-p-values/
	 * R commands in lib/rsteps.R
	 * @param args
	 */
	public static void main(final String[] args) throws Exception {
		final URI trainURI = new URI(args[0]);
		final Iterable<BurstMap> origSource = new TrivialReader(trainURI,',',null,false,null,false);
		final ArrayList<BurstMap> rows = new ArrayList<BurstMap>();
		final SortedMap<Integer,String> levelIndices = new TreeMap<Integer,String>();
		final SortedMap<String,Integer> checkableConditions = new TreeMap<String,Integer>();
		for(final BurstMap row: origSource) {
			final String cond = row.getAsString(CONDITIONCOL);
			final int dotIndex = cond.indexOf('.');
			if(dotIndex>0) {
				final String group = cond.substring(0,dotIndex);
				levelIndices.put(rows.size(),group);
			}
			//final double positives = row.getAsDouble(OUTCOMECOL);
			//final double total = row.getAsDouble(TOTALCOL);
			if(row.getAsString(cond)!=null) {
				checkableConditions.put(cond,rows.size());
			}
			//System.out.println("" + cond + "\t" + positives + "\t" + total);
			rows.add(row);
		}
		if(levelIndices.isEmpty()) {
			throw new Exception("expecting at least one level group (to carry total)");
		}
		// compute total
		final String firstGroup = levelIndices.get(levelIndices.firstKey());
		int totalCount = 0;
		for(final Entry<Integer, String> me: levelIndices.entrySet()) {
			if(firstGroup.equals(me.getValue())) {
				totalCount += rows.get(me.getKey()).getAsDouble(TOTALCOL);
			}
		}
		// build the vectors
		final int dim = rows.size() + 1;
		final ArrayList<BitSet> vecs = new ArrayList<BitSet>();
		final BitSet b = new BitSet(dim+1);
		while(true) {
			// find right most zero index
			int pos = dim-1;
			while((pos>=0)&&(b.get(pos))) {
				b.set(pos,false);
				--pos;
			}
			if(pos>=0) {
				b.set(pos,true);
			}
			boolean good = true;
			if(!levelIndices.isEmpty()) {
				// check exactly one set in level group (slow, but not doing this a lot)
				final Map<String,Integer> sum = new TreeMap<String,Integer>();
				for(final String group: levelIndices.values()) {
					sum.put(group,0);
				}
				for(final Map.Entry<Integer,String> me: levelIndices.entrySet()) {
					if(b.get(me.getKey())) {
						sum.put(me.getValue(),1+sum.get(me.getValue()));
					}
				}
				for(final Integer v: sum.values()) {
					good &= v==1;
				}
			}
			if(good) {
				vecs.add((BitSet)b.clone());
			}
			if(pos<0) {
				break;
			}
		}
		final int nWts = vecs.size();
		//System.out.println("possible data rows: " + nWts);
		// build the linear relations on data weights
		final int nChecks = 2 + checkableConditions.size();
		final boolean collarCounts = false;
		final int nSums = nChecks*rows.size();
		final int nRows = nSums + (collarCounts?nWts:0);
		final int nVars = nWts + (collarCounts?nWts:0);
		final ColtMatrix relnMat = ColtMatrix.factory.newMatrix(nRows,nVars,true);
		final double[] rhs = new double[nRows];
		for(int i=0;i<rows.size();++i) {
			final BurstMap row = rows.get(i);
			//final String cond = row.getAsString(CONDITIONCOL);
			final double total = row.getAsDouble(TOTALCOL);
			final double positives = row.getAsDouble(OUTCOMECOL);
			// rows in total and y*rows in total
			// own total
			rhs[nChecks*i] = total;
			// own times outcome total
			rhs[nChecks*i+1] = positives;
			// checkable conditions totals
			{
				int ci = 2;
				for(final Map.Entry<String,Integer> check: checkableConditions.entrySet()) {
					final String checkCond = check.getKey();
					rhs[nChecks*i + ci] = row.getAsDouble(checkCond);
					++ci;
				}
			}
			for(int j=0;j<vecs.size();++j) {
				final BitSet v = vecs.get(j);
				if(v.get(i)) {
					// own total
					relnMat.set(nChecks*i,j,1.0);
					// own times outcome total
					if(v.get(dim-1)) {
						relnMat.set(nChecks*i+1,j,1.0);
					}
					// checkable conditions totals
					{
						int ci = 2;
						for(final Map.Entry<String,Integer> check: checkableConditions.entrySet()) {
							final int checkIndex = check.getValue();
							if(v.get(checkIndex)) {
								relnMat.set(nChecks*i+ci,j,1.0);
							}
							++ci;
						}
					}					
				}
			}
		}
		// estimate expected independent counts (smoothed)
		final double[] expectation = new double[vecs.size()]; 
		for(int j=0;j<vecs.size();++j) {
			final BitSet v = vecs.get(j);
			for(int i=0;i<rows.size();++i) {
				if(v.get(i)) {
					expectation[j] += Math.log(1.0+rows.get(i).getAsDouble(TOTALCOL));
				} else {
					expectation[j] += Math.log(1.0+totalCount-rows.get(i).getAsDouble(TOTALCOL));
				}
			}
		}
		expNormalize(expectation,totalCount);
		if(collarCounts) {
			final double nDevs = 3.0;
			for(int i=0;i<nWts;++i) {
				relnMat.set(nSums+i,i,1.0);
				relnMat.set(nSums+i,i+nWts,-1.0); // slack
				rhs[nSums+i] = 0.5*expectation[i] - nDevs*(1+Math.sqrt(expectation[i])) - 10.0; // not tight enough to reliably set without user intervention
				// above is evidence (assuming no bug) that independence is not plausible for our data summaries
			}
		}
		// solve the LP to get the weights
		final Random rand = new Random(2524326L);
		final double[] c = new double[nVars];
		System.out.print("repgroup,repnum,wt,expectation");
		for(int i=0;i<rows.size();++i) {
			final BurstMap row = rows.get(i);
			final String cond = row.getAsString(CONDITIONCOL);
			System.out.print("," + cond);
		}
		System.out.println(",deceased");
		for(int repGroup=0;repGroup<10;++repGroup) {
			for(int i=0;i<nWts;++i) {
				c[i] = rand.nextGaussian();
			}
			final LPEQProb prob = new LPEQProb(relnMat.columnMatrix(),rhs,new DenseVec(c));
			//prob.print();
			final LPSolver solver = new RevisedSimplexSolver();
			final double tol = 1.0e-3;
			final LPSoln soln = solver.solve(prob, null, tol, 100000, JBlasMatrix.factory);
			prob.checkPrimFeas(soln.primalSolution,tol);
			final int nnz = soln.primalSolution.nIndices();
			for(int ii=0;ii<nnz;++ii) {
				final int i = soln.primalSolution.index(ii);
				if(i<nWts) {
					final double wi = soln.primalSolution.value(ii);
					final int intRowWt = (int)Math.round(wi);
					if(intRowWt>0) {
						final BitSet vec = vecs.get(i);
						for(int repnum=0;repnum<intRowWt;++repnum) {
							System.out.print(repGroup + "," + repnum + "," + wi + "," + expectation[i]);
							for(int j=0;j<dim;++j) {
								System.out.print("," + (vec.get(j)?"1":"0"));
							}
							System.out.println();
						}
					}
				}
			}
		}
	}

	private static void expNormalize(final double[] logExpectation, final double totalCount) {
		double max = Double.NEGATIVE_INFINITY;
		for(final double li: logExpectation) {
			max = Math.max(max,li);
		}
		double observed = 0.0;
		for(int i=0;i<logExpectation.length;++i) {
			final double ci = Math.exp(logExpectation[i]-max);
			logExpectation[i] = ci;
			observed += ci;
		}
		final double scale = totalCount/observed;
		for(int i=0;i<logExpectation.length;++i) {
			logExpectation[i] *= scale;
		}
	}

}
