
Generate synthetic data sets matching (to within +- a few individuals due to rounding) the published summaries of the data set:

“Association between muscular strength and mortality in men: prospective cohort study,” Ruiz et. al. BMJ 2008;337:a439 http://www.bmj.com/content/337/bmj.a439

This allows us to see a range of data sets that match the claimed summaries and explore a range of possible fit results you could experience from such data.  The point is to demonstrate the summary tables give are no replacement for having the actual data.  The paper's results look good, but there are data sets that match the given summaries that fail to support the results.  Most synthetic data sets generated do reproduce the published claims, so this is more about the desirability of sharing data than any actual complaint about the paper.

To run (assuming compiled source all jars in lib are on classpath)
  java com.winvector.ExperimentInspector file:lib/muscleData.csv > lib/syntheticData.csv

The synthetic data then has 10 different data sets that match the summaries from lib/muscleData.csv.  lib/rsteps.R shows how to fit the data and try to reproduce the paper's results.

This is related to the article:
  http://www.win-vector.com/blog/2013/04/worry-about-correctness-and-repeatability-not-p-values/

