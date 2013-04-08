
# read summaries
d <- read.table('muscleData.csv',header=T,sep=',')

# read synthetic data
d2 <- read.table('syntheticData.csv',header=T,sep=',')


# checks (with about 2 everywhere, had fractional weights)
dg = subset(d2,repgroup==2)
dim(dg)[[1]]
sapply(dg,sum)
sapply(subset(dg,deceased==1),sum)
sapply(subset(dg,MuscularStrength.lower==1),sum)
sapply(subset(dg,MuscularStrength.middle==1),sum)
sapply(subset(dg,MuscularStrength.upper==1),sum)

# show fits 
print(summary(lm(deceased~0+MuscularStrength.lower+MuscularStrength.middle+MuscularStrength.upper,data=dg)))
print(summary(lm(deceased~0+MuscularStrength.lower+MuscularStrength.middle+MuscularStrength.upper+sedentary+current.smoker+five.drinks.weekly+diabetes.millitus+ hypertension+hyercholeserolaemia +family.cardiovascular,data=dg)))

# show number of individuals in weighted version
dim(subset(dg,repnum==0))[[1]]


# how we found a bad group
# fit single variable model and all variable model over each group
for(g in unique(d2$repgroup)) {
   di = subset(d2,repgroup==g)
   print(paste('start group',g))
   print(summary(lm(deceased~0+MuscularStrength.lower+MuscularStrength.middle+MuscularStrength.upper,data=di)))
   print(summary(lm(deceased~0+MuscularStrength.lower+MuscularStrength.middle+MuscularStrength.upper+sedentary+current.smoker+five.drinks.weekly+diabetes.millitus+ hypertension+hyercholeserolaemia +family.cardiovascular,data=di)))
  print(paste('end group',g))
}


# confirm 2-ways all match
vars = c('sedentary', 'current.smoker', 'five.drinks.weekly',
  'diabetes.millitus', 'hypertension', 'hyercholeserolaemia',
  'family.cardiovascular')
v = vars[5]
for(g in unique(d2$repgroup)) {
   di = subset(d2,repgroup==g)
   print(paste('start group',g))
    f = as.formula(paste('deceased~0+MuscularStrength.lower+MuscularStrength.middle+MuscularStrength.upper',v,sep=' + '))
   print(summary(lm(f,data=di)))
  print(paste('end group',g))
}
for(v in vars) {
  f = as.formula(paste('deceased~0+MuscularStrength.lower+MuscularStrength.middle+MuscularStrength.upper',v,sep=' + '))
  print(summary(lm(f,data=di)))
}
