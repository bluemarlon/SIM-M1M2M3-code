
library(data.table)
library(nlme)

data.use <- fread("table 2 data.csv")
fm = lme(MissedDose ~ X1 +x2, random = ~0+ X2|PID, 
   correlation = corAR1(form = ~Day|PID), data = contsim, 
   control=lmeControl(maxIter = 8000,msMaxIter = 20))

summary(fm)