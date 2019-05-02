
library(data.table)
library(lme4)
library(nlme)
library(INLA)
####load the data set
data.use <- fread("table 1 data.csv")

#### lme4
summary(lmer(yjt ~ X1 +X2 + (0+X2|PID), data = data.test))

#### nlme
summary(lmer(yjt ~ X1 +X2 + (0+X2|PID), data = data.test))

summary(inla(yjt ~ X1 +X2+
                       f(PID, X2, model = "iid")
                      , data = data.test, 
                      control.predictor=list(compute=TRUE), control.compute = list(
                       dic = TRUE), control.inla = list(numint.maxfeval=80000000)))


