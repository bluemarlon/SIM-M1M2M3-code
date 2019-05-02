# Code for Simulating data from model M1 and 
# fitting model M3. Results from 100 replications 
# saved into M1M3.csv, which is called by 
# maketables.R to produce results in Table 4
# of the paper by Chen, Harel and Ravishanker.

library(nlme)
library(INLA)
library(msm)

# True parameters
n=50;nt=49
beta.0 = -1.8
beta.1 = 0.15
gamma.2 = 0.7
sigmasq.2 = 0.5
sigmasq.u = 0.25

set.seed(123457)

# Set up Fixed Effects true parameters
FE <- c(beta.0, beta.1, gamma.2)
# Generate predictor x1
x1.new <- rnorm(n * nt, 0, 1)
X1 = matrix(x1.new, nt)
# Generate predictor x2
x2.shape <- 1.5
x2.scale <- 1
x2.new <- rgamma(n * nt, shape = x2.shape, scale = x2.scale)
X2 = matrix(x2.new, nt)

para.tot1 = data.frame()

# kkkk: Index of the replications of the experiment
for(kkkk in 1:100){
# Start simulation of data from model M3
  para.tot = data.frame()
  
  yjt = c()
  random.save = c()
  residual.save = c()
  
  v2 = rnorm(n, 0, sqrt(sigmasq.2))
  random.save = c(random.save, v2)
  for(j in 1:n){
    u.jt = rnorm(nt, 0, sqrt(sigmasq.u))
    residual.save = c(residual.save, u.jt)
    for(tt in 1:nt){
      yjt = c(yjt,   u.jt[tt]+ beta.0 + beta.1 * X1[tt,j] + (gamma.2+v2[j])*X2[tt,j])
    }
  }
  
  data.test = cbind(rep(1:n, each = nt), rep(1:nt, n), yjt, x1.new, x2.new)
  colnames(data.test) = c("PID", "Day", "yjt","X1", "X2")
  data.test = data.frame(data.test)
# End of data simulation from model M3 
  
# Fit model M3 using R-INLA   
    phi = 0
      
      INLAdlm1 = inla(yjt ~ X1+X2+
                        f(Day, X2, model = "ar1", replicate = PID )
                      ,  data = data.test, quantiles = c(0.05, 0.5, 0.95),
                      control.predictor=list(compute=TRUE), control.compute = list(
                        dic = TRUE))
# Fixed Effects Estimates
    par.dlm <- c(
      t(cbind(
        INLAdlm1$summary.fixed[, c(3:5)],
        INLAdlm1$summary.fixed[, 5] - INLAdlm1$summary.fixed[, 3],
        INLAdlm1$summary.fixed[, 3] <= FE &
          FE <= INLAdlm1$summary.fixed[, 5],
        INLAdlm1$summary.fixed[, 4] - FE,
        INLAdlm1$summary.fixed[, 2]
      ))
    )
    marginaldm <- inla.tmarginal(fun = function(x){1/x},
                                 INLAdlm1$marginals.hyperpar$`Precision for Day`)
    sig2sd <- unlist(inla.zmarginal(marginaldm, silent = TRUE)[2])
    varsig2 <- inla.qmarginal(p = c(0.05, 0.5, 0.95), marginaldm)
    par.dlm <- c(par.dlm, c(varsig2,
                            varsig2[3] - varsig2[1],
                            varsig2[1] <= sigmasq.2/(1-phi^2) & sigmasq.2/(1-phi^2) <= varsig2[3],
                            varsig2[2] - sigmasq.2/(1-phi^2),
                            sig2sd,
                            varsig2[2]*(1- INLAdlm1$summary.hyperpar[3, 3]^2) - sigmasq.2))
    
    par.dlm <- c(par.dlm,
                 c(unlist(INLAdlm1$summary.hyperpar[3, 3:5]),
                   INLAdlm1$summary.hyperpar[3, 5] - INLAdlm1$summary.hyperpar[3, 3],
                   INLAdlm1$summary.hyperpar[3, 3] <= phi &
                     phi <= INLAdlm1$summary.hyperpar[3, 5],
                   INLAdlm1$summary.hyperpar[3, 4] - phi,
                   INLAdlm1$summary.hyperpar[3, 2]))
    marginaldm0 <- inla.tmarginal(fun = function(x){1/x},
                                  INLAdlm1$marginals.hyperpar$`Precision for the Gaussian observations`)
    sig0sd <-unlist(inla.zmarginal(marginaldm0, silent = TRUE)[2])
    varsig0 <- inla.qmarginal(p = c(0.05, 0.5, 0.95), marginaldm0)
    par.dlm <- c(par.dlm,
                 c(varsig0,
                   varsig0[3] - varsig0[1],
                   varsig0[1] <= sigmasq.u &
                     sigmasq.u <= varsig0[3],
                   varsig0[2] - sigmasq.u,
                   sig0sd)
    )
# RMSE
    par.dlm <- c(par.dlm,
                 mean((INLAdlm1$summary.fitted.values$mean-data.test$yjt)^2)
    )
    par.dlm <- data.frame(t(par.dlm))
    
    colnames(par.dlm) <-
      c(
        "intlower", "int", "inthigher", "intwidth", "intcoverage", "intbias","intsd",
        "x1lower", "X1", "x1higher", "x1width", "x1coverage","x1bias","x1sd",
        "x2lower", "X2", "x2higher", "x2width", "x2coverage","x2bias","x2sd",
        "Vx2lower", "Vx2", "Vx2higher", "Vx2width", "Vx2coverage", "Vx2bias","Vx2sd", "sig2bias",
        "philower", "phi", "phihigher", "phiwidth", "phicoverage", "phibias","phisd",
        "sig0lower", "sig0", "sig0higher", "sig0width", "sig0coverage", "sig0bias", "sig0sd",
        "RMSE"
      )
    
    
    para.tot1 <- rbind(para.tot1, par.dlm)
}
write.csv(as.matrix.data.frame(para.tot1), "M1M3.csv")
