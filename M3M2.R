# Code for Simulating data from model M3 and 
# fitting model M2. Results from 100 replications 
# saved into M3M2_1.csv, M3M2_2.csv, M3M2_3.csv, 
# M3M2_4.csv, which is called by 
# maketables.R to produce results in Table 6
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

for(iter in 1:4){
  # iter for four different phi cases
  if(iter==1){
    phi = -0.9
  }else if(iter==2){
    phi=-0.4
  }else if(iter==3){
    phi=0.4
  }else if(iter==4){
    phi=0.9
  }
  
  
  
  
  # Set up Fixed Effects True parameters
  FE <- c(beta.0, beta.1, gamma.2)
  # generate predictor x1
  x1.new <- rnorm(n * nt, 0, 1)
  X1 = matrix(x1.new, nt)
  # generate predictor x2
  x2.shape <- 1.5
  x2.scale <- 1
  x2.new <- rgamma(n * nt, shape = x2.shape, scale = x2.scale)
  X2 = matrix(x2.new, nt)
  
  para.tot1 = data.frame()
  
  # kkkk: Index of the replications of the experiment
  for(kkkk in 1:100){
    # start of data simulation from M3 model
    para.tot = data.frame()
    
    yjt = c()
    random.save = c()
    residual.save = c()
    
    
    for(j in 1:n){
      u.jt = rnorm(nt, 0, sqrt(sigmasq.u))
      residual.save = c(residual.save, u.jt)
      
      beta.2jt = arima.sim(n=nt+100, model=list(ar=c(phi)), sd = sqrt(sigmasq.2))[101:(nt+100)] #phi, sigma.w
      random.save = c(random.save, beta.2jt)
      
      for(tt in 1:nt){
        yjt = c(yjt,   u.jt[tt]+ beta.0 + beta.1 * X1[tt,j] + (gamma.2+beta.2jt[tt])*X2[tt,j])
      }
    }
    # end of data simulation from M3 model
    
    
    data.test = cbind(rep(1:n, each = nt), rep(1:nt, n), yjt, x1.new, x2.new)
    colnames(data.test) = c("PID", "Day", "yjt","X1", "X2")
    data.test = data.frame(data.test)
    
    ctrl <- lmeControl(opt='optim');
    
    fm.repeat <<- lme(yjt ~ X1 +X2, random = ~0+ X2|PID ,
                      correlation = corAR1(form = ~Day|PID), data = data.test)
    repeat.itv <<- intervals(fm.repeat, level = 0.95)
    
    ###1.FE estimate, width, intervals, FEbias, SE, sdbias
    par.rep <- c(
      t(cbind(repeat.itv[[1]],
              c(repeat.itv[[1]][, 3] - repeat.itv[[1]][, 1]),
              repeat.itv[[1]][, 1] <= FE &
                FE <= repeat.itv[[1]][, 3],
              fm.repeat$coefficients$fixed-FE,
              sqrt(diag(fm.repeat$varFix))
      ))
    )
    
    ##VC estimate
    par.rep <- c(par.rep, repeat.itv[[2]]$PID^2, diff(range(repeat.itv[[2]]$PID^2)))
    par.rep <- c(par.rep,
                 repeat.itv[[2]]$PID[1]^2<=sigmasq.2 &
                   sigmasq.2<=repeat.itv[[2]]$PID[3]^2)
    par.rep <- c(par.rep, as.numeric(VarCorr(fm.repeat))[1] - sigmasq.2)
    ###VC sd
    varrep <- fm.repeat$apVar
    parrep <- attr(varrep, "Pars")
    par.rep <- c(par.rep, deltamethod(~exp(x1)^2, parrep, varrep))
    
    par.rep <- c(par.rep, repeat.itv[[3]], diff(range(repeat.itv[[3]])))
    par.rep <- c(par.rep,
                 repeat.itv[[3]][1] <= phi &
                   phi <= repeat.itv[[3]][3])
    par.rep <- c(par.rep,
                 repeat.itv[[3]][2] - phi
    )
    par.rep <- c(par.rep, repeat.itv[[4]]^2, diff(range(repeat.itv[[4]]^2)))
    par.rep <- c(par.rep,
                 repeat.itv[[4]][1]^2 <= sigmasq.u/(1-phi^2) &
                   sigmasq.u/(1-phi^2) <= repeat.itv[[4]][3]^2 )
    par.rep <- c(par.rep,
                 as.numeric(VarCorr(fm.repeat))[2] - sigmasq.u/(1-phi^2),
                 deltamethod(~exp(x3)^2, parrep, varrep)
    )
    par.rep <- c(par.rep,
                 as.numeric(VarCorr(fm.repeat))[2] * (1-repeat.itv[[3]][2]^2)-sigmasq.u
    )
    ###RMSE
    par.rep <- c(par.rep,
                 mean((predict(fm.repeat)-data.test$yjt)^2)
    )
    par.rep <- data.frame(t(par.rep))
    
    colnames(par.rep) <-
      c(
        "intlower", "int", "inthigher", "intwidth", "intcoverage", "intbias","intsd",
        "x1lower", "X1", "x1higher", "x1width", "x1coverage","x1bias","x1sd",
        "x2lower", "X2", "x2higher", "x2width", "x2coverage","x2bias","x2sd",
        "sig2lower", "sig2", "sig2higher", "sig2width", "sig2coverage", "sig2bias", "sig2sd",
        "philower", "phi", "phihigher", "phiwidth", "phicoverage", "phibias",
        "reslower", "res", "reshigher", "reswidth", "rescoverage", "resbias", "resd",
        "sig0bias",
        "RMSE"
      )
    
    para.tot1 <- rbind(para.tot1, par.rep)
  }
  write.csv(as.matrix.data.frame(para.tot1), paste0("M3M2_", iter, ".csv", sep = ""))
}
