# Code for Simulating data from model M3 and 
# fitting model M1. Results from 100 replications 
# saved into M3M1_1.csv, M3M1_2.csv, M3M1_3.csv, 
# M3M1_4.csv, which is called by 
# maketables.R to produce results in Table 6
# of the paper by Chen, Harel and Ravishanker.
library(nlme)
library(INLA)
library(msm)

# True parameters
n=100;nt=98
beta.0 = -1.8
beta.1 = 0.15
gamma.2 = 0.7
sigmasq.2 = 0.5
sigmasq.u = 0.25
set.seed(123457)

for(iter in 1:4){
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

    data.test = cbind(rep(1:n, each = nt), rep(1:nt, n), yjt, x1.new, x2.new)
    colnames(data.test) = c("PID", "Day", "yjt","X1", "X2")
    data.test = data.frame(data.test)
    # end of data simulation from M3 model
    
    tryCatch({
      #fit M1 model
      fm.stat=lme(yjt ~ X1 +X2, random = ~0+X2|PID, data = data.test)
      #CI. for all parameters, FE and random effects
      static.itv <- intervals(fm.stat, level = 0.95)
      
      ###1.FE estimate, CI width, intervals, bias, SE
      # beta.0, beta.1, gamma.2
      # para.tot <- (fm.stat$coefficients$fixed)
      para.tot <- c(
        t(cbind(static.itv[[1]],
                c(static.itv[[1]][, 3] - static.itv[[1]][, 1]),
                static.itv[[1]][, 1] <= FE & 
                  FE <= static.itv[[1]][, 3],
                fm.stat$coefficients$fixed-FE,
                sqrt(diag(fm.stat$varFix))
                #,100*(fm.stat$coefficients$fixed-FE)/sqrt(diag(fm.stat$varFix))
        ))
      )
      fm.stat$sigma #residual SE
      ###2.VC estimate, CI width, intervals, bias, SE
      # sigmasq.2
      para.tot <- c(para.tot, static.itv[[2]]$PID^2, diff(range(static.itv[[2]]$PID^2)))
      para.tot <- c(para.tot, 
                    static.itv[[2]]$PID[1]^2 <= sigmasq.2 & 
                      sigmasq.2 <= static.itv[[2]]$PID[3]^2)
      para.tot <- c(para.tot, c(as.numeric(VarCorr(fm.stat))[1] - c(sigmasq.2)))
      # sigmasq.u
      varstat <- fm.stat$apVar
      parstat <- attr(varstat, "Pars")
      para.tot <- c(para.tot, deltamethod(~exp(x1)^2, parstat, varstat))
      
      para.tot <- c(para.tot, static.itv[[3]]^2, diff(range(static.itv[[3]]^2)))
      para.tot <- c(para.tot, 
                    static.itv[[3]][1]^2 <= sigmasq.u &
                      sigmasq.u <= static.itv[[3]][3]^2)
      para.tot <- c(para.tot, 
                    c(as.numeric(VarCorr(fm.stat))[2] - c(sigmasq.u)),
                    deltamethod(~exp(x2)^2, parstat, varstat)
      )
      ####RMSE
      para.tot <- c(para.tot, 
                    mean((predict(fm.stat)-data.test$yjt)^2)
      )
      para.tot <- data.frame(t(para.tot))
      colnames(para.tot) <- 
        c(
          "intlower", "int", "inthigher", "intwidth", "intcoverage", "intbias","intsd",
          "x1lower", "X1", "x1higher", "x1width", "x1coverage","x1bias","x1sd",
          "x2lower", "X2", "x2higher", "x2width", "x2coverage","x2bias","x2sd",
          "sig2lower", "sig2", "sig2higher", "sig2width","sig2coverage", "sig2bias", "sig2sd", 
          "sig0lower", "sig0", "sig0higher", "sig0width","sig0coverage", "sig0bias", "sig0sd",
          "RMSE"
        )
      
      
      para.tot1 <- rbind(para.tot1, para.tot)
    }, 
    error=function(e){
      print(paste0("wrong in seed:", kkkk, "; \n", e))
    })
    
  }
  write.csv(as.matrix.data.frame(para.tot1), paste0("M3M1_", iter, ".csv", sep = ""))
}