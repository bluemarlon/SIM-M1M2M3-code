library(dplyr)

# True parameters
n=50;nt=49
beta.0 = -1.8
beta.1 = 0.15
gamma.2 = 0.7
sigmasq.2 = 0.5
sigmasq.u = 0.25

alpha = 0.05
p=0.5
MM <- c()
f1 <- c()
ff1 <- c()

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
  
  MM1 <- c()
  
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
  
  # kkkk: Index of the replications of the experiment
  for(kkkk in 1:100){
    # start of data simulation from M3 model
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
    
    ###########################step 1###########################
    et <- c()
    for(k in unique(data.test$PID)){
      data.sub <- data.test %>% tbl_df() %>% filter(PID == k)
      lm.fit <- lm(yjt~X1+X2, data = data.sub)
      et <- rbind(et, cbind(lm.fit$residuals, k))
    }
    ###########################step 2###########################
    f <- c()
    for(k in unique(data.test$PID)){
      data.sub <- data.test %>% tbl_df() %>% filter(PID == k)
      data.sub$et2 <- et[et[,2]==k, 1]^2
      lm.fit2 <- lm(et2~X1+X2, data = data.sub)
      f <- c(f, any(summary(lm.fit2)$coefficients[-1, 4] < alpha) )
    }
    ###########################step 3###########################
    if (mean(f)> p){
      MM1 <- c(MM1, 3)
      f1 <- cbind(f1, f)
      ff1 <- cbind(ff1, NA)
      next
    }
    ###########################step 4###########################
    ff <- c()
    for(k in unique(data.test$PID)){
      ff <- c(ff, Box.test(et[et[,2]==k, 1])$p.value < alpha)
    }
    ff1 <- cbind(ff1, ff)
    if(mean(ff)> p){
      MM1 <- c(MM1, 2)
    }else{
      MM1 <- c(MM1, 1)
    }
  }
  
  MM <- cbind(MM, MM1)
}

apply(MM, 2, table)

# saveRDS(list(f1, ff1), "DMf.rds")
