


sdbias <- c()
rmse <- c()
width <- c()
coverage <- c()
##### if it is for the table 4, only use k=1, as there is only one csv file for M1 simulations
k <- 1
for(k in 1:4){
  #################################SM#######################################
  # dat1 <- read.csv(paste0("M3M1",k, ".csv", sep = ""))
  # dat1 <- read.csv(paste0("M2M1", k, ".csv", sep = ""))
  dat1 <- read.csv(paste0("M1M1", ".csv", sep = ""))
  
  

  need.list <- c("intbias", "intsd", "x1bias", "x1sd","x2bias", "x2sd", "sig2bias", "sig2sd", "sig0bias", "sig0sd")
  

  bias.list <- c("intbias", "x1bias", "x2bias", "sig2bias", "sig0bias")
  

  biasd.list <- c("intsd", "x1sd", "x2sd", "sig2sd", "sig0sd")
  

  cal1<- dat1[, need.list]

  
  bias1 <- cal1[, bias.list]
  biasd1 <- cbind(cal1[, biasd.list]) 
  #################################RM#######################################
  # dat2 <- read.csv(paste0("M3M2_", k, ".csv", sep = ""))
  
  # dat2 <- read.csv(paste0("M2M2_", k, ".csv", sep = ""))
  
  dat2 <- read.csv(paste0("M1M2", ".csv", sep = ""))
  
  need.list <- c("intbias", "intsd", "x1bias", "x1sd","x2bias", "x2sd", "sig2bias", "sig2sd", "phibias", "resbias", "resd") 
  
  bias.list <- c("intbias", "x1bias", "x2bias", "sig2bias", "resbias", "phibias")
  
  biasd.list <- c("intsd", "x1sd", "x2sd", "sig2sd", "resd")
  
  
  cal2<- dat2[, need.list]
  
  
  bias2 <- cal2[, bias.list]
  biasd2 <- cbind(cal2[, biasd.list],"phisd"=dat2[, "phiwidth"]/(2*1.96))
  
  colnames(biasd2)[ncol(biasd2)] <- "phisd"
  
  #################################SSM#######################################
  # dat3 <- read.csv(paste0("M3M3_", k, ".csv", sep = ""))
  # dat3 <- read.csv(paste0("M2M3_", k, ".csv", sep = ""))
  dat3 <- read.csv(paste0("M1M3", ".csv", sep = ""))
  
  need.list <- c("intbias", "intsd","x1bias", "x1sd","x2bias", "x2sd", "Vx2bias", "Vx2sd", "phibias", "phisd", "sig0bias", "sig0sd")
  bias.list <- c("intbias", "x1bias", "x2bias", "Vx2bias", "sig0bias", "phibias")
  biasd.list <- c("intsd", "x1sd", "x2sd", "Vx2sd", "sig0sd", "phisd")
  
  cal3<- dat3[, need.list]
  
  
  bias3 <- cal3[, bias.list]
  biasd3 <- cal3[, biasd.list]
  
  #################################make table#######################################
  #################################SM#######################################
  #sdbias
  sdbias1 <- round(colMeans(bias1/biasd1)*100, digits = 1)

  #RMSE
  rmse1 <- round(sqrt(colMeans(bias1^2)), digits = 2)
  wid.list <- c("intwidth","x1width", "x2width", "sig2width", "sig0width")
  
  coverage.list <- c("intcoverage","x1coverage", "x2coverage", "sig2coverage", "sig0coverage")
  
  #width
  width1<- round(colMeans(dat1[, wid.list]), digits = 2)

  #coverage
  coverage1 <- round(colMeans(dat1[, coverage.list])*100, digits = 1)
  #################################RM#######################################
  sdbias2 <- round(colMeans(bias2 /biasd2)*100, digits = 1)

  #RMSE
  print(rmse2 <- round(sqrt(colMeans(bias2^2)), digits = 2))

  wid.list <- c( "intwidth","x1width", "x2width", "sig2width", "reswidth", "phiwidth")
  coverage.list <- c("intcoverage","x1coverage", "x2coverage", "sig2coverage", "rescoverage", "phicoverage")
  #width
  print(width2<- round(colMeans(dat2[, wid.list]), digits = 2))

  #coverage
  print(coverage2 <- round(colMeans(dat2[, coverage.list])*100, digits = 1))
  #################################SSM#######################################
  #sdbias
  print(sdbias3 <- round(colMeans(bias3/biasd3)*100, digits = 1))
  #RMSE
  print(rmse3 <- round(sqrt(colMeans(bias3^2)), digits = 2))

  
  wid.list <- c("intwidth","x1width", "x2width", "Vx2width", "sig0width", "phiwidth")
  coverage.list <- c("intcoverage","x1coverage", "x2coverage", "Vx2coverage", "sig0coverage", "phicoverage")
  #width
  print(width3<- round(colMeans(dat3[, wid.list]), digits = 2))

  #coverage
  print(coverage3 <- round(colMeans(dat3[, coverage.list])*100, digits = 1))
  #################################combine#######################################
  print(sdbias <- cbind(sdbias,cbind(c(sdbias1, NaN), sdbias2, sdbias3)))
  print(rmse <- cbind(rmse, cbind(c(rmse1, NaN), rmse2, rmse3)))
  print(width <- cbind(width, cbind(c(width1, NaN), width2, width3)))
  print(coverage <- cbind(coverage, cbind(c(coverage1, NaN), coverage2, coverage3)))
}
apply(sdbias, 1, FUN = function(x) paste0(x, collapse = "&"))
apply(rmse, 1, FUN = function(x) paste0(x, collapse = "&"))
apply(coverage, 1, FUN = function(x) paste0(x, collapse = "&"))
apply(width, 1, FUN = function(x) paste0(x, collapse = "&"))

