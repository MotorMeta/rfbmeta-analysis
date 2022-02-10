# # # # # # # # # # # # # # # # # # # # # # # # # # #
#Meta-analysis of dependent effect sizes simulation #
# # # # # # # # # # # # # # # # # # # # # # # # # # #
library(MBESS)
library(metafor)
library(faux)
library(tidyverse)
set.seed(1000)



# individual meta-analysis simulation function

onemeta <- function(nSims, cxa,exa, cxb, exb, rho, sample_lb, sample_ub){
  esd_1 <-numeric(nSims) #empty container for all simulated ES
  esd_2 <-numeric(nSims) #empty container for all simulated ES
  SSn1 <-numeric(nSims) #empty container for random sample sizes group 1
  SSn2 <-numeric(nSims) #empty container for random sample sizes group 2
  for(i in 1:nSims){ #for each simulated experiment
    SampleSize<-sample(sample_lb:sample_ub, 1) #randomly draw a sample between lb and ub
    ge <-rnorm_multi(n = SampleSize, 2, 0, 1, r = rho, varnames = c("t1", "t2")) #sample from a multivariate normal
    gc <-rnorm_multi(n = SampleSize, 2, 0, 1, r = rho, varnames = c("t1", "t2")) #sample from a multivariate normal
    gc$t1 <- gc$t1 + cxa # add effect to time 1 for control
    ge$t1 <- ge$t1 + exa # add effect to time 1 for experimental
    gc$t2 <- gc$t2 + cxb # add effect to time 2 for experimental
    ge$t2 <- ge$t2 + exb # add effect to time 2 for experimental
    SSn1[i]<-SampleSize #save sample size group 1
    SSn2[i]<-SampleSize #save sample size group 2
    esd_1[i]<-smd(Mean.1= mean(ge$t1), Mean.2=mean(gc$t1), s.1=sd(ge$t1), s.2=sd(gc$t1), n.1=SampleSize, n.2=SampleSize, Unbiased=TRUE) #Use MBESS to calc Hedges g
    esd_2[i] <-smd(Mean.1= mean(ge$t2), Mean.2=mean(gc$t2), s.1=sd(ge$t2), s.2=sd(gc$t2), n.1=SampleSize, n.2=SampleSize, Unbiased=TRUE)}
  #Insert effect sizes and sample sizes
  n1<-c(SSn1)
  n2<-c(SSn2)
  J<-1-3/(4*( SSn1+ SSn2-2)-1) #correction for bias
  esdv_1 <-(((SSn1+SSn2)/(SSn1*SSn2))+(esd_1^2/(2*(SSn1+SSn2))))*J^2
  esdv_2 <-(((SSn1+SSn2)/(SSn1*SSn2))+(esd_2^2/(2*(SSn1+SSn2))))*J^2 
  dat <- as.data.frame(cbind(esd_1, esdv_1, esd_2, esdv_2, n1, n2))
  dat$id <- seq.int(nrow(dat))
  dat<-dat %>% 
    pivot_longer(
      cols = !c(n1,n2,id),
      names_to = c(".value", "num"),
      names_sep = "_")
  return(as.data.frame(dat))
 
  }



# simulation study ## multilevel models

iters <- 100 #number of simulated meta-analyses
mpval <- rep(NA, iters) #empty container for moderator p-values 
difb <- rep(NA, iters) #empty container for moderator beta values
nSims <- 32 #number of experiments per meta-analysis
cxa <- 0 #intervention effect for control group at first time point
exa <- 0 #intervention effect for experimental group at first time point
cxb <- 0 #intervention effect for control group at second time point
exb <- 0 #intervention effect for experimental group at second time point
rho <- .5 #correlation between time points in the population
sample_lb <- 10 #lower bound for number of participants per experiment
sample_ub <- 30 #upper bound for number of participants per experiment




for (i in 1:iters){
  dat <-onemeta(nSims = nSims, cxa = cxa, exa = exa, cxb = cxb, exb = exb, rho = rho, sample_lb = sample_lb, sample_ub = sample_ub)
  tryCatch(
    {
      res <- rma.mv(esd, esdv, mods = ~factor(num), random = ~1|id/num, data = dat) #fit multilevel model with time points nested in experiments
      mpval[i] <- res$QMp
      difb[i] <- res$b[2]},
    error=function(error_message) {
      message(error_message)
      return(NA) #when the model fails to converge, store as NA rather than stopping the loop
    }
  )
}

mean(mpval <= .05, na.rm = TRUE) # mean rejection rate with alpha = .05
mean(difb, na.rm = TRUE) # mean estimate of difference in intervention effect between time points






# simulation study ## cluster robust methods

iters <- 100 #number of simulated meta-analyses
mpval <- rep(NA, iters) #empty container for moderator p-values 
difb <- rep(NA, iters) #empty container for moderator beta values
nSims <- 32 #number of experiments per meta-analysis
cxa <- 0 #intervention effect for control group at first time point
exa <- 0 #intervention effect for experimental group at first time point
cxb <- 0 #intervention effect for control group at second time point
exb <- 0 #intervention effect for experimental group at second time point
rho <- .5 #correlation between time points in the population
sample_lb <- 10 #lower bound for number of participants per experiment
sample_ub <- 30 #upper bound for number of participants per experiment



for (i in 1:iters){
  dat <-onemeta(nSims = nSims, cxa = cxa, exa = exa, cxb = cxb, exb = exb, rho = rho, sample_lb = sample_lb, sample_ub = sample_ub)
  tryCatch(
    {
  res <- rma.mv(esd, esdv, mods = ~factor(num), random = ~1|id/num, data = dat)
  rob <-robust(res, cluster = dat$id) #use cluster robust inference methods
  mpval[i] <- rob$QMp
  difb[i] <- rob$b[2]},
  error=function(error_message) {
    message(error_message)
    return(NA)
  }
  )
}

mean(mpval <= .05, na.rm = TRUE) # mean rejection rate with alpha = .05
mean(difb, na.rm = TRUE) # mean estimate of difference in intervention effect between time points

