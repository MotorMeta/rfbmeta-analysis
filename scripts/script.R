library(dplyr)
library(metafor)
library(tidyverse)
library(compute.es)
library(zcurve)
library(clubSandwich)
library(RoBMA)





#read in dataset
dat  <- read_csv("data/dat_rob.csv")

# turn effect size data into escal object

dat <- escalc(measure = "SMD", yi = g, vi = var, n1i = n2, n2i = n100,  data = dat)





#-------------------------------------------------------------




#### create column of p-values for z-curve analysis

mes <- mes(m.1 = mean100, m.2 = mean2, sd.1 = sd100, sd.2 = sd2, n.1 = n100, n.2 = n2, dig = 9, id = effect, data = dat)

p1 <- dplyr::select(mes, pval.d)

mes2 <- mes(m.1 = mean100, m.2 = mean2, sd.1 = calsd100, sd.2 = calsd2, n.1 = n100, n.2 = n2, dig = 9, id = effect, data = dat)

p2 <- dplyr::select(mes2, pval.d)

tes <- tes(t = t, n.1 = n100, n.2 = n2, dig = 9, id = effect, data = dat)

p3 <- dplyr::select(tes, pval.d)

fes <- fes(f = F, n.1 = n100, n.2 = n2, dig = 9, id = effect, data = dat)

p4 <- dplyr::select(fes, pval.d)

afes <- a.fes(f = af, n.1 = n100, n.2 = n2, R = r, q= q, dig = 9, id = effect, data = dat)

p5 <- dplyr::select(afes, pval.d)

p <- coalesce(p1, p2, p3, p4, p5)

dat$pval <- coalesce(p1, p2, p3, p4, p5)



######------------------------------------------------------------


# select only primary analyses
prime <- filter(dat, dat$primary=="Yes")


## calculate total number of participants from studies in primary analysis

# create one row per study so participants are not double counted
 
a <- prime[!duplicated(prime$study),]

# add groups todether

primeN <- a$n100 + a$n2

# calculate the total N

sum(primeN, na.rm = TRUE)



#### univariate analyses at each time point for selection model fitting

# select only acquisition tests

acq <- filter(prime, prime$time=="Acquisition")

# create random effects meta-analysis object

acqres <- rma(yi, vi, test = "knha",  data = acq)



# fit selection model to acquisition data

selacq <- selmodel(acqres, type="stepfun", alternative="greater", steps=c(.025, 1))


summary(selacq)


# view random effects analysis

summary(acqres)
predict.rma(acqres)


# select only immediate retention tests

imm <- filter(prime, prime$time=="ImmRet")

# create random effects meta-analysis object

immres <- rma(yi, vi, data = imm)


# fit selection model to immediate retention data

immsel <- selmodel(immres, type = "stepfun", alternative = "two.sided", steps = c(.05, 1))

summary(immsel)


# view random effects analysis

summary(immres)


# select only delayed retention tests

ret <- filter(prime, prime$time=="DelayRet")

# create random effects meta-analysis object

retres <- rma(yi, vi, n1i = n2, n2i = n100, data = ret)


# fit selection model to delayed retention data

seldel <- selmodel(retres, type = "stepfun", alternative = "greater", steps = c(.025, 1))


summary(seldel)


## conduct robust baysian meta-analysis on delayed retention data


ret$se <- sqrt(ret$vi)

x <- ret$yi

y <- ret$se

x <- na.omit(x)

y <- na.omit(y)

## fit <- RoBMA(d = x, se = y, seed = 1)

## fits <- summary(fit)

## saveRDS(fits, file = "data/RoBMA.rds")

read_rds(file = "data/RoBMA.rds")


# create contour enhanced funnel plot of delayed retention data

funnel.rma(retres, refline = 0, level = c(90, 95, 99), shade = c("white", "gray55", 'gray75'), legend = TRUE)

#### fit z-curve to delayed retention data

retp <- filter(ret, ret$pval!="NA")

pp <- retp$pval

zedc <- zcurve(p = pp$pval.d)

summary(zedc)



# view random effects analysis

summary(retres)





#------------------------------------------------------------------
#
#
#
#  PRIMARY ANALYSES
#
#
#
#
#------------------------------------------------------------------


# fit multilevel model to acquisition, immediate retention, and delayed retention

# exclude transfer

mvm <- filter(prime, prime$time!="Transfer")

# rename delayed retention time point "retention" so time points are both chronological
# and alphabetical

mvm$time <- recode(mvm$time, DelayRet = "Retention")

# order variables for subsequent forest plot 

mvm <- mvm[order(mvm$yi),]
mvm <- mvm[order(mvm$time),]

test <- rma.mv(yi, vi, mods = ~factor(time)-1, random = ~ 1| study/time, 
               data = mvm)
test

robust(test, mvm$study)


predict(test)

#profile.rma.mv(test)



## calculate cooks distance with 3 level model

#cd <- cooks.distance(test, parallel = "snow", ncpus = 4)

#plot(cd)

#View(cd)

# remove influential case

mvm2 <- mvm[-c(34),]

# fit model with influential case removed

test2<- rma.mv(yi, vi, mods = ~factor(time)-1, random = ~ 1| study/time, 
               struct = "CS", data = mvm2)


test2

robust(test2, cluster = mvm2$study)

# create and save forest plot
# order by time point

mvm <- mvm[order(mvm$yi),]
mvm <- mvm[order(mvm$time),]

library(extrafont)
loadfonts()

metafor::forest(test, header = "Author, Year, and Time Point", cex= 0.7, cex.lab= .8,
                addfit= FALSE, slab = paste0(mvm$authors, "," , mvm$year, ":" , mvm$time), 
                xlim = c(-5.5,5), at= c(-3, -2, -1, 0, 1, 2, 3), rows=c(2:35,39:106,110:180), fonts = "Roboto Condensed")
text(2.5, -1, "Favours Reduced", cex = .8, font = 2, fonts = "Roboto Condensed")
text(-2.5, -1, "Favours 100", cex = .8, font = 2, fonts = "Roboto Condensed")
addpoly.default(x = .1932, ci.lb = -.1062, ci.ub = .4926, pi.lb = -1.0281, pi.ub = 1.4146,
                row = 0, cex = .8, efac = .7,  mlab = "Acquisition Estimate", addpred = TRUE, fonts = "Roboto Condensed")
addpoly.default(x=.0137, ci.lb = -.3006, ci.ub = .3281, pi.lb = -1.2113, pi.ub = 1.2388,
                row = 22, cex = .8, efac = .7, mlab = "Immediate Retention Estimate", addpred = TRUE, fonts = "Roboto Condensed")
addpoly.default(x=.1924, ci.lb = -.0453, ci.ub = .4302, pi.lb = -1.0153, pi.ub = 1.4002,
                row = 42, cex = .8, efac = .7, mlab = "Delayed Retention Estimate", addpred = TRUE, fonts = "Roboto Condensed")




######## test moderators#######---------------------------



###### test risk of bias moderators ######

# allocation concealment

alloc <- rma.mv(yi, vi, mods = ~ alloc, random = ~ 1 | study/time, 
              data = mvm)

alloc

# incomplete outcome reporting

incout <- rma.mv(yi, vi, mods = ~ incout, random = ~ 1 | study/time, 
                 data = mvm)

incout


# selective outcome reporting

select <- rma.mv(yi, vi, mods = ~ select, random = ~ 1 | study/time, 
                 data = mvm)


select


# sequence generation 

seqgen <- rma.mv(yi, vi, mods = ~ seqgen, random = ~ 1 | study/time, 
                 data = mvm)


seqgen

# blind assessment

blindass <- rma.mv(yi, vi, mods = ~ blindass, random = ~ 1 | study/time, 
                   data = mvm)

blindass

# double blind 

blindpers <- rma.mv(yi, vi, mods = ~ blindpers, random = ~ 1 | study/time, 
                    data = mvm)

blindpers


# other

other <- rma.mv(yi, vi, mods = ~ other, random = ~ 1 | study/time, 
                data = mvm)

other


###### risk of bias summary figure

data_rob1 <- read.csv("data/robdat.csv")

library(RColorBrewer)
robvis::rob_summary(data = data_rob1, tool = "ROB1", overall= FALSE, weighted= FALSE, colour = c("#1B9E77","#D95F02","#666666"))




#### prespecified moderators


# test age moderator

age <- rma.mv(yi, vi, mods = ~ time * age, random = ~ 1 | study/time, 
              data = mvm)

summary(age)



# sensitivity analysis with outlier removed

age2 <- rma.mv(yi, vi, mods = ~ time * age, random = ~ 1 | study/time, 
               data = mvm2)


summary(age2)

# test skill level moderator

skill <- rma.mv(yi, vi, mods = ~ time * skill, random = ~ 1 | study/time, 
                data = mvm)

summary(skill)

# sensitivity analysis with outlier removed

skill2 <- rma.mv(yi, vi, mods = ~ time * skill, random = ~ 1 | study/time, 
                 data = mvm2)


summary(skill2)

# test task moderator

task <- rma.mv(yi, vi, mods = ~ time * task, random = ~ 1 | study/time, 
               data = mvm)
task


## sensitivity analysis with outlier removed

mvm3 <- mvm[-c(41),]


task2 <- rma.mv(yi, vi, mods = ~ time * task, random = ~ 1 | study/time, 
                data = mvm3)

summary(task2)




# test faded moderator

faded <- rma.mv(yi, vi, mods = ~ time * faded, random = ~ 1 | study/time, 
                data = mvm)

summary(faded)


# sensitivity analysis with outlier removed

faded2 <- rma.mv(yi, vi, mods = ~ time * faded, random = ~ 1 | study/time, 
                 data = mvm2)

summary(faded2)


# test the yoked moderator

yoked <- rma.mv(yi, vi, mods = ~ time * yoked, random = ~ 1 | study/time, 
                data = mvm)

summary(yoked)

# sensitivity analysis with outlier removed


yoked2 <- rma.mv(yi, vi, mods = ~ time * yoked, random = ~ 1 | study/time, 
                 data = mvm2)
summary(yoked2)


# test the feedback moderator

feedback <- rma.mv(yi, vi, mods = ~ time*feedback, random = ~ 1 | study/time, 
                   data = mvm)
summary(feedback)


# sensitivity analysis with outlier removed

feedback2 <- rma.mv(yi, vi, mods = ~ time + feedback, random = ~ 1 | study/time, 
                    data = mvm2)


summary(feedback2)


# test the measure moderator

measure <- rma.mv(yi, vi, mods = ~ time * measure, random = ~ 1 | study/time, 
                  data = mvm)

summary(measure)


# sensitivity analysis with outlier removed


measure2 <- rma.mv(yi, vi, mods = ~ time * measure, random = ~ 1 | study/time, 
                   data = mvm2)

summary(measure2)

# sensitivity analysis with outlier removed and no test interaction


measure3 <- rma.mv(yi, vi, mods = ~ measure, random = ~ 1 | study/time, 
                   data = mvm2)

summary(measure3)


# test bandwidth moderator

bprime <- filter(dat, dat$bprime=="Yes" & dat$time!="Transfer")


band <- rbind(mvm, bprime)


band$id <- paste0(band$authors, band$year, band$experiment, sep = ".")

bandwidth <- rma.mv(yi, vi, mods = ~ bandwidth*time, random = ~ 1 | id/time, 
                    data = band)

summary(bandwidth)

# sensitivity analysis with outlier removed


band2 <- rbind(mvm2, bprime)

band2$id <- paste0(band2$authors, band2$year, band2$experiment, sep = ".")

bandwidth2 <- rma.mv(yi, vi, mods = ~ bandwidth*time, random = ~ 1 | id/time, 
                    data = band2)

summary(bandwidth2)


### continuous moderators

# center variables at minimum amounts

# find minimums

min(mvm$trials) #minimum number of trials in dataset was 5
min(mvm$days) # minimum number of days was 1
min(mvm$frequency, na.rm = TRUE) #minimum frequency was 10
min(mvm$irinterval, na.rm = TRUE) #minimum immediate retention interval was 1 min
min(mvm$drinterval, na.rm = TRUE) #minimum delayed retention interval was 1 day

#center variables

mvm$trials <- mvm$trials - 5
mvm$days <- mvm$days - 1
mvm$frequency <- mvm$frequency - 10
mvm$irinterval <- mvm$irinterval - 1
mvm$drinterval <- mvm$drinterval -1


# test trials moderator

trials <- rma.mv(yi, vi, mods = ~ trials:time, random = ~ 1 | study/time, 
                 data = mvm)

summary(trials)


# sensitivity analysis with outlier removed

trials2 <- rma.mv(yi, vi, mods = ~ trials:time, random = ~ 1 | study/time, 
                  data = mvm2)

summary(trials2)

# test days moderator


days <- rma.mv(yi, vi, mods = ~ days:time, random = ~ 1 | study/time, 
               data = mvm)

summary(days)


# sensitivity analysis with outlier removed


days2 <- rma.mv(yi, vi, mods = ~ time:days, random = ~ 1 | study/time, 
                data = mvm2)

summary(days2)


# test frequency moderator (of primary analysis, overall frequency tested in 
# sensitivity analyses below)

mvm$frequency <- as.numeric(mvm$frequency)

freq <- rma.mv(yi, vi, mods = ~ time:frequency, random = ~ 1 | study/time, 
               data = mvm)

summary(freq)


# sensitivity analysis with outlier removed

mvm2$frequency <- as.numeric(mvm2$frequency)

freq2 <- rma.mv(yi, vi, mods = ~ time:frequency, random = ~ 1 | study/time, 
                data = mvm2)



summary(freq2)


# test immediate retention interval moderator on imm and del retention


irint <- rma.mv(yi, vi, mods = ~ time:irinterval, random = ~ 1 | study/time, 
                data = mvm, btt = 3:4)


summary(irint)


#sensitivity analysis with outlier removed

irint2 <- rma.mv(yi, vi, mods = ~ time:irinterval, random = ~ 1 | study/time, 
                 data = mvm2, btt = 3:4)


summary(irint2)


# test delayed retention interval on delayed retention data

drint <- rma.mv(yi, vi, mods = ~ drinterval, random = ~ 1 | study,
                data = ret)

summary(drint)




####################### test time X frequency interactions ###-----------------


## acquisition vs delayed retention

# test for difference between time points (acquisition, delayed retention)


acq_ret <- rma.mv(yi, vi, mods = ~factor(time), random = ~ 1| study/time, 
                  data = mvm, btt = 2)


# view results

summary(acq_ret)

robust(acq_ret, cluster = mvm$study)

# conduct sensitivity analysis with outlier removed

acq_ret2 <- rma.mv(yi, vi, mods = ~factor(time), random = ~ 1| study/time, 
                  data = mvm2, btt = 2)

summary(acq_ret2)

robust(acq_ret2, mvm2$study)

## immediate retention vs delayed retention

# set immediate retention as reference level and test time (immediate vs. delayed)

imm_ret <- rma.mv(yi, vi, mods = ~relevel(factor(time), ref = "ImmRet"), random = ~ 1| study/time, 
                  data = mvm, btt = 3)


# view results

summary(imm_ret)

robust(imm_ret, mvm$study)

# conduct sensitivity analysis with outlier removed

imm_ret2 <- rma.mv(yi, vi, mods = ~relevel(factor(time), ref = "ImmRet"), random = ~ 1| study/time, 
                   data = mvm2, btt = 3)

# view results 

summary(imm_ret2)

robust(imm_ret2, mvm2$study)

#------------------------------------------------------------------------
#
#
#
#  sensitivity analyses
#
#
#
#------------------------------------------------------------------------

### sensitivity analyses


## cluster-robust inference method 

robust(test, cluster = mvm$study)


## construct approximate V matrix assuming r = .7 based on McKay & Ste-Marie, 2020.

V <- impute_covariance_matrix(mvm$vi, cluster=mvm$study, r=0.7)

## fit working model

rob1 <- rma.mv(yi, V, mods = ~ time-1,  random = ~factor(time) | study, struct = "CS", data = mvm)

## use cluster-robust inference methods

robust(rob1, cluster = mvm$study)

coef_test(rob1, vcov = "CR2", cluster=mvm$study)

conf_int(rob1, vcov = "CR2", cluster=mvm$study)

## test time moderator by including intercept

## fit working model (relevel time factor so that delayed retention is the reference level)

rob1 <- rma.mv(yi, V, mods = ~relevel(factor(time), ref="DelayRet"),  random = ~factor(time) | study, struct = "CS", data = mvm)

## use cluster-robust inference methods

robust(rob1, cluster = mvm$study)

coef_test(rob1, vcov = "CR2", cluster=mvm$study)

conf_int(rob1, vcov = "CR2", cluster=mvm$study)





## construct approximate V matrix assuming r = .5 

V <- impute_covariance_matrix(mvm$vi, cluster=mvm$study, r=0.5)

## fit working model

rob15 <- rma.mv(yi, V, mods = ~ time-1,  random = ~factor(time) | study, struct = "CS", data = mvm)

## use cluster-robust inference methods

robust(rob15, cluster = mvm$study)

coef_test(rob15, vcov = "CR2", cluster=mvm$study)

conf_int(rob15, vcov = "CR2", cluster=mvm$study)


## test time moderator assuming r = 0.5

rob15 <- rma.mv(yi, V, mods = ~ relevel(factor(time), ref = "DelayRet"),  random = ~factor(time) | study, struct = "CS", data = mvm)

## use cluster-robust inference methods

robust(rob15, cluster = mvm$study)

coef_test(rob15, vcov = "CR2", cluster=mvm$study)

conf_int(rob15, vcov = "CR2", cluster=mvm$study)



## construct approximate V matrix assuming r = .9 

V <- impute_covariance_matrix(mvm$vi, cluster=mvm$study, r=0.9)

## fit working model

rob19 <- rma.mv(yi, V, mods = ~ time-1,  random = ~factor(time) | study, struct = "CS", data = mvm)

## use cluster-robust inference methods

robust(rob19, cluster = mvm$study)

coef_test(rob19, vcov = "CR2", cluster=mvm$study)

conf_int(rob19, vcov = "CR2", cluster=mvm$study)



## test time moderator with r = .9

rob19 <- rma.mv(yi, V, mods = ~ relevel(factor(time), ref = "DelayRet"),  random = ~factor(time) | study, struct = "CS", data = mvm)

## use cluster-robust inference methods

robust(rob19, cluster = mvm$study)

coef_test(rob19, vcov = "CR2", cluster=mvm$study)

conf_int(rob19, vcov = "CR2", cluster=mvm$study)




## test time moderator with extreme case where approximate V matrix  r = -.9 

V <- impute_covariance_matrix(mvm$vi, cluster=mvm$study, r= -0.9)

## fit working model

robneg <- rma.mv(yi, V, mods = ~ relevel(factor(time), ref = "DelayRet"),  random = ~factor(time) | study, struct = "CS", data = mvm)

## use cluster-robust inference methods

robust(robneg, cluster = mvm$study)

coef_test(robneg, vcov = "CR2", cluster=mvm$study)

conf_int(robneg, vcov = "CR2", cluster=mvm$study)




#
# Results not sensitive to approximate range of values for V matrix
#


## cluster-robust inference method with influential case removed

robust(test2, cluster = mvm2$study)


## construct approximate V matrix assuming r = .7 based on McKay & Ste-Marie, 2020.

V <- impute_covariance_matrix(mvm2$vi, cluster=mvm2$study, r=0.7)

# fit working model with influential case removed

rob <- rma.mv(yi, V, mods = ~ time-1, random = ~factor(time)-1 | study, data = mvm2)

# use cluster robust inference methods

robust(rob, cluster = mvm2$study)

coef_test(rob, vcov = "CR2", cluster=mvm2$study)

conf_int(rob, vcov = "CR2", cluster=mvm2$study)



###### Testing the spatial error subset as requested by reviewer

# recode all spatial error outcomes into one category

mvm$rmeasure <- recode(mvm$measure, AE = "serror", RMSE = "serror", ACE = "serror", 
                        E = "serror")

# select only spatial error

serror <- filter(mvm, mvm$rmeasure=="serror")

# select only delayed retention

sret <- filter(serror, serror$time=="Retention")

# test moderators

# age

rma(yi, vi, mods = ~factor(age), data = sret)

# skill

rma(yi, vi, mods = ~factor(skill), data = sret)

# task

rma(yi, vi, mods = ~factor(task), data = sret)

# faded

rma(yi, vi, mods = ~factor(faded), data = sret)

# yoked

rma(yi, vi, mods = ~factor(yoked), data = sret)

# feedback

rma(yi, vi, mods = ~factor(feedback), data = sret)

# trials

rma(yi, vi, mods = ~ trials, data = sret)

# days

rma(yi, vi, mods = ~days, data = sret)

# frequency 

rma(yi, vi, mods = ~frequency, data = sret)



######  with measure nested in time nested in experiment


# remove transfer data

ntra <- filter(dat, dat$time!="Transfer")

# create study variable

ntra$author.year.study <- paste(ntra$authors, ntra$year, ntra$experiment, sep = ".")


# fit four level model

sens <- rma.mv(yi, vi, mods = ~ factor(time)-1,  random = ~ 1|author.year.study/measure/time, data = ntra)

sens


# run profile analysis to check for overparamterization

profile(sens)

# check for influential cases with cooks distance

cdtest <- cooks.distance.rma.mv(sens)

View(cdtest)

plot(cdtest)

# remove outlier

ntra2 <- ntra[-c(34),]


# refit 4 level model with outlier removed

sens2 <- rma.mv(yi, vi, mods = ~ factor(time)-1,  random = ~ 1|author.year.study/measure/time, data = ntra2)


sens2

##### moderator analyses with 4 level model


# age

sensage <- rma.mv(yi, vi, mods = ~ time*age,  random = ~ 1|author.year.study/measure/time, data = ntra)

summary(sensage)

# with influential case removed

sensage2 <- rma.mv(yi, vi, mods = ~ time*age,  random = ~ 1|author.year.study/measure/time, data = ntra2)

summary(sensage2)

# skill


sensskill <- rma.mv(yi, vi, mods = ~ time*skill,  random = ~ 1|author.year.study/measure/time, data = ntra)

summary(sensskill)

# with influential case removed


sensskill2 <- rma.mv(yi, vi, mods = ~ time*skill,  random = ~ 1|author.year.study/measure/time, data = ntra2)

summary(sensskill2)


# task


senstask <- rma.mv(yi, vi, mods = ~ time*task,  random = ~ 1|author.year.study/measure/time, data = ntra)

summary(senstask)

# with influential case removed


senstask2 <- rma.mv(yi, vi, mods = ~ time*task,  random = ~ 1|author.year.study/measure/time, data = ntra2)

summary(senstask2)


# with Drews et al. -2.14 effect removed

ntra3 <- ntra[-c(41),]

senstask3 <- rma.mv(yi, vi, mods = ~ time*task,  random = ~ 1|author.year.study/measure/time, data = ntra3)

summary(senstask3)


# bandwidth

sensband <- rma.mv(yi, vi, mods = ~ time*bandwidth,  random = ~ 1|author.year.study/measure/time, data = ntra)

summary(sensband)


# with influential case removed


sensband2 <- rma.mv(yi, vi, mods = ~ time*bandwidth,  random = ~ 1|author.year.study/measure/time, data = ntra2)

summary(sensband2)


# faded

sensfaded <- rma.mv(yi, vi, mods = ~ time*faded,  random = ~ 1|author.year.study/measure/time, data = ntra)

summary(sensfaded)

# with influential case removed

sensfaded2 <- rma.mv(yi, vi, mods = ~ time*faded,  random = ~ 1|author.year.study/measure/time, data = ntra2)

summary(sensfaded2)

# yoked

sensyoked <- rma.mv(yi, vi, mods = ~ time*yoked,  random = ~ 1|author.year.study/measure/time, data = ntra)

summary(sensyoked)

# with influential case removed

sensyoked2 <- rma.mv(yi, vi, mods = ~ time*yoked,  random = ~ 1|author.year.study/measure/time, data = ntra2)

summary(sensyoked2)


# feedback

sensfeedback <- rma.mv(yi, vi, mods = ~ time*feedback,  random = ~ 1|author.year.study/measure/time, data = ntra)

summary(sensfeedback)

# with influential case removed

sensfeedback2 <- rma.mv(yi, vi, mods = ~ time*feedback,  random = ~ 1|author.year.study/measure/time, data = ntra2)

summary(sensfeedback2)

# feedback no interaction

sensfeedback3 <- rma.mv(yi, vi, mods = ~ feedback,  random = ~ 1|author.year.study/measure/time, data = ntra)

summary(sensfeedback3)

# measure

sensmeasure <- rma.mv(yi, vi, mods = ~ time*measure,  random = ~ 1|author.year.study/measure/time, data = ntra)

summary(sensmeasure)

# collapse measurements into broader categories as suggested by reviewer

ntra$rmeasure <- recode(ntra$measure, AE = "serror", RMSE = "serror", ACE = "serror", 
                         E = "serror")

ntra$rmeasure <- recode(ntra$rmeasure, "ABS time" = "terror", "REL time" = "terror")

ntra$rmeasure <- recode(ntra$rmeasure, Other = "OTHER")

sensr <- rma.mv(yi, vi, mods = ~ time*rmeasure,  random = ~ 1|author.year.study/rmeasure/time, data = ntra)

sensr


# with influential case removed

sensmeasure2 <- rma.mv(yi, vi, mods = ~ time*measure,  random = ~ 1|author.year.study/measure/time, data = ntra2)

summary(sensmeasure2)

# collapse measurements and with outlier removed

ntra2$rmeasure <- recode(ntra2$measure, AE = "serror", RMSE = "serror", ACE = "serror", 
                         E = "serror")

ntra2$rmeasure <- recode(ntra2$rmeasure, "ABS time" = "terror", "REL time" = "terror")

ntra2$rmeasure <- recode(ntra2$rmeasure, Other = "OTHER")


sens2r <- rma.mv(yi, vi, mods = ~ time*rmeasure,  random = ~ 1|author.year.study/rmeasure/time, data = ntra2)

sens2r

# measure with interaction removed

sensmeasure3 <- rma.mv(yi, vi, mods = ~ measure,  random = ~ 1|author.year.study/measure/time, data = ntra)

summary(sensmeasure3)

# with measures collapsed and interaction removed

sensr3 <- rma.mv(yi, vi, mods = ~ factor(rmeasure),  random = ~ 1|author.year.study/rmeasure/time, data = ntra)

sensr3


# just at delayed retention

rret <- filter(ntra, ntra$time=="DelayRet")

rma.mv(yi, vi, mods = ~ factor(rmeasure),  random = ~ 1|author.year.study/rmeasure, data = ntra)


# test whether primary coded measures differ from secondary measures

sensprimary <- rma.mv(yi, vi, mods = ~ time*primary,  random = ~ 1|author.year.study/measure/time, data = ntra)

summary(sensprimary)



# with influential case removed

sensprimary2 <- rma.mv(yi, vi, mods = ~ time*primary,  random = ~ 1|author.year.study/measure/time, data = ntra2)

summary(sensprimary2)

### meta-regressions

# trials

senstrials <- rma.mv(yi, vi, mods = ~ trials*time,  random = ~ 1|author.year.study/measure/time, data = ntra)

summary(senstrials)

# with influential case removed

senstrials2 <- rma.mv(yi, vi, mods = ~ trials*time,  random = ~ 1|author.year.study/measure/time, data = ntra2)

summary(senstrials2)

# days

sensdays <- rma.mv(yi, vi, mods = ~ days*time,  random = ~ 1|author.year.study/measure/time, data = ntra)

summary(sensdays)

# with influential case removed


sensdays2 <- rma.mv(yi, vi, mods = ~ days*time,  random = ~ 1|author.year.study/measure/time, data = ntra2)

summary(sensdays2)

######   frequency; this is main test of frequency moderator overall as it accounts for all frequency groups

####### frequency moderator analysis with recalculated effects to account for shared 100% comparison group

# select only effects from studies that compare multiple frequencies

freqgrp <- filter(dat, dat$freqgrps=="Yes")

# divide the 100% group sample size by number of comparisons made to other frequencies

freqgrp$n100 <- freqgrp$n100/freqgrp$frgrn

# calculate effect sizes and variances with adjusted sample sizes

gs <- compute.es::mes(mean100, mean2, sd100, sd2, n100, n2, id = effect, data = freqgrp)

# create a vector of hedges g values

g <- dplyr::pull(gs, g)

# create a vector of variances

var <- dplyr::pull(gs, var.g)

# change the g and variance values in the freqgrp df to the newly computed vectors

freqgrp$g <- g

freqgrp$var <- var

# create dataframe of effects not from studies that compared frequencies

nfreqgrp <- filter(dat, dat$freqgrps!="Yes")

# combine the effects from studies comparing frequencies with those that do not

freqcomp <- rbind(nfreqgrp, freqgrp)

# create author.year.study variable

freqcomp$author.year.study <- paste(freqcomp$authors, freqcomp$year, freqcomp$experiment, sep = ".")

# fit multilevel model with frequency X time moderator

frequencycomp <- rma.mv(g, var, mods = ~ frequency*time,  random = ~ 1|author.year.study/measure/time, data = freqcomp)

summary(frequencycomp)





#-----------------------------------------------------------------------------

# additional sensitivity analyses

#### attempting to account for experiments nested in articles produces overparameterized model 

# create article.study variable

mvm$article.study <- paste(mvm$article, mvm$study, sep = ".")

# estimate effect of reduced feedback frequency at each time point

pres <- rma.mv(yi, vi, mods = ~factor(time)-1, random = list(~ 1 | article/study, ~ time | article.study), 
               struct = "CS", data = mvm)


profile.rma.mv(pres)



#----------------------------------------------------------
#
#  supplemental analysis - transfer data not intended for main
#  analyses because transfer tests differ from acquisition in
#  idiosyncratic ways. data were collected for future research
#----------------------------------------------------------

# select only delayed transfer tests

tra <- filter(prime, prime$time=="Transfer")

# create random effects meta-analysis object
trares <- rma(g, var, data = tra)

# compute standardized residuals

rstudent(trares)

# create influence statistic object
inf4 <- influence(trares)

# plot influence stats
plot(inf4, dfb = TRUE)

# remove one outlier

tra2 <- tra[-c(8),]

# create new random effects meta-analysis object

trares2 <- rma(g,var, data = tra2)

# view first and second random effects analyses

summary(trares)
summary(trares2)


 