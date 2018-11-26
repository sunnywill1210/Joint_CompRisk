rm(list=ls(all=TRUE))
library(MASS)
library(survival)
library(cmprsk)
#setwd('~/Dropbox/Dissertation/first project/real data/Hodgkin')

hd_all <- read.table("~/Dropbox/Dissertation/first project/realdata/Hodgkin/hodgkin.txt", sep = ",", header = TRUE)


## event of interest: hd$maltime is time to second malignancy
## competing risks: hd$dftime is time to relapse or death
## censor: end of followup

# dfcens: Censoring variable:0=Censored, 1=Failure
# mcens: Second malignancy: 0= No, 1= Yes

# time to the first event
hd_all$alltime=apply(cbind(hd_all$dftime,hd_all$maltime),1,min)

m1cens=(hd_all$alltime==hd_all$maltime & hd_all$mcens==1)+0 # 0: 
cr1cens=(hd_all$alltime==hd_all$dftime & hd_all$dfcens==1)+0

# cause of event 0=end of followup, 1=second malignancy, 2= death or relapse
hd_all$cens=m1cens+2*cr1cens 

# cause of event 0=end of followup or relapse, 1=second malignancy, 2= death
hd_all$cens2=hd_all$mcens+2*(hd_all$mcens==0 & hd_all$stat==1)

# event of interest, second malignancy 1= yes  0= no
hd_all$m1cens=(hd_all$alltime==hd_all$maltime & hd_all$mcens==1)+0 
#hd$cr1cens=(alltime==hd$dftime & hd$dfcens==1)+0

# group variable
hd_all$x=(hd_all$age>=30)+0

#Definition of Model A-D are in pintilie's book P 156
crr(hd_all$alltime,hd_all$cens,hd_all$x) # Model A death or relapse as competing risk
coxph(Surv(hd_all$alltime,hd_all$m1cens) ~ hd_all$x) # Model B death or relapse as competing risk

#crr(hd_all$maltime,hd_all$cens2,hd_all$x)  # Model C (the model used in paper for testing CI) time to death is the only competing risk
#coxph(Surv(hd_all$maltime,hd_all$mcens) ~ hd_all$x) # Model D (the model used in paper for testing CSH) time to death is the only competing risk

### Joint Test for model C and D
#time=hd_all$maltime
#event=hd_all$mcens
#group=hd_all$x
#cens=hd_all$dftime
#type=hd_all$cens2

### Joint Test for model A and B
time=hd_all$alltime
event=hd_all$m1cens
group=hd_all$x
cens=hd_all$dftime
type=hd_all$cens


### plot of cumulative cause specific function between two age groups

fit1 <- survfit(Surv(time,event==1) ~ group, data = hd_all) 

#quartz(height=4,width=6)
plot(seq(0,36, by=6),seq(0,0.6,by=0.1),type="n",xlab="Time to second malignancy",ylab="",axe=F)
lines(fit1, type="s", mark.time=FALSE, lty = 1:2, fun="cumhaz")
title("(a)")
axis(side=1)
axis(side=2, at=seq(0,0.6,length=7),las=1)
mtext("Cumulative cause specific ", side=2, line=3, cex=.8)
mtext("hazard", side=2, line=2,cex=0.8)
legend(list(x=0,y=0.5),cex=0.9,c("Age < 30 years","Age >= 30 years"),
       lty = c(1,2))
box()



fit2 <- cuminc(time, type, group)
#quartz(height=4,width=6)
plot(seq(0,36, by=6),seq(0,0.6,by=0.1),type="n",xlab="Time to second malignancy",ylab="",axe=F)
lines(fit2$'0 1'$time,fit2$'0 1'$est,type='l',lty=1)
lines(fit2$'1 1'$time,fit2$'1 1'$est,type='l',lty=2)
title("(b)")
axis(side=1)
axis(side=2, at=seq(0,0.6,length=7),las=1)
#mtext("Time to second malignancy", side=1, line=2.7, cex=.8)
mtext("Cumulative incidence of ", side=2, line=3, cex=.8)
mtext("Second malignancy", side=2, line=2,cex=0.8)
legend(list(x=0,y=0.5),cex=0.9,c("Age < 30 years","Age >= 30 years"),
       lty = c(1,2))
box()



