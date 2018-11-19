
########################################################################
### Example: Hodgkin disease
### Condcut two sample Maximum and Chisquare tests for CSH and ACH
### Produce confidence region figure
########################################################################

rm(list=ls(all=TRUE))
library(MASS)
library(survival)
library(cmprsk)
library(fields)
library(ellipse)

###Load dataset hodgkin.txt
hd_all <- read.table("~/Dropbox/Dissertation/first project/realdata/Hodgkin/hodgkin.txt", sep = ",", header = TRUE)
#hd_all <- read.table(".\\hodgkin.txt", sep = ",", header = TRUE)

###Cleaning dataset and create necessary variables

# dfcens: Censoring variable:0=Censored, 1=Failure
# mcens: Second malignancy: 0= No, 1= Yes

# time to the first event
hd_all$alltime=apply(cbind(hd_all$dftime,hd_all$maltime),1,min)

m1cens=(hd_all$alltime==hd_all$maltime & hd_all$mcens==1)+0 
cr1cens=(hd_all$alltime==hd_all$dftime & hd_all$dfcens==1)+0

## event of interest: hd$maltime is time to second malignancy
## competing risks: hd$dftime is time to relapse or death
## censor: end of followup

# cause of event 0=end of followup, 1=second malignancy, 2= death or relapse
hd_all$cens=m1cens+2*cr1cens 

# event of interest, second malignancy 1= yes  0= no
hd_all$m1cens=(hd_all$alltime==hd_all$maltime & hd_all$mcens==1)+0 
#hd$cr1cens=(alltime==hd$dftime & hd$dfcens==1)+0

# Define eoi and failure variable to fit CSH model and ACH model
hd_all$eoi=hd_all$cens
hd_all$failure=(hd_all$cens==1 | hd_all$cens==2)

# group variable
hd_all$x=(hd_all$age>=30)+0

####################################################
###  Conduct joint test#############################
####################################################

teststat=matrix(NA,nrow=1,ncol=5)
pvalue=matrix(NA,nrow=1,ncol=5)

  ### Test 1: Wald test for group coefficient in cox model for cause specific hazard for event of interest### 
  sdfcox1 <- coxph(Surv(hd_all$alltime,hd_all$eoi==1) ~ hd_all$x) # CSH model
  teststat[1]=(sdfcox1$coef/sqrt(sdfcox1$var))^2  
  pvalue[1]=1-pchisq(teststat[1],1)

  ### Test 2: Wald test for group coefficient in cox model for all cause hazard### 
  sdfcox_all <- coxph(Surv(hd_all$alltime,hd_all$failure==1) ~ hd_all$x) # ACH model
  teststat[2]=(sdfcox_all$coef/sqrt(sdfcox_all$var))^2  
  pvalue[2]=1-pchisq(teststat[2],1)

  ### Test 3: bonferroni adjustment for test 1 and 2 ###
  pvalue[3]=min(pvalue[1],pvalue[2])*2

  ### Test 4: Joint logrank test statistics ###

  var_1=sdfcox1$var[1,1] 		#variance of beta_j
  var_2=sdfcox_all$var[1,1]		#variance of beta_\dot
  var_12=var_2

  corr<-var_12/sqrt(var_1*var_2)

  W_1=sdfcox1$coef/sqrt(var_1)
  W_2=sdfcox_all$coef/sqrt(var_2)

  W=matrix(c(W_1,W_2),nrow=1)
  cov=matrix(c(1,corr,corr,1),nrow=2) #covariance matrix
  if (det(cov)==0) {browser()}
  teststat[4]=W%*%solve(cov)%*%t(W)   #chisquare test stat
  pvalue[4]=1-pchisq(teststat[4],2)

  ### maximum test ###
  teststat[5]=max(abs(sdfcox1$coef[1]/sqrt(sdfcox1$var[1,1])), abs(sdfcox_all$coef[1]/sqrt(sdfcox_all$var[1,1])))
  nbootstrap=100000000
  sample=mvrnorm(n = nbootstrap, c(0,0), Sigma=matrix(c(1,corr,corr,1),nrow=2), tol=1e-4)
  max=matrix(NA,nrow=nbootstrap,ncol=1)
  for (j in 1:nbootstrap){max[j]=max(abs(sample[j,]))}
  index<-seq(1,nbootstrap)
  index.a<-index[sort(max)> teststat[5]]
  pvalue[5]=1-index.a[1]/nbootstrap
  colnames(pvalue)=c("CSH","ACH","Bonf","Chi","Max")
  pvalue

####################################################
### produce confidence region plots for Hodgkin data 
####################################################

mean1=sdfcox1$coef[1]
se1=sqrt(sdfcox1$var[1,1])

mean2=sdfcox_all$coef[1]
se2=sqrt(sdfcox_all$var[1,1])
corr=corr

max_quantile=quantile(max,0.95)

plot(c(exp(mean1-3*se1),exp(mean1+3*se1)),c(exp(mean2-3*se2),exp(mean2+3*se2)), type='n',xlab="",ylab="",axes=F)
lines(exp(ellipse(corr,scale=c(se1,se2), centre=c(mean1,mean2),level=0.95)),type='l',lty=4)
rect(exp(mean1-2.24*se1),exp(mean2-2.24*se2),exp(mean1+2.24*se1),exp(mean2+2.24*se2), lty=2)
rect(exp(mean1-max_quantile*se1),exp(mean2-max_quantile*se2),exp(mean1+max_quantile*se1),exp(mean2+max_quantile*se2),lty=1)
points (exp(mean1), exp(mean2),pch=18)
text(exp(mean1), exp(mean2)+0.05,"(1.55, 1.67)")

axis(side=1)
axis(side=2)
#title(main="(b)")
legend(list(x=2.0,y=1.5),c("Chisq", "Bonf","Max"), cex=0.8, lty = c(4,2,1))
mtext("CSH ratio", side=1, line=2.4, cex=1)
mtext("ACH ratio", side=2, line=2.5, cex=1)



