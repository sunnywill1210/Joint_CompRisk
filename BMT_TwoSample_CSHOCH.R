########################################################################
### Example: BMT
### Condcut two sample Maximum and Chisquare tests for CSH and OCH
### Produce confidence region figure
########################################################################

rm(list=ls(all=TRUE))
library(MASS)
library(survival)
library(cmprsk)
library(ellipse)

#load data
source(".\\bmt_datacleaning.R")

#Specify which group variable to test
group=mini

teststat_CSHOCH=matrix(NA,nrow=1,ncol=5)
pvalue_CSHOCH=matrix(NA,nrow=1,ncol=5)

### Test 1: Wald test for group coefficient in cox model for cause specific hazard for event of interest### 
sdfcox1 <- coxph(Surv(time,type==1) ~ group) # CSH model
teststat_CSHOCH[1]=(sdfcox1$coef/sqrt(sdfcox1$var))^2  
pvalue_CSHOCH[1]=1-pchisq(teststat_CSHOCH[1],1)

### Test 2: Wald test for group coefficient in cox model for other cause hazard### 
sdfcox2 <- coxph(Surv(time,type==2) ~ group) # OCH model
teststat_CSHOCH[2]=(sdfcox2$coef/sqrt(sdfcox2$var))^2  
pvalue_CSHOCH[2]=1-pchisq(teststat_CSHOCH[2],1)

### Test 3: bonferroni adjustment for test 1 and 2 ###
pvalue_CSHOCH[3]=min(pvalue_CSHOCH[1],pvalue_CSHOCH[2])*2

### Test 4: Joint logrank test statistics ###

var_1=sdfcox1$var[1,1]   	#variance of beta_j
var_2=sdfcox2$var[1,1]		#variance of beta_\dot
var_12=0

corr<-var_12/sqrt(var_1*var_2)

W_1=sdfcox1$coef/sqrt(var_1)
W_2=sdfcox2$coef/sqrt(var_2)

W=matrix(c(W_1,W_2),nrow=1)
cov=matrix(c(1,corr,corr,1),nrow=2) #covariance matrix
if (det(cov)==0) {browser()}
teststat_CSHOCH[4]=W%*%solve(cov)%*%t(W)   #chisquare test stat
pvalue_CSHOCH[4]=1-pchisq(teststat_CSHOCH[4],2)

### maximum test ###
teststat_CSHOCH[5]=max(abs(W_1), abs(W_2))
nbootstrap=1000000
sample=mvrnorm(n = nbootstrap, c(0,0), Sigma=matrix(c(1,corr,corr,1),nrow=2), tol=1e-4)
max=matrix(NA,nrow=nbootstrap,ncol=1)
for (j in 1:nbootstrap){max[j]=max(abs(sample[j,]))}
pvalue_CSHOCH[5]=sum(max>teststat_CSHOCH[5])/nbootstrap

colnames(pvalue_CSHOCH)=c("CSH","OCH","Bonf","Chi","Max")

pvalue_CSHOCH

#####################################
### producing confidence region graph
####################################
mean1=sdfcox1$coef[1]
se1=sqrt(sdfcox1$var[1,1])

mean2=sdfcox2$coef[1]
se2=sqrt(sdfcox2$var[1,1])
corr=0

max_quantile=quantile(max,0.95)

plot(c(exp(mean1-3*se1),exp(mean1+3*se1)),c(exp(mean2-3*se2),exp(mean2+3*se2)), type='n',xlab="",ylab="",axes=F)
lines(exp(ellipse(corr,scale=c(se1,se2), centre=c(mean1,mean2),level=0.95)),type='l',lty=4)
rect(exp(mean1-2.24*se1),exp(mean2-2.24*se2),exp(mean1+2.24*se1),exp(mean2+2.24*se2), lty=2)
rect(exp(mean1-max_quantile*se1),exp(mean2-max_quantile*se2),exp(mean1+max_quantile*se1),exp(mean2+max_quantile*se2),lty=1)
points (exp(mean1), exp(mean2),pch=18)
text(exp(mean1), exp(mean2)+0.1, "(1.56, 0.51)")

axis(side=1)
axis(side=2)
#title(main="(c)")
legend(list(x=exp(mean1+se1),y=1.2),c("Chisq", "Bonf","Max"), cex=1, lty = c(4,2,1))
mtext("CSH ratio", side=1, line=2.4, cex=1)
mtext("OCH ratio", side=2, line=2.5, cex=1)


