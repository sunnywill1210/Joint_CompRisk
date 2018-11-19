########################################################################
### Example: BMT
### Condcut two sample Maximum and Chisquare tests for CSH and CIF
### Produce confidence region figure
########################################################################

rm(list=ls(all=TRUE))
library(MASS)
library(survival)
library(cmprsk)
library(ellipse)

#load data
source(".\\bmt_datacleaning.R")

#load the required program to calculate correlation between CSH and CIF test statistics
source(".\\Correlation_CSHCIF.R")

#Specify which group variable to test
group=mini

teststat_CSHCIF=matrix(NA,nrow=1,ncol=5)
pvalue_CSHCIF=matrix(NA,nrow=1,ncol=5)

### Test 1: Wald test for group coefficient in cox model for cause specific hazard for event of interest### 
sdfcox1 <- coxph(Surv(time,type==1) ~ group)
teststat_CSHCIF[1]=(sdfcox1$coef/sqrt(sdfcox1$var))^2  
pvalue_CSHCIF[1]=1-pchisq(teststat_CSHCIF[1],1)

### Test 2: Wald test for group coefficient in cox model for other cause hazard### 
ci<-crr(time,type,group,failcode=1,cencode=0)
teststat_CSHCIF[2]=(ci$coef/sqrt(ci$var))^2  
pvalue_CSHCIF[2]=1-pchisq(teststat_CSHCIF[2],1)

### Test 3: bonferroni adjustment for test 1 and 2 ###
pvalue_CSHCIF[3]=min(pvalue_CSHCIF[1],pvalue_CSHCIF[2])*2

### Calculate correlation
Sigma_pq=Sigma_pq(time=time,type=type,group=group)

corr_beta=Sigma_pq$corr_beta
cov_beta=Sigma_pq$cov_beta

### Test 4: Chisquare Test ########
var_1=sdfcox1$var[1,1]     
var_2=ci$var
var_12=cov_beta
W_1=sdfcox1$coef[1]
W_2=ci$coef[1]
W=matrix(c(W_1,W_2),nrow=1)
cov=matrix(c(var_1,var_12,var_12,var_2),nrow=2) #covariance matrix
if (det(cov)==0) {browser()}
chisq=W%*%solve(cov)%*%t(W)        
pvalue_CSHCIF[4]=pchisq(chisq,df=2,lower.tail=FALSE)


### Test 5: Maximum Test ########
maxtest=max(abs(sdfcox1$coef[1]/sqrt(sdfcox1$var[1,1])), abs(ci$coef[1]/sqrt(ci$var[1,1])))
nbootstrap=50000
sample=mvrnorm(n = nbootstrap, c(0,0), Sigma=matrix(c(1,corr_beta,corr_beta,1),nrow=2), tol = 1e-6)
max=matrix(NA,nrow=nbootstrap,ncol=1)
for (j in 1:nbootstrap){max[j]=max(abs(sample[j,]))}
pvalue_CSHCIF[5] = sum(max>maxtest)/nbootstrap

colnames(pvalue_CSHCIF)=c("CSH","CIF","Bonf","Chi","Max")
pvalue_CSHCIF
