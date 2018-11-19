########################################################################
### Example: Hodgkin disease
### Condcut two sample Maximum and Chisquare tests for CSH and CIF
########################################################################

rm(list=ls(all=TRUE))
library(MASS)
library(survival)
library(cmprsk)


# load program required to calculate correlation between CSH and CIF test statistics
source("Correlation_CSHCIF.R")

# load data
#hd_all <- read.table("~/Dropbox/Dissertation/first project/realdata/Hodgkin/hodgkin.txt", sep = ",", header = TRUE)
hd_all <- read.table(".\\hodgkin.txt", sep = ",", header = TRUE)

###Cleaning dataset and create necessary variables
# dfcens: Censoring variable:0=Censored, 1=Failure
# mcens: Second malignancy: 0= No, 1= Yes

# time to the first event
hd_all$alltime=apply(cbind(hd_all$dftime,hd_all$maltime),1,min)

m1cens=(hd_all$alltime==hd_all$maltime & hd_all$mcens==1)+0 
cr1cens=(hd_all$alltime==hd_all$dftime & hd_all$dfcens==1)+0

hd_all$cens=m1cens+2*cr1cens 

# event of interest, second malignancy 1= yes  0= no
hd_all$m1cens=(hd_all$alltime==hd_all$maltime & hd_all$mcens==1)+0 
#hd$cr1cens=(alltime==hd$dftime & hd$dfcens==1)+0

# group variable
hd_all$x=(hd_all$age>=30)+0

time=hd_all$alltime  #X
type=hd_all$cens 
group=hd_all$x

data<-data.frame(time,type,group)

####################################################
###  Conduct joint tests############################
####################################################
teststat_CSHCIF=matrix(NA,nrow=1,ncol=5)
pvalue_CSHCIF=matrix(NA,nrow=1,ncol=5)

### Test 1: Wald test for group coefficient in cox model for cause specific hazard for event of interest### 
sdfcox1 <- coxph(Surv(time,type==1) ~ group, data = data)
teststat_CSHCIF[1]=(sdfcox1$coef/sqrt(sdfcox1$var))^2  
pvalue_CSHCIF[1]=1-pchisq(teststat_CSHCIF[1],1)

### Test 2: Wald test for group coefficient in cox model for other cause hazard### 
ci<-crr(time,type,group,failcode=1,cencode=0)
teststat_CSHCIF[2]=(ci$coef/sqrt(ci$var))^2  
pvalue_CSHCIF[2]=1-pchisq(teststat_CSHCIF[2],1)

### Test 3: bonferroni adjustment for test 1 and 2 ###
pvalue_CSHCIF[3]=min(pvalue_CSHCIF[1],pvalue_CSHCIF[2])*2


### calculate correlation between two test statistics
Sigma_pq=Sigma_pq(data$time,data$type,data$group)

corr_beta=Sigma_pq$corr_beta
cov_beta=Sigma_pq$cov_beta

### Test 4:  Chisquare Test ########
var_1=sdfcox1$var[1,1]     
var_2=ci$var
var_12=cov_beta
W_1=sdfcox1$coef[1]
W_2=ci$coef[1]
W=matrix(c(W_1,W_2),nrow=1)
cov=matrix(c(var_1,var_12,var_12,var_2),nrow=2) #covariance matrix
if (det(cov)==0) {browser()}
teststat_CSHCIF[4]=W%*%solve(cov)%*%t(W)         #chisquare test stat
pvalue_CSHCIF[4]=pchisq(teststat_CSHCIF[4],df=2,lower.tail=FALSE)

### Test 5: Maximum Test ########
maxtest=max(abs(sdfcox1$coef[1]/sqrt(sdfcox1$var[1,1])), abs(ci$score/sqrt(ci$var[1,1])))
nbootstrap=100000
sample=mvrnorm(n = nbootstrap, c(0,0), Sigma=matrix(c(1,corr_beta,corr_beta,1),nrow=2), tol = 1e-6)
max=matrix(NA,nrow=nbootstrap,ncol=1)
for (j in 1:nbootstrap){max[j]=max(abs(sample[j,]))}
sortmax<-sort(max)
pvalue_CSHCIF[5]=1-(seq(1,length(max))[maxtest<sortmax][1]-1)/nbootstrap

colnames(pvalue_CSHCIF)=c("CSH","CIF","Bonf","Chi","Max")
pvalue_CSHCIF
