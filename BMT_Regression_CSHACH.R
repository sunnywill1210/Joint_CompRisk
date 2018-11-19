########################################################################
### Example: BMT
### Condcut Maximum and Chisquare tests for CSH and ACH based on regression
### Produce confidence region figure
########################################################################

rm(list=ls(all=TRUE))
library(MASS)
library(survival)
library(cmprsk)
library(ellipse)

#load data
source(".\\bmt_datacleaning.R")

ntest=5
teststat=matrix(NA,nrow=1,ncol=ntest)
count=matrix(0,nrow=1,ncol=ntest)
p.val=matrix(NA,nrow=1,ncol=ntest)


### Test 1: Cox model for cause specific hazard for event of interest ###

sdfcox1 <- coxph(Surv(time,type==1) ~ mini+adv, data = data)
teststat[1]=sdfcox1$coef[1]/sqrt(sdfcox1$var[1,1])
p.val[1]=pnorm(abs(teststat[1]),lower.tail=FALSE)*2
p.val[1]


### Test 2: Cox Model for cause-specific hazard ###
sdfcox_all <- coxph(Surv(time,type != 0) ~ mini+adv, data = data)
teststat[2]=sdfcox_all $coef[1]/sqrt(sdfcox_all$var[1,1])
p.val[2]=pnorm(abs(teststat[2]),lower.tail=FALSE)*2
p.val[2]

### Test 3: bonferroni adjustment for test 1 and 2 ###
p.val[3]=min(p.val[1],p.val[2])*2
p.val[3]

### Test 4: Joint logrank test statistics ###
var_1=sdfcox1$var[1,1] 		#variance of beta_j
var_2=sdfcox_all$var[1,1]		#variance of beta_\dot
var_12=var_2
### Chisquare test ###
corr<-var_12/sqrt(var_1*var_2)

W_1=sdfcox1$coef/sqrt(var_1)
W_2=sdfcox_all$coef/sqrt(var_2)

W=matrix(c(W_1[1],W_2[1]),nrow=1)
cov=matrix(c(1,corr,corr,1),nrow=2) #covariance matrix
if (det(cov)==0) {browser()}
teststat[4]=W%*%solve(cov)%*%t(W)   #chisquare test stat
p.val[4]=1-pchisq(teststat[4],2)

### maximum test ###
teststat[5]=max(abs(sdfcox1$coef[1]/sqrt(sdfcox1$var[1,1])), abs(sdfcox_all$coef[1]/sqrt(sdfcox_all$var[1,1])))
nbootstrap=1000000
sample=mvrnorm(n = nbootstrap, c(0,0), Sigma=matrix(c(1,corr,corr,1),nrow=2), tol=1e-4)
max=matrix(NA,nrow=nbootstrap,ncol=1)
for (j in 1:nbootstrap){max[j]=max(abs(sample[j,]))}
p.val[5] = sum(max>teststat[5])/nbootstrap

colnames(p.val)=c("CSH","ACH","Bonf","Chi","Max")
p.val


#####################################
### producing confidence region graph
#####################################
mean1=sdfcox1$coef[2]
se1=sqrt(sdfcox1$var[2,2])

mean2=sdfcox_all$coef[2]
se2=sqrt(sdfcox_all$var[2,2])
corr=corr

max_quantile=quantile(max,0.95)


plot(c(exp(mean1-3*se1),exp(mean1+3*se1)),c(exp(mean2-3*se2),exp(mean2+3*se2)), type='n',xlab="",ylab="",axes=F)
lines(exp(ellipse(corr,scale=c(se1,se2), centre=c(mean1,mean2),level=0.95)),type='l',lty=4)
rect(exp(mean1-2.24*se1),exp(mean2-2.24*se2),exp(mean1+2.24*se1),exp(mean2+2.24*se2), lty=2)
rect(exp(mean1-max_quantile*se1),exp(mean2-max_quantile*se2),exp(mean1+max_quantile*se1),exp(mean2+max_quantile*se2),lty=1)
points (exp(mean1), exp(mean2),pch=18)
text(exp(mean1), exp(mean2)+0.1,"(1.46,0.9)")

axis(side=1)
axis(side=2)
legend(list(x=exp(mean1+1.2*se1),y=2.1),c("Chisq", "Bonf","Max"), cex=0.8, lty = c(4,2,1))
mtext("CSH ratio", side=1, line=2.4, cex=1)
mtext("ACH ratio", side=2, line=2.5, cex=1)





