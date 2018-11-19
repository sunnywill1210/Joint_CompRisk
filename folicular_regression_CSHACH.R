########################################################################
### Example: Follicular
### Condcut Maximum and Chisquare tests for CSH and ACH based on regression
### Produce confidence region figure
########################################################################

rm(list=ls(all=TRUE))
library(MASS)
library(survival)
library(fields)
library(ellipse)

#load data
survival <- read.table("~/Dropbox/Dissertation/first project/realdata/follicular/follicular.txt", sep = ",", header = TRUE)

#clean data
evcens <- as.numeric(survival$resp == "NR" | survival$relsite != "")
crcens <- as.numeric(survival$resp == "CR" & survival$relsite == "" & survival$stat == 1)
survival$cause <- ifelse(evcens == 1, 1, ifelse(crcens == 1, 2, 0))
table(survival$cause)

survival$stage <- as.numeric(survival$clinstg == 2)
chemo <- as.numeric(survival$ch == "Y")
times1 <- sort(unique(survival$dftime[survival$cause == 1]))

survival$time=survival$dftime
survival$event<- ifelse(survival$cause==0, 0, 1)
survival$type<- ifelse(survival$cause==1, 1, 0)
survival$group<-chemo
head(survival)

####################################################
###  Conduct joint test         ###################
####################################################
ntest=5
teststat=matrix(NA,nrow=1,ncol=ntest)
count=matrix(0,nrow=1,ncol=ntest)
p.val=matrix(NA,nrow=1,ncol=ntest)

### Test 1: Cox model for overall survival ###
sdfcox1 <- coxph(Surv(time,event==1) ~ group+age+stage+hgb, data = survival)
teststat[1]=sdfcox1$coef[1]/sqrt(sdfcox1$var[1,1])
p.val[1]=pnorm(abs(teststat[1]),lower.tail=FALSE)

### Test 2: Cox Model for cause-specific hazard ###
sdfcox2 <- coxph(Surv(time, type==1) ~ group+age+stage+hgb, data = survival)
teststat[2]=sdfcox2$coef[1]/sqrt(sdfcox2$var[1,1])
p.val[2]=pnorm(abs(teststat[2]),lower.tail=FALSE)

### Test 3: bonferroni adjustment for test 1 and 2 ###
p.val[3]=min(p.val[1],p.val[2])*2

### Test 4: Chi-square test ###
var_1=sdfcox2$var[1,1]     #variance of S_j
var_12=sdfcox1$var[1,1]		#variance of S
var_2=sdfcox1$var[1,1]	
corr=var_12/(sqrt(var_1)*sqrt(var_2))

S_1=sdfcox2$coef[1] #S_j
S_12=sdfcox1$coef[1]	#S
S=matrix(c(S_12,S_1),nrow=1)
cov=matrix(c(var_12,var_12,var_12,var_1),nrow=2) #covariance matrix
if (det(cov)==0) {browser()}
teststat[4]=S%*%solve(cov)%*%t(S)		     #chisquare test stat
p.val[4]=pchisq(teststat[4],df=2,lower.tail=FALSE)

###Test 5: maximum test (one sided)
nbootstrap=50000
sample=mvrnorm(n = nbootstrap, c(0,0), Sigma=matrix(c(1,corr,corr,1),nrow=2), tol=1e-4)
max=matrix(NA,nrow=nbootstrap,ncol=1)
teststat[5]=max(-teststat[1], -teststat[2]) # negative sign because the origninal values for teststat are negative
for (j in 1:nbootstrap){max[j]=max(-sample[j,1],-sample[j,2])}
index<-seq(1,nbootstrap)
index.a<-index[sort(max)>teststat[5]]
p.val[5]=1-index.a[1]/nbootstrap

names(p.val) = c("OS", "CSH","Bonf","Chisq","max") 
round(p.val,3)
p.val

##############################
## Confidence Region for HR 
##############################

mean1=sdfcox2$coef[1]
se1=sqrt(sdfcox2$var[1,1])

mean2=sdfcox1$coef[1]
se2=sqrt(sdfcox1$var[1,1])

corr=corr
max_quantile=quantile(max,0.95)

plot(c(exp(mean1-3.5*se1),exp(mean1+3*se1)),c(exp(mean2-3.5*se2),exp(mean2+3*se2)), type='n',xlab="",ylab="",axes=F)
lines(exp(ellipse(corr,scale=c(se1,se2), centre=c(mean1,mean2),level=0.95)),type='l',lty=4)
rect(exp(mean1-2.24*se1),exp(mean2-2.24*se2),exp(mean1+2.24*se1),exp(mean2+2.24*se2), lty=2)
rect(exp(mean1-max_quantile*se1),exp(mean2-max_quantile*se2),exp(mean1+max_quantile*se1),exp(mean2+max_quantile*se2),lty=1)
points (exp(mean1), exp(mean2),pch=18)
text(exp(mean1), exp(mean2)+0.05,"(0.74,0.76)")

axis(side=1)
axis(side=2)
#title(main="Figure 2: Confidence Region for Hazard Ratio")
legend(list(x=0.85,y=0.7),c("Chisq", "Bonf","Max"), cex=0.8, lty = c(4,2,1))
mtext("CSH ratio", side=1, line=2.4, cex=1)
mtext("ACH ratio", side=2, line=2.5, cex=1)


