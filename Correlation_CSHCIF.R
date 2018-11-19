###############################################################################################
##Function name: Sigma_pq
##Purpose: Calculate the covariance matrix in equation (20) in JASA paper
##Input: dataset with time (X=min (T,C)), Type (0=censor, 1=event of interest, 2=competing risk)
      # and group variable
##Reference:Li & Yang JASA, Fine & Gray 1999
##Note: only works for one covariate in the model
##Author: Qing Yang
##Last updated: April 2018
###############################################################################################


Sigma_pq<-function (time,type,group){

  time=time
  type=type 
  group=group
  nsubject=length(time)
  
###### Testing Cause specific hazard #######
sdfcox1 <- coxph(Surv(time,type==1) ~ group)
#sdfcox1

utime=length(coxph.detail(sdfcox1)$time)
emean=coxph.detail(sdfcox1)$means
riskmat=coxph.detail(sdfcox1,riskmat=TRUE)$riskmat
hazard=coxph.detail(sdfcox1)$hazard
nrisk=coxph.detail(sdfcox1)$nrisk

dN=matrix(0,nrow=nsubject,ncol=utime)  
cox.time =  coxph.detail(sdfcox1)$time 
sur.time =  time 
for (i in 1:nsubject){
  dN[i,cox.time==sur.time[i]] = 1
} 

# only works for one covariate
mean_x=mean(group)
hh=matrix(NA,nrow=nsubject,ncol=utime)
ll=matrix(NA,nrow=nsubject,ncol=utime)
l=matrix(NA,nrow=nsubject,ncol=1)
for (i in 1:nsubject){
  for (j in 1:utime){
    hh[i,j]=(group[i]-emean[j])^2*riskmat[i,j]*(hazard[j]*exp((group[i]-mean_x)*sdfcox1$coef))
    ll[i,j]=(group[i]-emean[j])*(dN[i,j]-riskmat[i,j]*hazard[j]*exp((group[i]-mean_x)*sdfcox1$coef))
  }
  l[i]=sum(ll[i,])
}

omega_pp=1/sdfcox1$var
#print(sum(hh))
#print(sum(l^2))
#print(sum(coxph.detail(sdfcox1)$imat*t(coxph.detail(sdfcox1)$nrisk*coxph.detail(sdfcox1)$hazard)))
#print(sum(coxph.detail(sdfcox1)$imat))
#print(1/sdfcox1$var)

##### Testing Cumulative incidence hazard ######
ci<-crr(time,type,group,failcode=1,cencode=0)
utime_ci=length(ci$uftime) # get the unique failure time points of event of interest

# create Kaplan-Meier Estimator for censoring distribution
censfit=survfit(Surv(time,type == 0) ~ 1) # fit a survival curve for censoring variable, a MAJOR change on april 2018
censfit2=crr(time,type,group,failcode=0,cencode=1)  # get the unique censoring time points
uctime=length(censfit2$uftime)    # number of distinct censoring time

rY_j=matrix(NA,nrow=nsubject,ncol=utime_ci) #\tilde{Y}
rN_j=matrix(NA,nrow=nsubject,ncol=utime_ci) #\tilde{N}
rdN_j=matrix(NA,nrow=nsubject,ncol=utime_ci) # d\tilde{N} in paper
r=matrix(NA,nrow=nsubject,ncol=utime_ci) #indicator I(C_i>=T_i^t_j)
w=matrix(NA,nrow=nsubject,ncol=utime_ci)  #omega in paper, IPCW


for (i in 1:nsubject){
  for (j in 1:utime_ci){
    if (time[i]<ci$uftime[j] & type[i]==0) {r[i,j]=0}
    else {r[i,j]=1}
    w[i,j]=r[i,j]*summary(censfit,times=ci$uftime[j])$surv/summary(censfit,times=min(ci$uftime[j],time[i]))$surv
    
    if (r[i,j]==1) {
      if ((time[i]==ci$uftime[j]) & type[i]==1) {rdN_j[i,j]=1} 
      else {rdN_j[i,j]=0}
      if ((time[i]<=ci$uftime[j]) & type[i]==1) {rN_j[i,j]=1} 
      else {rN_j[i,j]=0}
    }
    else {
      rdN_j[i,j]=0
      rN_j[i,j]=0    
    }
    if (j==1) {rY_j[i,j]=1} #always at risk at time 1
    else {rY_j[i,j]=1-rN_j[i,j-1]}  
  }
}

emean_ci=matrix(NA,nrow=utime_ci,ncol=1)
### S_x(beta,u) in Jason Fine's paper
### first and second moments score function at each different time point
S_2=matrix(NA,nrow=utime_ci,ncol=1)
S_1=matrix(NA,nrow=utime_ci,ncol=1)
S_0=matrix(NA,nrow=utime_ci,ncol=1)

for (j in 1:utime_ci){
  S_1[j]=mean(w[,j]*rY_j[,j]*group[]*exp(ci$coef*group[]))
  S_0[j]=mean(w[,j]*rY_j[,j]*exp(ci$coef*group[]))
  S_2[j]=mean(w[,j]*rY_j[,j]*group[]^2*exp(ci$coef*group[]))
  emean_ci[j]=S_1[j]/S_0[j]
}

### calculate \hat{\ita _i}
dm_1=matrix(NA,nrow=nsubject,ncol=utime_ci)
ita=matrix(NA,nrow=nsubject,ncol=1)
for (i in 1:nsubject){
  for (j in 1:utime_ci){
    dm_1[i,j]=rdN_j[i,j]-rY_j[i,j]*exp(ci$coef*group[i])*ci$bfitj[j]
  }
  ita[i]=sum((group[i]-emean_ci[1:utime_ci])*w[i,1:utime_ci]*dm_1[i,1:utime_ci])
}

## censoring distribution
Y_c=matrix(NA,nrow=nsubject,ncol=uctime)  # I(X_i>t)
N_c=matrix(NA,nrow=nsubject,ncol=uctime)  # I(X_i<=t, delta_i=0)
dN_c=matrix(NA,nrow=nsubject,ncol=uctime) # I(X_i=t, delta_i=0)

for (i in 1:nsubject){
  for (j in 1:uctime){
    if ((time[i]<=censfit2$uftime[j]) & type[i]==0) {N_c[i,j]=1} 
    else {N_c[i,j]=0}  
    if ((time[i]==censfit2$uftime[j]) & type[i]==0) {dN_c[i,j]=1} 
    else {dN_c[i,j]=0}
    
    if (time[i]>=censfit2$uftime[j]) {Y_c[i,j]=1}
    else {Y_c[i,j]=0}
  }
}

### Calculate \hat{\pi_i}
## Calculate \hat{M}^c_i(t)
dlambda_c=matrix(NA,nrow=uctime,ncol=1)
dm_c=matrix(NA,nrow=nsubject, ncol=uctime)
for (i in 1:nsubject){
  for (j in 1:uctime){
    dlambda_c[j]=sum(dN_c[,j])/sum(Y_c[,j])
    dm_c[i,j]=dN_c[i,j]-Y_c[i,j]*dlambda_c[j]
  }
}

### calculate q
q=array(NA,dim=c(nsubject,uctime,utime_ci))
q_u=matrix(NA,nrow=uctime)
indi=array(0,dim=c(nsubject,uctime,utime_ci))
phi=matrix(NA,nrow=nsubject,ncol=uctime)
phi_i=matrix(NA,nrow=nsubject,ncol=1)
## calculate \hat{q}_i
for (i in 1:nsubject){
  for (u in 1:uctime){
    for (s in 1:utime_ci){
      if (ci$uftime[s]>=censfit2$uftime[u] & censfit2$uftime[u]>time[i]) {indi[i,u,s]=1}
      else {indi[i,u,s]=0}
      q[i,u,s]=-(group[i]-emean_ci[s])*w[i,s]*dm_1[i,s]*indi[i,u,s]
    }
  }
}

for (u in 1:uctime){
  q_u[u]=sum(q[,u,])/nsubject
  for (i in 1:nsubject){
    phi[i,u]=q_u[u]/(sum(Y_c[,u])/nsubject)*dm_c[i,u]
  }
}

for (i in 1:nsubject){ phi_i[i]=sum(phi[i,])}

#sigma is the variance-covariance maxtrix for score test stat with respect to gamma_j
#sigma=sum((ita[]+phi_i[])^2)/nsubject
omega_qq_star=sum((ita[]+phi_i[])^2)/nsubject

#omega is the partial derivative matrix of the score function
#type_indi=matrix(NA,nrow=nsubject,ncol=1)
omega_qq=matrix(NA,nrow=nsubject,ncol=utime_ci)
for (i in 1:nsubject){
  for (j in 1:utime_ci){
    omega_qq[i,j]=(S_2[j]/S_0[j]-emean_ci[j]^2)*rdN_j[i,j]/nsubject
  }
}

Var_test1=sum(omega_qq[,])^{-2}*omega_qq_star/nsubject
Var_test2=ci$var

###### covariance #######
ll_joint=matrix(NA,nrow=nsubject,ncol=1)
for (i in 1:nsubject){  
  ll_joint[i]=l[i]*ita[i]+l[i]*phi[i]
}

omega_pq=sum(ll_joint)/nsubject

cov_beta=sdfcox1$var*omega_pq*sum(omega_qq)^{-1}
corr_beta=cov_beta/(sqrt(ci$var)*sqrt(sdfcox1$var))
corr_beta

return(list(corr_beta=corr_beta,cov_beta=cov_beta,Var_test1=Var_test1,Var_test2=Var_test2))
}
