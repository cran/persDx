###
#0. empirical auc
###
np_lpd_survival_auc=function(D,Y,Stime,predict.time,span){
  AUC=survivalROC::survivalROC.C(Stime=Stime,status=D,marker=Y,predict.time=predict.time,span=span)$AUC
  return(AUC)
}

###
#1. subgroup auc with 1 cov
###
np_lpd_cov1_survival_subauc=function(D,YA,YB,Xq,alpha0,alpha1,Stime,predict.time,span){
  #1. tau.hat
  eqL=alpha1*Xq
  eqR=alpha1*alpha0

  tau=NA
  tau[eqL>=eqR]="A"
  tau[eqL<eqR]="B"

  #2. subAUC
  StimeA=Stime[tau=="A"]
  DA=D[tau=="A"]
  YAA=YA[tau=="A"]

  StimeB=Stime[tau=="B"]
  DB=D[tau=="B"]
  YBB=YB[tau=="B"]

  Stime.comb=c(StimeA,StimeB)
  D.comb=c(DA,DB)
  YAB.comb=c(YAA,YBB)

  AUC=np_lpd_survival_auc(D.comb,YAB.comb,Stime.comb,predict.time,span)

  #3. sample size
  #n1.A=sum(DA==1)
  #n0.A=sum(DA==0)
  #n1.B=sum(DB==1)
  #n0.B=sum(DB==0)

  return(AUC=AUC)
}

###
#1. subgroup auc with 2 cov
###
np_lpd_cov2_survival_subauc=function(D,YA,YB,Xq,Xr,alpha0q,alpha0r,beta0,beta1,beta2,Stime,predict.time,span){ #alpha0q is not needed
  #1. tau.hat
  if(is.infinite(beta2)){
    eqL=beta1*Xr
    eqR=beta1*alpha0r
  }else{ #it can handle when beta2==0
    eqL=beta1*Xq
    eqR=beta1*(beta2*Xr+beta0)
  }

  tau=NA
  tau[eqL>=eqR]="A"
  tau[eqL<eqR]="B"

  #2. sub.sample
  StimeA=Stime[tau=="A"]
  DA=D[tau=="A"]
  YAA=YA[tau=="A"]

  StimeB=Stime[tau=="B"]
  DB=D[tau=="B"]
  YBB=YB[tau=="B"]

  Stime.comb=c(StimeA,StimeB)
  D.comb=c(DA,DB)
  YAB.comb=c(YAA,YBB)

  AUC=np_lpd_survival_auc(D.comb,YAB.comb,Stime.comb,predict.time,span)

  #3. sample size
  #n1.A=sum(DA==1)
  #n0.A=sum(DA==0)
  #n1.B=sum(DB==1)
  #n0.B=sum(DB==0)

  return(AUC=AUC)
}
