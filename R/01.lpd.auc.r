###
#0. empirical auc
###
lpd.auc=function(D,Y){
  Y1=Y[D==1]
  Y0=Y[D==0]
  n1=length(Y1)
  n0=length(Y0)
  AUC=sum(outer(Y0,Y1,"<")+0.5*outer(Y0,Y1,"=="))/n1/n0

  return(AUC)
}

###
#1. subgroup auc with 1 cov
###
lpd.cov1.subauc=function(D,YA,YB,Xq,alpha0,alpha1){
  #1. tau.hat
  eqL=alpha1*Xq
  eqR=alpha1*alpha0

  tau=NA
  tau[eqL>=eqR]="A"
  tau[eqL<eqR]="B"

  #2. subAUC
  DA=D[tau=="A"]
  YAA=YA[tau=="A"]

  DB=D[tau=="B"]
  YBB=YB[tau=="B"]

  D.comb=c(DA,DB)
  YAB.comb=c(YAA,YBB)

  AUC=lpd.auc(D.comb,YAB.comb)

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
lpd.cov2.subauc=function(D,YA,YB,Xq,Xr,alpha0q,alpha0r,beta0,beta1,beta2){ #alpha0q is not needed
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
  DA=D[tau=="A"]
  YAA=YA[tau=="A"]

  DB=D[tau=="B"]
  YBB=YB[tau=="B"]

  D.comb=c(DA,DB)
  YAB.comb=c(YAA,YBB)

  AUC=lpd.auc(D.comb,YAB.comb)

  #3. sample size
  #n1.A=sum(DA==1)
  #n0.A=sum(DA==0)
  #n1.B=sum(DB==1)
  #n0.B=sum(DB==0)

  return(AUC=AUC)
}
