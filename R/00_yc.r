YC1.ft=function(YA,YB,Xq,alpha0,alpha1, n){ #combined (or personalized) Y based on alpha
  YC=tau=rep(NA,n)

  eqL=alpha1*Xq
  eqR=alpha1*alpha0

  tau[eqL>=eqR]="A"
  tau[eqL<eqR]="B"

  YC[tau=="A"]=YA[tau=="A"]
  YC[tau=="B"]=YB[tau=="B"]

  return(YC)
}

YC2.ft=function(YA,YB,Xq,Xr,alpha0r,beta0,beta1,beta2, n){ #combined (or personalized) Y based on beta
                                                     #alpha0q is not needed
  YC=tau=rep(NA,n)

  if(is.infinite(beta2)){
    eqL=beta1*Xr
    eqR=beta1*alpha0r
  }else{ #it can handle when beta2==0
    eqL=beta1*Xq
    eqR=beta1*(beta2*Xr+beta0)
  }

  tau[eqL>=eqR]="A"
  tau[eqL<eqR]="B"

  YC[tau=="A"]=YA[tau=="A"]
  YC[tau=="B"]=YB[tau=="B"]

  return(YC)
}

YC3.ft=function(df,theta0,theta1,n,p, num.out){

  YC=tau=rep(NA,n)

  #data
  D=df$D
  YA=df$YA
  YB=df$YB
  X=df[,(1:p)+num.out] #num of outomce: 3 for D,YA,YB; 4 for Stime,D,YA,YB

  LP=as.matrix(X,nrow=n,ncol=p)%*%as.matrix(theta1,nrow=p,ncol=1)

  tau[which(LP>=theta0)]="A"
  tau[which(LP<theta0)]="B"

  #5.2. subAUC
  YC[tau=="A"]=YA[tau=="A"]
  YC[tau=="B"]=YB[tau=="B"]

  df$YC=YC
  df$tau=tau


  return(df=df)
}

