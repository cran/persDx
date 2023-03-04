np_lpd_survival=function(Stime,D,YA,YB,X,dirA="<",dirB="<",predict.time,span,eps=0.01, PLOT=TRUE){
  ###
  #1. direction
  ###
  if(dirA==">") YA=-YA
  if(dirB==">") YB=-YB

  ###
  #2. initial setting
  ###
  #2.1. sample size
  n1=sum(D==1)
  n0=sum(D==0)
  n=n1+n0

  #2.2. subggroup & x
  sub=rep(NA,n)
  df=data.frame(Stime,D,YA,YB,X,sub)

  p=ncol(df)-5 #in case X is a vector.
  if(p==1)
    X=as.matrix(X)

  ###
  #3. global AUCs
  ###
  AUCA=np_lpd_survival_auc(D,YA,Stime,predict.time,span)
  AUCB=np_lpd_survival_auc(D,YB,Stime,predict.time,span)

  ###
  #4. forward grid rotation
  ###
  #4.1. grid search
  STEP=1
  sel.p=AUC.p=NA  #which cov is selected & corresponding AUC
  cov.idx=1:p     #candidate cov idex
  theta=rep(0,p); theta0=NA
  res1=list()
  for(q in cov.idx){
    Xq=X[,q]
    res1[[q]]=np_lpd_survival_cov1(D,YA,YB,Xq,Stime,predict.time,span)
  }
  res1.mat=do.call("rbind",res1)
  if(sum(!is.na(res1.mat[,1]))==0){ #no improvements
    if(AUCA>=AUCB){; AUC=AUCA; df$sub="A"; theta=c(rep(0,p)); theta0=-Inf
    }else{;          AUC=AUCB; df$sub="B"; theta=c(rep(0,p)); theta0=Inf
    }
    theta.hat=c(-theta0,theta)
    names(theta.hat)=paste0("theta",0:p)
    return(list(df=df,AUCA=AUCA,AUCB=AUCB,AUC=AUC,theta=theta.hat))
  }

  sel.p[STEP]=which.max(unlist(res1.mat[,3])) #column 5 is AUC (alpha0, alpha1, AUC.AA, AUC.BB, AUC)
  AUC.p[STEP]=unlist(res1.mat[sel.p[STEP],3])
  Xq=X[,sel.p[STEP]]
  res1q=res1[[sel.p[STEP]]]
  alpha1q=res1q$alpha1     #direction
  alpha0q=res1q$alpha0     #threshold

  theta[sel.p[STEP]]=alpha1q
  theta0=alpha1q*alpha0q

  #4.2. forwarded grid rotation
  p.useful=sum(!is.na(unlist(res1.mat[,3]))) #NA appears when max(AUC) has no improvement
  if(p.useful<=1)
    return(list(df=df,AUCA=AUCA,AUCB=AUCB,AUC=AUC,theta=theta.hat))

  beta0=beta1=beta2=NA
  for(STEP in 2:p.useful){
    res2=list()
    for(r in cov.idx){
      res2[[r]]=NA
      Xr=X[,r]
      alpha0r=res1[[r]]$alpha0
      if(!(r%in%sel.p))     #rth variable is not selected;
        if(!is.na(alpha0r)) #rth variable is useful (NA appears when max(AUC) has no improvement)
          res2[[r]]=np_lpd_survival_cov2(D,YA,YB,Xq,Xr,alpha0q,alpha0r,Stime,predict.time,span)
    }
    res2.mat=do.call("rbind",res2)

    if(sum(!is.na(res2.mat[,1]))==0){ #no more improvement
      STEP=STEP-1
      break
    }

    #4.5.2. sel.p: which cov is choosed; AUC.p: correspoinding AUC
    sel.p[STEP]=which.max(unlist(res2.mat[,4])) #column 4 for AUC
    AUC.p[STEP]=unlist(res2.mat[sel.p[STEP],4])
    AUC.inc=AUC.p[STEP]-AUC.p[STEP-1] #AUC increment
    if(AUC.inc<eps){ #stop
      #sel.p=sel.p[-STEP]
      #AUC.p=AUC.p[-STEP]
      STEP=STEP-1
      break
    }

    Xr=X[,sel.p[STEP]]

    res2q=res2[[sel.p[STEP]]]
    beta0[STEP]=beta0qr=res2q$beta0
    beta1[STEP]=beta1qr=res2q$beta1
    beta2[STEP]=beta2qr=res2q$beta2

    Xq=Xq-beta2qr*Xr#lin comb
    alpha0q=beta0qr #threshold
  }

  #reparameterization
  if(STEP>=2){
    theta=rep(0,p); theta0=NA
    theta[sel.p[1]]=beta1[STEP]    #last direction parameter is only used
    theta0=beta1[STEP]*beta0[STEP] #last direction parameter & intercept are only used
    if(STEP>2){
      for(j in STEP:2){         #for theta0 & theta, multiplication of beta1[STEP] for adjusting the direction
        theta[sel.p[j]]=-beta1[STEP]*beta2[j]
      }
    }
  }

  ###
  #5. subgroup AUCs
  ###
  theta.hat=matrix(c(-theta0,theta),ncol=1) #equivalent to X1%*%theta>theta0
  X1=as.matrix(cbind(1,X))
  LP=X1%*%theta.hat #linear predicttor for personalized diagnostics rule

  subA=which(LP>=0)
  subB=which(LP<0)

  df$sub[subA]="A"
  df$sub[subB]="B"

  StimeA=df$Stime[subA]
  DA=df$D[subA]
  YAA=df$YA[subA]

  StimeB=df$D[subB]
  DB=df$D[subB]
  YBB=df$YB[subB]

  #5.2 subAUC
  Stime.comb=c(StimeA,StimeB)
  D.comb=c(DA,DB)
  YAB.comb=c(YAA,YBB)

  AUC=np_lpd_survival_auc(D.comb,YAB.comb,Stime.comb,predict.time,span)

  ###
  #7. plot
  ###
  if(PLOT==TRUE){
    roc.A=survivalROC::survivalROC.C(Stime=Stime,status=D,marker=YA,predict.time=predict.time,span=span)
    roc.B=survivalROC::survivalROC.C(Stime=Stime,status=D,marker=YB,predict.time=predict.time,span=span)

    roc.AB=survivalROC::survivalROC.C(Stime=Stime.comb,status=D.comb,marker=YAB.comb,predict.time=predict.time,span=span)

    tp.A=roc.A$TP;   fp.A=roc.A$FP
    tp.B=roc.B$TP;   fp.B=roc.B$FP
    tp.AB=roc.AB$TP; fp.AB=roc.AB$FP

    plot(tp.AB~fp.AB,type='l',ylab="Sensitivity",xlab="1-Specificity",col=1,lwd=2,main="ROC")
    points(tp.A~fp.A,type='l',col=2,lwd=1)
    points(tp.B~fp.B,type='l',col=4,lwd=1)
    abline(a=0,b=1,col="darkgray",lwd=1)
    legend("bottomright",c("AUC","AUC_A","AUC_B"),col=c(1,2,4),lwd=c(2,1,1),bty='n')
  }

  df$YC=YAB.comb
  if(dirA==">") df$YA=-YA
  if(dirB==">") df$YB=-YB

  theta.hat=c(theta.hat)
  names(theta.hat)=paste0("theta",0:p)
  colnames(df)=c("D","YA","YB","X","tau","YC")
  return(list(df=df,
              AUCA=AUCA,AUCB=AUCB,     #global AUC
              AUC=AUC,
              theta=theta.hat #decision rule. -theta0 is t1X1+...tpXp>t0 => -t0 + vs t1X1+...tpXp>0.
              ))
}
