lin.persDx=function(D,YA,YB,X,dirA="<",dirB="<",eps=0.01, PLOT=TRUE){
  ###
  #1. direction
  ###
  if(dirA==">")
    YA=-YA
  if(dirB==">")
    YB=-YB

  ###
  #2. iniitial setting
  ###
  #2.1. sample size
  n1=sum(D==1)
  n0=sum(D==0)
  n=n1+n0

  #2.2. subggroup & x
  sub=rep(NA,n)
  df=data.frame(D,YA,YB,X,sub)

  p=ncol(df)-4 #in case X is a vector.
  if(p==1)
    X=as.matrix(X)

  ###
  #3. global AUCs
  ###
  AUC.A=lpd.auc(D,YA)
  AUC.B=lpd.auc(D,YB)

  ###
  #4. forward grid rotation
  ###
  #4.1. grid search
  STEP=1
  sel.p=AUC.p=c() #which cov is selected & corresponding AUC
  cov.idx=1:p     #candidate cov idex
  theta=rep(0,p); theta0=NA
  res1=list()
  for(q in cov.idx){
    Xq=X[,q]
    res1[[q]]=lpd.cov1(D,YA,YB,Xq)
  }
  res1.mat=do.call("rbind",res1)
  if(sum(!is.na(res1.mat[,1]))==0){ #no improvements
    theta=c(rep(0,p))
    if(AUC.A>=AUC.B){
      AUC=AUC.A
      df$sub="A"
      theta=c(rep(0,p))
      theta0=-Inf
    }else{
      AUC=AUC.B
      df$sub="B"
      theta=c(rep(0,p))
      theta0=Inf
    }

    theta.hat=c(-theta0,theta)
    theta.hat=c(theta.hat)
    names(theta.hat)=paste0("theta",0:p)

    return(list(df=df,
                AUC.A=AUC.A,AUC.B=AUC.B,     #global AUC
                AUC=AUC, #subgroup AUC
                theta=theta.hat #decision rule. -theta0 is t1X1+...tpXp>t0 => -t0 + vs t1X1+...tpXp>0.
    ))
  }

  sel.p[STEP]=which.max(unlist(res1.mat[,3])) #column 5 is AUC (alpha0, alpha1, AUC.AA, AUC.BB, AUC)
  AUC.p[STEP]=unlist(res1.mat[sel.p[STEP],3])
  Xq=X[,sel.p[STEP]]
  res1q=res1[[sel.p[STEP]]]
  alpha1q=res1q$alpha1 #direction
  alpha0q=res1q$alpha0 #threshold

  #4.2. forwarded grid rotation
  p.useful=sum(!is.na(unlist(res1.mat[,3]))) #NA appears when max(AUC) has no improvement
  if(p.useful>1){
    beta0=beta1=beta2=c()
    for(STEP in 2:p.useful){
      res2=list()
      for(r in cov.idx){
        res2[[r]]=NA
        Xr=X[,r]
        alpha0r=res1[[r]]$alpha0
        if(!(r%in%sel.p))     #rth variable is not selected;
          if(!is.na(alpha0r)) #rth variable is useful (NA appears when max(AUC) has no improvement)
            res2[[r]]=lpd.cov2(D,YA,YB,Xq,Xr,alpha0q,alpha0r)
      }
      res2.mat=do.call("rbind",res2)

      if(sum(!is.na(res2.mat[,1]))==0){ #no more improvement
        STEP=STEP-1
        break
      }else{
        #4.5.2. sel.p: which cov is choosed; AUC.p: correspoinding AUC
        sel.p[STEP]=which.max(unlist(res2.mat[,4])) #column 4 for AUC
        AUC.p[STEP]=unlist(res2.mat[sel.p[STEP],4])
        AUC.inc=AUC.p[STEP]-AUC.p[STEP-1] #AUC increment

        if(AUC.inc<eps){ #stop
          sel.p=sel.p[-STEP]
          AUC.p=AUC.p[-STEP]
          STEP=STEP-1
          break
        }else{ #selected cov & beta & Xq (i.e. X.hat) with threshold
          Xr=X[,sel.p[STEP]]

          res2q=res2[[sel.p[STEP]]]
          beta0[STEP]=beta0qr=res2q$beta0
          beta1[STEP]=beta1qr=res2q$beta1
          beta2[STEP]=beta2qr=res2q$beta2

          Xq=Xq-beta2qr*Xr#lin comb
          alpha0q=beta0qr #threshold
        }
      }
    }
  }

  #reparameterization
  if(STEP==1){
    theta[sel.p[1]]=alpha1q #re-parameterization
    theta0=alpha1q*alpha0q  #multiplication of alpha1q for adjusting the direction
  }else if(STEP>=2){
    theta[sel.p[STEP]]=beta1[STEP] #last direction parameter is only used
    theta0=beta1[STEP]*beta0[STEP] #last direction parameter & intercept are only used
    for(j in 2:STEP)               #for theta0 & theta, multiplication of beta1[STEP] for adjusting the direction
      theta[sel.p[j]]=-beta1[STEP]*beta2[j]
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

  DA=df$D[subA]
  YAA=df$YA[subA]

  DB=df$D[subB]
  YBB=df$YB[subB]

  #5.2 subAUC
  D.comb=c(DA,DB)
  YAB.comb=c(YAA,YBB)

  AUC=lpd.auc(D.comb,YAB.comb)

  ###
  #6. plot
  ###
  if(PLOT==TRUE){
    roc.A=pROC::roc(D~YA,direction="<",levels=c(0,1))
    roc.B=pROC::roc(D~YB,direction="<",levels=c(0,1))
    roc.AB=pROC::roc(D.comb~YAB.comb,direction="<",levels=c(0,1))

    tp.A=roc.A$sensitivities;   fp.A=1-roc.A$specificities
    tp.B=roc.B$sensitivities;   fp.B=1-roc.B$specificities
    tp.AB=roc.AB$sensitivities; fp.AB=1-roc.AB$specificities

    plot(tp.AB~fp.AB,type='l',ylab="Sensitivity",xlab="1-Specificity",col=1,lwd=2,main="ROC")
    points(tp.A~fp.A,type='l',col=2,lwd=1)
    points(tp.B~fp.B,type='l',col=4,lwd=1)
    abline(a=0,b=1,col="darkgray",lwd=1)
    legend("bottomright",c("AUC","AUC.A","AUC.B"),col=c(1,2,4),lwd=c(2,1,1),bty='n')
  }

  if(dirA==">")
    df$YA=-YA
  if(dirB==">")
    df$YB=-YB

  theta.hat=c(theta.hat)
  names(theta.hat)=paste0("theta",0:p)
  return(list(df=df,
              AUC.A=AUC.A,AUC.B=AUC.B,     #global AUC
              AUC=AUC,
              theta=theta.hat #decision rule. -theta0 is t1X1+...tpXp>t0 => -t0 + vs t1X1+...tpXp>0.
              ))
}
