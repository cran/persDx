np_lpd_cov1=function(D,YA,YB,X,cov.idx,p,
                     STEP,AUC.p,sel.p,theta0,theta1,
                     AUCA,AUCB, n){

  #1. grid search for each Xq
  res1=data.frame(alpha0=rep(NA,p),alpha1=NA,AUC=NA)
  for(q in cov.idx){
    Xq=X[,q]
    res1[q,]=np_lpd_cov1_comp(D,YA,YB,Xq, n)
  }

  if(sum(!is.na(res1$alpha0))==0){ #no improvements
    return(list(theta0=NA,theta1=NA,
                sel.p=NA,AUC.p=pmax(AUCA,AUCB),res1=NA))
  }

  #2. select p that maximizes AUC
  sel.p[STEP]=which.max(res1$AUC) #column 5 is AUC (alpha0, alpha1, AUC.AA, AUC.BB, AUC)
  AUC.p[STEP]=res1$AUC[sel.p[[STEP]]]
  #Xq=X[,sel.p[STEP]]

  res1q=res1[sel.p[STEP],]
  alpha1q=res1q$alpha1     #direction
  alpha0q=res1q$alpha0     #threshold

  #3. reparameterization
  theta0[STEP]=alpha0q
  theta1[STEP,sel.p[STEP]]=alpha1q

  return(list(theta0=theta0,theta1=theta1,
              sel.p=sel.p,AUC.p=AUC.p,res1=res1))
}

#estimate decision rule based on emprical AUC with univariate covariate (X)
np_lpd_cov1_comp=function(D,YA,YB,Xq, n){
  X.unq=sort(unique(Xq))
  n.unq=length(X.unq)

  alpha1.pt=+1; AUC.pt=rep(NA,n.unq) #positive
  alpha1.nt=-1; AUC.nt=rep(NA,n.unq) #negative

  #1. grid search
  for(j in 1:n.unq){
    alpha0=X.unq[j]

    YC=YC1.ft(YA,YB,Xq,alpha0,alpha1.pt, n)
    AUC.pt[j]=np_lpd_auc(D,YC)

    YC=YC1.ft(YA,YB,Xq,alpha0,alpha1.nt, n)
    AUC.nt[j]=np_lpd_auc(D,YC)
  }

  #2. argmax
  Midx.pt=is.finite(AUC.pt)
  Midx.nt=is.finite(AUC.nt)
  if(sum(Midx.pt)==0 & sum(Midx.nt)==0){ #no grid search points that pass the restrictions.
    alpha0.hat=NA;
    alpha1.hat=NA;
    AUC.max=NA
  }else{
    if(sum(Midx.pt)==0){ #for YA: no grid search points that pass the restrictions.
      pt.max=0 #negative
    }else if(sum(Midx.nt)==0){ #for YB: no grid search points that pass the restrictions.
      pt.max=1 #positive
    }else{
      pt.max=max(AUC.pt[Midx.pt],na.rm=TRUE)>max(AUC.nt[Midx.nt],na.rm=TRUE)
            #0 for negative
            #1 for positive
    }
    if(pt.max){ #max when alpha1.pt
      alpha1.hat=alpha1.pt
      AUC=AUC.pt
      Midx=Midx.pt
    }else{ #max when alpha1.nt
      alpha1.hat=alpha1.nt
      AUC=AUC.nt
      Midx=Midx.nt
    }
    AUC.max=max(AUC,na.rm=TRUE)
    sel.max=which(AUC==AUC.max)[1]
    alpha0.hat=X.unq[sel.max]
  }
  return(data.frame(alpha0=alpha0.hat,alpha1=alpha1.hat,AUC=AUC.max))
}
