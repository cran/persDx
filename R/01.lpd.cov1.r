#estimate decision rule based on emprical AUC with univariate covariate (X)
lpd.cov1=function(D,YA,YB,Xq){
  X.unq=sort(unique(Xq))
  n.unq=length(X.unq)

  AUC.pt=rep(NA,n.unq) #subgroup AUC #a1 is positive
  AUC.nt=rep(NA,n.unq)               #a0 is negative

  #1. grid search
  for(j in 1:n.unq){
    alpha0=X.unq[j]

    alpha1.pt=1
    AUC.pt[j]=lpd.cov1.subauc(D,YA,YB,Xq,alpha0=alpha0,alpha1=alpha1.pt)  #alpha=1

    alpha1.nt=-1
    AUC.nt[j]=lpd.cov1.subauc(D,YA,YB,Xq,alpha0=alpha0,alpha1=alpha1.nt) #alpha=-1
  }

  Midx.pt=is.finite(AUC.pt)
  Midx.nt=is.finite(AUC.nt)

  #2. argmax
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
    }
    if(pt.max){ #max when alpha1.pt
      alpha1.hat=alpha1.pt
      AUC=AUC.pt
      Midx=Midx.pt
      pch.hat=16
    }else{ #max when alpha1.nt
      alpha1.hat=alpha1.nt
      AUC=AUC.nt
      Midx=Midx.nt
      pch.hat=17
    }
    sel.max=which.max(AUC[Midx])
    AUC.max=(AUC[Midx])[sel.max]
    alpha0.hat=(X.unq[Midx])[sel.max]
  }
  return(list(alpha0=alpha0.hat,alpha1=alpha1.hat,AUC=AUC.max))
}
