#estimate decision rule based on emprical AUC with univariate covariate (X)
np_lpd_cov2=function(D,YA,YB,Xq,Xr,alpha0q,alpha0r){

  beta2.unq=sort(unique(c((Xq-alpha0q)/(Xr-alpha0r)))) #0 and (Inf or -Inf) are always included
  beta0.unq=alpha0q-beta2.unq*alpha0r
  n.unq=length(beta2.unq)

  AUC.pt=rep(NA,n.unq)
  AUC.nt=rep(NA,n.unq)

  for(j in 1:n.unq){ #compute subgroup AUC at each j with two direction
    beta2=beta2.unq[j]
    beta0=beta0.unq[j]

    beta1.pt=1 #when beta2=1
    AUC.pt[j]=np_lpd_cov2_subauc(D,YA,YB,Xq,Xr,alpha0q,alpha0r,beta0,beta1.pt,beta2)

    beta1.nt=-1 #when beta2=-1
    AUC.nt[j]=np_lpd_cov2_subauc(D,YA,YB,Xq,Xr,alpha0q,alpha0r,beta0,beta1.nt,beta2)
  }
  Midx.pt=is.finite(AUC.pt)
  Midx.nt=is.finite(AUC.nt)

  #2. argmax
  if(sum(Midx.pt)==0 & sum(Midx.nt)==0){ #no grid search points that pass the restrictions.
    beta0.hat=NA;
    beta1.hat=NA;
    beta2.hat=NA;
    AUC.max=NA
  }else{
    if(sum(Midx.pt)==0){ #for YA: no grid search points that pass the restrictions.
      pt.max=0
    }else if(sum(Midx.nt)==0){ #for YB: no grid search points that pass the restrictions.
      pt.max=1
    }else{
      pt.max=max(AUC.pt[Midx.pt],na.rm=TRUE)>max(AUC.nt[Midx.nt],na.rm=TRUE)
    }
    if(pt.max){
      beta1.hat=1
      AUC=AUC.pt
      Midx=Midx.pt
      pch.hat=16
    }else{
      beta1.hat=-1
      AUC=AUC.nt
      Midx=Midx.nt
      pch.hat=17
    }
    sel.max=which.max(AUC[Midx])
    AUC.max=(AUC[Midx])[sel.max]
    beta2.hat=(beta2.unq[Midx])[sel.max]
    beta0.hat=(beta0.unq[Midx])[sel.max]

    #beta2.angle=atan(beta2.unq)*180/pi
    #beta2.angle.hat=atan(beta2.hat)*180/pi
  }

  return(list(beta0=beta0.hat,beta1=beta1.hat,beta2=beta2.hat,AUC=AUC.max))
}

