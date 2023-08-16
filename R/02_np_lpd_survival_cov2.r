np_lpd_survival_cov2=function(Stime,D,YA,YB,Xq,X,cov.idx,p, predict.time,span,
                     STEP,sel.p,AUC.p,theta0,theta1,
                     Q,A,THETA0,THETA1,ALPHA0, n){

  #1. grid rotation for each Xq and Xr
  res2=list()
  k=0
  for(r in setdiff(cov.idx,sel.p)){
    if(!is.na(ALPHA0[[r]][1])){ #rth variable is useful (NA appears when max(AUC) has no improvement)
      for(q in 1:Q){
        for(a in 1:(A+1)){
          k=k+1
          res2[[k]]=data.frame(r=r,q=q,a=a,
                               np_lpd_survival_cov2_comp(Stime,D,YA,YB,Xq[[q]],X[,r],THETA0[q],ALPHA0[[r]][a],predict.time,span, n))
        }
      }
    }
  }

  res2.mat=do.call(rbind,res2)
  res2.max=res2.mat[which.max(res2.mat$AUC),]

  #2. sel.p: which cov is choosen; AUC.p: correspoinding AUC
  r.hat=res2.max$r #sel.p
  q.hat=res2.max$q
  a.hat=res2.max$a

  sel.p[STEP]=r.hat
  AUC.p[STEP]=res2.max$AUC

  Xq=Xq[[q.hat]]
  Xr=X[,r.hat]
  THETA1=(THETA1[q.hat][[1]])
  THETA0=THETA0[q.hat]

  beta0=res2.max$beta0
  beta1=res2.max$beta1
  beta2=res2.max$beta2

  #3. reparameterization
  theta0[STEP]=beta0
  theta1[STEP,]=c(THETA1)
  theta1[STEP,sel.p[1]]=beta1
  theta1[STEP,r.hat]=-beta2

  return(list(theta0=theta0,theta1=theta1,
              sel.p=sel.p,AUC.p=AUC.p))
}

#estimate decision rule based on emprical AUC with univariate covariate (X)
np_lpd_survival_cov2_comp=function(Stime,D,YA,YB,Xq,Xr,alpha0q,alpha0r,predict.time,span,n){
  YC=tau=rep(NA,n)

  beta2.unq=sort(unique(c((Xq-alpha0q)/(Xr-alpha0r)))) #0 and (Inf or -Inf) are always included
  beta0.unq=alpha0q-beta2.unq*alpha0r
  n.unq=length(beta2.unq)

  beta1.pt=+1; AUC.pt=rep(NA,n.unq) #positive
  beta1.nt=-1; AUC.nt=rep(NA,n.unq) #negative

  for(j in 1:n.unq){ #compute subgroup AUC at each j with two direction
    beta2=beta2.unq[j]
    beta0=beta0.unq[j]

    YC=YC2.ft(YA,YB,Xq,Xr,alpha0r,beta0,beta1.pt,beta2, n)
    AUC.pt[j]=survivalROC::survivalROC.C(Stime=Stime,status=D,marker=YC,predict.time=predict.time,span=span)$AUC

    YC=YC2.ft(YA,YB,Xq,Xr,alpha0r,beta0,beta1.nt,beta2, n)
    AUC.nt[j]=survivalROC::survivalROC.C(Stime=Stime,status=D,marker=YC,predict.time=predict.time,span=span)$AUC
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
      beta1.hat=beta1.pt
      AUC=AUC.pt
      Midx=Midx.pt
    }else{
      beta1.hat=beta1.nt
      AUC=AUC.nt
      Midx=Midx.nt
    }
    AUC.max=max(AUC,na.rm=TRUE)
    sel.max=which(AUC==AUC.max)[1]
    beta2.hat=(beta2.unq)[sel.max]
    beta0.hat=(beta0.unq)[sel.max]

    #beta2.angle=atan(beta2.unq)*180/pi
    #beta2.angle.hat=atan(beta2.hat)*180/pi
  }

  return(list(beta0=beta0.hat,beta1=beta1.hat,beta2=beta2.hat,AUC=AUC.max))
}

