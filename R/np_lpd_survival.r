np_lpd_survival=function(Stime,D,YA,YB,X,dirA="<",dirB="<",
                         predict.time,span,
                         eps=0.01, plot=TRUE,
                         A=0,B=0,c=2,d=2){

  num.out=4 #number of outcome: Stime,D,YA,YB
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
  tau=YC=rep(NA,n)
  df=data.frame(Stime,D,YA,YB,X,tau,YC)

  p=ncol(df)-(num.out+2) #in case X is a vector, where 4 for Stime,D,YA,YB and 2 for tau, YC
  if(p==1)
    X=as.matrix(X)

  ###
  #3. global AUCs
  ###
  AUCA=survivalROC::survivalROC.C(Stime=Stime,status=D,marker=YA,predict.time=predict.time,span=span)$AUC
  AUCB=survivalROC::survivalROC.C(Stime=Stime,status=D,marker=YB,predict.time=predict.time,span=span)$AUC

  ###
  #4. forward grid rotation
  ###
  #4.1. grid search
  STEP=1

  sel.p=AUC.p=rep(NA,p)  #which cov is selected & corresponding AUC
  cov.idx=1:p     #candidate cov idex

  theta0=rep(0,p)
  theta1=data.frame(matrix(0,nrow=p,ncol=p))
  colnames(theta1)=1:p

  fit1=np_lpd_survival_cov1(Stime=Stime,D=D,YA=YA,YB=YB,X=X,cov.idx=cov.idx,p=p, predict.time=predict.time,span=span,
                             STEP=STEP,AUC.p=AUC.p,sel.p=sel.p,theta0=theta0,theta1=theta1,
                             AUCA=AUCA,AUCB=AUCB, n=n)
  res1=fit1$res1
  AUC.p=fit1$AUC.p
  sel.p=fit1$sel.p
  theta0=fit1$theta0
  theta1=fit1$theta1

  #4.2. grid rotation
  AUC.inc=fit1$AUC.p[1]-max(AUCA,AUCB)
  if(AUC.inc>=eps){
    if(p>=2){
      STEP=2

      #4.2.1. threshold modification for Xr
      SD=apply(X,2,sd)
      MIN=apply(X,2,min)
      MAX=apply(X,2,max)
      LW0=pmax(MIN,res1$alpha0-c*SD) #block
      UP0=pmin(MAX,res1$alpha0+c*SD)
      ALPHA0=list()
      for(j in 1:p)
        ALPHA0[[j]]=c(res1$alpha0[j],runif(A,LW0[j],UP0[j])) #ALPHA0

      Xq=THETA1=list()
      Q=1
      Xq[[Q]]=X[,sel.p[STEP-1]]
      THETA0=unlist(ALPHA0[[sel.p[STEP-1]]])
      THETA1[[Q]]=theta1[,1]

      fit2=np_lpd_survival_cov2(Stime=Stime, D=D,YA=YA,YB=YB,Xq=Xq,X=X,cov.idx=cov.idx,p=p, predict.time=predict.time,span=span,
                                 STEP=STEP,sel.p=sel.p,AUC.p=AUC.p,theta0=theta0,theta1=theta1,
                                 Q=Q,A=A,THETA0=THETA0,THETA1=THETA1,ALPHA0=ALPHA0, n=n)
      if(is.na(fit2$AUC.p[STEP]))
        STEP=1

      AUC.inc=fit2$AUC.p[STEP]-AUC.p[STEP-1]
      if(AUC.inc<eps) #no improvement
        STEP=1 #Will be stopped

      theta0=fit2$theta0
      theta1=fit2$theta1
      sel.p=fit2$sel.p
      AUC.p=fit2$AUC.p
    }
    if(p>=3 & STEP>=2){
      for(STEP in 3:p){
        Xqa=THETA0a=THETA1a=list() #each B & A
        Xq =THETA0= THETA1 =list() #combine B & A to Q

        #4.2.2. slope modification
        STEPb=STEP-1 #previous step
        theta0b=theta0[STEPb]
        theta1b=unlist(theta1[STEPb,])

        theta1a.mat=matrix(0,B+1,p)
        theta1a.mat[1,]=theta1b  #first is without modification
        theta1a.mat[,sel.p[1]]=1 #ignore direction

        if(B>0){
          sel.pb=sel.p[2:STEPb]
          for(k in sel.pb){
            lwk=theta1b[k]-d*abs(theta1b[k])
            upk=theta1b[k]+d*abs(theta1b[k])
            theta1a.mat[2:(B+1),k]=runif(B,lwk,upk)
          }
        }

        #4.2.3. threshold modification for Xq
        mat.X=as.matrix(X)
        for(b in 1:(B+1)){
          THETA1a[[b]]=matrix(theta1a.mat[b,],ncol=1)
          Xqa[[b]]=mat.X %*% THETA1a[[b]]
        }

        SD =apply(Xqa[[1]],2,sd)
        MIN=apply(Xqa[[1]],2,min)
        MAX=apply(Xqa[[1]],2,max)
        LW1=pmax(MIN,theta0b-c*SD) #block
        UP1=pmin(MAX,theta0b+c*SD)
        THETA0a=c(theta0b,runif(A,LW1,UP1)) #ALPHA0

        Q=0
        for(b in 1:(B+1)){
          for(a in 1:(A+1)){
            Q=Q+1
            Xq[[Q]]=Xqa[[b]]
            THETA1[Q]=THETA1a[b]
            THETA0[Q]=THETA0a[a]
          }
        }
        THETA0=unlist(THETA0)

        #4.2.4. grid rotation
        fit2=np_lpd_survival_cov2(Stime=Stime, D=D,YA=YA,YB=YB,Xq=Xq,X=X,cov.idx=cov.idx,p=p, predict.time=predict.time,span=span,
                                  STEP=STEP,sel.p=sel.p,AUC.p=AUC.p,theta0=theta0,theta1=theta1,
                                  Q=Q,A=A,THETA0=THETA0,THETA1=THETA1,ALPHA0=ALPHA0, n=n)

        if(is.na(fit2$AUC.p[STEP])){
          STEP=STEP-1
          break
        }
        AUC.inc=fit2$AUC.p[STEP]-AUC.p[STEP-1]
        if(AUC.inc<eps){
          STEP=STEP-1
          break
        }

        theta0=fit2$theta0
        theta1=fit2$theta1
        sel.p=fit2$sel.p
        AUC.p=fit2$AUC.p
      }
    }
  }

  #4.2.5. estimated parameters
  theta1=as.numeric(theta1[STEP,])
  theta0=as.numeric(theta0[STEP])

  if(theta1[sel.p[1]]==-1){ #reparameterization
    theta1=theta1*(-1)
    theta1[sel.p[1]]=-1
    theta0=theta0*(-1)
  }
  names(theta1)=paste0("theta",1:p)

  ###
  #5. subgroup AUCs
  ###
  df=YC3.ft(df=df,theta0=theta0,theta1=theta1,n=n,p=p,num.out=num.out)
  roc.C=survivalROC::survivalROC.C(Stime=df$Stime,status=df$D,marker=df$YC,predict.time=predict.time,span=span)
    AUC.C=roc.C$AUC
    cut.C=roc.C$cut.values
    tp.C=roc.C$TP;
    fp.C=roc.C$FP
    tpfp.C=data.frame(cutoff=cut.C,tp=tp.C,fp=fp.C)

  ###
  #6. plot
  ###
  if(plot==TRUE){
    roc.A=survivalROC::survivalROC.C(Stime=df$Stime,status=df$D,marker=df$YA,predict.time=predict.time,span=span)
    roc.B=survivalROC::survivalROC.C(Stime=df$Stime,status=df$D,marker=df$YB,predict.time=predict.time,span=span)

    tp.A=roc.A$TP;   fp.A=roc.A$FP
    tp.B=roc.B$TP;   fp.B=roc.B$FP

    plot  (tp.C~fp.C,type='l',ylab="Sensitivity",xlab="1-Specificity",col=1,lwd=1,main="ROC")
    points(tp.A~fp.A,type='l',col=2,lwd=1)
    points(tp.B~fp.B,type='l',col=4,lwd=1)
    #abline(a=0,b=1,col="darkgray",lwd=1)
    legend("bottomright",c("AUC","AUCA","AUCB"),col=c(1,2,4),lwd=c(2,1,1),bty='n')
  }

  if(dirA==">") df$YA=-df$YA
  if(dirB==">") df$YB=-df$YB

  return(list(df=df,AUCA=AUCA,AUCB=AUCB,
              AUC=AUC.C,tpfp=tpfp.C, theta=theta1, theta0=theta0))
}
