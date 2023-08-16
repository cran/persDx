###
#0. empirical auc
###
phi=function(x,y)
  return(I(x<y)+0.5*I(x==y))

np_lpd_auc=function(D,Y){
  Y1=Y[D==1]
  Y0=Y[D==0]
  n1=length(Y1)
  n0=length(Y0)
  #AUC=sum(outer(Y0,Y1,"<")+0.5*outer(Y0,Y1,"=="))/n1/n0
  AUC=sum(outer(Y0, Y1, phi))/n1/n0

  return(AUC)
}
