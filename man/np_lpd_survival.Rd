\name{np_lpd_survival}
\alias{np_lpd_survival}
\title{Nonparametric estimation of linear personalized diagnostic rules with right-censored survival outcome.}
\description{Nonparametric estimation of personalized diagnostics rule to find subgroup-specific biomarkers according to linear combination of predictors.}
\usage{
  np_lpd_survival(Stime,D,YA,YB,X,dirA,dirB,predict.time,span,eps,plot,A,B,c,d)
}
\arguments{
  \item{Stime}{Event time or censoring time for subjects (n x 1 vector).}
  \item{D}{Indicator of status, where D=1 if death or event, and
  D=0 otherwise (n X 1 vector).}
  \item{YA}{Biomarker A, measured on a continuous or orinal scale (n X 1 vector).}
  \item{YB}{Biomarker B, measured on a continuous or ordinal scale (n X 1 vector).}
  \item{X}{Predictors (n x p matrix).}
  \item{predict.time}{Time point to evaluate YA and YB.}
  \item{span}{Span for the nearest neighbor estimation (NNE).}
  \item{dirA}{Direction of YA to D, where dirA="<" (or dirA=">") indicates higher (or lower) YA is assoicated with Pr(D=1)). Default is dirA="<".}
  \item{dirB}{Direction of YB to D, where dirB="<" (or dirB=">") indicates higher (or lower) YB is assoicated with Pr(D=1)). Default is dirB="<".}
  \item{eps}{Tuning parameter for predictor selections. Default is eps=0.01.}
  \item{plot}{plot=TRUE (or FALSE) shows (or does not show) the receiver operating charactriestics (ROC) curve.}
  \item{A}{Grid search parameter (Discrete). Default is A=0}
  \item{B}{Grid search parameter (Discrete). Default is B=0}
  \item{c}{Grid search parameter. Default is c=2}
  \item{d}{Grid search parameter. Default is d=2}
}
\details{The np.lin.survival.persDx function estimates the personalized diagnostics rule \eqn{\tau(X)}, where \eqn{\tau(X)}=A recommends \eqn{YA} if \eqn{\theta_1 X_1+...+\theta_p X_p > \theta_0} or \eqn{\tau(X)}=B recommends \eqn{YB} otherwise by maximizing (empirical) survival area under the ROC curve (AUC) at the predict.time using the Nearest
Neighbor Estimation. Here, the survival AUC is computed based on \eqn{YC} with the direction of "<", i.e. higher \eqn{YC} is associated with Pr(D=1), where \eqn{YC=YA} if \eqn{\tau(X)}=A and dirA="<", or \eqn{YC=YB} if \eqn{\tau(X)}=B and dirB="<". If dirA=">" (or dirB=">"), negative YA (or YB) is used.

A forward grid rotation algorithm (FGR) is used to estimate \eqn{\theta_0,\theta_1,...,\theta_p} by sequentially adding each of the predictors to \eqn{\tau(X)} that increases the AUC the most. The stopping criteria is AUC increasement is less than or equal to eps. The eps controls the model complexity. The cross-validation techniques can be used to find the optimal eps.

The FGR results in a suboptimal solution. The accuracy is improved by setting higher A, B, c, d, but it increases computational costs. We thus recond this function when p is small or around 10.
}

\value{A list of class np.lin.persDx:
\item{df}{Data frame with Stime, D, YA, YB, X, tau, YC, where tau=A or B for recommending YA or YB, respectively.}
\item{AUCA}{Survival AUC for YA at predict.time.}
\item{AUCB}{Survival AUC for YB at predict.time.}
\item{AUC}{Survival AUC for YC.}
\item{tpfp}{Data frame with cutoff, tp, fp, where tp and fp are true and false positive positives at the cutoff values of YC.}
\item{theta}{Estimated regression parameters.}
\item{theta0}{Estimated threshold parameter.}
\item{PLOT}{TRUE or FALSE to show survival ROC curves.}
}

\references{Yaliang Zhang and Yunro Chung, Nonparametric estimation of linear personalized diagnostics rules via efficient grid algorithm (submitted)}
\author{Yunro Chung [aut, cre]}

\examples{
#simulate data
set.seed(1)
n=100
X=abs(rnorm(n,1,1))
C=abs(rnorm(n,1,1))
Stime=pmin(X,C)
D=as.numeric(X<=C)

X1=runif(n,0,1)
X2=runif(n,0,1)
X3=runif(n,0,1)
X=data.frame(X1,X2,X3)

tau=rep("B",n)
tau[X1+X2>=1]="A"

D2=rep(0,n) #event by time 2
D2[which(Stime<=3 & D==1)]=1

YA=D2*(rnorm(n,2,1)*(tau=="A")+rnorm(n,0,1)*(tau=="B"))+
   (1-D2)*rnorm(n,0,1)
YB=D2*(rnorm(n,1,1)*(tau=="B")+rnorm(n,0,1)*(tau=="A"))+
   (1-D2)*rnorm(n,0,1)

#run
span=0.1
fit=np_lpd_survival(Stime, D, YA, YB, X, predict.time=1, span=span)
fit
}
