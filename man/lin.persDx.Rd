\name{lin.persDx}
\alias{lin.persDx}
\title{Estimating linear personalized diagnostic rules.}
\description{Estimate personalized diagnostics rule that recommends biomarker (or test) A versus B by maximizing area under the receiver operating (ROC) curve (AUC) based on linear combinatino of predictors.}
\usage{
  lin.persDx(D,YA,YB,X,dirA,dirB,eps,PLOT)
}
\arguments{
  \item{D}{Binary outcome with D=1 for disease (or case) and D=0 for non-diseased (or control) (n X 1 vector).}
  \item{YA}{Biomarker A, measured on a continuous or orinal scale (n X 1 vector).}
  \item{YB}{Biomarker B, measured on a continuous or ordinal scale (n X 1 vector).}
  \item{X}{Predictors (n x p matrix).}
  \item{dirA}{Direction of YA to D, where dirA="<" (or dirA=">") indicates higher (or lower) YA is assoicated with Pr(D=1)). Default is dirA="<".}
  \item{dirB}{Direction of YB to D, where dirB="<" (or dirB=">") indicates higher (or lower) YB is assoicated with Pr(D=1)). Default is dirB="<".}
  \item{eps}{Tuning parameter for predictor selections. Default is eps=0.01.}
  \item{PLOT}{TRUE of FALSE for showing plot.}
}
\details{The lin.persDx function estimates the personalized diagnostics rule \eqn{\tau(X)}, where \eqn{\tau(X)}=A recommends YA if \eqn{\theta_0+\theta_1 X_1+\cdots+\theta_p X_p > 0} or \eqn{\tau(X)}=B recommends YB otherwise by maximizing AUC. It guarantess to increase AUC compared to the one-size-fits all strategh. That is, AUC is always greater than maximum of AUC.A and AUC.B, where AUC.A is an empirical AUC value using (D,YA), AUC.B is an empirical AUC value using (D,YB), AUC is an empirical AUC value using (D,YC), where YC=YA if \eqn{\tau(X)}=A and YC=YB if \eqn{\tau(X)}=B. Since the empirical AUC is not a continuous function on (\eqn{\theta_0,\theta_1,\ldots,\theta_p}), a grid-based search algorithm is used. Specifically, for a univariate predictor \eqn{X=X_1}, a grid search algorithm is used to estimate (\eqn{\theta_0,\theta_1}), and for multiple predictors \eqn{X=(X_1,X_2,\ldots,X_p)}, a forward selection algorithm is further implemented that sequentially adds each of the predictors to \eqn{\tau(X)} that increases the AUC most. The stopping criteria is AUC increasement \eqn{\le} eps. It allows eps=0, which searches all \eqn{X}, while some of \eqn{\theta} can be estimated to zeros when corresponding predictors have no contribution to increase the AUC.
}

\value{A list of class lin.persDx:
\item{df}{Data frame with YA, YB, X, and estimated subgroup (A or B).}
\item{AUC.A}{AUC for YA.}
\item{AUC.B}{AUC for YB.}
\item{AUC}{AUC for YA or YB based on \eqn{\tau(X)}.}
\item{theta}{\eqn{\theta_0, \theta_1, \ldots, \theta_p}.}
}

\references{Estimation of linear decision rules for personalized diagnostics using area under the receiver operating characteristics curve (in progress)}
\author{Yunro Chung [aut, cre]}

\examples{
set.seed(1)
n1=n0=100
n=n1+n0
D=c(rep(1,n1),rep(0,n0))
X1=runif(n,0,1)
X2=runif(n,0,1)
X3=runif(n,0,1)
X=data.frame(X1,X2,X3)
muA=2
muB=1
mu0=0
sigma=1

tau=rep("B",n)
tau[X1+X2>=1]="A"

YA=D*(rnorm(n,muA,sigma)*(tau=="A")+rnorm(n,mu0,sigma)*(tau=="B"))+(1-D)*rnorm(n,mu0,sigma)
YB=D*(rnorm(n,muB,sigma)*(tau=="B")+rnorm(n,mu0,sigma)*(tau=="A"))+(1-D)*rnorm(n,mu0,sigma)

fit=lin.persDx(D, YA, YB, X)
fit
}
