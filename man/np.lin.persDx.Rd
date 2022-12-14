\name{np.lin.persDx}
\alias{np.lin.persDx}
\title{Nonparametric estimation of linear personalized diagnostic rules.}
\description{Nonparametric estimation of personalized diagnostics rule to find subgroup-specific biomarkers according to linear combination of predictors.}
\usage{
  np.lin.persDx(D,YA,YB,X,dirA,dirB,eps,PLOT)
}
\arguments{
  \item{D}{Binary outcome with D=1 for disease (or case) and D=0 for non-diseased (or control) (n X 1 vector).}
  \item{YA}{Biomarker A, measured on a continuous or orinal scale (n X 1 vector).}
  \item{YB}{Biomarker B, measured on a continuous or ordinal scale (n X 1 vector).}
  \item{X}{Predictors (n x p matrix).}
  \item{dirA}{Direction of YA to D, where dirA="<" (or dirA=">") indicates higher (or lower) YA is assoicated with Pr(D=1)). Default is dirA="<".}
  \item{dirB}{Direction of YB to D, where dirB="<" (or dirB=">") indicates higher (or lower) YB is assoicated with Pr(D=1)). Default is dirB="<".}
  \item{eps}{Tuning parameter for predictor selections. Default is eps=0.01.}
  \item{PLOT}{PLOT=TRUE (or FALSE) shows (or does not show) the receiver operating charactriestics (ROC) curve.}

}
\details{The np.lin.persDx function estimates the personalized diagnostics rule \eqn{\tau(X)}, where \eqn{\tau(X)}=A recommends \eqn{YA} if \eqn{\theta_0+\theta_1 X_1+...+\theta_p X_p > 0} or \eqn{\tau(X)}=B recommends \eqn{YB} otherwise by maximizing (empirical) area under the ROC curve (AUC). Here, AUC is computed based on \eqn{YC} with the direction of "<", i.e. higher \eqn{YC} is associated with Pr(D=1), where \eqn{YC=YA} if \eqn{\tau(X)}=A and dirA="<", or \eqn{YC=YB} if \eqn{\tau(X)}=B and dirB="<". If dirA=">" (or dirB=">"), negative YA (or YB) is used. A forward grid-based search algorithm is used to estimate \eqn{\theta=(\theta_0,\theta_1,...,\theta_p)} by sequentially adding each of the predictors to \eqn{\tau(X)} that increases AUC the most. The stopping criteria is AUC increasement is less than or equal to eps. Cross-validation technique can be used to chose the optimal eps.
}

\value{A list of class np.lin.persDx:
\item{df}{Data frame with YA, YB, X, estimated subgroup (A or B), YC}
\item{AUC.A}{AUC for YA.}
\item{AUC.B}{AUC for YB.}
\item{AUC}{AUC for YA or YB based on \eqn{\tau(X)}.}
\item{theta}{Estimated regression parameters.}
\item{PLOT}{TRUE or FALSE to show ROC curve.}
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

fit=np.lin.persDx(D, YA, YB, X)
fit
}
