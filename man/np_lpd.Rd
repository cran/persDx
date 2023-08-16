\name{np_lpd}
\alias{np_lpd}
\title{Nonparametric estimation of linear personalized diagnostic rules.}
\description{Nonparametric estimation of the personalized diagnostics rule to find subgroup-specific biomarkers according to linear combination of predictors.}
\usage{
  np_lpd(D,YA,YB,X,dirA,dirB,eps,plot,A,B,c,d)
}
\arguments{
  \item{D}{Binary outcome with D=1 for disease (or case) and D=0 for non-diseased (or control) (n X 1 vector).}
  \item{YA}{Biomarker A, measured on a continuous scale (n X 1 vector).}
  \item{YB}{Biomarker B, measured on a continuous scale (n X 1 vector).}
  \item{X}{Predictors (n x p matrix).}
  \item{dirA}{Direction of YA to D, where dirA="<" (or dirA=">") indicates higher (or lower) YA is assoicated with Pr(D=1)). Default is dirA="<".}
  \item{dirB}{Direction of YB to D, where dirB="<" (or dirB=">") indicates higher (or lower) YB is assoicated with Pr(D=1)). Default is dirB="<".}
  \item{eps}{Tuning parameter for predictor selections. Default is eps=0.01.}
  \item{plot}{plot=TRUE (or FALSE) shows (or does not show) the receiver operating charactriestics (ROC) curve.}
  \item{A}{Grid search parameter (Discrete). Default is A=0}
  \item{B}{Grid search parameter (Discrete). Default is B=0}
  \item{c}{Grid search parameter. Default is c=2}
  \item{d}{Grid search parameter. Default is d=2}
}
\details{The np.lin.persDx function estimates the personalized diagnostics rule \eqn{\tau(X)}, where \eqn{\tau(X)}=A recommends \eqn{YA} if \eqn{\theta_1 X_1+...+\theta_p X_p > \theta_0} or \eqn{\tau(X)}=B recommends \eqn{YB} otherwise by maximizing (empirical) area under the ROC curve (AUC). Here, the AUC is computed based on \eqn{YC} with the direction of "<", i.e. higher \eqn{YC} is associated with Pr(D=1), where \eqn{YC=YA} if \eqn{\tau(X)}=A and dirA="<", or \eqn{YC=YB} if \eqn{\tau(X)}=B and dirB="<". If dirA=">" (or dirB=">"), negative YA (or YB) is used.

A forward grid rotation algorithm (FGR) is used to estimate \eqn{\theta_0,\theta_1,...,\theta_p} by sequentially adding each of the predictors to \eqn{\tau(X)} that increases the AUC the most. The stopping criteria is AUC increasement is less than or equal to eps. The eps controls the model complexity. The cross-validation techniques can be used to find the optimal eps.

The FGR results in a suboptimal solution. The accuracy is improved by setting higher A, B, c, d, but it increases increase computational costs, or vice versa. We thus recond this function when p is small or around 10.
}

\value{A list of class np.lin.persDx:
\item{df}{Data frame with D, YA, YB, X, tau, YC, where tau=A or B for recommending YA or YB, respectively.}
\item{AUCA}{AUC for YA.}
\item{AUCB}{AUC for YB.}
\item{AUC}{AUC for YC.}
\item{tpfp}{Data frame with cutoff, tp, fp, where tp and fp are true and false positive positives at the cutoff values of YC.}
\item{theta}{Estimated regression parameters.}
\item{theta0}{Estimated threshold parameter.}
\item{PLOT}{TRUE or FALSE to show ROC curves.}
}

\references{Yaliang Zhang and Yunro Chung, Nonparametric estimation of linear personalized diagnostics rules via efficient grid algorithm (submitted)}
\author{Yunro Chung [aut, cre]}

\examples{
#simulate data
set.seed(1)
n=100
D=c(rep(1,n/2),rep(0,n/2))

X1=runif(n,0,1)
X2=runif(n,0,1)
X3=runif(n,0,1)
X=data.frame(X1,X2,X3)

tau=rep("B",n)
tau[X1+X2>=1]="A"

YA=D*(rnorm(n,2,1)*(tau=="A")+rnorm(n,0,1)*(tau=="B"))+
   (1-D)*rnorm(n,0,1)
YB=D*(rnorm(n,1,1)*(tau=="B")+rnorm(n,0,1)*(tau=="A"))+
   (1-D)*rnorm(n,0,1)

#run
fit=np_lpd(D, YA, YB, X)
fit
}
