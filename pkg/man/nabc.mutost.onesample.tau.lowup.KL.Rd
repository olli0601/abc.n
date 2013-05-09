\name{nabc.mutost.onesample.tau.lowup.KL}
\alias{nabc.mutost.onesample.tau.lowup.KL}
\title{Calibrate the number of simulated summary values and the equivalence region for the test of location equivalence by minimising the Kullback-Leibler divergence between the power function and the summary likelihood.}
\usage{
  nabc.mutost.onesample.tau.lowup.KL(n.of.x, s.of.x,
    n.of.y, s.of.y, mx.pw, alpha, tau.u.ub, debug = 0,
    plot = F)
}
\arguments{
  \item{plot}{if \code{TRUE}, the summary likelihood and
  abc approximation are plotted}

  \item{max.it}{this algorithm stops prematurely when the
  number of iterations to calibrate the number of simulated
  data points exceeds 'max.it'}

  \item{debug}{Flag if C implementation is used.}

  \item{n.of.x}{number of observed summary values}

  \item{s.of.x}{standard deviation of observed summary
  values}

  \item{n.of.y}{number of simulated summary values}

  \item{s.of.y}{standard deviation of simulated summary
  values}

  \item{mx.pw}{maximum power at the point of reference
  (rho.star=0) (only when \code{calibrate.tau.u==TRUE}).}

  \item{alpha}{level of the equivalence test}
}
\value{
  vector of length 4 \item{n.of.y}{number of simulated
  summary values} \item{KL_div}{the Kullback Leibler
  divergence} \item{tau.u}{upper tolerance of the
  equivalence region} \item{pw.cmx}{actual maximum power
  associated with the equivalence region}
}
\description{
  Calibrate the number of simulated summary values and the
  equivalence region for the test of location equivalence
  by minimising the Kullback-Leibler divergence between the
  power function and the summary likelihood.
}
\examples{
prior.u<- 2; prior.l<- -prior.u; tau.u<- 0.75; xn<- yn<- 60; xmu<- 0.5;
xsigma2<- ysigma2<- 2; alpha<- 0.01
rho<- seq(prior.l,prior.u,length.out=1e3)
#summary likelihood
y<-dnorm(rho,0,sqrt(xsigma2/xn))
y<- y / diff(pnorm(c(prior.l,prior.u),0,sqrt(xsigma2/xn)))
#abc approximation to summary likelihood based on equivalence test
tmp <- nabc.mutost.onesample.n.of.y.KL(xn, sqrt(xsigma2/xn), yn,
sqrt(ysigma2),0.9, alpha, tau.u.ub=2*tau.u,plot=T)
yn <- tmp["n.of.y"]
tau.u	<- tmp["tau.u"]
y2<- nabc.mutost.pow(rho, yn-1, tau.u, sqrt(ysigma2/yn), alpha)
rho2<- rho[which(y2!=0)]
y2<- y2[which(y2!=0)]
y2<- y2/sum(diff(rho2)*y2[-1])
#plot summary likelihood and abc approximation thereof
plot(1,1,type='n',xlim=range(rho),ylim=range(c(y,y2)),xlab=expression(rho))
lines(rho,y,col="red")
lines(rho2,y2,col="blue")
abline(v=0,col="red")
}

