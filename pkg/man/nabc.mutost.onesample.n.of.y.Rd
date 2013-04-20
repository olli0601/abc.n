\name{nabc.mutost.onesample.n.of.y}
\alias{nabc.mutost.onesample.n.of.y}
\title{Calibrate the number of simulated summary values and the equivalence region for the test of location equivalence}
\usage{
  nabc.mutost.onesample.n.of.y(n.of.x, s.of.Sx, mx.pw,
    s.of.y, alpha, tau.u.ub = 2, tol = 1e-05, max.it = 100,
    debug = 0)
}
\arguments{
  \item{n.of.x}{number of observed summary values}

  \item{s.of.Sx}{standard deviation in the observed summary
  likelihood}

  \item{mx.pw}{maximum power at the point of reference
  (rho.star).}

  \item{s.of.y}{standard deviation in the simulated summary
  values}

  \item{alpha}{level of the equivalence test}

  \item{tau.up.ub}{guess on an upper bound on the upper
  tolerance of the equivalence region}

  \item{tol}{this algorithm stops when the actual variation
  in the ABC approximation to the summary likelihood is
  less than 'tol' from 's.of.Sx*s.of.Sx'}

  \item{max.it}{this algorithm stops prematurely when the
  number of iterations to calibrate the number of simulated
  data points exceeds 'max.it'}

  \item{debug}{Flag if C implementation is used.}
}
\value{
  vector of length 6 \item{1}{number of simulated summary
  values} \item{2}{lower tolerance of the equivalence
  region} \item{3}{upper tolerance of the equivalence
  region} \item{4}{actual variation of the power}
  \item{5}{actual maximum power associated with the
  equivalence region} \item{6}{error ie abs(actual
  variation - variation in the observed summary
  likelihood)}
}
\description{
  Calibrate the number of simulated summary values and the
  equivalence region for the test of location equivalence
}
\examples{
prior.u<- 2; prior.l<- -prior.u; tau.u<- 0.75; xn<- yn<- 60; xmu<- 0.5; xsigma2<- ysigma2<- 2; alpha<- 0.01
	rho<- seq(prior.l,prior.u,length.out=1e3)
	#summary likelihood
	y<-	dnorm(rho,0,sqrt(xsigma2/xn))
	y<- y / diff(pnorm(c(prior.l,prior.u),0,sqrt(xsigma2/xn)))
	#abc approximation to summary likelihood based on equivalence test
	tmp	<- nabc.mutost.onesample.n.of.y(xn, sqrt(xsigma2/xn), 0.9, sqrt(ysigma2), alpha, tau.u.ub=2*tau.u )
	yn	<- tmp[1]
	tau.u	<- tmp[3]
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

