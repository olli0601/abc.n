\name{nabc.chisqstretch.n.of.y}
\alias{nabc.chisqstretch.n.of.y}
\title{Calibrate the number of simulated summary values and the equivalence region for the test of dispersion equivalence}
\usage{
  nabc.chisqstretch.n.of.y(n.of.x, s.of.Sx, mx.pw, alpha,
    tau.u.ub = 2, tol = 1e-05, max.it = 100)
}
\arguments{
  \item{n.of.x}{number of observed summary values}

  \item{s.of.Sx}{standard deviation in the observed summary
  likelihood}

  \item{mx.pw}{maximum power at the point of reference
  (rho.star).}

  \item{alpha}{level of the equivalence test}

  \item{tau.up.ub}{guess on an upper bound on the upper
  tolerance of the equivalence region}

  \item{tol}{this algorithm stops when the actual variation
  in the ABC approximation to the summary likelihood is
  less than 'tol' from 's.of.Sx*s.of.Sx'}

  \item{max.it}{this algorithm stops prematurely when the
  number of iterations to calibrate the number of simulated
  data points exceeds 'max.it'}
}
\value{
  vector of length 8 \item{1}{number of simulated summary
  values} \item{2}{lower tolerance of the equivalence
  region} \item{3}{upper tolerance of the equivalence
  region} \item{4}{lower ABC tolerance c^-} \item{5}{upper
  ABC tolerance c^+} \item{6}{actual variation of the
  power} \item{7}{actual maximum power associated with the
  equivalence region} \item{8}{error ie abs(actual
  variation - variation in the observed summary
  likelihood)}
}
\description{
  Calibrate the number of simulated summary values and the
  equivalence region for the test of dispersion equivalence
}
\examples{
xn<- 60; alpha	<- 0.01; prior.u	<- 3; prior.l	<- 1/3; tau.u<- 2.5; xsig2	<- 1
	#summary likelihood of sigma2 given sample mean and sum of squares
	th		<- seq(prior.l,prior.u,length.out=1e3)
	shape		<- (xn-2)/2
	scale		<- xsig2*xn*xn/(xn-1)/2
	y		<- densigamma(th, shape, scale)
	var.Sx	<- scale*scale/((shape-1)*(shape-1)*(shape-2))
	#abc approximation to summary likelihood
	nabc.chisqstretch.n.of.y(xn, sqrt(var.Sx), 0.9, alpha, tau.u.ub=tau.u)
	yn		<- tmp[1]
	tau.l	<- tmp[2]
	tau.u	<- tmp[3]
	c.l		<- tmp[4]
	c.u		<- tmp[5]
	y2		<- nabc.chisqstretch.pow(th, yn-1, c.l, c.u)
#plot the summary likelihood and the abc approximation
plot(th,y/mean(y),ylim=range(c(y/mean(y),y2/mean(y2))),type='l')
lines(th,y2/mean(y2),col="blue")
}

