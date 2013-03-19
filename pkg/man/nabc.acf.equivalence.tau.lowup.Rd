\name{nabc.acf.equivalence.tau.lowup}
\alias{nabc.acf.equivalence.tau.lowup}
\title{Calibrate the equivalence region of the asymptotic equivalence test for autocorrelations at lag 1 for given maximum power}
\usage{
  nabc.acf.equivalence.tau.lowup(mx.pw, tau.up.ub, n,
    alpha, rho.star = 0, tol = 1e-05, max.it = 100)
}
\arguments{
  \item{mx.pw}{maximum power at the point of reference
  (rho.star).}

  \item{tau.up.ub}{guess on an upper bound on the upper
  tolerance of the equivalence region}

  \item{n}{number of pairs (x_i,x_i-1) after thinning of
  the time series x_1, x_2, ...}

  \item{alpha}{level of the equivalence test}

  \item{rho.star}{point of reference. Defaults to the point
  of equality rho.star=0.}

  \item{tol}{this algorithm stops when the actual maximum
  power is less than 'tol' from 'mx.pw'}

  \item{max.it}{this algorithm stops prematurely when the
  number of iterations to find the equivalence region
  exceeds 'max.it'}
}
\value{
  vector of length 4 \item{1}{lower tolerance of the
  equivalence region} \item{2}{upper tolerance of the
  equivalence region} \item{3}{actual maximum power
  associated with the equivalence region} \item{4}{error ie
  abs(actual power - mx.pw)}
}
\description{
  Calibrate the equivalence region of the asymptotic
  equivalence test for autocorrelations at lag 1 for given
  maximum power
}
\examples{
tau.u<- 0.09
	tau.l<- -tau.u
	sim.n<-	5e3
	leave.out<- 2
	nabc.acf.equivalence.tau.lowup(0.9, 2, floor(sim.n / (1+leave.out)), 0.01)
}

