\name{corrz.calibrate.tolerances}
\alias{corrz.calibrate.tolerances}
\title{Calibrate the equivalence region of the asymptotic equivalence test for autocorrelations at lag 1 for given maximum power}
\usage{
corrz.calibrate.tolerances(mx.pw, tau.up.ub, sigma, alpha = 0.01,
  rho.star = 0, tol = 1e-05, max.it = 100, pow_scale = 2, verbose = 0)
}
\arguments{
  \item{mx.pw}{maximum power at the point of reference
  (rho.star).}

  \item{pow_scale}{Used to set the support of the pdf
  associated to the power function. The power is truncated
  between \code{[tau.l/pow_scale,tau.u*pow_scale]} and then
  standardized.}

  \item{tau.up.ub}{guess on an upper bound on the upper
  tolerance of the equivalence region}

  \item{rho.star}{point of reference. Defaults to the point
  of equality rho.star=0.}

  \item{tol}{this algorithm stops when the actual maximum
  power is less than 'tol' from 'mx.pw'}

  \item{max.it}{this algorithm stops prematurely when the
  number of iterations to find the equivalence region
  exceeds 'max.it'}

  \item{verbose}{Logical. If \code{verbose==TRUE}, detailed
  calibration output is printed to the console}

  \item{sigma}{Standard deviation of the test statistic.}

  \item{alpha}{Level of the equivalence test.}
}
\value{
vector of length 6 \item{tau.l}{lower tolerance of the
equivalence region} \item{tau.u}{upper tolerance of the
equivalence region} \item{curr.mx.pw}{actual maximum power
associated with the equivalence region} \item{error}{actual
error between the maximum power and mx.pw} \item{cl}{lower
point of the critical region, i.e. lower standard ABC
tolerance} \item{cu}{upper point of the critical region,
i.e. upper standard ABC tolerance}
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
	corrz.calibrate.tolerances(0.9, 2, 1/sqrt(floor(sim.n / (1+leave.out))-3), 0.01)
}
\seealso{
\code{\link{corrz.calibrate}}
}

