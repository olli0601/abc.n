\name{corrz.pow}
\alias{corrz.pow}
\title{Power of the asymptotic equivalence test for autocorrelations at lag 1}
\usage{
corrz.pow(rho, tau.u, alpha, sigma, norm = 1, support = c(-Inf, Inf),
  log = FALSE)
}
\arguments{
  \item{rho}{Auxiliary error parameter (vector).}

  \item{tau.u}{Upper tolerance of the equivalence region.}

  \item{alpha}{Level of the equivalence test.}

  \item{sigma}{Standard deviation of the test statistic.}

  \item{norm}{Normalization constant for the truncated
  power function.}

  \item{support}{Vector of dimension 2. Support of the
  truncated power function.}

  \item{log}{Logical; if \code{TRUE}, densities d are given
  as log(d).}
}
\value{
Power of the test for the values of \code{rho}.
}
\description{
This function computes the power function of the two-sided,
one-sample scaled TOSZ test for testing the equivalence of
autocorrelations for normal summary values. The power
function is not exact because the TOSZ is test is
asymptotic in the number of summary values.
}
\note{
The summary likelihood can be truncated to \code{support}
and then standardized with \code{norm}. For computational
efficiency, both \code{norm} and \code{support} must be
provided although each one can be derived (numerically)
from the other.
}
\examples{
#
#	illustrate power function of the autocorrelation test for normal summary values
#
alpha	<- 0.01
tau.u	<- 0.09
tau.l	<- -tau.u
sim.n	<-	5e3
rho		<- seq(tau.l,tau.u,0.001)
pw		<- corrz.pow(rho, tau.u, alpha, 1/sqrt(floor(sim.n/3)-3))
plot(rho,pw,col="red", ylim=c(0,1), type='l')
#	flat power function when the number of simulated summary values is much lower	
sim.n	<-	5e2
pw		<- corrz.pow(rho, tau.u, alpha, 1/sqrt(floor(sim.n/3)-3))
lines(rho,pw,col="blue")
#	flat power function when the number of simulated summary values is much higher	
sim.n	<-	2e4
pw		<- corrz.pow(rho, tau.u, alpha, 1/sqrt(floor(sim.n/3)-3))
lines(rho,pw,col="green")
}
\seealso{
\code{\link{chisqstretch.calibrate}},
\code{\link{corrz.pow.norm}}
}

