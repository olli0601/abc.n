\name{corrz.criticalregion}
\alias{corrz.criticalregion}
\title{Compute the ABC tolerances of the asymptotic equivalence test for autocorrelations at lag 1}
\usage{
corrz.criticalregion(tau.u, sigma, alpha = 0.01)
}
\arguments{
  \item{tau.u}{Upper tolerance of the equivalence region.}

  \item{sigma}{Standard deviation of the test statistic.}

  \item{alpha}{Level of the equivalence test.}
}
\value{
vector of length 2 \item{cl}{lower point of the critical
region, i.e. lower standard ABC tolerance} \item{cu}{upper
point of the critical region, i.e. upper standard ABC
tolerance}
}
\description{
Compute the ABC tolerances of the asymptotic equivalence
test for autocorrelations at lag 1
}
\examples{
tau.u<- 0.09
	tau.l<- -tau.u
	sim.n<-	5e3
leave.out<- 2
corrz.criticalregion(tau.u, floor(sim.n / (1+leave.out)), alpha=0.01)
}
\seealso{
\code{\link{corrz.calibrate}}
}

