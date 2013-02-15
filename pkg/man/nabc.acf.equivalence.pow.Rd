\name{nabc.acf.equivalence.pow}
\alias{nabc.acf.equivalence.pow}
\title{Compute power of the asymptotic equivalence test for autocorrelations at lag 1}
\usage{
  nabc.acf.equivalence.pow(rho, tau.u, alpha, s)
}
\arguments{
  \item{rho}{true difference in simulated and observed
  autocorrelation at lag 1}

  \item{tau.u}{upper tolerance of the equivalence region}

  \item{alpha}{level of the equivalence test}

  \item{s}{standard deviation of the test statistic}
}
\value{
  power of the asymptotic test. this is approximate because
  the test is asymptotic
}
\description{
  Compute power of the asymptotic equivalence test for
  autocorrelations at lag 1
}
\examples{
tau.u<- 0.09
	tau.l<- -tau.u
	sim.n<-	5e3
	rho<- seq(tau.l,tau.u,0.001)
	pw<- nabc.acf.equivalence.pow(rho, tau.u, alpha, 1/sqrt(floor(sim.n/3)-3))
}

