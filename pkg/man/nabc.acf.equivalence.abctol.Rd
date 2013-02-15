\name{nabc.acf.equivalence.abctol}
\alias{nabc.acf.equivalence.abctol}
\title{Compute the ABC tolerances of the asymptotic equivalence test for autocorrelations at lag 1}
\usage{
  nabc.acf.equivalence.abctol(tau.l, tau.u, n, alpha)
}
\arguments{
  \item{tau.l}{lower tolerance of the equivalence region}

  \item{tau.u}{upper tolerance of the equivalence region}

  \item{n}{number of pairs (x_i,x_i-1) after thinning of
  the time series x_1, x_2, ...}

  \item{alpha}{level of the equivalence test}
}
\value{
  vector of length 2, first entry is lower ABC tolerance,
  second entry is upper ABC tolerance
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
nabc.acf.equivalence.abctol(tau.l, tau.u, floor(sim.n / (1+leave.out)), 0.01)
}

