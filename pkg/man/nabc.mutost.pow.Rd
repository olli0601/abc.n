\name{nabc.mutost.pow}
\alias{nabc.mutost.pow}
\title{Compute power of the equivalence test for population means of normal summary values}
\usage{
  nabc.mutost.pow(rho, df, tau.u, s.of.T, alpha,
    rtn.fun = FALSE)
}
\arguments{
  \item{rho}{true difference in simulated and observed
  population means}

  \item{df}{degrees of freedom of the simulated summary
  values}

  \item{tau.u}{upper tolerance of the equivalence region}

  \item{s.of.T}{standard deviation of the test statistic}

  \item{alpha}{level of the equivalence test}

  \item{rtn.fun}{indicator if a function to compute the
  power should be returned. Defaults to 0.}
}
\value{
  approximate power of the exact test. this is approximate
  because the standard deviation of the normal model for
  the simulated summary values is not known.
}
\description{
  Compute power of the equivalence test for population
  means of normal summary values
}
\examples{
prior.u<- 5; prior.l<- -prior.u; tau.u	<- 0.75; yn<- 60; ysigma2<- 1; alpha<- 0.01
rho	<- seq(prior.l,prior.u,length.out=1e3)
nabc.mutost.pow(rho, yn-1, tau.u, sqrt(ysigma2/yn), alpha)
}

