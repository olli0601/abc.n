\name{nabc.chisqstretch.tau.low}
\alias{nabc.chisqstretch.tau.low}
\title{Calibrate the lower tolerance interval of the equivalence region for the test of dispersion equivalence}
\usage{
  nabc.chisqstretch.tau.low(tau.up, df, alpha,
    rho.star = 1, tol = 1e-05, max.it = 100, for.mle = 0)
}
\arguments{
  \item{tau.up}{upper tolerance of the equivalence region}

  \item{df}{degrees of freedom}

  \item{alpha}{level of the equivalence test}

  \item{rho.star}{point of reference. Defaults to the point
  of equality rho.star=1}

  \item{tol}{this algorithm stops when the actual point of
  reference is less than 'tol' from 'rho.star'}

  \item{max.it}{this algorithm stops prematurely when the
  number of iterations to find the equivalence region
  exceeds 'max.it'}

  \item{for.mle}{calibrate so that the mode of the power is
  at the MLE}
}
\value{
  tau.low, lower tolerance of the equivalence region
}
\description{
  Calibrate the lower tolerance interval of the equivalence
  region for the test of dispersion equivalence
}
\examples{
tau.u<- 2.2
 yn<- 60
	tau.l<- nabc.chisqstretch.tau.low(tau.u, yn-1, 0.01)
}

