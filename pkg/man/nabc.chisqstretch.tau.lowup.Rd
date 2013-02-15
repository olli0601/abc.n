\name{nabc.chisqstretch.tau.lowup}
\alias{nabc.chisqstretch.tau.lowup}
\title{Calibrate the equivalence region for the test of dispersion equivalence for given maximum power}
\usage{
  nabc.chisqstretch.tau.lowup(mx.pw, tau.up.ub, df, alpha,
    rho.star = 1, tol = 1e-05, max.it = 100)
}
\arguments{
  \item{mx.pw}{maximum power at the point of reference
  (rho.star).}

  \item{tau.up.ub}{guess on an upper bound on the upper
  tolerance of the equivalence region}

  \item{df}{degrees of freedom}

  \item{alpha}{level of the equivalence test}

  \item{rho.star}{point of reference. Defaults to the point
  of equality rho.star=1.}

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
  Calibrate the equivalence region for the test of
  dispersion equivalence for given maximum power
}
\examples{
yn<- 60
	nabc.chisqstretch.tau.lowup(0.9, 2.5, yn-1, 0.01)
}

