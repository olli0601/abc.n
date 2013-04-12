\name{nabc.mutost.onesample.tau.lowup.var}
\alias{nabc.mutost.onesample.tau.lowup.var}
\title{Calibrate the equivalence region for the test of location equivalence for given variance of the summary likelihood}
\usage{
  nabc.mutost.onesample.tau.lowup.var(s.of.Sx, df, s.of.T,
    tau.up.ub, alpha, rho.star = 0, tol = 1e-05,
    max.it = 100, debug = 0)
}
\arguments{
  \item{s.of.Sx}{standard deviation of the summary
  likelihood}

  \item{df}{degrees of freedom}

  \item{s.of.T}{standard deviation of the test statistic}

  \item{tau.up.ub}{guess on an upper bound on the upper
  tolerance of the equivalence region}

  \item{alpha}{level of the equivalence test}

  \item{rho.star}{point of reference. Defaults to the point
  of equality rho.star=0.}

  \item{tol}{this algorithm stops when the actual maximum
  power is less than 'tol' from 'mx.pw'}

  \item{max.it}{this algorithm stops prematurely when the
  number of iterations to find the equivalence region
  exceeds 'max.it'}

  \item{debug}{Flag if C implementation is used.}
}
\value{
  vector of length 4 \item{1}{lower tolerance of the
  equivalence region} \item{2}{upper tolerance of the
  equivalence region} \item{3}{actual variance associated
  with the power} \item{4}{error ie abs(actual var(power) -
  var(summary likelihood))}
}
\description{
  Calibrate the equivalence region for the test of location
  equivalence for given variance of the summary likelihood
}
\examples{
yn<- 60; ysigma2<- 1; alpha<- 0.01
	nabc.mutost.onesample.tau.lowup.var(0.002, yn-1, sqrt(ysigma2/yn), 2, alpha )
}

