\name{mutost.calibrate.tolerances}
\alias{mutost.calibrate.tolerances}
\title{Calibrate the equivalence region for muTOST}
\usage{
mutost.calibrate.tolerances(mx.pw, df, s.of.T, tau.up.ub, alpha, rho.star = 0,
  tol = 1e-05, max.it = 100, debug = 0)
}
\arguments{
  \item{mx.pw}{maximum power at the point of reference
  (rho.star).}

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
equivalence region} \item{3}{actual maximum power
associated with the equivalence region} \item{4}{error ie
abs(actual power - mx.pw)}
}
\description{
Calibrate the equivalence region for the test of location
equivalence for given maximum power
}
\examples{
yn<- 60; ysigma2<- 1; alpha<- 0.01
	mutost.calibrate.tolerances(0.9, yn-1, sqrt(ysigma2/yn), 2, alpha )
}

