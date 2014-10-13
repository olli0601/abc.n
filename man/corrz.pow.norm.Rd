\name{corrz.pow.norm}
\alias{corrz.pow.norm}
\title{Area under the power of the asymptotic equivalence test for autocorrelations at lag 1.}
\usage{
corrz.pow.norm(tau.u, sigma, alpha = 0.01, support = c(-Inf, Inf))
}
\arguments{
  \item{tau.u}{Upper tolerance of the equivalence region.}

  \item{sigma}{Standard deviation of the test statistic.}

  \item{alpha}{Level of the equivalence test.}

  \item{support}{Vector of dimension 2. Support of the
  truncated power function.}
}
\value{
Integral of the power function
}
\description{
Area under the power of the asymptotic equivalence test for
autocorrelations at lag 1.
}
\seealso{
\code{\link{corrz.pow}}
}

