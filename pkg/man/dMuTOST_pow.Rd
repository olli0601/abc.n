\name{dMuTOST_pow}
\alias{dMuTOST_pow}
\title{Compute the density of the (possible truncated) power of the equivalence test for population means of normal summary values}
\usage{
  dMuTOST_pow(rho, df, s.of.T, tau.u, alpha, norm = 1,
    support = c(-Inf, Inf), log = FALSE)
}
\arguments{
  \item{rho}{vector of quantile}

  \item{norm}{normalization constant for the truncated
  power function.}

  \item{support}{vector of dimension 2. Support of the
  truncated power function.}

  \item{log}{logical; if \code{TRUE}, densities d are given
  as log(d).}

  \item{tau.u}{upper tolerance of the equivalence region.
  If \code{calibrate.tau.u==TRUE}, it's just a guess on an
  upper bound on the upper tolerance of the equivalence
  region.}

  \item{alpha}{level of the equivalence test}

  \item{df}{degrees of freedom of the simulated summary
  values}

  \item{s.of.T}{standard deviation of the test statistic}
}
\description{
  Compute the density of the (possible truncated) power of
  the equivalence test for population means of normal
  summary values
}
\note{
  The summary likelihood can be truncated to \code{support}
  and then standardized with \code{norm}. For computational
  efficiency, both \code{norm} and \code{support} must be
  provided although each one can be derived (numerically)
  from the other.
}

