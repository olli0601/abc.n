\name{nabc.chisqstretch.pow}
\alias{nabc.chisqstretch.pow}
\title{Compute the density of the (possible truncated) power of the equivalence test for population dispersion of normal summary values}
\usage{
  nabc.chisqstretch.pow(rho, scale, df, c.l, c.u, norm = 1,
    trafo = 1, support = c(0, Inf), log = FALSE)
}
\arguments{
  \item{rho}{vector of quantile}

  \item{scale}{scaling of T apart from rho, either n-1 for
  unbiased ABC or n for exact MAP}

  \item{df}{degrees of freedom}

  \item{c.l,c.u}{lower and upper ABC tolerance}

  \item{norm}{normalization constant for the truncated
  power function.}

  \item{support}{vector of dimension 2. Support of the
  truncated power function.}

  \item{log}{logical; if \code{TRUE}, densities d are given
  as log(d).}
}
\description{
  Compute the density of the (possible truncated) power of
  the equivalence test for population dispersion of normal
  summary values
}
\note{
  The summary likelihood can be truncated to \code{support}
  and then standardized with \code{norm}. For computational
  efficiency, both \code{norm} and \code{support} must be
  provided although each one can be derived (numerically)
  from the other.
}

