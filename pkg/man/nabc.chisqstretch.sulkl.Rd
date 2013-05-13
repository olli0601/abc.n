\name{nabc.chisqstretch.sulkl}
\alias{nabc.chisqstretch.sulkl}
\title{Compute the density of the (possibly truncated) summary likelihood for population dispersion of normal summary values}
\usage{
  nabc.chisqstretch.sulkl(rho, n.of.x, s.of.x, norm = 1,
    support = c(0, Inf), log = FALSE)
}
\arguments{
  \item{rho}{vector of quantile}

  \item{norm}{scalar, 0<\code{norm}<=1, normalization
  constant for the truncated summary likelihood.}

  \item{support}{vector of dimension 2, support of the
  truncated summary likelihood.}

  \item{log}{logical; if \code{TRUE}, densities d are given
  as log(d).}

  \item{n.of.x}{number of observed summary values}

  \item{s.of.x}{standard deviation of observed summary
  values}
}
\description{
  Compute the density of the (possibly truncated) summary
  likelihood for population dispersion of normal summary
  values
}
\note{
  The summary likelihood can be truncated to \code{support}
  and then standardized with \code{norm}. For computational
  efficiency, both \code{norm} and \code{support} must be
  provided although each one can be derived from the other.
  \code{support=qigamma(c(1-norm,1+norm)/2,(n.of.x-2)/2,s.of.x^2*(n.of.x-1)/2)}
  and
  \code{norm=diff(pigamma(support,(n.of.x-2)/2,s.of.x^2*(n.of.x-1)/2)}.
}

