\name{corrz.sulkl}
\alias{corrz.sulkl}
\title{Density of the summary likelihood of the autocorrelation error for normal summary values}
\usage{
corrz.sulkl(rho, sigma, norm = 1, support = c(-Inf, Inf), log = FALSE)
}
\arguments{
  \item{rho}{Auxiliary error parameter (vector).}

  \item{sigma}{Standard deviation of the test statistic.}

  \item{norm}{scalar, 0<\code{norm}<=1, normalization
  constant for the truncated summary likelihood.}

  \item{support}{vector of dimension 2, support of the
  truncated summary likelihood.}

  \item{log}{logical; if \code{TRUE}, densities d are given
  as log(d).}
}
\value{
Summary likelihood for the error parameter \code{rho}
}
\description{
The density of the summary likelihood of the
autocorrelation error is normal.
}
\note{
The summary likelihood can be truncated to \code{support}
and then standardized with \code{norm}. For computational
efficiency, both \code{norm} and \code{support} must be
provided although each one can be derived from the other.
}
\seealso{
\code{\link{chisqstretch.calibrate}},
\code{\link{corrz.sulkl.norm}}
}

