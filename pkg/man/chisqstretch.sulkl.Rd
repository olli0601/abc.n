\name{chisqstretch.sulkl}
\alias{chisqstretch.sulkl}
\title{Density of the summary likelihood of the variance for normal summary values}
\usage{
chisqstretch.sulkl(rho, n.of.x, s.of.x, trafo = (n.of.x - 1)/n.of.x * s.of.x *
  s.of.x, norm = 1, support = c(0, Inf), log = FALSE)
}
\arguments{
  \item{rho}{Auxiliary error parameter}

  \item{n.of.x}{Number of observed summary values}

  \item{s.of.x}{Standard deviation of the summary values}

  \item{trafo}{Parameter transformation to compute the
  summary likelihood on the \code{rho} error space.}

  \item{norm}{scalar, 0<\code{norm}<=1, normalization
  constant for the truncated summary likelihood.}

  \item{support}{vector of dimension 2, support of the
  truncated summary likelihood.}

  \item{log}{logical; if \code{TRUE}, densities d are given
  as log(d).}
}
\description{
The density of the summary likelihood of the variance for
normal summary values is scaled inverse chi square.
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
\seealso{
\code{\link{chisqstretch.calibrate}},
\code{\link{chisqstretch.calibrate.tolerances.getkl}} for
an example.
}

