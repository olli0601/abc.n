\name{chisqstretch.su.lkl.norm}
\alias{chisqstretch.su.lkl.norm}
\title{Area under the summary likelihood of the variance for normal summary values}
\usage{
chisqstretch.su.lkl.norm(n.of.x, s.of.x, trafo = 1, support = c(0, Inf))
}
\arguments{
  \item{n.of.x}{Number of observed summary values}

  \item{s.of.x}{Standard deviation of the summary values}

  \item{trafo}{Parameter transformation to compute the
  summary likelihood on the \code{rho} error space.}

  \item{support}{vector of dimension 2, support of the
  truncated summary likelihood.}
}
\description{
This function computes the area under the summary
likelihood \code{chisqstretch.sulkl}.
}
\seealso{
\code{\link{chisqstretch.calibrate}}
}

