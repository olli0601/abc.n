\name{chisqstretch.pow.norm}
\alias{chisqstretch.pow.norm}
\title{Area under the power function of the variance test for normal summary values}
\usage{
chisqstretch.pow.norm(scale, df, c.l, c.u, trafo = 1, support = c(0, Inf))
}
\arguments{
  \item{scale}{scaling parameter; typically \code{n}, the
  number of observed summary values}

  \item{df}{degrees of freedom; typically \code{m-1} where
  \code{m} is the number of simulated summary values}

  \item{c.l}{lower point of the critical region of the test
  (equivalent to the lower standard ABC tolerance
  \code{epsilon^-})}

  \item{c.u}{upper point of the critical region of the test
  (equivalent to the upper standard ABC tolerance
  \code{epsilon^+})}

  \item{trafo}{Parameter transformation to translate the
  power to \code{nu} space.}

  \item{support}{vector of dimension 2. Support of the
  truncated power function.}
}
\description{
This function computes the area under the power function
\code{chisqstretch.pow.norm}.
}
\seealso{
\code{\link{chisqstretch.pow}},
\code{\link{chisqstretch.calibrate}}
}

