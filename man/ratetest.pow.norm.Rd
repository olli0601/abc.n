% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ratetest.R
\name{ratetest.pow.norm}
\alias{ratetest.pow.norm}
\title{Area under the \code{ratetest} power function}
\usage{
ratetest.pow.norm(c.l, c.u, m, support = c(0, Inf))
}
\arguments{
\item{c.l}{Lower boundary point of the critical region (equivalent to the lower ABC tolerance \code{epsilon^-})}

\item{c.u}{Upper boundary point of the critical region (equivalent to the upper ABC tolerance \code{epsilon^+})}

\item{m}{Number of simulated values}

\item{support}{Support of the truncated power function}
}
\description{
This function computes the area under the power function \code{ratetest.pow}.
}
\seealso{
\code{\link{ratetest.pow}}, \code{\link{ratetest.calibrate}}
}

