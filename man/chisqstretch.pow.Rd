\name{chisqstretch.pow}
\alias{chisqstretch.pow}
\title{Power of the variance test for normal summary values}
\usage{
chisqstretch.pow(rho, scale, df, c.l, c.u, norm = 1, trafo = 1,
  support = c(0, Inf), log = FALSE)
}
\arguments{
  \item{rho}{Auxiliary error parameter (vector).}

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

  \item{norm}{normalization constant for the truncated
  power function.}

  \item{support}{vector of dimension 2. Support of the
  truncated power function.}

  \item{log}{l ogical; if \code{TRUE}, densities d are
  given as log(d).}
}
\description{
This function computes the power function of the two-sided,
one-sample scaled chi sqare test for testing the
equivalence of variances for normal summary values.
}
\note{
The summary likelihood can be truncated to \code{support}
and then standardized with \code{norm}. For computational
efficiency, both \code{norm} and \code{support} must be
provided although each one can be derived (numerically)
from the other.
}
\examples{
#
#	illustrate power function of the variance test for normal summary values
#
alpha		<- 0.01
n.of.x		<- 60 
n.of.y		<- 60
df			<- n.of.y-1
rho			<- seq(0.1, 5, len=1024)
#	compute tolerances
tmp			<- .Call("abcScaledChiSq",	c(n.of.x, df, 1/2.2, 2.2, alpha, 1e-10, 100, 0.05)	)
#	compute power function for given tolerances
pw			<- chisqstretch.pow(rho, n.of.x, df, tmp[1], tmp[2])
plot(rho,pw,col="red", type='l')
lines(rho,pw,col="red")
#	repeat for more simulated summary values
n.of.y		<- 110
df			<- n.of.y-1
tmp			<- .Call("abcScaledChiSq",	c(n.of.x, df, 1/2.2, 2.2, alpha, 1e-10, 100, 0.05)	)
pw			<- chisqstretch.pow(rho, n.of.x, df, tmp[1], tmp[2])
plot(rho,pw,col="red", type='l')
lines(rho,pw,col="red")
}
\seealso{
\code{\link{chisqstretch.pow.norm}},
\code{\link{chisqstretch.calibrate}}
}

