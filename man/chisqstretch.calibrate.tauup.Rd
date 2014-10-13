\name{chisqstretch.calibrate.tauup}
\alias{chisqstretch.calibrate.tauup}
\title{Calibrate the upper tolerance interval of the equivalence region for \code{chisqstretch}}
\usage{
chisqstretch.calibrate.tauup(mx.pw, tau.up.ub, scale, df, alpha = 0.01,
  rho.star = 1, tol = 1e-05, max.it = 100, pow.scale = 1.5,
  verbose = 0)
}
\arguments{
  \item{mx.pw}{maximum power at the point of reference
  (rho.star).}

  \item{tau.up.ub}{guess on an upper bound on the upper
  tolerance of the equivalence region}

  \item{rho.star}{point of reference. Defaults to the point
  of equality rho.star=1.}

  \item{tol}{this algorithm stops when the actual maximum
  power is less than 'tol' from 'mx.pw'}

  \item{max.it}{this algorithm stops prematurely when the
  number of iterations to find the equivalence region
  exceeds 'max.it'}

  \item{pow.scale}{Used to set the support of the power
  function. The power is truncated between
  \code{[tau.l/pow_scale,tau.u*pow_scale]} and then
  standardized.}

  \item{scale}{scaling parameter; typically \code{n}, the
  number of observed summary values}

  \item{df}{degrees of freedom; typically \code{m-1} where
  \code{m} is the number of simulated summary values}

  \item{alpha}{size of the test}

  \item{verbose}{Logical. If \code{verbose==TRUE}, details
  of the calibration are printed to the console.}
}
\value{
vector of length 6 \item{tau.low}{lower tolerance of the
equivalence region} \item{tau.up}{upper tolerance of the
equivalence region} \item{curr.mx.pw}{actual maximum power
associated with the equivalence region} \item{error}{actual
error between the power at rho.star and mx.pw}
\item{cl}{lower point of the critical region, i.e. lower
standard ABC tolerance} \item{cu}{upper point of the
critical region, i.e. upper standard ABC tolerance}
}
\description{
This function calibrates the upper tolerance interval of
the equivalence region for the \code{chisqstretch}
equivalence test so that the mode of the power function is
at \code{rho.star=1} and so that the power function at
\code{rho.star} euqals \code{mx.pw}. This involves
recursive recursive calls to re-calibrate the lower
tolerance region.
}
\examples{
#
#	illustrate calibration of tau.l and tau.u for given target maximum power at rho.star
#
alpha		<- 0.01
scale		<- n.of.x<- 60
n.of.y		<- 299
df			<- n.of.y-1
x			<- rnorm(n.of.x,0,1)
s.of.x		<- sd(x)
tmp			<- chisqstretch.calibrate.tauup(0.9, 3*s.of.x, scale, df, alpha, rho.star=1, tol= 1e-5, max.it=100, verbose=1)
rho			<- seq(0.4, 2, len=1024)
pw			<- chisqstretch.pow(rho, n.of.x, df, tmp["cl"], tmp["cu"])		
plot(rho,pw,col="black", type='l')
}
\seealso{
\code{\link{chisqstretch.calibrate}}
}

