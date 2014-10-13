\name{chisqstretch.calibrate.taulow}
\alias{chisqstretch.calibrate.taulow}
\title{Calibrate the lower tolerance interval of the equivalence region for \code{chisqstretch}}
\usage{
chisqstretch.calibrate.taulow(tau.up, scale, df, alpha = 0.01, rho.star = 1,
  tol = 1e-05, max.it = 100, pow_scale = 1.5, verbose = 0)
}
\arguments{
  \item{alpha}{size of the test}

  \item{tau.up}{Upper tolerance of the equivalence region}

  \item{rho.star}{point of reference. Defaults to the point
  of equality \code{rho.star=1}.}

  \item{tol}{this algorithm stops when the actual point of
  reference is less than 'tol' from 'rho.star'}

  \item{max.it}{this algorithm stops prematurely when the
  number of iterations to find the equivalence region
  exceeds 'max.it'}

  \item{pow_scale}{Used to set the support of the power
  function. The power is truncated between
  \code{[tau.l/pow_scale,tau.u*pow_scale]} and then
  standardized.}

  \item{verbose}{Logical. If \code{verbose==TRUE}, details
  of the calibration are printed to the console.}

  \item{scale}{scaling parameter; typically \code{n}, the
  number of observed summary values}

  \item{df}{degrees of freedom; typically \code{m-1} where
  \code{m} is the number of simulated summary values}
}
\value{
vector of length 4 \item{tau.low}{lower tolerance of the
equivalence region} \item{cl}{lower point of the critical
region, i.e. lower standard ABC tolerance} \item{cu}{upper
point of the critical region, i.e. upper standard ABC
tolerance} \item{error}{actual error between the mode of
the power function and rho.star}
}
\description{
This function calibrates the lower tolerance interval of
the equivalence region for the \code{chisqstretch}
equivalence test so that the mode of the power function is
at \code{rho.star=1}.
}
\examples{
#
#	illustrate calibration of tau.l for given tau.u such that the mode of the
#	power function is at rho.star
#
alpha		<- 0.01
scale		<- n.of.x<- 60
df			<- 299
tmp			<- chisqstretch.calibrate.taulow(2, scale, df, alpha, rho.star=1, tol= 1e-5, max.it=100, verbose=1)
rho			<- seq(0.1, 5, len=1024)
pw			<- chisqstretch.pow(rho, n.of.x, df, tmp["cl"], tmp["cu"])
plot(rho,pw,col="black", type='l')
tmp			<- chisqstretch.calibrate.taulow(1.8, scale, df, alpha, rho.star=1, tol= 1e-5, max.it=100, verbose=1)
pw			<- chisqstretch.pow(rho, n.of.x, df, tmp["cl"], tmp["cu"])
lines(rho,pw,col="blue")
tmp			<- chisqstretch.calibrate.taulow(1.5, scale, df, alpha, rho.star=1, tol= 1e-5, max.it=100, verbose=1)
pw			<- chisqstretch.pow(rho, n.of.x, df, tmp["cl"], tmp["cu"])
lines(rho,pw,col="red")
tmp			<- chisqstretch.calibrate.taulow(1.3, scale, df, alpha, rho.star=1, tol= 1e-5, max.it=100, verbose=1)
pw			<- chisqstretch.pow(rho, n.of.x, df, tmp["cl"], tmp["cu"])
lines(rho,pw,col="green")		
}
\seealso{
\code{\link{chisqstretch.calibrate}}
}

