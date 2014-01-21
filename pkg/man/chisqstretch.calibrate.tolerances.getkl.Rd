\name{chisqstretch.calibrate.tolerances.getkl}
\alias{chisqstretch.calibrate.tolerances.getkl}
\title{KL divergence between the summary likelihood and the power function of \code{chisqstretch}}
\usage{
chisqstretch.calibrate.tolerances.getkl(n.of.x, s.of.x, scale, df, tau.u,
  mx.pw = 0.9, alpha = 0.01, pow_scale = 1.5, calibrate.tau.u = T,
  plot = F)
}
\arguments{
  \item{mx.pw}{Power at the point of reference rho.star=1
  (only used when \code{calibrate.tau.u==TRUE}).}

  \item{alpha}{Level of the equivalence test}

  \item{calibrate.tau.u}{If \code{TRUE} the upper tolerance
  of the equivalence region (\code{tau.u}) is calibrated so
  that power at the point of reference is equal to
  \code{mx.pw}}

  \item{tau.u}{Upper tolerance of the equivalence region.
  If \code{calibrate.tau.u==TRUE}, \code{tau.u} is just a
  guess on an upper bound on the upper tolerance of the
  equivalence region to speed up calibration.}

  \item{pow_scale}{Used to set the support of the pdf
  associated to the power function. The power is truncated
  between \code{[tau.l/pow_scale,tau.u*pow_scale]} and then
  standardized.}

  \item{plot}{Logical. If \code{plot==TRUE}, the power of
  the calibrated test is plotted along with the summary
  likelihood.}

  \item{n.of.x}{Number of observed summary values}

  \item{s.of.x}{Standard deviation of the summary values}

  \item{scale}{scaling parameter; typically \code{n}, the
  number of observed summary values}

  \item{df}{degrees of freedom; typically \code{m-1} where
  \code{m} is the number of simulated summary values}
}
\value{
vector of length 6 \item{KL_div}{the Kullback Leibler
divergence} \item{tau.l}{lower tolerance of the equivalence
region} \item{tau.u}{upper tolerance of the equivalence
region} \item{c.l}{lower point of the critical region, i.e.
lower standard ABC tolerance} \item{c.u}{upper point of the
critical region, i.e. upper standard ABC tolerance}
\item{pw.cmx}{actual maximum power at the point of
equality}
}
\description{
Compute the Kullback-Leibler divergence between the summary
likelihood and the power function of the
\code{chisqstretch} equivalence test. The KL divergence is
required to calibrate the number of simulated data points
of the test.
}
\note{
Whatever the value of \code{calibrate.tau.u}, the lower
tolerance of the equivalence region (\code{tau.l}) is
always numerically calibrated so that the mode of the power
function is at the point of equality rho.star.
}
\examples{
#
#	illustrate computation of the KL divergence
#
n.of.x		<- 60 
n.of.y		<- 80
x			<- rnorm(n.of.x,0,1)
df			<- n.of.y-1
#	upper tolerance is not calibrated
#	the power function plateaus around 1, which is not desirable 
abc.param	<- chisqstretch.calibrate.tolerances.getkl(length(x), sd(x), length(x), df, 3*sd(x), mx.pw=0.9, alpha=0.01, pow_scale=1.5, calibrate.tau.u=FALSE, plot=TRUE)
abc.param
#	upper tolerance is calibrated
#	the power function has a nice peak at 1
abc.param	<- chisqstretch.calibrate.tolerances.getkl(length(x), sd(x), length(x), df, 3*sd(x), mx.pw=0.9, alpha=0.01, pow_scale=1.5, calibrate.tau.u=TRUE, plot=TRUE)
abc.param
}

