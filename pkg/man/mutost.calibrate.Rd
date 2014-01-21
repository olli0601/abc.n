\name{mutost.calibrate}
\alias{mutost.calibrate}
\title{Calibration for muTOST}
\usage{
mutost.calibrate(KL_args, max.it = 100, debug = 0, plot = FALSE,
  plot_debug = FALSE, verbose = FALSE)
}
\arguments{
  \item{KL_args}{list of arguments needed to run
  \code{\link{mutost.calibrate.tolerances.getkl}} (see note
  below)}

  \item{max.it}{this algorithm stops prematurely when the
  number of iterations to find the equivalence region
  exceeds 'max.it'}

  \item{debug}{flag if C implementation is used}

  \item{plot}{logical; if \code{TRUE} the summary
  likelihood and calibrated power are ploted (only when
  \code{debug==TRUE}).}

  \item{plot_debug}{logical; if \code{TRUE} the summary
  likelihood and power are ploted at each calibration step
  (only when \code{debug==TRUE}).}

  \item{verbose}{flag if detailed information on the
  computations should be printed to standard out}
}
\value{
vector of length 4 \item{KL_div}{the Kullback Leibler
divergence} \item{n.of.y}{number of simulated summary
values} \item{tau.u}{upper tolerance of the equivalence
region} \item{pw.cmx}{actual maximum power associated with
the equivalence region}
}
\description{
Calibrate the equivalence region for the test of location
equivalence for given maximum power
}
\note{
\code{KL_args} is a named list with the following elements:
\code{n.of.x}, \code{s.of.x}, \code{n.of.y}, \code{s.of.y},
\code{mx.pw}, \code{alpha}, \code{tau.u} and
\code{pow_scale}.
}
\examples{
xn <- 60
yn <- xn
xmean <- 1
ymean <- 1
xsigma <- 1

obs <- rnorm(xn, xmean, xsigma)
obs <- (obs - mean(obs))/sd(obs) * xsigma + xmean

mx.pw <- 0.9
tau.u.ub <- 0.5
alpha <- 0.01

# Example 1: Same mean but sd(sim)>sd(obs): the power function is not "too tight" for yn=xn. 
# We increase ny while calibrating tau.u and mx.pw in order to minimize KL.
ysigma <- 1.5
sim <- rnorm(yn, ymean, ysigma)

KL_args <- list(n.of.x=length(obs), s.of.x= sd(obs), n.of.y=length(sim), s.of.y=sd(sim), mx.pw=mx.pw, alpha=alpha, tau.u=tau.u.ub, pow_scale=1.5, debug=0)

tmp <- mutost.calibrate(KL_args, max.it = 100, debug = TRUE, plot = TRUE, plot_debug = FALSE, verbose=TRUE) 

# Example 2: Same mean but sd(obs)>sd(sim): the power function is "too tight" for yn=xn.
# In this case, the only way to minimize KL is to use yn<xn, which is not what we like.
# Instead, we give up on mx.pw and calibrate tau.u.
ysigma <- 0.5
sim <- rnorm(yn, ymean, ysigma)
KL_args$s.of.y <- sd(sim)
tmp <- mutost.calibrate(KL_args, max.it = 100, debug = TRUE, plot = TRUE, plot_debug = FALSE, verbose=TRUE) 
}
\seealso{
\code{\link{mutost.calibrate.tolerances.getkl}}
}

