\name{nabc.calibrate.m.and.tau.yesmxpw.yesKL}
\alias{nabc.calibrate.m.and.tau.yesmxpw.yesKL}
\title{Calibrate the number of simulated summary values and the equivalence region for a specified test of equivalence by minimising the Kullback-Leibler divergence from the standardized power function to the summary likelihood.}
\usage{
  nabc.calibrate.m.and.tau.yesmxpw.yesKL(test_name,
    KL_args, max.it = 100, debug = 0, plot = FALSE,
    plot_debug = FALSE)
}
\arguments{
  \item{KL_divergence}{character, name of the function to
  compute the KL divergence for the test of equivalence.
  This function must have at least two arguments:
  \code{n.of.y}, the number of simulated summaries, and
  \code{plot}, which flag some plotting options.}

  \item{KL_args}{list of argument to be passed to
  \code{KL_divergence}.}

  \item{max.it}{this algorithm stops prematurely when the
  number of iterations to calibrate the number of simulated
  data points exceeds 'max.it'}

  \item{debug}{flag if C implementation is used.}

  \item{plot}{if \code{TRUE}, plot the result of
  minimization. Only if implemented in
  \code{KL_divergence}.}
}
\value{
  vector of variable length. The first value is
  \code{n.of.y}: the adjusted number of simulated
  summaries. The rest is returned from
  \code{KL_divergence}.
}
\description{
  Calibrate the number of simulated summary values and the
  equivalence region for a specified test of equivalence by
  minimising the Kullback-Leibler divergence from the
  standardized power function to the summary likelihood.
}
\examples{
xn <- 60
yn <- xn
xmean <- 1
xsigma <- 1
ymean <- 1
ysigma <- 2

ymean <- xmean <- 0
ysigma <- 1

obs <- rnorm(xn, xmean, xsigma)
obs <- (obs - mean(obs))/sd(obs) * xsigma + xmean
sim <- rnorm(yn, ymean, ysigma)

n.of.x <- xn
s.of.x <- sd(obs)
n.of.y <- yn
s.of.y <- sd(sim)
mx.pw <- 0.9
alpha <- 0.01
tau.u.ub <- 2

## example 1: test of location equivalence for normally distributed variables (muTOST)
#compute the Kullback-Leibler divergence between the summary likelihood and the standardized power; and plot. 
nabc.mutost.kl(n.of.x, s.of.x, n.of.y, s.of.y, mx.pw, alpha, calibrate.tau.u = T, tau.u = tau.u.ub, plot = T)

#adjust n.of.y to minimize the Kullback-Leibler divergence, and plot result.
nabc.calibrate.m.and.tau.yesmxpw.yesKL("nabc.mutost.kl", args = list(n.of.x = n.of.x, s.of.x = s.of.x, n.of.y = n.of.y, s.of.y = s.of.y, 
	mx.pw = mx.pw, alpha = alpha, calibrate.tau.u = T, tau.u = tau.u.ub), plot = T)

## example 2: test of dispersion equivalence for normally distributed variables (chisqstretch)
#compute the Kullback-Leibler divergence between the summary likelihood and the standardized power; and plot. 
nabc.chisqstretch.kl(n.of.x, s.of.x, n.of.y, s.of.y, mx.pw, alpha, calibrate.tau.u = T, tau.u = tau.u.ub, plot = T)

#adjust n.of.y to minimize the Kullback-Leibler divergence, and plot result.
nabc.calibrate.m.and.tau.yesmxpw.yesKL("nabc.chisqstretch.kl", args = list(n.of.x = n.of.x, s.of.x = s.of.x, n.of.y = n.of.y, 
	s.of.y = s.of.y, mx.pw = mx.pw, alpha = alpha, calibrate.tau.u = T, tau.u = tau.u.ub), plot = T)



}

