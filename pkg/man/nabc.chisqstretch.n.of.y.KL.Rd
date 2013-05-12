\name{nabc.chisqstretch.n.of.y.KL}
\alias{nabc.chisqstretch.n.of.y.KL}
\title{Calibrate the number of simulated summary values and the equivalence region for the test of location equivalence by minimising the Kullback-Leibler divergence between the power function and the summary likelihood.}
\usage{
  nabc.chisqstretch.n.of.y.KL(n.of.x, s.of.x, n.of.y,
    s.of.y, mx.pw, alpha, tau.u.ub = 2, for.mle = 0,
    max.it = 100, debug = 0, plot = F)
}
\arguments{
  \item{tau.u.ub}{guess on an upper bound on the upper
  tolerance of the equivalence region}

  \item{plot}{if \code{TRUE}, the summary likelihood and
  abc approximation are plotted}

  \item{max.it}{this algorithm stops prematurely when the
  number of iterations to calibrate the number of simulated
  data points exceeds 'max.it'}

  \item{debug}{Flag if C implementation is used.}

  \item{n.of.x}{number of observed summary values}

  \item{s.of.x}{standard deviation of observed summary
  values}

  \item{n.of.y}{number of simulated summary values}

  \item{s.of.y}{standard deviation of simulated summary
  values}

  \item{mx.pw}{maximum power at the point of reference
  (rho.star=0) (only when \code{calibrate.tau.u==TRUE}).}

  \item{alpha}{level of the equivalence test}
}
\value{
  vector of length 4 \item{n.of.y}{number of simulated
  summary values} \item{KL_div}{the Kullback Leibler
  divergence} \item{tau.u}{upper tolerance of the
  equivalence region} \item{pw.cmx}{actual maximum power
  associated with the equivalence region}
}
\description{
  Calibrate the number of simulated summary values and the
  equivalence region for the test of location equivalence
  by minimising the Kullback-Leibler divergence between the
  power function and the summary likelihood.
}
\examples{
xn <- 60
yn <- xn
xmean <- 2
xsigma <- 2
ymean <- 2
ysigma <- 3

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

#compute the Kullback-Leibler divergence between the summary likelihood and the power; and plot. 
KL_divergence_mutost(n.of.x, s.of.x, n.of.y, s.of.y, mx.pw, alpha, calibrate.tau.u = T, tau.u = tau.u.ub, 
	plot = T)

#adjust n.of.y to minimize the Kullback-Leibler divergence, and plot result.
nabc.mutost.onesample.n.of.y.KL(n.of.x, s.of.x, n.of.y, s.of.y, mx.pw, alpha, plot = T)
}

