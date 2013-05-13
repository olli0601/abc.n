\name{nabc.adjust.tau.lowup.KLdiv}
\alias{nabc.adjust.tau.lowup.KLdiv}
\title{Calibrate the equivalence region for a specified test of equivalence by minimising the Kullback-Leibler divergence between the power function and the summary likelihood.}
\usage{
  nabc.adjust.tau.lowup.KLdiv(KL_divergence, args,
    tau.u.lb = 1, max.it = 100, debug = 0, plot = F)
}
\arguments{
  \item{KL_divergence}{character, name of the function to
  compute the KL divergence for the test of equivalence.
  This function must have at least two arguments:
  \code{tau.u}, the upper tolerance of the equivalence
  region, and \code{plot}, which flag some plotting
  options.}

  \item{tau.u.lb}{guess on an lower bound on the upper
  tolerance of the equivalence region}

  \item{args}{list of argument to be passed to
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
  \code{n.of.y}: the number of simulated summaries. The
  rest is returned from \code{KL_divergence}.
}
\description{
  Calibrate the equivalence region for a specified test of
  equivalence by minimising the Kullback-Leibler divergence
  between the power function and the summary likelihood.
}
\note{
  This procedure is not relevant for the test of dispersion
  equivalence for normally distributed variables
  (chisqstretch).
}
\examples{
xn <- 60
yn <- xn
xmean <- 2
xsigma <- 1
ymean <- 2
ysigma <- 0.5

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

## test of location equivalence for normally distributed variables (muTOST)
#compute the Kullback-Leibler divergence between the summary likelihood and the standardized power; and plot. 
tmp <- KL_divergence_mutost(n.of.x, s.of.x, n.of.y, s.of.y, mx.pw, alpha, calibrate.tau.u = T, tau.u = tau.u.ub, 
	plot = T)
print(tmp)

#adjust tau.u to minimize the Kullback-Leibler divergence, and plot result.
nabc.adjust.tau.lowup.KLdiv("KL_divergence_mutost", args = list(n.of.x = n.of.x, s.of.x = s.of.x, n.of.y = n.of.y, 
	s.of.y = s.of.y, mx.pw = mx.pw, alpha = alpha), tau.u.lb = tmp["tau.u"], plot = T)
}

