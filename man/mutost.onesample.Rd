\name{mutost.onesample}
\alias{mutost.onesample}
\title{Exact TOST for location equivalence.}
\usage{
mutost.onesample(sim, obs, args = NA, verbose = FALSE, tau.u = 0,
  tau.l = -tau.u, alpha = 0, mx.pw = 0.9, annealing = 1,
  sd.tolerance = 0.05, normal.test = "sf.test", plot = 0,
  legend.txt = "")
}
\arguments{
  \item{sim}{simulated summary values}

  \item{obs}{observed summary values}

  \item{args}{argument that contains the equivalence region
  and the level of the test (see Examples). This is the
  preferred method for specifying arguments and overwrites
  the dummy default values}

  \item{verbose}{flag if detailed information on the
  computations should be printed to standard out}

  \item{tau.u}{upper tolerance of the equivalence region}

  \item{tau.l}{lower tolerance of the equivalence region}

  \item{alpha}{level of the equivalence test}

  \item{mx.pw}{maximum power at the point of equality}

  \item{annealing}{inflation factor of tolerances of the
  equivalence region}

  \item{sd.tolerance}{tolerance on the difference between
  the standart deviation of the simulated data and the one
  used in power calibrations. If this difference if too
  high, the actual maximum of the power might be
  considerably different from mx.pw.}

  \item{normal.test}{name of function with which normality
  of the summary values is tested}

  \item{plot}{logical; if \code{TRUE} the summary
  likelihood and calibrated power are ploted.}

  \item{legend.txt}{character; title of the plot (only when
  \code{plot==TRUE}).}
}
\value{
vector of length 16 containing \item{lkl}{likelihood of T
(test statistic for equality) under the TOST}
\item{error}{test statistic, here p-value of TOST}
\item{pval}{rescaled p-value of TOST that is expected to
follow U(0,1) under the point null hypothesis}
\item{link.mc.obs}{mean(obs)} \item{link.mc.sim}{mean(sim)}
\item{rho.mc}{mean(sim) - mean(obs)} \item{cil}{lower ABC
tolerance, here 0} \item{cir}{upper ABC tolerance, here
alpha} \item{tl}{lower tolerance of the equivalence region}
\item{tr}{upper tolerance of the equivalence region}
\item{al}{lower ABC tolerance, here 0} \item{ar}{upper ABC
tolerance, here alpha} \item{pfam.pval}{p-value of the
normal test for the simulated summary values}
\item{nsim}{number of simulated summary values}
\item{mx.pow}{maximum power at the point of equality}
\item{rho.pow}{power at the point rho.mc}
}
\description{
Perform the exact TOST for location equivalence when the
summary values are normally distributed
}
\examples{
xn <- 60
xmean <- 1
xsigma <- 1

obs <- rnorm(xn, xmean, xsigma)
obs <- (obs - mean(obs))/sd(obs) * xsigma + xmean

mx.pw <- 0.9
annealing <- 1
tau.u.ub <- 0.5
alpha <- 0.01

args <- paste("mutost",annealing,tau.u.ub,alpha,sep='/')

# Example 1: Same mean but sd(sim)>sd(obs): the power function is not "too tight" for yn=xn. 
# We increase ny while calibrating tau.u and mx.pw in order to minimize KL.
yn <- xn*3
ymean <- 1
ysigma <- 1.5
sim <- rnorm(yn, ymean, ysigma)
tmp <- mutost.onesample(sim, obs, args= args, verbose= TRUE, tau.u= 1, mx.pw=mx.pw, sd.tolerance=0.05, normal.test= "sf.test", plot=1, legend.txt="test")

# Example 2: If we don't have enough simulations, use ny and calibrate tau.u and mx.pw in order to minimize KL. 
# The KL is greater than in Example 1.
sim <- sim[1:(xn+10)]
tmp <- mutost.onesample(sim, obs, args= args, verbose= TRUE, tau.u= 1, mx.pw=mx.pw, sd.tolerance=0.05, normal.test= "sf.test", plot=1, legend.txt="test")

# Example 3: Same mean but sd(obs)>sd(sim): the power function is "too tight" for yn=xn.
# In this case, the only way to minimize KL is to use yn<xn, which is not what we like.
# Instead, we give up on mx.pw and calibrate tau.u.
ymean <- 1
ysigma <- 0.5
sim <- rnorm(yn, ymean, ysigma)
tmp <- mutost.onesample(sim, obs, args= args, verbose= TRUE, tau.u= 1, mx.pw=mx.pw, sd.tolerance=0.05, normal.test= "sf.test", plot=1, legend.txt="test")
}

