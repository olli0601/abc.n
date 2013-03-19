\name{nabc.chisqstretch}
\alias{nabc.chisqstretch}
\title{Perform the exact test for dispersion equivalence when the summary values are normally distributed}
\usage{
  nabc.chisqstretch(sim, obs.mc, args = NA,
    verbose = FALSE, tau.l = 1, tau.u = 1, guess.tau.l = 0,
    alpha = 0, normal.test = "sf.test", for.mle = 0)
}
\arguments{
  \item{sim}{simulated summary values}

  \item{obs.mc}{variance of the observed summary values}

  \item{args}{argument that contains the equivalence region
  and the level of the test (see Examples). This is the
  preferred method for specifying arguments and overwrites
  the dummy default values}

  \item{verbose}{flag if detailed information on the
  computations should be printed to standard out}

  \item{tau.l}{lower tolerance of the equivalence region}

  \item{tau.u}{upper tolerance of the equivalence region}

  \item{guess.tau.l}{guess on the lower tolerance of the
  equivalence region. Used when the tolerances are annealed
  and calibration is numerically unstable.}

  \item{alpha}{level of the equivalence test}

  \item{leave.out}{thinning, how many values in the pair
  sequence (x_i,x_i-1) should be left out. Defaults to
  zero.}

  \item{normal.test}{name of function with which normality
  of the summary values is tested}

  \item{for.mle}{calibrate so that the mode of the power is
  at the MLE}
}
\value{
  vector containing \item{error}{test statistic, here
  var(sim)/obs.mc} \item{cil}{lower ABC tolerance c^-}
  \item{cir}{upper ABC tolerance c^+} \item{mx.pw}{Maximum
  power at the point of equality}
  \item{rho.mc}{log(var(sim) / obs.mc)}
}
\description{
  Perform the exact test for dispersion equivalence when
  the summary values are normally distributed
}
\examples{
alpha<- 0.01; xn<- yn<- 60; xsigma2<- 1; tau.u<- 2.2
	tau.l<- nabc.chisqstretch.tau.low(tau.u, yn-1, alpha)
	args<- paste("chisqstretch",tau.l,tau.u,alpha,sep='/')
	x<- rnorm(xn,0,sd=sqrt(xsigma2))
	y<- rnorm(yn,0,sd=sqrt(xsigma2))
	nabc.chisqstretch(y, var(x), args=args, verbose= 0)
}

