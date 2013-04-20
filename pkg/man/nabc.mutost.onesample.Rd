\name{nabc.mutost.onesample}
\alias{nabc.mutost.onesample}
\title{Perform the exact TOST for location equivalence when the summary values are normally distributed}
\usage{
  nabc.mutost.onesample(sim, obs, obs.n = NA, obs.sd = NA,
    args = NA, verbose = FALSE, tau.u = 0, tau.l = -tau.u,
    alpha = 0, mx.pw = 0.9, annealing = 1,
    normal.test = "sf.test", plot = 0, legend.txt = "")
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

  \item{s.of.x}{standard deviation of the observed summary
  values}

  \item{tau.u}{upper tolerance of the equivalence region}

  \item{tau.l}{lower tolerance of the equivalence region}

  \item{alpha}{level of the equivalence test}

  \item{mx.pw}{maximum power at the point of equality}

  \item{annealing}{inflation factor of tolerances of the
  equivalence region}

  \item{normal.test}{name of function with which normality
  of the summary values is tested}
}
\value{
  vector containing \item{error}{test statistic, here
  p-value of TOST} \item{cil}{lower ABC tolerance, here 0}
  \item{cir}{upper ABC tolerance, here alpha}
  \item{mx.pw}{Maximum power at the point of equality}
  \item{rho.mc}{mean(sim) - obs.mean}
}
\description{
  Perform the exact TOST for location equivalence when the
  summary values are normally distributed
}
\examples{
tau.u<- 0.5; tau.l<- -tau.u; alpha<- 0.01; xn<- yn<- 60; xmu<- ymu<- 0.5; xsigma2<- ysigma2<- 2
	args<- paste("mutost",1,tau.u,alpha,sep='/')
	x<- rnorm(xn,xmu,sd=sqrt(xsigma2))
	y<- rnorm(yn,ymu,sd=sqrt(ysigma2))
	nabc.mutost.onesample(y, x, args= args, verbose= 0)
}

