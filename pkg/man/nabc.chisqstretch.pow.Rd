\name{nabc.chisqstretch.pow}
\alias{nabc.chisqstretch.pow}
\title{Compute power of the exact equivalence test for dispersion}
\usage{
  nabc.chisqstretch.pow(rho, df, cl, cu)
}
\arguments{
  \item{rho}{true ratio in simulated variance / observed
  variance}

  \item{df}{degrees of freedom}

  \item{cl}{lower ABC tolerance}

  \item{cu}{upper ABC tolerance}
}
\value{
  power of the exact test. this is exact.
}
\description{
  Compute power of the exact equivalence test for
  dispersion
}
\examples{
alpha<- 0.01
	tau.up<- 1.09
	yn<- 5e3
	tau.low<- nabc.chisqstretch.tau.low(tau.up, yn-1, alpha)
	rej<- .Call("abcScaledChiSq",	c(yn-1,tau.low,tau.up,alpha,1e-10,100,0.05)	)
	rho<- seq(tau.low,tau.up,by=0.001)
	nabc.chisqstretch.pow(rho,yn-1,rej[1],rej[2])
}

