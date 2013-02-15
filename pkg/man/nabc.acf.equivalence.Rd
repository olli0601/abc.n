\name{nabc.acf.equivalence}
\alias{nabc.acf.equivalence}
\title{Perform the asymptotic equivalence test for autocorrelations at lag 1}
\usage{
  nabc.acf.equivalence(sim, obs, args = NA,
    verbose = FALSE, alpha = 0, leave.out = 0,
    normal.test = "sf.test")
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

  \item{alpha}{level of the equivalence test}

  \item{leave.out}{thinning, how many values in the pair
  sequence (x_i,x_i-1) should be left out. Defaults to
  zero.}

  \item{normal.test}{name of function with which normality
  of the summary values is tested}
}
\value{
  vector containing \item{error}{test statistic. Here,
  instead of T we return the p-value of the TOST.}
  \item{cil}{lower ABC tolerance. Here, instead of c^- we
  return 0.} \item{cir}{upper ABC tolerance. Here, instead
  of c^+ we return 'alpha'.} \item{al}{free entry. Here set
  to c^-.} \item{ar}{free entry. Here set to c^+.}
  \item{mx.pw}{Maximum power at the point of equality}
  \item{rho.mc}{sample estimate of 'rho'}
}
\description{
  Perform the asymptotic equivalence test for
  autocorrelations at lag 1
}
\examples{
leave.out<- 2
tau.u<- 0.09
alpha<- 0.01
n<- 5e3
sigma<- 1
a<- 0.1
args<- paste("acfequiv",leave.out,tau.u,alpha,sep='/')
x<-	rnorm(n+1,0,sigma)
x<- x[-1] + x[-(n+1)]*a
y<-	rnorm(n+1,0,sigma)
y<- y[-1] + y[-(n+1)]*a
nabc.acf.equivalence(y,x,args)
}

