\name{mutost.pow}
\alias{mutost.pow}
\title{Density of the power for muTOST}
\usage{
mutost.pow(rho, df, s.of.T, tau.u, alpha, norm = 1, support = c(-Inf, Inf),
  log = FALSE)
}
\arguments{
  \item{rho}{vector of quantile}

  \item{norm}{normalization constant for the truncated
  power function.}

  \item{support}{vector of dimension 2. Support of the
  truncated power function.}

  \item{log}{logical; if \code{TRUE}, densities d are given
  as log(d).}

  \item{tau.u}{upper tolerance of the equivalence region.
  If \code{calibrate.tau.u==TRUE}, it's just a guess on an
  upper bound on the upper tolerance of the equivalence
  region.}

  \item{alpha}{level of the equivalence test}
}
\description{
Compute the density of the (possible truncated) power of
the equivalence test for population means of normal summary
values
}
\note{
The summary likelihood can be truncated to \code{support}
and then standardized with \code{norm}. For computational
efficiency, both \code{norm} and \code{support} must be
provided although each one can be derived (numerically)
from the other.
}
\examples{
require(plyr)
require(ggplot2)

rho <- seq(-2,2,0.01)
n.of.y <- 30
s.of.T <- 0.2
alpha <- 0.01

# power increases with tau.u and becomes flat
tmp <- ldply(seq(0.2,1.6,0.2),function(tau.u){
	return(data.frame(tau.u=as.factor(tau.u),rho=rho,power=mutost.pow(rho, df=n.of.y-1, s.of.T, tau.u=tau.u, alpha)))
})

p <- ggplot(tmp,aes(x=rho,y=power,colour=tau.u))
p <- p+geom_line()
print(p)


# power increases with n.of.y and becomes flat
tau.u <- 1.5
tmp <- ldply(c(seq(1,5,1),seq(20,100,20)),function(n.of.y){
	return(data.frame(n.of.y=as.factor(n.of.y),rho=rho,power=mutost.pow(rho, df=n.of.y-1, s.of.T, tau.u, alpha)))
})

p <- ggplot(tmp,aes(x=rho,y=power,colour=n.of.y))
p <- p+geom_line()
print(p)
}

