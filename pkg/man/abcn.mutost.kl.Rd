\name{abcn.mutost.kl}
\alias{abcn.mutost.kl}
\title{Compute Kullback-Leibler divergence between the summary likelihood and the power function of mutost}
\usage{
  abcn.mutost.kl(n.of.x, s.of.x, n.of.y, s.of.y, mx.pw,
    alpha, calibrate.tau.u = F, tau.u = 1, pow_scale = 1.5,
    debug = 0, plot = F)
}
\arguments{
  \item{n.of.x}{number of observed summary values}

  \item{s.of.x}{standard deviation of observed summary
  values}

  \item{n.of.y}{number of simulated summary values}

  \item{s.of.y}{standard deviation of simulated summary
  values}

  \item{mx.pw}{maximum power at the point of reference
  (rho.star=0) (only when \code{calibrate.tau.u==TRUE}).}

  \item{alpha}{level of the equivalence test}

  \item{calibrate.tau.u}{if \code{TRUE} the upper tolerance
  of the equivalence region (\code{tau.u}) is calibrated so
  that power at the point of reference is equal to
  \code{mx.pw}}

  \item{tau.u}{upper tolerance of the equivalence region.
  If \code{calibrate.tau.u==TRUE}, it's just a guess on an
  upper bound on the upper tolerance of the equivalence
  region.}

  \item{pow_scale}{scale for the support of the
  standardized power. The power is truncated between
  \code{pow_scale*[-tau.u,tau.u]} and then standardized.}

  \item{debug}{flag if C implementation is used}

  \item{plot}{whether to plot the two distributions}
}
\value{
  vector of length 3 \item{KL_div}{the Kullback Leibler
  divergence} \item{tau.u}{upper tolerance of the
  equivalence region} \item{pw.cmx}{actual maximum power
  associated with the equivalence region}
}
\description{
  Compute Kullback-Leibler divergence between the summary
  likelihood and the power function of mutost
}
\examples{
abcn.mutost.kl(n.of.x=60,s.of.x=0.1,n.of.y=60,s.of.y=0.3, mx.pw=0.9,
alpha=0.01, calibrate.tau.u=T, tau.u=1, plot=T)
}

