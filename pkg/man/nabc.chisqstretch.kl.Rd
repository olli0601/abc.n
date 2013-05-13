\name{nabc.chisqstretch.kl}
\alias{nabc.chisqstretch.kl}
\title{Compute Kullback-Leibler divergence between the summary likelihood and the power function of chisqstretch}
\usage{
  nabc.chisqstretch.kl(n.of.x, s.of.x, n.of.y, s.of.y,
    mx.pw, alpha, calibrate.tau.u = F, tau.u = 1,
    for.mle = 0, pow_scale = 1.5, debug = 0, plot = F)
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
  If \code{calibrate.tau.u==TRUE}, \code{tau.u} is just a
  guess on an upper bound on the upper tolerance of the
  equivalence region.}

  \item{for.mle}{calibrate so that the mode of the power is
  at the MLE}

  \item{pow_scale}{scale for the support of the
  standardized power. The power is truncated between
  \code{[tau.l/pow_scale,tau.u*pow_scale]} and then
  standardized.}

  \item{debug}{flag if C implementation is used}

  \item{plot}{whether to plot the two distributions}
}
\value{
  vector of length 6 \item{KL_div}{the Kullback Leibler
  divergence} \item{tau.l}{lower tolerance of the
  equivalence region} \item{tau.u}{upper tolerance of the
  equivalence region} \item{c.l}{lower ABC tolerance}
  \item{c.u}{upper ABC tolerance} \item{pw.cmx}{actual
  maximum power associated with the equivalence region}
}
\description{
  Compute Kullback-Leibler divergence between the summary
  likelihood and the power function of chisqstretch
}
\note{
  Whatever the value of \code{calibrate.tau.u}, the lower
  tolerance of the equivalence region (\code{tau.l}) is
  always numerically calibrated using
  \link{nabc.chisqstretch.tau.low}.
}
\examples{
nabc.chisqstretch.kl(n.of.x=60,s.of.x=0.1,n.of.y=60,s.of.y=0.3, mx.pw=0.9,
alpha=0.01, calibrate.tau.u=T, tau.u=1, plot=T)
}

