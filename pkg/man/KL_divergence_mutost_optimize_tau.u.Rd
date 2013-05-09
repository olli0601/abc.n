\name{KL_divergence_mutost_optimize_tau.u}
\alias{KL_divergence_mutost_optimize_tau.u}
\title{A wrapper to minimize \code{KL_divergence_mutost} over tau.u using the function \code{optimize}}
\usage{
  KL_divergence_mutost_optimize_tau.u(tau.u, n.of.x,
    s.of.x, n.of.y, s.of.y, mx.pw, calibrate.tau.u = F,
    alpha, lkl_norm, lkl_support, debug = 0, plot = F)
}
\arguments{
  \item{tau.u}{upper tolerance of the equivalence region.
  If \code{calibrate.tau.u==TRUE}, it's just a guess on an
  upper bound on the upper tolerance of the equivalence
  region.}

  \item{n.of.x}{number of observed summary values}

  \item{s.of.x}{standard deviation of observed summary
  values}

  \item{n.of.y}{number of simulated summary values}

  \item{s.of.y}{standard deviation of simulated summary
  values}

  \item{mx.pw}{maximum power at the point of reference
  (rho.star=0) (only when \code{calibrate.tau.u==TRUE}).}

  \item{calibrate.tau.u}{if \code{TRUE} the upper tolerance
  of the equivalence region (\code{tau.u}) is calibrated so
  that power at the point of reference is equal to
  \code{mx.pw}}

  \item{alpha}{level of the equivalence test}

  \item{debug}{flag if C implementation is used}

  \item{plot}{whether to plot the two distributions}
}
\description{
  A wrapper to minimize \code{KL_divergence_mutost} over
  tau.u using the function \code{optimize}
}

