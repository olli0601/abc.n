\name{KL_divergence_mutost_optimize_n.of.y}
\alias{KL_divergence_mutost_optimize_n.of.y}
\title{A wrapper to minimize \code{KL_divergence_mutost} over n.of.y using the function \code{optimize}}
\usage{
  KL_divergence_mutost_optimize_n.of.y(n.of.y, n.of.x,
    s.of.x, s.of.y, mx.pw, calibrate.tau.u, tau.u, alpha,
    lkl_norm, lkl_support, debug = 0, plot = F)
}
\arguments{
  \item{n.of.y}{number of simulated summary values}

  \item{n.of.x}{number of observed summary values}

  \item{s.of.x}{standard deviation of observed summary
  values}

  \item{s.of.y}{standard deviation of simulated summary
  values}

  \item{mx.pw}{maximum power at the point of reference
  (rho.star=0) (only when \code{calibrate.tau.u==TRUE}).}

  \item{calibrate.tau.u}{if \code{TRUE} the upper tolerance
  of the equivalence region (\code{tau.u}) is calibrated so
  that power at the point of reference is equal to
  \code{mx.pw}}

  \item{tau.u}{upper tolerance of the equivalence region.
  If \code{calibrate.tau.u==TRUE}, it's just a guess on an
  upper bound on the upper tolerance of the equivalence
  region.}

  \item{alpha}{level of the equivalence test}

  \item{lkl_norm}{normalization constant for truncated
  summary likelihood. Default=1}

  \item{lkl_support}{support of truncated summary
  likelihood. Default=c(-Inf,Inf).}

  \item{debug}{flag if C implementation is used}

  \item{plot}{whether to plot the two distributions}
}
\description{
  A wrapper to minimize \code{KL_divergence_mutost} over
  n.of.y using the function \code{optimize}
}

