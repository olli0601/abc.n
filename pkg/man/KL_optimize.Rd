\name{KL_optimize}
\alias{KL_optimize}
\title{A wrapper to minimize \code{KL_divergence_mutost_n} over n.of.y using the function \code{optimize}}
\usage{
  KL_optimize(n.of.y, n.of.x, s.of.x, s.of.y, mx.pw,
    tau.u.ub, alpha, debug = 0, plot = F)
}
\arguments{
  \item{n.of.y}{number of simulated summary values}

  \item{n.of.x}{number of observed summary values}

  \item{s.of.x}{standard deviation of observed summary
  values}

  \item{s.of.y}{standard deviation of simulated summary
  values}

  \item{mx.pw}{maximum power at the point of reference
  (rho.star=0)}

  \item{tau.u.ub}{guess on an upper bound on the upper
  tolerance of the equivalence region}

  \item{alpha}{level of the equivalence test}

  \item{debug}{flag if C implementation is used}

  \item{plot}{whether to plot the two distributions}
}
\description{
  A wrapper to minimize \code{KL_divergence_mutost_n} over
  n.of.y using the function \code{optimize}
}

