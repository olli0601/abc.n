\name{integrand_KL_divergence_1D_mutost}
\alias{integrand_KL_divergence_1D_mutost}
\title{Integrand of the Kullback-Leibler divergence for mutost}
\usage{
  integrand_KL_divergence_1D_mutost(rho, n.of.x, s.of.x,
    n.of.y, s.of.y, tau.u, alpha, lkl_norm = 1,
    pow_norm = 1, lkl_support, pow_support)
}
\arguments{
  \item{rho}{value at which integrand is evaluated. Can be
  a vector.}

  \item{tau.u}{upper tolerance of the equivalence region.}

  \item{pow_norm}{constant for normalization for the
  truncated power. Default=1.}
}
\description{
  Integrand of the Kullback-Leibler divergence for mutost
}

