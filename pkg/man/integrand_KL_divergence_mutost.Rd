\name{integrand_KL_divergence_mutost}
\alias{integrand_KL_divergence_mutost}
\title{Integrand of the Kullback-Leibler divergence for mutost}
\usage{
  integrand_KL_divergence_mutost(rho, n.of.x, n.of.y,
    tau.u, s.of.y, alpha, lkl_norm = 1, pow_norm = 1)
}
\arguments{
  \item{rho}{value at which integrand is evaluated. Can be
  a vector.}

  \item{n.of.y}{number of simulated summary values}

  \item{tau.u}{upper tolerance of the equivalence region}

  \item{lkl_norm}{constant for normalization for the
  truncated summary likelihood (default to 1)}

  \item{pow_norm}{constant for normalization for the
  truncated power (default to 1)}

  \item{n.of.x}{number of observed summary values}

  \item{s.of.y}{standard deviation in the simulated summary
  values}

  \item{alpha}{level of the equivalence test}
}
\description{
  Integrand of the Kullback-Leibler divergence for mutost
}

