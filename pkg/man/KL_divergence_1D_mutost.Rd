\name{KL_divergence_1D_mutost}
\alias{KL_divergence_1D_mutost}
\title{Compute the Kullback-Leibler divergence D(P||Q) where P is the summary likelihood and Q is the mutost power.}
\usage{
  KL_divergence_1D_mutost(n.of.x, n.of.y, tau.u, s.of.y,
    alpha, lkl_norm = 1, pow_norm = 1, lower, upper)
}
\arguments{
  \item{lower}{lower bound for integration. Can be
  infinite.}

  \item{upper}{upper bound for integration. Can be
  infinite.}

  \item{n.of.x}{number of observed summary values}

  \item{n.of.y}{number of simulated summary values}

  \item{tau.u}{upper tolerance of the equivalence region}

  \item{s.of.y}{standard deviation in the simulated summary
  values}

  \item{alpha}{level of the equivalence test}

  \item{lkl_norm}{constant for normalization for the
  truncated summary likelihood (default to 1)}

  \item{pow_norm}{constant for normalization for the
  truncated power (default to 1)}
}
\description{
  Compute the Kullback-Leibler divergence D(P||Q) where P
  is the summary likelihood and Q is the mutost power.
}

