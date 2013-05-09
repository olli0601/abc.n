\name{KL_divergence_1D_mutost}
\alias{KL_divergence_1D_mutost}
\title{Compute the Kullback-Leibler divergence D(P||Q) where P is the summary likelihood and Q is the mutost power.}
\usage{
  KL_divergence_1D_mutost(n.of.x, s.of.x, n.of.y, s.of.y,
    tau.u, alpha, lkl_norm = 1, pow_norm = 1, lkl_support,
    pow_support, lower, upper)
}
\arguments{
  \item{tau.u}{upper tolerance of the equivalence region}

  \item{lkl_norm}{constant for normalization for the
  truncated summary likelihood (default to 1)}

  \item{pow_norm}{constant for normalization for the
  truncated power (default to 1)}

  \item{lower}{lower bound for integration. Can be
  infinite.}

  \item{upper}{upper bound for integration. Can be
  infinite.}

  \item{n.of.x}{number of observed summary values}

  \item{s.of.x}{standard deviation of observed summary
  values}

  \item{n.of.y}{number of simulated summary values}

  \item{s.of.y}{standard deviation of simulated summary
  values}

  \item{alpha}{level of the equivalence test}
}
\description{
  Compute the Kullback-Leibler divergence D(P||Q) where P
  is the summary likelihood and Q is the mutost power.
}

