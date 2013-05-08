\name{KL_divergence_1D}
\alias{KL_divergence_1D}
\title{Compute the Kullback-Leibler divergence D(P||Q) from the 1D-distribution P and Q, when the density is available.}
\usage{
  KL_divergence_1D(P_dist, Q_dist, P_arg, Q_arg, lower,
    upper)
}
\arguments{
  \item{lower}{lower bound for integration. Can be
  infinite.}

  \item{upper}{upper bound for integration. Can be
  infinite.}

  \item{P_dist}{name of the function that compute the
  density of P}

  \item{Q_dist}{name of the function that compute the
  density of Q}

  \item{P_arg}{list of arguments for \code{P_dist}}

  \item{Q_arg}{list of arguments for \code{Q_dist}}
}
\description{
  Compute the Kullback-Leibler divergence D(P||Q) from the
  1D-distribution P and Q, when the density is available.
}
\examples{
KL_divergence_1D(P_dist="dexp",Q_dist="dexp",P_arg=list(rate=0.2),Q_arg=list(rate=0.4),lower=0,upper=Inf)
#theoretical divergence = log(0.2/0.4)+(0.4-0.2)-1 = 1-log(2) = 0.307
}

