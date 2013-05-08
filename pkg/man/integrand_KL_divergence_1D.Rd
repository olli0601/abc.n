\name{integrand_KL_divergence_1D}
\alias{integrand_KL_divergence_1D}
\title{Integrand of the Kullback-Leibler divergence for generic distribution}
\usage{
  integrand_KL_divergence_1D(x, P_dist, Q_dist, P_arg,
    Q_arg)
}
\arguments{
  \item{x}{value at which integrand is evaluated. Can be a
  vector.}

  \item{P_dist}{name of the function that compute the
  density of P}

  \item{Q_dist}{name of the function that compute the
  density of Q}

  \item{P_arg}{list of arguments for \code{P_dist}}

  \item{Q_arg}{list of arguments for \code{Q_dist}}
}
\description{
  Integrand of the Kullback-Leibler divergence for generic
  distribution
}

