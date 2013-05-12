\name{integrand_KL_divergence_1D}
\alias{integrand_KL_divergence_1D}
\title{Integrand of the Kullback-Leibler divergence D(P||Q) for generic distributions}
\usage{
  integrand_KL_divergence_1D(x, dP, dQ, P_arg, Q_arg)
}
\arguments{
  \item{x}{value at which integrand is evaluated. Can be a
  vector.}

  \item{dP,dQ}{name of the functions that compute the
  density of P and Q}

  \item{P_arg,Q_arg}{list of arguments for \code{dP} and
  \code{dQ}}
}
\description{
  Integrand of the Kullback-Leibler divergence D(P||Q) for
  generic distributions
}

