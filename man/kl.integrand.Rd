\name{kl.integrand}
\alias{kl.integrand}
\title{Integrand of the Kullback-Leibler divergence D(P||Q) for generic distributions}
\usage{
kl.integrand(x, dP, dQ, P_arg, Q_arg)
}
\arguments{
  \item{x}{value at which integrand is evaluated. Can be a
  vector.}

  \item{dP}{name of the functions that compute the density
  of P}

  \item{dQ}{name of the functions that compute the density
  of Q}

  \item{P_arg}{list of arguments for \code{dP}}

  \item{Q_arg}{list of arguments for \code{dQ}}
}
\value{
KL divergence
}
\description{
Integrand of the Kullback-Leibler divergence D(P||Q) for
generic distributions
}

