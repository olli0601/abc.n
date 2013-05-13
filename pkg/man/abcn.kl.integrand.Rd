\name{abcn.kl.integrand}
\alias{abcn.kl.integrand}
\title{Integrand of the Kullback-Leibler divergence D(P||Q) for generic distributions}
\usage{
  abcn.kl.integrand(x, dP, dQ, P_arg, Q_arg)
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

