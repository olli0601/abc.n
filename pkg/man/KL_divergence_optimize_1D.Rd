\name{KL_divergence_optimize_1D}
\alias{KL_divergence_optimize_1D}
\title{A wrapper to minimize \code{KL_divergence} over the parameter \code{x_name} using the function \link{optimize}}
\usage{
  KL_divergence_optimize_1D(x_value, x_name,
    is_integer = FALSE, KL_divergence, args,
    verbose = FALSE)
}
\arguments{
  \item{x_value}{value tested.}

  \item{x_name}{name of the parameter over which
  minimization is performed.}

  \item{is_integer}{if \code{TRUE}, \code{x_value} is
  rounded to the closest integer. See \link{round}.}

  \item{KL_divergence}{character, name of the function that
  computes the KL divergence.}

  \item{args}{additional arguments to be passed to
  \code{KL_divergence}.}

  \item{verbose}{logical, if \code{TRUE}, print warnings.}
}
\value{
  the minimized Kullback-Leibler divergence (scalar).
}
\description{
  A wrapper to minimize \code{KL_divergence} over the
  parameter \code{x_name} using the function
  \link{optimize}
}

