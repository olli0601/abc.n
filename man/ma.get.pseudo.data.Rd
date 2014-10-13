\name{ma.get.pseudo.data}
\alias{ma.get.pseudo.data}
\title{Generate an MA(1) pseudo data set}
\usage{
ma.get.pseudo.data(n, mu, a, sd, leave.out.a = 2, leave.out.s2 = 1,
  verbose = 0, tol = 0.002, return_eps_0 = FALSE)
}
\arguments{
  \item{n}{length of time series}

  \item{mu}{mean of time series}

  \item{a}{\code{a} parameter}

  \item{sd}{\code{sigma2} parameter}

  \item{leave.out.a}{thinning of summary values for the
  autocorrelation test}

  \item{leave.out.s2}{thinning of summary values for the
  variance test}

  \item{verbose}{print verbose output to the console}

  \item{tol}{tolerance with which expected MAP is matched
  in pseudo data}

  \item{return_eps_0}{logical. If
  \code{return_eps_0==TRUE}, u0 used to create the time
  series is returned}
}
\value{
MA(1) time series
}
\description{
This function generates an MA(1) pseudo data set subject to
constraints on the sample autocorrelation, the sample
variance and the MLE of the time series. All these
statistics of the time series are within a certain
tolerance value of the corresponding population level
statistics.
}

