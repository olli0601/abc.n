\name{ma.cor}
\alias{ma.cor}
\title{Compute the autocorrelation in a time series along with some other info}
\usage{
ma.cor(x, leave.out = 0, len = ceiling(length(x)/(1 + leave.out)))
}
\arguments{
  \item{x}{numerical vector of time series values}

  \item{leave.out}{thinning, how many values in the pair
  sequence (x_i,x_i-1) should be left out. Defaults to
  zero.}

  \item{len}{number of (x_i,x_i-1) pairs}
}
\value{
vector of length 3 \item{cor}{autocorrelation in the
thinned sequence} \item{z}{Z-transformation of the
autocorrelation (this is atanh of "cor")} \item{n}{Number
of pairs (x_i,x_i-1) after thinning}
}
\description{
Compute the autocorrelation in a time series along with
some other info
}
\examples{
ma.cor(rnorm(100,0,1), leave.out=2)
}

