\name{nabc.acf.equivalence.cor}
\alias{nabc.acf.equivalence.cor}
\title{Compute the autocorrelation in a time series along with some other info}
\usage{
  nabc.acf.equivalence.cor(x, leave.out = 0)
}
\arguments{
  \item{x}{time series (simply vector)}

  \item{leave.out}{thinning, how many values in the pair
  sequence (x_i,x_i-1) should be left out. Defaults to
  zero.}
}
\value{
  \item{cor}{autocorrelation in the thinned sequence}
  \item{z}{Z-transformation of the autocorrelation (this is
  atanh of "cor")} \item{n}{Number of pairs (x_i,x_i-1)
  after thinning}
}
\description{
  Compute the autocorrelation in a time series along with
  some other info
}
\examples{
nabc.acf.equivalence.cor(rnorm(100,0,1), leave.out=2)
}

