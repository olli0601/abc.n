\name{corrz.sulkl.norm}
\alias{corrz.sulkl.norm}
\title{Area under the density of the summary likelihood of the autocorrelation error for normal summary values}
\usage{
corrz.sulkl.norm(Tsd, support = c(-Inf, Inf))
}
\arguments{
  \item{Tsd}{standard deviation of the test statistic}

  \item{support}{vector of dimension 2, support of the
  truncated summary likelihood.}
}
\value{
Integral of the summary likelihood on the support
}
\description{
Area under the density of the summary likelihood of the
autocorrelation error for normal summary values
}
\seealso{
\code{\link{corrz.sulkl}}
}

