\name{kl.2D}
\alias{kl.2D}
\title{Compute the KL divergence for two dimensional densities P and Q}
\usage{
kl.2D(df1, df2, nbin = 100)
}
\arguments{
  \item{df1}{data.table with two columes th1 and th2,
  random draws from P}

  \item{df2}{data.table with two columes th1 and th2,
  random draws from Q}

  \item{nbin}{number of bins for both densities}
}
\value{
KL divergence
}
\description{
Compute the KL divergence for two dimensional densities P
and Q
}

