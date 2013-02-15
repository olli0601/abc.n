\name{nabc.get.pfam.pval}
\alias{nabc.get.pfam.pval}
\title{Test if summary values are normally distributed}
\usage{
  nabc.get.pfam.pval(x, normal.test)
}
\arguments{
  \item{x}{summary values}

  \item{norma.test}{name of function with which normality
  of the summary values is tested}
}
\value{
  p value of the test
}
\description{
  Test if summary values are normally distributed
}
\examples{
nabc.get.pfam.pval(rnorm(1e4),"shapiro.test")
}

