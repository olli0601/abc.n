% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abcstar_util.R
\name{abccheck.normal}
\alias{abccheck.normal}
\title{Test if summary values are normally distributed}
\usage{
abccheck.normal(x, normal.test)
}
\arguments{
\item{x}{summary values}

\item{normal.test}{name of function with which normality of the summary values is tested}
}
\value{
p value of the test
}
\examples{
abccheck.normal(rnorm(1e4),"shapiro.test")
}

