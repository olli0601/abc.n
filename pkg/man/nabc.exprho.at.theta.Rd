\name{nabc.exprho.at.theta}
\alias{nabc.exprho.at.theta}
\title{Estimate summary parameter errors rho from unbiased Monte Carlo estimates rho.mc for all proposed theta including rejections}
\usage{
  nabc.exprho.at.theta(df, theta.names, rho.names,
    thin = 1)
}
\arguments{
  \item{df}{data frame with all proposed theta and
  corresponding rho.mc for each summary of interest}

  \item{theta.names}{vector of theta names (columns in df)}

  \item{rho.names}{vector of rho names (columns in df)}

  \item{thin}{thinning factor in case there are many rows
  in df}
}
\value{
  matrix containing the estimated rho (per column). The ith
  row corresponds to the ith theta in df.
}
\description{
  Estimate summary parameter errors rho from unbiased Monte
  Carlo estimates rho.mc for all proposed theta including
  rejections
}

