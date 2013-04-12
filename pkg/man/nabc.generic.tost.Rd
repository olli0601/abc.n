\name{nabc.generic.tost}
\alias{nabc.generic.tost}
\title{Perform a generic two one sided test. This is an internal function.}
\usage{
  nabc.generic.tost(tost.args, tau.l, tau.u, alpha,
    tost.distr = "t")
}
\arguments{
  \item{tost.args}{vector of arguments for generic TOST}

  \item{tau.l}{lower tolerance of equivalence region}

  \item{tau.u}{upper tolerance of equivalence region}

  \item{alpha}{level of equivalence test}

  \item{tost.distr}{name of distribution of tost}
}
\value{
  vector of length 7
}
\description{
  Perform a generic two one sided test. This is an internal
  function.
}

