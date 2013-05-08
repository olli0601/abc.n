\name{KL_divergence_mutost_n}
\alias{KL_divergence_mutost_n}
\title{Compute Kullback-Leibler divergence between the summary likelihood and the power function of mutost after calibration of tau.u}
\usage{
  KL_divergence_mutost_n(n.of.x, n.of.y, mx.pw, s.of.y,
    tau.u.ub, alpha, debug = 0)
}
\arguments{
  \item{tau.u.ub}{guess on an upper bound on the upper
  tolerance of the equivalence region}

  \item{debug}{Flag if C implementation is used.}

  \item{n.of.x}{number of observed summary values}

  \item{n.of.y}{number of simulated summary values}

  \item{s.of.y}{standard deviation in the simulated summary
  values}

  \item{alpha}{level of the equivalence test}
}
\description{
  Compute Kullback-Leibler divergence between the summary
  likelihood and the power function of mutost after
  calibration of tau.u
}
\examples{
KL_divergence_mutost_n(n.of.x=20,n.of.y=30, mx.pw=0.9, s.of.y=5, tau.u.ub=6, alpha=0.01, debug = 0)
}

