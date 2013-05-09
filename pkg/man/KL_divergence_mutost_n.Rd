\name{KL_divergence_mutost_n}
\alias{KL_divergence_mutost_n}
\title{Compute Kullback-Leibler divergence between the summary likelihood and the power function of mutost after calibration of tau.u}
\usage{
  KL_divergence_mutost_n(n.of.x, s.of.x, n.of.y = n.of.x,
    s.of.y, mx.pw, tau.u.ub, alpha, debug = 0, plot = F)
}
\arguments{
  \item{tau.u.ub}{guess on an upper bound on the upper
  tolerance of the equivalence region}

  \item{plot}{whether to plot the two distributions.}

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
KL_divergence_mutost_n(n.of.x=20,s.of.x=0.5,n.of.y=20, mx.pw=0.9, s.of.y=0.5, tau.u.ub=2, alpha=0.01, debug = 0)
}

