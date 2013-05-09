\name{KL_divergence_mutost_tau.u}
\alias{KL_divergence_mutost_tau.u}
\title{Compute Kullback-Leibler divergence between the summary likelihood and the power function of mutost after calibration of tau.u}
\usage{
  KL_divergence_mutost_tau.u(n.of.x, s.of.x, n.of.y,
    s.of.y, mx.pw, tau.u.ub, alpha, debug = 0, plot = F)
}
\arguments{
  \item{n.of.x}{number of observed summary values}

  \item{s.of.x}{standard deviation of observed summary
  values}

  \item{n.of.y}{number of simulated summary values}

  \item{s.of.y}{standard deviation of simulated summary
  values}

  \item{mx.pw}{maximum power at the point of reference
  (rho.star=0)}

  \item{tau.u.ub}{guess on an upper bound on the upper
  tolerance of the equivalence region}

  \item{alpha}{level of the equivalence test}

  \item{debug}{flag if C implementation is used}

  \item{plot}{whether to plot the two distributions}
}
\value{
  vector of length 3 \item{KL_div}{the Kullback Leibler
  divergence} \item{tau.u}{upper tolerance of the
  equivalence region} \item{pw.cmx}{actual maximum power
  associated with the equivalence region}
}
\description{
  Compute Kullback-Leibler divergence between the summary
  likelihood and the power function of mutost after
  calibration of tau.u
}
\examples{
KL_divergence_mutost_n(n.of.x=60,s.of.x=0.1,n.of.y=60,s.of.y=0.3, mx.pw=0.9, tau.u.ub=1, alpha=0.01,plot=T)
}

