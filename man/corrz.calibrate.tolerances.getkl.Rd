\name{corrz.calibrate.tolerances.getkl}
\alias{corrz.calibrate.tolerances.getkl}
\title{KL divergence between the summary likelihood and the power function of \code{corrz}}
\usage{
corrz.calibrate.tolerances.getkl(s.of.x, s.of.T, tau.u, mx.pw = 0.9,
  alpha = 0.01, pow_scale = 2, calibrate.tau.u = T, plot = F,
  verbose = 0)
}
\arguments{
  \item{s.of.x}{standard deviation of the test statistic
  for the observed summary values; used to construct the
  summary likelihood}

  \item{s.of.T}{standard deviation of the test statistic
  for the simulted summary values; used to construct the
  power function}

  \item{calibrate.tau.u}{If \code{calibrate.tau.u==TRUE}
  the upper tolerance of the equivalence region
  (\code{tau.u}) is calibrated so that power at the point
  of reference is equal to \code{mx.pw}}

  \item{plot}{Logical. If \code{plot==TRUE}, the power of
  the calibrated test is plotted along with the summary
  likelihood.}

  \item{tau.u}{Upper tolerance of the equivalence region.}

  \item{alpha}{Level of the equivalence test.}

  \item{mx.pw}{maximum power at the point of reference
  (rho.star).}

  \item{pow_scale}{Used to set the support of the pdf
  associated to the power function. The power is truncated
  between \code{[tau.l/pow_scale,tau.u*pow_scale]} and then
  standardized.}

  \item{verbose}{Logical. If \code{verbose==TRUE}, detailed
  calibration output is printed to the console}
}
\value{
vector of length 6 \item{KL_div}{KL divergence between the
power and the summary likelihood} \item{tau.l}{lower
tolerance of the equivalence region} \item{tau.u}{upper
tolerance of the equivalence region} \item{pw.cmx}{actual
maximum power associated with the equivalence region}
\item{cl}{lower point of the critical region, i.e. lower
standard ABC tolerance} \item{cu}{upper point of the
critical region, i.e. upper standard ABC tolerance}
}
\description{
Compute the Kullback-Leibler divergence between the summary
likelihood and the power function of the \code{corrz}
equivalence test. The KL divergence is required to
calibrate the number of simulated data points of the test.
}
\seealso{
\code{\link{corrz.calibrate}}
}

