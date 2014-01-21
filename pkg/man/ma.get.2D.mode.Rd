\name{ma.get.2D.mode}
\alias{ma.get.2D.mode}
\title{Compute the 2d mode of a 2D density}
\usage{
ma.get.2D.mode(x, y, xlim = NA, ylim = NA, xlab = "x", ylab = "y",
  n.hists = 5, nbin = 2, nlevels = 5, width.infl = 0.25,
  gridsize = c(100, 100), method = "kde", plot = 0,
  contour.col = "black", cols = head(rev(gray(seq(0, 0.95, len = trunc(50 *
  1.4)))), 50), ...)
}
\arguments{
  \item{x}{first dimension of random variable}

  \item{y}{second dimension of random variable}

  \item{xlim}{range of first dimension of random variable}

  \item{ylim}{range of second dimension of random variable}

  \item{xlab}{name of first dimension of random variable}

  \item{ylab}{name of second dimension of random variable}

  \item{n.hists}{number of average shifted histograms}

  \item{nbin}{number of bins for histogram}

  \item{nlevels}{number of contour levels}

  \item{width.infl}{bandwidth inflation of kernel density
  estimate}

  \item{gridsize}{gridsize of kernel density estimate}

  \item{method}{either 'kde' or 'ash'}

  \item{plot}{Logical}

  \item{contour.col}{color of contours}

  \item{cols}{color gradient of density estimate}

  \item{...}{further options}
}
\description{
This function computes the 2d mode of a 2D density based on
average shifted histograms
}

