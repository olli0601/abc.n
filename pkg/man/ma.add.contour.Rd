\name{ma.add.contour}
\alias{ma.add.contour}
\title{Compute the 2d mode of a 2D density}
\usage{
ma.add.contour(x, y, xlim = NA, ylim = NA, nlevels = 5,
  width.infl = 0.25, gridsize = c(100, 100), contour.col = "black", ...)
}
\arguments{
  \item{x}{first dimension of random variable}

  \item{y}{second dimension of random variable}

  \item{xlim}{range of first dimension of random variable}

  \item{ylim}{range of second dimension of random variable}

  \item{nlevels}{number of contour levels}

  \item{width.infl}{bandwidth inflation of kernel density
  estimate}

  \item{gridsize}{gridsize of kernel density estimate}

  \item{contour.col}{color of contours}

  \item{...}{further options}
}
\description{
This function computes the 2d mode of a 2D density based on
average shifted histograms
}

