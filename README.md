# Calibration procedures for accurate ABC

The corresponding paper is now on arxiv http://arxiv.org/abs/1305.4283	

## Authors:

* Oliver Ratmann <oliver.ratmann@ic.ac.uk>
* Anton Camacho <ntncmch@gmail.com>	
* Sen Hu <ethansen.hu@gmail.com>
* Caroline Colijn <c.colijn@imperial.ac.uk>

## Contributors:

If you'd like to add calibration routines for different equivalence tests to this library,
please let us know and we are happy to give you write access to the github repository

# Installation

The easiest way to install `abc.star` is to use the `devtools` package:

```r
# install.packages("devtools")
library(devtools)
install_github("olli0601/abc.star")
```

To work with this `R` package:
fire up `R`, and type 

```r
library(help=abc.star)
```


# Content:

* Calibration routines for dispersion equivalence of normally distributed summaries start with `chisqstretch.*`
* Calibration routines for location equivalence of normally distributed summaries start with `nabc.mutost.onesample.*`
* Calibration routines for asymptotic equivalence of autocorrelation start with `ma.*`
* Various helper functions to reconstruct and display link functions and level sets

