# Calibration procedures for accurate ABC

Approximate Bayesian Computations (ABC) are a Monte Carlo technique to perform approximate parameter inference when the likelihood term cannot be easily evaluated. 
ABC proceeds by summarising the data, simulating from the model, comparing simulated summaries to observed summaries with a distance function, and accepting the simulated summaries if they do not differ from the observed summaries by more than a user-defined tolerance parameter. These steps are repeated through many Monte Carlo iterations to obtain an approximation to the true posterior density of the model parameters. The process by which precise ABC tolerances and ABC distance functions can be obtained is often referred to as ABC calibrations. These calibrations make use of decision theoretic arguments to construct the ABC accept/reject step, so that the ABC accept/reject step enjoys certain desirable properties. 

The `abc.star` R package implements ABC calibration routines for the most commonly occurring scenarios.

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

ABC calibration routines for particular testing scenarios. These can be combined like lego blocks to build your calibrated ABC algorithm.

* Calibration routines to test if the location parameter of simulated and observed summaries is similar. These start with `mutost.*`
* Calibration routines to test if the variance parameter of simulated and observed summaries is similar. These start with `vartest.*`
* Calibration routines to test if the rate parameter of simulated and observed summaries is similar. These start with `ratetest.*`
* Calibration routines for the Z-test. This is an asymptotic test that can be used for several testing problems, for example testing if autocorrelations are similar. This test starts with `ztest.*`

