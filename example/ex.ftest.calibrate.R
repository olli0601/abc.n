# calibrating the F-test

#set number of variables (i.e. summary statistics)
p <- 3
#set number of simulations
n <- 100

# Example 1: calibrate critical region and power of ABC accept/reject step
# this requires to specify alpha, mx.pw (calibration parameters), n (the number of 
# simulated data sets) and p (the number of dimensions)

ftest.calibrate(n, p, what = 'MXPW', mx.pw = 0.9, alpha = 0.01, plot = TRUE)

# Example 2: calculate ABC false positive rate for given ABC tolerance
# this requires to specify c, tau2 (ad-hoc ABC parameters), n (the number of 
# simulated data sets) and p (the number of dimensions)
# note: can be useful to compute the ABC false positive rate for uncalibrated ABC routines

ftest.calibrate(n, p, what = 'ALPHA', c = 6.5, tau2 = 0.21, plot = TRUE)

# Example 3: calibrate critical region for given ABC false positive rate and equivalence region
# this requires to specify alpha (calibration parameters), tau2 (ad-hoc ABC parameter), 
# n (the number of simulated data sets) and p (the number of dimensions)
# note: this is just an intermediate calibration and may result in unsuitable power properties

ftest.calibrate(n, p, what = 'CR', tau2 = 0.21, alpha = 0.01, plot = TRUE)
