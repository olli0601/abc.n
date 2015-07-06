# calibrating the F-test

# set number of variables (i.e. summary statistics)
p 		<- 3
# set number of observations
n.of.x 	<- 100
# set number of simulations
n.of.y 	<- 100
# set T2 calculated on the observed data
t2.x	<- 0.25

# Example 1: calculate ABC false positive rate for given ABC tolerance
# this requires to specify c, tau2 (ad-hoc ABC parameters), n (the number of 
# simulated data sets) and p (the number of dimensions)
# note: can be useful to compute the ABC false positive rate for uncalibrated ABC routines
ftest.calibrate(n.of.y=n.of.y, p=p, what = 'ALPHA', c = 6.5, tau = 0.21, plot = TRUE)

# Example 2: calibrate critical region for given ABC false positive rate and equivalence region
# this requires to specify alpha (calibration parameters), tau2 (ad-hoc ABC parameter), 
# n (the number of simulated data sets) and p (the number of dimensions)
# note: this is just an intermediate calibration and may result in unsuitable power properties
ftest.calibrate(n.of.y=n.of.y, p=p, what = 'CR', tau = 0.21, alpha = 0.01, plot = TRUE)

# Example 3: calibrate critical region and power of ABC accept/reject step
# this requires to specify alpha, mx.pw (calibration parameters), n (the number of 
# simulated data sets) and p (the number of dimensions)
# note: this is just an intermediate calibration and may result in unsuitable power properties
ftest.calibrate(n.of.y=n.of.y, p=p, what = 'MXPW', mx.pw = 0.9, alpha = 0.01, use.R= TRUE, plot = TRUE)
ftest.calibrate(n.of.y=n.of.y, p=p, what = 'MXPW', mx.pw = 0.9, alpha = 0.01, use.R= FALSE, plot = TRUE)

# Example 4: calibrate critical region, power of ABC accept/reject step, and #simulated data points
# this requires to specify alpha, mx.pw (calibration parameters), and n.of.x, t2.of.x (summary parameters)
# note: these is the default calibration
##	ftest.calibrate(n.of.x=n.of.x, t2.x=t2.x, p=p, what='KL', mx.pw=0.9, alpha=0.01, use.R= TRUE, plot=TRUE, verbose=FALSE)
ftest.calibrate(n.of.x=n.of.x, t2.x=t2.x, p=p, what='KL', mx.pw=0.9, alpha=0.01, use.R= TRUE, debug=TRUE, plot=TRUE, verbose=FALSE)
ftest.calibrate(n.of.x=n.of.x, t2.x=t2.x, p=p, what='KL', mx.pw=0.9, alpha=0.01, use.R= FALSE, plot=FALSE, verbose=FALSE)
