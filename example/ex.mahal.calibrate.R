n.of.x	<- 60
n.of.y  <- 60
p <- 3

# Example 1: calibrate critical region and power of ABC accept/reject step (default)
# this requires to specify alpha, mx.pw (calibration parameters), 
# p (dimensionality of MVN),
# n.of.y (number of simulations), 
# tau.u.ub (for numerical optimization)
# note: this is the default calibration because it specifies all ABC parameters
# for sensible calibration parameters, and only requires minimal information on 
# the observed summary values. If a set of observed summary values is available,
# use the calibrations in Example 4.

mahaltest.calibrate(n.of.x = n.of.x, p = p, n.of.y = n.of.y, tau.u.ub = (p - 2) + 1, what = 'MXPW', mx.pw = 0.9, alpha = 0.01, plot = TRUE, verbose = FALSE)

# Example 2: calculate ABC false positive rate for given ABC tolerance
# this requires to specify c.l, c.u, tau.l, tau.u (ad-hoc ABC parameters), n.of.x, p, n.of.y (summary parameters)
# note: can be useful to compute the ABC false positive rate for uncalibrated ABC routines

#mahaltest.calibrate(n.of.x = n.of.x, p = p, n.of.y = n.of.y, c.l = 2/3, c.u = 3/2, tau.l = 1/2, tau.u = 2, what = 'ALPHA', alpha = 0.01, plot = TRUE, verbose = FALSE)

# Example 3: calibrate critical region for given ABC false positive rate and equivalence region
# this requires to specify alpha (calibration parameters), tau.l, tau.u (ad-hoc ABC parameter), n.of.x, n.of.y (summary parameters)
# note: this is just an intermediate calibration and may result in unsuitable power properties

#vartest.calibrate(n.of.x=n.of.x, n.of.y=n.of.y, tau.l=1/2, tau.u=2, what='CR', alpha=0.01, plot=TRUE, verbose=FALSE)

# Example 4: calibrate critical region, power of ABC accept/reject step, and #simulated data points
# this requires to specify alpha, mx.pw (calibration parameters), n.of.x (summary parameters),
# and p (dimensionality of MVN)
# note: the advantage here is that the KL divergence is also minimised, but we also
# need to have a set of observed summary values

n.of.y <- NA	#this is now a calibration output
mahaltest.calibrate(p = 3, what = 'KL', mx.pw = 0.9, alpha = 0.01, plot = TRUE, verbose = FALSE)
