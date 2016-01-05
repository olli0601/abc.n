
# Example 1: calculate ABC false positive rate (alpha) for given ABC tolerance
# this requires to specify c.l, c.u, tau.l, tau.u (ad-hoc ABC parameters), n.of.y, what='ALPHA'
# note: can be useful to compute the ABC false positive rate for uncalibrated ABC routines

ratetest.calibrate(n.of.y=40, c.l=0.8, c.u=1.2, tau.l=0.5, tau.u=1.5, what='ALPHA', plot=TRUE)
ratetest.calibrate(n.of.y=40, c.l=0.8, c.u=1.2, tau.l=0.8, tau.u=1.2, what='ALPHA', plot=TRUE)

# Example 2: calibrate critical region for given ABC false positive rate and equivalence region
# this requires to specify alpha (default is 0.01), tau.l, tau.u, n.of.y, what='CR'
# note: this is just an intermediate calibration and may result in unsuitable power properties

ratetest.calibrate(n.of.y=40, tau.l=1/1.4, tau.u=1.4, alpha=0.01, what='CR', plot=TRUE)


# Example 3: calibrate critical region and power of ABC accept/reject step (default method when not specifying 'what')
# this requires to specify alpha (default is 0.01), mx.pw (desired maximum power), n.of.y, tau.u.ub (for numerical optimization)
# note: this is the default calibration because it specifies all ABC parameters
# for sensible calibration parameters, and only requires minimal information on 
# the observed summary values. If a set of observed summary values is available,
# use the calibrations in Example 4.

ratetest.calibrate(n.of.y=40, mx.pw=0.9, tau.u.ub=1.4, what='MXPW', plot=TRUE)


# Example 4: calibrate critical region, power of ABC accept/reject step, and number of simulated data points
# this requires to specify alpha (default is 0.01), mx.pw (desired maximum power, default is 0.9), and n.of.x, mean.x
# and a starting value for n.of.y, which must be equal to n.of.x
# note: the advantage here is that the KL divergence is also minimised, but we also
# need to have a set of observed summary values

ratetest.calibrate(n.of.x=30, mean.x=2, n.of.y=30, what='KL', plot=TRUE)

