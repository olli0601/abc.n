# calibrating the Z-test,
# to test the equivalence of autocorrelations of normal summary values in the MA(1) model

n2s		<- function(n){ 1/sqrt(floor(n)-3) }	#need formula to convert n.of.y into s.of.T, depends on application
s2n		<- function(s){ (1/s)^2+3 }				#need formula to convert s.of.T into n.of.y, depends on application
n.of.y	<- 1500		#number of independent pairs of the form (y_t, y_{t+1}) for simulated summary values y_t, t=1, ...

# Example 1: calibrate critical region and power of ABC accept/reject step
# this requires to specify alpha, mx.pw (calibration parameters), and n.of.y (summary parameters)
# note: this is **not** the default ABC calibration because it is for this test always possible to 
# minimize the KL divergence with respect to the summary likelihood, see Example 4

ztest.calibrate(n.of.y=n.of.y, n2s=n2s, s2n=s2n, mx.pw=0.9, alpha=0.01, what='MXPW', plot=TRUE)

# Example 2: calculate ABC false positive rate for given ABC tolerance
# this requires to specify c.u, tau.u (ad-hoc ABC parameters), n.of.y (summary parameters)
# note: can be useful to compute the ABC false positive rate for uncalibrated ABC routines

ztest.calibrate(n.of.y=n.of.y, n2s=n2s, s2n=s2n, c.u=1.5427, tau.u=0.1, what='ALPHA')

# Example 3: calibrate critical region for given ABC false positive rate and equivalence region
# this requires to specify alpha (calibration parameters), tau.u (ad-hoc ABC parameter), n.of.y  (summary parameters)
# note: this is just an intermediate calibration and may result in unsuitable power properties

ztest.calibrate(n.of.y=n.of.y, n2s=n2s, s2n=s2n, tau.u=0.1, alpha=0.01, what='CR', plot=TRUE)

# Example 4: calibrate critical region, power of ABC accept/reject step, and #simulated data points (default)
# this requires to specify alpha, mx.pw (calibration parameters), n.of.x (summary parameters)
# note: this is recommended to calibrate the ABC accept/reject step

n.of.x	<- 1500	
ztest.calibrate(n.of.x=n.of.x, n2s=n2s, s2n=s2n, mx.pw=0.9, alpha=0.01, what='KL', plot=TRUE)


