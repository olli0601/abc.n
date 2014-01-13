#
#	illustrate power function of the autocorrelation test for normal summary values
#
xa				<- 0.1
xsigma2			<- 1	
xn				<- 150	
#	NOTE set this to 5e-3; only low to pass R CMD check
tol				<- 5e-2
#	generate an MA(1) pseudo data set such that the sample autocorrelation and variance are very close to the true values
x				<- ma.get.pseudo.data(xn, 0, xa, xsigma2, tol=tol, verbose=0)
#	calibrate the power function for the subset (x1,x2), (x4,x5), (x7,x8), ...
zx				<- ma.cor(x, leave.out=2)
abc.param.a		<- corrz.calibrate(zx["n"], mx.pw=0.9, alpha=0.01, max.it=100, pow_scale=2, debug=FALSE, plot=TRUE)
abc.param.a
#	calibrate the power function for the subset (x3,x4), (x6,x7), (x9,x10), ...
#	only almost the same as for the subset above because the number of pairs is 49
zx				<- ma.cor(x[-c(1,2)], leave.out=2)
abc.param.a		<- corrz.calibrate(zx["n"], mx.pw=0.9, alpha=0.01, max.it=100, pow_scale=2, debug=FALSE, plot=TRUE)
abc.param.a
