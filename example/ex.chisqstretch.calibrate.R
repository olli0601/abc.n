#
#	illustrate calibration of tau.l, tau.u and m 
#
n.of.x		<- 60
x			<- rnorm(n.of.x,0,1)
s.of.x		<- sd(x)
abc.param	<- chisqstretch.calibrate(n.of.x, s.of.x, scale=n.of.x, n.of.y=n.of.x, mx.pw=0.9, alpha=0.01, max.it=100, debug=FALSE, plot=TRUE)
abc.param
