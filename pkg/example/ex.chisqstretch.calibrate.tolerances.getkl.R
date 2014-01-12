#
#	illustrate computation of the KL divergence
#
n.of.x		<- 60 
n.of.y		<- 80
x			<- rnorm(n.of.x,0,1)
df			<- n.of.y-1
#	upper tolerance is not calibrated
#	the power function plateaus around 1, which is not desirable 
abc.param	<- chisqstretch.calibrate.tolerances.getkl(length(x), sd(x), length(x), df, 3*sd(x), mx.pw=0.9, alpha=0.01, pow_scale=1.5, calibrate.tau.u=FALSE, plot=TRUE)
abc.param
#	upper tolerance is calibrated
#	the power function has a nice peak at 1
abc.param	<- chisqstretch.calibrate.tolerances.getkl(length(x), sd(x), length(x), df, 3*sd(x), mx.pw=0.9, alpha=0.01, pow_scale=1.5, calibrate.tau.u=TRUE, plot=TRUE)
abc.param
