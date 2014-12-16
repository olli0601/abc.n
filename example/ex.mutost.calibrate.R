alpha		<- 0.01
# Example 1A: calibrate critical region and power of ABC accept/reject step
yn			<- 60
ymean 		<- 1.3
ysigma 		<- 0.8
tau.u.ub 	<- 0.5
sim 		<- rnorm(yn, ymean, ysigma)
sim 		<- (sim - mean(sim))/sd(sim) * ysigma + ymean
ans 		<- mutost.calibrate(n.of.y=length(sim), s.of.y=sd(sim), 
	tau.u.ub=tau.u.ub, what='MXPW', alpha=alpha, plot=TRUE)
						
# Example 1B: calibrate critical region and power of ABC accept/reject step
#	critical region depends on ysigma
ysigma 		<- 1.3
tau.u.ub 	<- 0.5
sim 		<- rnorm(yn, ymean, ysigma)
sim 		<- (sim - mean(sim))/sd(sim) * ysigma + ymean
ans 		<- mutost.calibrate(n.of.y=length(sim), s.of.y=sd(sim), 
	tau.u.ub=tau.u.ub, what='MXPW', alpha=alpha, plot=TRUE)
						
# Example 2A: calculate ABC false positive rate for given ABC tolerance
cus			<- c(0.1, 0.2, 0.3, 0.4, 0.5)
tau.u		<- 0.6817	
ans 		<- sapply(cus, function(c.u)	mutost.calibrate(n.of.y=length(sim), s.of.y=sd(sim), c.u=c.u, tau.u=tau.u, what='ALPHA')	) 

# Example 2B: calibrate critical region for given ABC false positive rate and equivalence region
# very low power
ans 		<- mutost.calibrate(n.of.y=length(sim), s.of.y=sd(sim), tau.u=1.2,
	what='CR', alpha=alpha, plot=TRUE)

# Example 3A: calibrate critical region, power of ABC accept/reject step, and #simulated data points
# need set of observed summary values
# case sd(sim)>sd(obs): the power function is not "too tight" for yn=xn.
xn 			<- 60
xmean 		<- 1
xsigma 		<- 1
obs 		<- rnorm(xn, xmean, xsigma)
obs 		<- (obs - mean(obs))/sd(obs) * xsigma + xmean
ans 		<- mutost.calibrate(n.of.x=length(obs), s.of.x= sd(obs), n.of.y=length(sim), s.of.y=sd(sim), tau.u.ub=5,
	what='KL', alpha=alpha, plot=TRUE, debug=TRUE)
# Example 3B: calibrate critical region, power of ABC accept/reject step, and #simulated data points
# need set of observed summary values
# case sd(sim)<sd(obs): the power function is "too tight" for yn=xn.
ysigma 		<- 0.7
sim 		<- rnorm(yn, ymean, ysigma)
sim 		<- (sim - mean(sim))/sd(sim) * ysigma + ymean
ans 		<- mutost.calibrate(n.of.x=length(obs), s.of.x= sd(obs), n.of.y=length(sim), s.of.y=sd(sim), 
	tau.u.ub=3, what='KL', alpha=alpha, plot=TRUE, debug=TRUE)

