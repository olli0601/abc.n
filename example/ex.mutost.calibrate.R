# observed summary values
xn 			<- 60
xmean 		<- 1
xsigma 		<- 1
obs 		<- rnorm(xn, xmean, xsigma)
obs 		<- (obs - mean(obs))/sd(obs) * xsigma + xmean

# Example 1: Same mean but sd(sim)>sd(obs): the power function is not "too tight" for yn=xn. 
# We increase ny while calibrating tau.u and mx.pw in order to minimize KL.
yn 			<- xn
ymean 		<- 1
ysigma 		<- 1.5
tau.u.ub 	<- 0.5
sim 		<- rnorm(yn, ymean, ysigma)
ans 		<- mutost.calibrate(n.of.x=length(obs), s.of.x= sd(obs), n.of.y=length(sim), s.of.y=sd(sim), tau.u=tau.u.ub, mx.pw=0.9, alpha=0.01, debug=TRUE, plot=TRUE, verbose=0) 

# Example 2: Same mean but sd(obs)>sd(sim): the power function is "too tight" for yn=xn.
# In this case, the only way to minimize KL is to use yn<xn, which is not what we like.
# Instead, we give up on mx.pw and calibrate tau.u.
ysigma 	<- 0.5
sim 	<- rnorm(yn, ymean, ysigma)
ans 	<- mutost.calibrate(n.of.x=length(obs), s.of.x= sd(obs), n.of.y=length(sim), s.of.y=sd(sim), tau.u=tau.u.ub, mx.pw=0.9, alpha=0.01, plot = TRUE, verbose=TRUE) 
