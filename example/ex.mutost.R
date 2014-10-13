xn <- 60
xmean <- 1
xsigma <- 1

obs <- rnorm(xn, xmean, xsigma)
obs <- (obs - mean(obs))/sd(obs) * xsigma + xmean

mx.pw <- 0.9
annealing <- 1
tau.u.ub <- 0.5
alpha <- 0.01

args <- paste("mutost",annealing,tau.u.ub,alpha,sep='/')

# Example 1: Same mean but sd(sim)>sd(obs): the power function is not "too tight" for yn=xn. 
# We increase ny while calibrating tau.u and mx.pw in order to minimize KL.
yn <- xn*3
ymean <- 1
ysigma <- 1.5
sim <- rnorm(yn, ymean, ysigma)
tmp <- mutost.onesample(sim, obs, args= args, verbose= TRUE, tau.u= 1, mx.pw=mx.pw, sd.tolerance=0.05, normal.test= "sf.test", plot=1, legend.txt="test")

# Example 2: If we don't have enough simulations, use ny and calibrate tau.u and mx.pw in order to minimize KL. 
# The KL is greater than in Example 1.
sim <- sim[1:(xn+10)]
tmp <- mutost.onesample(sim, obs, args= args, verbose= TRUE, tau.u= 1, mx.pw=mx.pw, sd.tolerance=0.05, normal.test= "sf.test", plot=1, legend.txt="test")

# Example 3: Same mean but sd(obs)>sd(sim): the power function is "too tight" for yn=xn.
# In this case, the only way to minimize KL is to use yn<xn, which is not what we like.
# Instead, we give up on mx.pw and calibrate tau.u.
ymean <- 1
ysigma <- 0.5
sim <- rnorm(yn, ymean, ysigma)
tmp <- mutost.onesample(sim, obs, args= args, verbose= TRUE, tau.u= 1, mx.pw=mx.pw, sd.tolerance=0.05, normal.test= "sf.test", plot=1, legend.txt="test")
