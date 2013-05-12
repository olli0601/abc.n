xn <- 60
yn <- xn
xmean <- 2
xsigma <- 2
ymean <- 2
ysigma <- 1

obs <- rnorm(xn, xmean, xsigma)
obs <- (obs - mean(obs))/sd(obs) * xsigma + xmean
sim <- rnorm(yn, ymean, ysigma)

n.of.x <- xn
s.of.x <- sd(obs)
n.of.y <- yn
s.of.y <- sd(sim)
mx.pw <- 0.9
alpha <- 0.01
tau.u.ub <- 2

#compute the Kullback-Leibler divergence between the summary likelihood and the power; and plot. 
tmp <- KL_divergence_mutost(n.of.x, s.of.x, n.of.y, s.of.y, mx.pw, alpha, calibrate.tau.u = T, tau.u = tau.u.ub, 
	plot = T)

#adjust tau.u to minimize the Kullback-Leibler divergence, and plot result.
nabc.mutost.onesample.tau.lowup.KL(n.of.x, s.of.x, n.of.y, s.of.y, mx.pw, alpha, tau.u.lb = tmp["tau.u"], 
	plot = T)
