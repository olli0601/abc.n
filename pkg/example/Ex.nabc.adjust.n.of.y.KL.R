xn <- 60
yn <- xn
xmean <- 1
xsigma <- 1
ymean <- 1
ysigma <- 2

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

## example 1: test of location equivalence for normally distributed variables (muTOST)
#compute the Kullback-Leibler divergence between the summary likelihood and the standardized power; and plot. 
abcn.mutost.kl(n.of.x, s.of.x, n.of.y, s.of.y, mx.pw, alpha, calibrate.tau.u = T, tau.u = tau.u.ub, plot = T)

#adjust n.of.y to minimize the Kullback-Leibler divergence, and plot result.
nabc.calibrate.m.and.tau.yesmxpw.yesKL("abcn.mutost.kl", args = list(n.of.x = n.of.x, s.of.x = s.of.x, n.of.y = n.of.y, s.of.y = s.of.y, 
	mx.pw = mx.pw, alpha = alpha, calibrate.tau.u = T, tau.u = tau.u.ub), plot = T)

## example 2: test of dispersion equivalence for normally distributed variables (chisqstretch)
#compute the Kullback-Leibler divergence between the summary likelihood and the standardized power; and plot. 
nabc.chisqstretch.kl(n.of.x, s.of.x, n.of.y, s.of.y, mx.pw, alpha, calibrate.tau.u = T, tau.u = tau.u.ub, plot = T)

#adjust n.of.y to minimize the Kullback-Leibler divergence, and plot result.
nabc.calibrate.m.and.tau.yesmxpw.yesKL("nabc.chisqstretch.kl", args = list(n.of.x = n.of.x, s.of.x = s.of.x, n.of.y = n.of.y, 
	s.of.y = s.of.y, mx.pw = mx.pw, alpha = alpha, calibrate.tau.u = T, tau.u = tau.u.ub), plot = T)



