xn <- 60
yn <- xn
xmean <- 2
xsigma <- 1
ymean <- 2
ysigma <- 0.5

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

## test of location equivalence for normally distributed variables (muTOST)
#compute the Kullback-Leibler divergence between the summary likelihood and the standardized power; and plot. 
tmp <- nabc.mutost.kl(n.of.x, s.of.x, n.of.y, s.of.y, mx.pw, alpha, calibrate.tau.u =TRUE, tau.u = tau.u.ub, 
	plot =TRUE)
print(tmp)

#adjust tau.u to minimize the Kullback-Leibler divergence, and plot result.
nabc.calibrate.tau.nomxpw.yesKL("nabc.mutost.kl", args = list(n.of.x = n.of.x, s.of.x = s.of.x, n.of.y = n.of.y, 
	s.of.y = s.of.y, mx.pw = mx.pw, alpha = alpha), tau.u.lb = tmp["tau.u"], plot =TRUE)
