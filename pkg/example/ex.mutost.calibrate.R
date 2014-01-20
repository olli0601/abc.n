xn <- 60
yn <- xn
xmean <- 1
ymean <- 1
xsigma <- 1

obs <- rnorm(xn, xmean, xsigma)
obs <- (obs - mean(obs))/sd(obs) * xsigma + xmean

mx.pw <- 0.9
tau.u.ub <- 0.5
alpha <- 0.01

# Example 1: Same mean but sd(sim)>sd(obs): the power function is not "too tight" for yn=xn. 
# We increase ny while calibrating tau.u and mx.pw in order to minimize KL.
ysigma <- 1.5
sim <- rnorm(yn, ymean, ysigma)

KL_args <- list(n.of.x=length(obs), s.of.x= sd(obs), n.of.y=length(sim), s.of.y=sd(sim), mx.pw=mx.pw, alpha=alpha, tau.u=tau.u.ub, pow_scale=1.5, debug=0)

tmp <- mutost.calibrate(KL_args, max.it = 100, debug = TRUE, plot = TRUE, plot_debug = FALSE, verbose=TRUE) 

# Example 2: Same mean but sd(obs)>sd(sim): the power function is "too tight" for yn=xn.
# In this case, the only way to minimize KL is to use yn<xn, which is not what we like.
# Instead, we give up on mx.pw and calibrate tau.u.
ysigma <- 0.5
sim <- rnorm(yn, ymean, ysigma)
KL_args$s.of.y <- sd(sim)
tmp <- mutost.calibrate(KL_args, max.it = 100, debug = TRUE, plot = TRUE, plot_debug = FALSE, verbose=TRUE) 
