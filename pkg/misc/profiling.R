#library(profr)

DIR_PKG<-"~/Documents/GitProjects/nABC/git_abc.n/pkg/"

library(ggplot2)
library(reshape2)
library(pscl)
library(nortest)
library(stats)
library(plyr)

xn <- 60
yn <- 20
xmean <- 1
xsigma <- 1
ymean <- 1
ysigma <- 2

ymean <- xmean <- 0
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

yn<- 60; ysigma2<- 1; alpha<- 0.01

if(0){
system.time(nabc.mutost.kl(n.of.x, s.of.x, n.of.y, s.of.y, mx.pw, alpha, calibrate.tau.u = T, tau.u = tau.u.ub,plot=T,debug=1))
system.time(nabc.mutost.kl(n.of.x, s.of.x, n.of.y, s.of.y, mx.pw, alpha, calibrate.tau.u = T, tau.u = tau.u.ub,plot=F,debug=0))





#test C code
rho<-seq(-1,1,length.out=100000)
system.time(ans<-nabc.mutost.sulkl(rho,n.of.x,s.of.x,debug=1))
system.time(ans<-nabc.mutost.sulkl(rho,n.of.x,s.of.x,debug=0))

lower=-1
upper=1
system.time(tmp<- .Call("abcMuTOST_sulkl_integrate_qng",lower,upper,.Machine$double.eps^0.25,.Machine$double.eps^0.25,as.double(n.of.x),s.of.x,1,1))
system.time(tmp2<-integrate(nabc.mutost.sulkl,lower,upper,n.of.x,s.of.x,log=T))
print(tmp)
print(tmp2)

lower=-1
upper=1
tau.u=5
system.time(tmp<- .Call("abcMuTOST_pow_integrate_qng",lower,upper,.Machine$double.eps^0.25,.Machine$double.eps^0.25,as.double(n.of.y-1),s.of.y/sqrt(n.of.y),tau.u,alpha,1,0))
system.time(tmp2<-integrate(dMuTOST_pow,lower,upper,n.of.y-1,s.of.y/sqrt(n.of.y),tau.u,alpha,log=F))
print(tmp)
print(tmp2)
}

#test of
if(0){
## example 1: test of location equivalence for normally distributed variables (muTOST)
interval=0.02
Rprof(interval=interval)

#nabc.mutost.onesample.tau.lowup.pw(0.9, yn-1, sqrt(ysigma2/yn), 2, alpha, debug=0)

#compute the Kullback-Leibler divergence between the summary likelihood and the standardized power; and plot. 
nabc.mutost.kl(n.of.x, s.of.x, n.of.y, s.of.y, mx.pw, alpha, calibrate.tau.u = T, tau.u = tau.u.ub,plot=F,debug=1)
nabc.mutost.kl(n.of.x, s.of.x, n.of.y, s.of.y, mx.pw, alpha, calibrate.tau.u = T, tau.u = tau.u.ub,plot=F,debug=0)

#adjust n.of.y to minimize the Kullback-Leibler divergence, and plot result.
#nabc.calibrate.m.and.tau.yesmxpw.yesKL("nabc.mutost.kl", args = list(n.of.x = n.of.x, s.of.x = s.of.x, n.of.y = n.of.y, s.of.y = s.of.y, mx.pw = mx.pw, alpha = alpha, calibrate.tau.u = T, tau.u = tau.u.ub), plot = F)
nabc.calibrate.m.and.tau.yesmxpw.yesKL("nabc.chisqstretch.kl", args = list(n.of.x = n.of.x, s.of.x = s.of.x, n.of.y = n.of.y, 
	s.of.y = s.of.y, mx.pw = mx.pw, alpha = alpha, calibrate.tau.u = T, tau.u = tau.u.ub), plot = F)
Rprof(NULL)
df_prof<-parse_rprof("Rprof.out", interval= interval)
head(df_prof)

#adjust n.of.y to minimize the Kullback-Leibler divergence, and plot result.
system.time(ans<-nabc.calibrate.m.and.tau.yesmxpw.yesKL("nabc.mutost.kl", args = list(n.of.x = n.of.x, s.of.x = s.of.x, n.of.y = n.of.y, s.of.y = s.of.y, 
	mx.pw = mx.pw, alpha = alpha, calibrate.tau.u = T, tau.u = tau.u.ub, debug=1), plot = F))
system.time(ans<-nabc.calibrate.m.and.tau.yesmxpw.yesKL("nabc.mutost.kl", args = list(n.of.x = n.of.x, s.of.x = s.of.x, n.of.y = n.of.y, s.of.y = s.of.y, 
	mx.pw = mx.pw, alpha = alpha, calibrate.tau.u = T, tau.u = tau.u.ub, debug=0), plot = F))

## example 2: test of dispersion equivalence for normally distributed variables (chisqstretch)
#compute the Kullback-Leibler divergence between the summary likelihood and the standardized power; and plot. 
nabc.chisqstretch.kl(n.of.x, s.of.x, n.of.y, s.of.y, mx.pw, alpha, calibrate.tau.u = T, tau.u = tau.u.ub, plot = T)

#adjust n.of.y to minimize the Kullback-Leibler divergence, and plot result.
system.time(nabc.calibrate.m.and.tau.yesmxpw.yesKL("nabc.chisqstretch.kl", args = list(n.of.x = n.of.x, s.of.x = s.of.x, n.of.y = n.of.y, 
	s.of.y = s.of.y, mx.pw = mx.pw, alpha = alpha, calibrate.tau.u = T, tau.u = tau.u.ub), plot = F))

}


test_nomxpw.yesK<-function(){


KL_args <- list(n.of.x= n.of.x, s.of.x= s.of.x, n.of.y=n.of.y, s.of.y=s.of.y, mx.pw=mx.pw, alpha=alpha, pow_scale=1.5)
tau.u.lb <- 0.01
max.it <- 100
test_name <- "mutost"

KL_args$debug=0
system.time(ans<-nabc.calibrate.tau.nomxpw.yesKL(test_name, KL_args, tau.u.lb, max.it,debug=0))
print(ans)

KL_args2<-KL_args
KL_args2$tau.u=ans["tau.u"]
KL_args2$calibrate.tau.u=F
KL_args2$debug=1
do.call("nabc.mutost.kl", KL_args2)

flush.console()
print("ref")
KL_args$debug=1
system.time(ans<-nabc.calibrate.tau.nomxpw.yesKL(test_name, KL_args, tau.u.lb, max.it,debug=1))
print(ans)

}


test_mxpw.yesK<-function(){

KL_args <- list(n.of.x= n.of.x, s.of.x= s.of.x, n.of.y=n.of.y, s.of.y=s.of.y, mx.pw=mx.pw, alpha=alpha, pow_scale=1.5)
KL_args$tau.u <- 0.01
max.it <- 100
test_name <- "mutost"

print("R")
KL_args$debug=1
system.time(ans<-nabc.calibrate.m.and.tau.yesmxpw.yesKL(test_name, KL_args, max.it,debug=1,plot_debug=F))
print(ans)
flush.console()
print("C")
KL_args$debug=0
system.time(ans<-nabc.calibrate.m.and.tau.yesmxpw.yesKL(test_name, KL_args, max.it,debug=0))
print(ans)


}


main<-function(){
	#setwd(DIR_PKG)
	#dev_mode()
	#load_all(recompile=T)
	#load_all()
	test_mxpw.yesK()
	
}

main()
