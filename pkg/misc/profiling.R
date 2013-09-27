#library(profr)

DIR_PKG<-"~/Documents/GitProjects/nABC/git_abc.n/pkg/"
DIR_PKG<- "/Users/Oliver/git/abc.n/pkg"



if(0)
{
	system.time(nabc.mutost.kl(n.of.x, s.of.x, n.of.y, s.of.y, mx.pw, alpha, calibrate.tau.u = T, tau.u = tau.u.ub,plot=T,debug=1))
	system.time(nabc.mutost.kl(n.of.x, s.of.x, n.of.y, s.of.y, mx.pw, alpha, calibrate.tau.u = T, tau.u = tau.u.ub,plot=F,debug=0))
	#test C code
	rho	<- seq(-1,1,length.out=100000)
	system.time(ans<-nabc.mutost.sulkl(rho,n.of.x,s.of.x,debug=1))
	system.time(ans<-nabc.mutost.sulkl(rho,n.of.x,s.of.x,debug=0))

	lower=-1
	upper=1
	system.time(tmp<- .Call("abc_mutost_integrate_sulkl",lower,upper,.Machine$double.eps^0.25,.Machine$double.eps^0.25,as.double(n.of.x),s.of.x,1,1))
	system.time(tmp2<-integrate(nabc.mutost.sulkl,lower,upper,n.of.x,s.of.x,log=T))
	print(tmp)
	print(tmp2)
	
	lower=-1
	upper=1
	tau.u=5
	system.time(tmp<- .Call("abc_mutost_integrate_pow",lower,upper,.Machine$double.eps^0.25,.Machine$double.eps^0.25,as.double(n.of.y-1),s.of.y/sqrt(n.of.y),tau.u,alpha,1,0))
	system.time(tmp2<-integrate(dMuTOST_pow,lower,upper,n.of.y-1,s.of.y/sqrt(n.of.y),tau.u,alpha,log=F))
	print(tmp)
	print(tmp2)
}

#test of
if(0)
{
	require(profr)
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

#main()
