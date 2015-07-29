#------------------------------------------------------------------------------------------------------------------------
nabc.test.chi2stretch.calibrate<- function()
{
	#require(devtools)
	#setwd(DIR_PKG)
	#dev_mode()
	#load_all(recompile=T)
	#load_all()
	library(abc.star)
	library(nortest)
	library(stats)
	library(plyr)
	
	
	verbose	<- 1
	#simulate data
	n.of.x	<- n.of.y	<- 60
	ymu		<- xmu		<- 0
	xsigma2	<- 1
	ysigma2	<- 1.2
	x		<- rnorm(n.of.x,xmu,sd=sqrt(xsigma2))
	x 		<- (x - mean(x))/sd(x) * sqrt(xsigma2) + xmu
	if(verbose)	cat(paste("n of x=",length(x),"mean of x=",mean(x),"sd of x=",sd(x)))
	y		<- rnorm(n.of.y,ymu,sd=sqrt(ysigma2))
	y 		<- (y - mean(y))/sd(y) * sqrt(ysigma2) + ymu
	if(verbose)	cat(paste("n of y=",length(y),"mean of y=",mean(y),"sd of y=",sd(y)))
	
	#set up test
	for.mle		<- 1
	s.of.x		<- sd(x)
	s.of.y		<- sd(y)	
	alpha		<- 0.01		
	tau.u		<- 2.2 		
	#tau.l		<- chisqstretch.calibrate.taulow(tau.u, df, alpha, for.mle=for.mle)
	#tau.h		<- 0.65
	calibrate.tau.u<- F
	mx.pw		<- 0.9
	pow_scale	<- 1.5
	
	scale	<- n.of.x
	#chisqstretch.calibrate.taulow(tau.u, scale, df, alpha, rho.star=1, tol= 1e-5, max.it=100)
	if(1)
	{
		df			<- 299
		tmp			<- .Call("abcScaledChiSq",	c(scale,df, 1/2.2, 2.2, alpha,1e-10,100,0.05)	)
		rho			<- seq(0.1, 5, len=1024)
		pw			<- chisqstretch.pow(rho, length(x), df, tmp[1], tmp[2])
		plot(rho,pw,col="red", type='l')
		lines(rho,pw,col="red")
	}
	if(1)
	{
		df			<- 299
		tmp			<- .Call("abcScaledChiSq",	c(scale,df, 1/2.2, 2.2, alpha,1e-10,100,0.05)	)
		rho			<- seq(0.1, 5, len=1024)
		pw			<- chisqstretch.pow(rho, length(x), df, tmp[1], tmp[2])
		plot(rho,pw,col="black", type='l')
		tmp			<- .Call("abcScaledChiSq",	c(scale,df, 0.8, 2.2, alpha,1e-10,100,0.05)	)
		pw			<- chisqstretch.pow(rho, length(x), df, tmp[1], tmp[2])
		lines(rho,pw,col="red")
		tmp			<- .Call("abcScaledChiSq",	c(scale,df, 0.99, 2.2, alpha,1e-10,100,0.05)	)
		pw			<- chisqstretch.pow(rho, length(x), df, tmp[1], tmp[2])
		lines(rho,pw,col="green")
		tmp			<- .Call("abcScaledChiSq",	c(scale,df, 1.3, 2.2, alpha,1e-10,100,0.05)	)
		pw			<- chisqstretch.pow(rho, length(x), df, tmp[1], tmp[2])
		lines(rho,pw,col="yellow")		
	}
	if(1)
	{
		df			<- 299
		tmp			<- chisqstretch.calibrate.taulow(2.2, scale, df, alpha, rho.star=1, tol= 1e-5, max.it=100, verbose=1)
		rho			<- seq(0.1, 5, len=1024)
		pw			<- chisqstretch.pow(rho, length(x), df, tmp["cl"], tmp["cu"])
		plot(rho,pw,col="black", type='l')		
	}
	if(1)
	{
		df			<- 299
		tmp			<- chisqstretch.calibrate.taulow(2, scale, df, alpha, rho.star=1, tol= 1e-5, max.it=100, verbose=1)
		rho			<- seq(0.1, 5, len=1024)
		pw			<- chisqstretch.pow(rho, length(x), df, tmp["cl"], tmp["cu"])
		plot(rho,pw,col="black", type='l')
		tmp			<- chisqstretch.calibrate.taulow(1.8, scale, df, alpha, rho.star=1, tol= 1e-5, max.it=100, verbose=1)
		pw			<- chisqstretch.pow(rho, length(x), df, tmp["cl"], tmp["cu"])
		lines(rho,pw,col="blue")
		tmp			<- chisqstretch.calibrate.taulow(1.5, scale, df, alpha, rho.star=1, tol= 1e-5, max.it=100, verbose=1)
		pw			<- chisqstretch.pow(rho, length(x), df, tmp["cl"], tmp["cu"])
		lines(rho,pw,col="red")
		tmp			<- chisqstretch.calibrate.taulow(1.3, scale, df, alpha, rho.star=1, tol= 1e-5, max.it=100, verbose=1)
		pw			<- chisqstretch.pow(rho, length(x), df, tmp["cl"], tmp["cu"])
		lines(rho,pw,col="green")		
		#tau.up should be within 1.3,1.5  resulting rej region is a bit larger than [4.6655219802 5.3111693847]
		
		tmp			<- chisqstretch.calibrate.tauup(mx.pw, 3*s.of.y, scale, df, alpha, rho.star=1, tol= 1e-5, max.it=100, verbose=1)
		rho			<- seq(0.1, 5, len=1024)
		pw			<- chisqstretch.pow(rho, length(x), df, tmp["cl"], tmp["cu"])		
		lines(rho,pw,col="black",lty=2)		
	}
	if(1)
	{
		tmp			<- chisqstretch.calibrate.tauup(mx.pw, 3*s.of.y, scale, 299, alpha, rho.star=1, tol= 1e-5, max.it=100, verbose=0)
		rho			<- seq(0.1, 5, len=1024)
		pw			<- chisqstretch.pow(rho, length(x), df, tmp["cl"], tmp["cu"])
		plot(rho,pw,col="black", type='l')
		tmp			<- chisqstretch.calibrate.tauup(mx.pw, 3*s.of.y, scale, 149, alpha, rho.star=1, tol= 1e-5, max.it=100, verbose=0)
		rho			<- seq(0.1, 5, len=1024)
		pw			<- chisqstretch.pow(rho, length(x), 149, tmp["cl"], tmp["cu"])
		lines(rho,pw,col="blue")
		tmp			<- chisqstretch.calibrate.tauup(mx.pw, 3*s.of.y, scale, 99, alpha, rho.star=1, tol= 1e-5, max.it=100, verbose=0)
		rho			<- seq(0.1, 5, len=1024)
		pw			<- chisqstretch.pow(rho, length(x), 99, tmp["cl"], tmp["cu"])
		lines(rho,pw,col="red")
		
		tmp			<- chisqstretch.calibrate(n.of.x, s.of.x, scale=n.of.x, n.of.y=n.of.x, mx.pw=0.9, alpha=0.01, max.it=100, debug=F, plot=T)	
	}
	if(1)
	{
		n.of.x		<- 100
		s.of.x		<- 1
		tmp			<- chisqstretch.calibrate(n.of.x, s.of.x, scale=n.of.x, n.of.y=n.of.x, mx.pw=0.9, alpha=0.01, max.it=100, debug=F, plot=T)
	}
	if(1)
	{
		n.of.x		<- 150
		s.of.x		<- 1
		tmp			<- chisqstretch.calibrate(n.of.x, s.of.x, scale=n.of.x, n.of.y=n.of.x, mx.pw=0.9, alpha=0.01, max.it=100, debug=F, plot=T)
	}
	if(1)
	{
		n.of.x		<- 200
		s.of.x		<- 1
		tmp			<- chisqstretch.calibrate(n.of.x, s.of.x, scale=n.of.x, n.of.y=n.of.x, mx.pw=0.9, alpha=0.01, max.it=100, debug=F, plot=T)
	}
	if(1)
	{
		n.of.x		<- 150
		s.of.x		<- sqrt(1.01)
		tmp			<- chisqstretch.calibrate(n.of.x, s.of.x, scale=n.of.x, n.of.y=n.of.x, mx.pw=0.9, alpha=0.01, max.it=100, debug=F, plot=T)
	}	
	stop()
}
#------------------------------------------------------------------------------------------------------------------------
nabc.test.chi2stretch.illustrate.chi2<- function()
{
	yn		<- 60
	yn		<- 70
	tau.up	<- 2.2
	df		<- yn-1
	scale	<- 60
	alpha	<- 0.01
	tau.low	<- chisqstretch.tau.low(tau.up, df, alpha) 		
	verbose	<- 1
	
	rej		<- .Call("abcScaledChiSq",	c(df,df,tau.low,tau.up,alpha,1e-10,100,0.05)	)
	#print(rej)
	
	cols<- c("gray60","gray80","black","black","black")
	xu	<- seq(1/1.5,3,0.001)
	xl	<- seq(0,1.5,0.001)		
	yu	<- dchisq(xu/tau.up*scale,df)
	#yu<- yu / mean(yu)
	yl	<- dchisq(xl/tau.low*scale,df)
	#yl<- yl / mean(yl)
	
	f.name<- paste(dir.name,"/nABC.Chisq_example.pdf",sep='')
	#pdf(f.name,version="1.4",width=pdf.width,height=pdf.height)
	par(mar=c(5,4,0.5,0.5))
	plot(1,1,type='n',bty='n',xlim=range(c(xu,xl)),ylim=range(c(yu,yl,-0.001)),ylab="scaled density",xlab="T")
	cl		<- seq(rej[1],1,0.001)
	cu		<- seq(1,rej[2],0.001)
	polygon( c(cl,cu,rej[2:1]), c(dchisq(cl/tau.low*scale,df),dchisq(cu/tau.up*scale,df),0,0), col=cols[1], border=NA)
	polygon( c(tau.low,tau.up,tau.up,tau.low), c(-0.0005,-0.0005,-1,-1), col=cols[2], border=NA )		
	lines(xu,yu, col=cols[3], lty=4)		
	lines(xl,yl, col=cols[4], lty=1)
	#abline(v=1,lty=2)
	
	legend<- expression( f[m-1](T/tau^'-'), f[m-1](T/tau^'+'),'['*c^'-'*", "*c^'+'*"]",'['*tau^'-'*", "*tau^'+'*"]","")						
	legend(x=1.38,y=0.01,lty=c(1,4,NA,NA,NA),fill=c(cols[3],cols[4],cols[1],cols[2],"transparent"),legend=legend, bty= 'n', border=NA)
	#dev.off()
	
	if(verbose)	cat(paste("\ncase: rho.max at 1\ntau.low is",tau.low,"\ttau.up is",tau.up,"\tcl is",rej[1],"\tcu is",rej[2]))
}
#------------------------------------------------------------------------------------------------------------------------
nabc.test.chi2stretch.montecarlo.calibrated.tau.and.increasing.m<- function()		#check MLE, yn>xn
{
	package.mkdir(DATA,"nABC.chi2stretch")
	dir.name	<- paste(DATA,"nABC.chi2stretch",sep='/')
	pdf.width	<- 4
	pdf.height	<- 5	
	resume		<- 1
	m			<- NA	
	xn			<- yn	<- 60
	df			<- yn-1
	alpha		<- 0.01		
	tau.u		<- 2.2 		
	pw.cmx		<- KL_div	<- NA
	tau.h		<- 0.65
	
	ymu			<- xmu	<- 0
	xsigma2		<- 1
	prior.u		<- 4
	prior.l		<- 0.2
	N			<- 1e6
	
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									m= return(as.numeric(substr(arg,3,nchar(arg)))),NA)	}))
		if(length(tmp)>0) m<- tmp[1]
	}
	
	
	simu.chi2stretch.fix.x.uprior.ysig2<- function(N, prior.l, prior.u, x, yn, ymu)		
	{		
		if(prior.u<1)	stop("project.nABC.StretchedChi2.fix.x.uprior.ysig2: error at 1c")
		if(prior.l>1)	stop("project.nABC.StretchedChi2.fix.x.uprior.ysig2: error at 1d")
		ans						<- vector("list",3)
		names(ans)				<- c("x","xsigma2","data")
		ans[["x"]]				<- x
		ans[["xsigma2"]]		<- (length(x)-1)*var(x)/length(x)			#MLE of sig2
		ans[["data"]]			<- sapply(1:N,function(i)
				{					
					ysigma2		<- runif(1,prior.l,prior.u)
					y			<- rnorm(yn,ymu,sd=sqrt(ysigma2))
					tmp			<- c( ysigma2, (var(y)*(yn-1))/(var(x)*(length(x)-1)), log( var(y)/var(x) ), var(y)-var(x) )									
					tmp					
				})								
		rownames(ans[["data"]])	<- c("ysigma2","T","rho.mc","sy2-sx2")
		ans
	}
	
	
	if(!is.na(m))
	{		
		f.name<- paste(dir.name,"/nABC.Chisq_mle_incrm_",N,"_",xn,"_",prior.u,"_",prior.l,"_m",m,".R",sep='')
		cat(paste("\nnABC.Chisq: compute ",f.name))
		options(show.error.messages = FALSE, warn=1)		
		readAttempt<-try(suppressWarnings(load(f.name)))						
		options(show.error.messages = TRUE)						
		if(!resume || inherits(readAttempt, "try-error"))
		{
			x		<- rnorm(xn,xmu,sd=sqrt(xsigma2))
			x 		<- (x - mean(x))/sd(x) * sqrt(xsigma2) + xmu			
			
			yn		<- m										
			ans.m	<- simu.chi2stretch.fix.x.uprior.ysig2(N, prior.l, prior.u, x, yn, ymu)
			cat(paste("\nnABC.Chisq: save ",f.name))
			save(ans.m,file=f.name)										
		}
		else		
		{
			cat(paste("\nnABC.MA: resumed ",f.name))
		}
	}
	else
	{
		f.name	<- list.files(dir.name, pattern=paste("^nABC.Chisq_mle_incrm_",sep=''), full.names = TRUE)
		tmp		<- sapply( strsplit(f.name,'_',fixed=1), function(x) tail(x,1) )
		tmp		<- as.numeric( sapply( tmp, function(x) substr(x[1], 2, nchar(x[1])-2) ) )
		tmp		<- sort(tmp, index.return=1)
		yns		<- tmp$x
		f.name	<- f.name[ tmp$ix ]
		
		cat(paste("\nfound data, n=", length(f.name)))
		ans		<- sapply(seq_along(f.name),function(j)
				{
					#j<- 46
					#j<- 3; plot<- 1
					load(f.name[j])
					if(verbose)	cat(paste("\nprocess", f.name[j]))
					
					yn												<- yns[j]
					x												<- ans.m[["x"]]
					abc.param										<- chisqstretch.calibrate.tauup(0.9, 3*sd(x), length(x), yn-1, alpha=alpha)					
					#compute (theoretical KL and) abc parameters
					abc.param										<- chisqstretch.calibrate.tolerances.getkl(length(x), sd(x), length(x), yn-1, 3*sd(x), mx.pw=0.9, alpha=0.01, pow_scale=1.5, debug = 0, calibrate.tau.u=T, plot = F)
					#
					#	clean up empirical abc posterior density
					#
					acc												<- which( ans.m[["data"]]["T",]>=abc.param["c.l"]  &  ans.m[["data"]]["T",]<=abc.param["c.u"] )					
					tmp												<- range(ans.m[["data"]]["ysigma2",acc])
					acc.h											<- project.nABC.movingavg.gethist(ans.m[["data"]]["ysigma2",acc], ans.m[["xsigma2"]], breaks= seq(tmp[1],tmp[2],len=50), width= 0.5, plot=plot, rtn.dens=1)
					acc.h$dens$y[acc.h$dens$y<1e-3]					<- 0	
					tmp												<- range(which(acc.h$dens$y!=0))
					acc.h$dens$x									<- acc.h$dens$x[seq.int(tmp[1],tmp[2])]
					acc.h$dens$y									<- acc.h$dens$y[seq.int(tmp[1],tmp[2])]
					#
					#	need summary likelihood to compute KL
					#
					lkl_norm										<- integrate(chisqstretch.sulkl, lower= min(acc.h$dens$x), upper= max(acc.h$dens$x),  n.of.x=length(x), s.of.x=sd(x), trafo=1, support=range(acc.h$dens$x), log=FALSE)$value
					if(plot)
					{
						lkl											<- chisqstretch.sulkl(acc.h$dens$x, length(x), sd(x), trafo=1, norm=lkl_norm, support= range(acc.h$dens$x), log=FALSE)
						lines(acc.h$dens$x,lkl,col="green")	
					}
					lkl_arg											<- list(n.of.x= length(x), s.of.x= sd(x), norm = lkl_norm, trafo=1, support=range(acc.h$dens$x))					
					#
					#	compute theoretical KL
					#
					pow_norm 										<- integrate(chisqstretch.pow, lower=min(acc.h$dens$x), upper=max(acc.h$dens$x), scale=length(x), df=yn-1, c.l=abc.param["c.l"], c.u=abc.param["c.u"], norm=1, trafo=(length(x)-1)/length(x)*sd(x)*sd(x), support= range(acc.h$dens$x), log=FALSE)$value
					if(plot)
					{
						pw												<- chisqstretch.pow(seq(min(acc.h$dens$x),max(acc.h$dens$x),len=1000), length(x), yn-1, abc.param["c.l"], abc.param["c.u"], norm=pow_norm, trafo=(length(x)-1)/length(x)*sd(x)*sd(x) )
						lines(seq(min(acc.h$dens$x),max(acc.h$dens$x),len=1000),pw, col="red")						
					}
					pow_arg											<- list(scale=length(x), df=yn-1, c.l=abc.param["c.l"], c.u=abc.param["c.u"], norm=pow_norm, support=range(acc.h$dens$x), trafo=(length(x)-1)/length(x)*sd(x)*sd(x))
					suppressWarnings({ #suppress numerical inaccuracy warnings if any								
								KL.div.th									<- integrate(kl.integrand, lower=min(acc.h$dens$x), upper=max(acc.h$dens$x), dP=chisqstretch.sulkl, dQ=chisqstretch.pow, P_arg=lkl_arg, Q_arg=pow_arg)
							})
					#
					#	compute empirical KL between summary lkl (uniform prior so this is posterior on prior support) and ABC posterior
					#
					acc.mc.dens										<- approxfun(x= acc.h$dens$x, y= acc.h$dens$y, method="linear", yleft=0, yright=0, rule=2 )
					tmp												<- min(acc.h$dens$y[acc.h$dens$y!=0])
					acc.h$dens$y[acc.h$dens$y==0]					<- tmp
					acc.mc.dens.log									<- approxfun(x= acc.h$dens$x, y= log(acc.h$dens$y), method="linear", yleft=-Inf, yright=-Inf, rule=2 )
					tmp												<- function(x, log=T){ if(log){ acc.mc.dens.log(x)	}else{	acc.mc.dens(x)	}	}										
					suppressWarnings({ #suppress numerical inaccuracy warnings
								KL.div.mc 									<- integrate(kl.integrand, lower=min(acc.h$dens$x), upper=max(acc.h$dens$x), dP=chisqstretch.sulkl, dQ=tmp, P_arg=lkl_arg, Q_arg=list())						
							})
					
					c(yn=yn, KL.div.th=KL.div.th$value, KL.div.mc=KL.div.mc$value)					
				})
		
		f.name<- paste(dir.name,"/nABC.Chisq_mle_incrm_",N,"_",xn,"_",prior.u,"_",prior.l,"_",prior.u,".pdf",sep='')
		print(f.name)
		pdf(f.name,version="1.4",width=4,height=5)
		par(mar=c(5,5,0.5,0.5))				
		plot(1,1,type="n",bty='n',xlim=range(ans["yn",]),ylim=range(ans[2:3,]),xlab="m",ylab="KL divergence")		
		cols		<- rev(c(my.fade.col("black",0.6),my.fade.col("black",0.2)))
		pch			<- 15:16
		cex			<- 0.6
		points(ans["yn",], ans["KL.div.th",], pch=pch[1], col=cols[1], cex=cex)
		points(ans["yn",], ans["KL.div.mc",], pch=pch[2], col=cols[2], cex=cex)		
		legend("topleft",bty='n',pch=pch,legend=c("value used for calibration","Monte Carlo estimate\nafter ABC run"),col=cols)
		dev.off()			
	}
}
#------------------------------------------------------------------------------------------------------------------------
ms.vartest.montecarlo.ABCii.plugin.MLE<- function()		#check MLE, yn>xn
{
	require(pscl)
		
	dir.name	<- paste(DATA,"nABC.vt",sep='/')
	resume		<- 1	
	xn			<- yn	<- 60	
	alpha		<- 0.01		 		
	pw.cmx		<- KL_div	<- NA	
	
	ymu			<- xmu	<- 0
	m			<- 1
	xsigma2		<- 1
	prior.u		<- 4
	prior.l		<- 0.2
	N			<- 1e6
	
	f.names		<- list.files(dir.name, pattern='m[0-9]+\\.R$')
	dfa			<- as.data.table(t(sapply(f.names, function(f.name)
			{
				fname		<- paste(dir.name,"/",f.name, sep='')
				cat(paste("\nnABC.Chisq: load ",fname))
				options(show.error.messages = FALSE, warn=1)		
				readAttempt<-try(suppressWarnings(load(fname)))						
				options(show.error.messages = TRUE)						
				
				x				<- ans[["x"]]
				df				<- as.data.table(t(ans[["data"]]))
				df[, dMLE:= dgamma(ans[["xsigma2"]], shape=xn/2, scale=2*df[, MLE]/yn )]				
				#df[, dMLE:= densigamma(df[, MLE], xn/2-1, sum(x^2)/2 )]
				set(df, NULL, 'aMLE', df[, dMLE/max(dMLE)])
				df[, U:= runif(nrow(df))]	
				acc				<- df[, which(U<=aMLE)]
				#	get KDE from ABC output			
				acc.h			<- histo2(df[acc,ysigma2], ans[["xsigma2"]], nbreaks=70, width= 0.5, plot=F, rtn.dens=1)
				acc.h$dens$y[acc.h$dens$y<1e-3]		<- 0	
				tmp									<- which(acc.h$dens$y!=0)
				acc.h$dens$x						<- acc.h$dens$x[tmp]
				acc.h$dens$y						<- acc.h$dens$y[tmp]
				#	plot against exact posterior
				if(0)
				{
					plot(acc.h$dens$x, acc.h$dens$y, type='l', ylim=c(0,2.5))
					lines(acc.h$dens$x, densigamma(acc.h$dens$x, (xn-2)/2, sum(x^2)/2), col='red')							
					abline(v=var(x)*(xn-1)/xn)			
				}
				#	compute KDE of ABC posterior on sigma2
				acc.mc.dens							<- approxfun(x= acc.h$dens$x, y= acc.h$dens$y, method="linear", yleft=0, yright=0, rule=2 )
				tmp									<- min(acc.h$dens$y[acc.h$dens$y!=0])
				acc.h$dens$y[acc.h$dens$y==0]		<- tmp
				acc.mc.dens.log						<- approxfun(x= acc.h$dens$x, y= log(acc.h$dens$y), method="linear", yleft=-Inf, yright=-Inf, rule=2 )
				tmp									<- function(x, log=T){ if(log){ acc.mc.dens.log(x)	}else{	acc.mc.dens(x)	}	}
				#	get exact posterior
				lkl.norm							<- diff(pigamma(range(acc.h$dens$x), xn/2-1, sum(x^2)/2 ))
				tmp2								<- function(x, log=T)
				{
					if(log)
					{
						ans	<- log(densigamma(x, xn/2-1, sum(x^2)/2))-log(lkl.norm)
						ans[ which(x<prior.l | x>prior.u)]	<- -Inf
					}													
					if(!log)
					{
						ans	<- densigamma(x, xn/2-1, sum(x^2)/2)/lkl.norm
						ans[ which(x<prior.l | x>prior.u)]	<- 0
					}
					ans
				} 
				#	compute empirical KL between summary lkl (uniform prior so this is posterior on prior support) and ABC posterior
				suppressWarnings({ #suppress numerical inaccuracy warnings
							KL.div.mc 				<- integrate(kl.integrand, lower=min(acc.h$dens$x), upper=max(acc.h$dens$x), dP=tmp2, dQ=tmp, P_arg=list(), Q_arg=list())						
						})
				#	compute ABC MAP - MAP 
				c('MAP.diff'=acc.h$hmode - var(x)*(xn-1)/xn, 'KL.div.mc'=KL.div.mc$value, 'acc.prob'=length(acc)/nrow(df))				
			})))
	fname		<- paste(dir.name,"/",gsub('m[0-9]+\\.R','accurate\\.R',f.names[1]), sep='')
	cat(paste('\nsave file to', fname))
	save(dfa, file=fname)
	
	fname		<- list.files(dir.name, pattern='accurate\\.R$')
	load(paste(dir.name, fname, sep='/'))
	dfa[, mean(MAP.diff)]
	dfa[, mean(KL.div.mc)]
	dfa[, mean(acc.prob)]
}
#------------------------------------------------------------------------------------------------------------------------
ms.vartest.montecarlo.precompute<- function()		#check MLE, yn>xn
{
	require(pscl)
	package.mkdir(DATA,"nABC.vt")
	dir.name	<- paste(DATA,"nABC.vt",sep='/')
	resume		<- 1	
	xn			<- yn	<- 60	
	alpha		<- 0.01		 		
	pw.cmx		<- KL_div	<- NA	
	
	ymu			<- xmu	<- 0
	xsigma2		<- 1
	prior.u		<- 4
	prior.l		<- 0.2
	N			<- 1e6
	
	simu.chi2stretch.fix.x.uprior.ysig2<- function(N, prior.l, prior.u, x, yn, ymu)		
	{		
		if(prior.u<1)	stop("project.nABC.StretchedChi2.fix.x.uprior.ysig2: error at 1c")
		if(prior.l>1)	stop("project.nABC.StretchedChi2.fix.x.uprior.ysig2: error at 1d")
		ans						<- vector("list",3)
		names(ans)				<- c("x","xsigma2","data")
		ans[["x"]]				<- x
		ans[["xsigma2"]]		<- (length(x)-1)*var(x)/length(x)			#MLE of sig2
		ans[["data"]]			<- sapply(1:N,function(i)
				{					
					ysigma2		<- runif(1,prior.l,prior.u)
					y			<- rnorm(yn,ymu,sd=sqrt(ysigma2))
					tmp			<- c( ysigma2, (var(y)*(yn-1))/(var(x)*(length(x)-1) ), var(y)*(yn-1)/yn, var(y)-var(x)  )														
					tmp					
				})								
		rownames(ans[["data"]])	<- c("ysigma2", "T", "MLE", "sy2-sx2")
		ans
	}
	
	for(m in 58:1000)
	{		
		f.name	<- paste(dir.name,"/nABC.vartest_yneqxn_",N,"_",xn,"_",prior.u,"_",prior.l,"_m",m,".R",sep='')
		cat(paste("\nnABC.Chisq: compute ",f.name))
		options(show.error.messages = FALSE, warn=1)		
		readAttempt<-try(suppressWarnings(load(f.name)))						
		options(show.error.messages = TRUE)						
		if(!resume || inherits(readAttempt, "try-error"))
		{
			x		<- rnorm(xn, xmu,sd=sqrt(xsigma2))
			x 		<- (x - mean(x))/sd(x) * sqrt(xsigma2) + xmu			
			#	not calibrated
			yn		<- length(x)
			ans		<- simu.chi2stretch.fix.x.uprior.ysig2(N, prior.l, prior.u, x, yn, ymu)
			cat(paste("\nnABC.Chisq: save ",f.name))
			save(ans,file=f.name)				
			ans		<- NULL
		}
		f.name	<- paste(dir.name,"/nABC.vartest_yncali_",N,"_",xn,"_",prior.u,"_",prior.l,"_m",m,".R",sep='')
		cat(paste("\nnABC.Chisq: compute ",f.name))
		options(show.error.messages = FALSE, warn=1)		
		readAttempt<-try(suppressWarnings(load(f.name)))						
		options(show.error.messages = TRUE)						
		if(0 && (!resume || inherits(readAttempt, "try-error")))
		{
			x		<- rnorm(xn, xmu,sd=sqrt(xsigma2))
			x 		<- (x - mean(x))/sd(x) * sqrt(xsigma2) + xmu			
			#	calibrated	
			cali	<- vartest.calibrate(n.of.x=xn, s.of.x=sd(x), df=0, what='KL', mx.pw=0.9, alpha=0.01, plot=FALSE, verbose=FALSE)
			yn		<- cali['n.of.y']
			ans		<- simu.chi2stretch.fix.x.uprior.ysig2(N, prior.l, prior.u, x, yn, ymu)
			cat(paste("\nnABC.Chisq: save ",f.name))
			save(ans,file=f.name)				
			ans		<- NULL
		}
		gc()
	}
}
#------------------------------------------------------------------------------------------------------------------------
nabc.test.chi2stretch.montecarlo.calibrated.tau.and.m<- function()		#check MLE, yn>xn
{
	package.mkdir(DATA,"nABC.chi2stretch")
	dir.name	<- paste(DATA,"nABC.chi2stretch",sep='/')
	pdf.width	<- 4
	pdf.height	<- 5	
	resume		<- 1
	m			<- 3	
	xn			<- yn	<- 60
	df			<- yn-1
	alpha		<- 0.01		
	tau.u		<- 2.2 		
	pw.cmx		<- KL_div	<- NA
	tau.h		<- 0.65
	
	ymu			<- xmu	<- 0
	xsigma2		<- 1
	prior.u		<- 4
	prior.l		<- 0.2
	N			<- 1e6
	
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									m= return(as.numeric(substr(arg,3,nchar(arg)))),NA)	}))
		if(length(tmp)>0) m<- tmp[1]
	}
	
	
	simu.chi2stretch.fix.x.uprior.ysig2<- function(N, prior.l, prior.u, x, yn, ymu)		
	{		
		if(prior.u<1)	stop("project.nABC.StretchedChi2.fix.x.uprior.ysig2: error at 1c")
		if(prior.l>1)	stop("project.nABC.StretchedChi2.fix.x.uprior.ysig2: error at 1d")
		ans						<- vector("list",3)
		names(ans)				<- c("x","xsigma2","data")
		ans[["x"]]				<- x
		ans[["xsigma2"]]		<- (length(x)-1)*var(x)/length(x)			#MLE of sig2
		ans[["data"]]			<- sapply(1:N,function(i)
				{					
					ysigma2		<- runif(1,prior.l,prior.u)
					y			<- rnorm(yn,ymu,sd=sqrt(ysigma2))
					tmp			<- c( ysigma2, var(y)/var(x), (var(y)*(yn-1))/(var(x)*(length(x)-1)/length(x)), log( var(y)/var(x) ), var(y)-var(x) )									
					tmp					
				})								
		rownames(ans[["data"]])	<- c("ysigma2","vy/vx","T","rho.mc","sy2-sx2")
		ans
	}
	
	
	if(!is.na(m))
	{		
		f.name<- paste(dir.name,"/nABC.Chisq_mle_yncalibrated_",N,"_",xn,"_",prior.u,"_",prior.l,"_m",m,".R",sep='')
		cat(paste("\nnABC.Chisq: compute ",f.name))
		options(show.error.messages = FALSE, warn=1)		
		readAttempt<-try(suppressWarnings(load(f.name)))						
		options(show.error.messages = TRUE)						
		if(!resume || inherits(readAttempt, "try-error"))
		{
			x		<- rnorm(xn,xmu,sd=sqrt(xsigma2))
			x 		<- (x - mean(x))/sd(x) * sqrt(xsigma2) + xmu			
			
			ans		<- simu.chi2stretch.fix.x.uprior.ysig2(N, prior.l, prior.u, x, yn, ymu)
			g(yn, tau.l, tau.u, c.l, c.u, pw.cmx, KL_div)	%<-% chisqstretch.calibrate(length(x), sd(x), mx.pw=0.9, alpha=alpha, max.it=100, debug=F, plot=F)
			cat(paste("\nn.of.y=",yn,"tau.l=", tau.l,"tau.u=", tau.u,"c.l=", c.l,"c.u=", c.u,"pw.cmx=", pw.cmx,"KL_div=", KL_div))
			f.name	<- paste(dir.name,"/nABC.Chisq_mle_yncalibrated_",N,"_",xn,"_",prior.u,"_",prior.l,"_m",m,".R",sep='')
			ans.ok	<- simu.chi2stretch.fix.x.uprior.ysig2(N, prior.l, prior.u, x, yn, ymu)
			cat(paste("\nnABC.Chisq: save ",f.name))
			save(ans.ok,file=f.name)				
			ans.ok	<- NULL				
			
			yn		<- 300			
			f.name	<- paste(dir.name,"/nABC.Chisq_mle_yntoolarge_",N,"_",xn,"_",prior.u,"_",prior.l,"_m",m,".R",sep='')				
			ans.too	<- simu.chi2stretch.fix.x.uprior.ysig2(N, prior.l, prior.u, x, yn, ymu)
			cat(paste("\nnABC.Chisq: save ",f.name))
			save(ans.too,file=f.name)				
			ans.too	<- NULL
			
			yn		<- length(x)
			g(tau.l, tau.u, curr.mx.pw,	error, cl, cu)		%<-% chisqstretch.calibrate.tauup(0.9, 3*sd(x), length(x), yn-1, alpha=alpha)
			f.name	<- paste(dir.name,"/nABC.Chisq_mle_yneqxn_",N,"_",xn,"_",prior.u,"_",prior.l,"_m",m,".R",sep='')				
			ans.eq	<- simu.chi2stretch.fix.x.uprior.ysig2(N, prior.l, prior.u, x, yn, ymu)
			cat(paste("\nnABC.Chisq: save ",f.name))
			save(ans.eq,file=f.name)				
			ans.eq	<- NULL			
		}
		else		
		{
			cat(paste("\nnABC.MA: resumed ",f.name))
			if(0)		#check if calibration OK
			{			
				x												<- 	ans.ok[["x"]]
				g(yn, tau.l, tau.u, c.l, c.u, pw.cmx, KL_div)	%<-% chisqstretch.calibrate(length(x), sd(x), mx.pw=0.9, alpha=alpha, max.it=100, debug=F, plot=F)
				tstat											<- ans.ok[["data"]]["T",] / length(x)
				acc.ok											<- which( tstat>=c.l  &  tstat<=c.u )
				acc.h.ok	<- project.nABC.movingavg.gethist(ans.ok[["data"]]["ysigma2",acc.ok], ans.ok[["xsigma2"]], nbreaks= 50, width= 0.5, plot=1, ylim=c(0,2.25))
				
				sig2		<- seq(min(acc.h.ok$breaks),max(acc.h.ok$breaks),len=1000)
				su.lkl		<- chisqstretch.sulkl(sig2, length(x), sd(x), trafo=1, norm = 1, support= c(0,Inf), log=FALSE)
				lines(sig2,su.lkl,col="blue")
				abline(v=sig2[which.max(su.lkl)],col="blue",lty=2)
			}
			if(1)	#produce Fig for paper
			{
				x												<- ans.ok[["x"]]
				abc.param.ok									<- chisqstretch.calibrate(length(x), sd(x), mx.pw=0.9, alpha=alpha, max.it=100, debug=F, plot=F)
				tmp												<- abc.param.ok
				tstat											<- ans.ok[["data"]]["T",] / length(x)
				acc.ok											<- which( tstat>=tmp["cl"]  &  tstat<=tmp["cu"] )
				acc.h.ok										<- project.nABC.movingavg.gethist(ans.ok[["data"]]["ysigma2",acc.ok], ans.ok[["xsigma2"]], nbreaks= 50, width= 0.5, plot=0, ylim=c(0,2.25))
				
				#read also 
				f.name<- paste(dir.name,"/nABC.Chisq_mle_yntoolarge_",N,"_",xn,"_",prior.u,"_",prior.l,"_m",m,".R",sep='')
				cat(paste("\nnABC.Chisq: resume ",f.name))						
				readAttempt<-try(suppressWarnings(load(f.name)))
				yn												<- 300
				x												<- ans.too[["x"]]
				abc.param.toolarge								<- chisqstretch.calibrate.tauup(0.9, 3*sd(x), length(x), yn-1, alpha=alpha)
				tmp												<- abc.param.toolarge
				tstat											<- ans.too[["data"]]["T",] / length(x)
				acc.too											<- which( tstat>=tmp["cl"]  &  tstat<=tmp["cu"] )
				acc.h.too										<- project.nABC.movingavg.gethist(ans.too[["data"]]["ysigma2",acc.too], ans.too[["xsigma2"]], nbreaks= 50, width= 0.5, plot=0, ylim=c(0,3.3))
				
				#read also 
				f.name<- paste(dir.name,"/nABC.Chisq_mle_yneqxn_",N,"_",xn,"_",prior.u,"_",prior.l,"_m",m,".R",sep='')												
				cat(paste("\nnABC.Chisq: resume ",f.name))						
				readAttempt<-try(suppressWarnings(load(f.name)))	
				abc.param.yneqxn								<- chisqstretch.calibrate.tauup(0.9, 3*sd(x), xn, xn-1, alpha, rho.star=1, tol= 1e-5, max.it=100, verbose=0)
				tmp												<- abc.param.yneqxn
				tstat											<- ans.eq[["data"]]["T",] / length(x)
				acc.yneqxn										<- which( tstat>=tmp["cl"]  &  tstat<=tmp["cu"] )
				acc.h.yneqxn									<- project.nABC.movingavg.gethist(ans.eq[["data"]]["ysigma2",acc.yneqxn], ans.eq[["xsigma2"]], nbreaks= 50, width= 0.5, plot=0, ylim=c(0,3.3))
				
				#get KL of standard ABC
				tstat		<- ans.eq[["data"]]["sy2-sx2",]
				#tstat		<- ans.eq[["data"]]["rho.mc",]				
				s2			<- seq( min( ans.eq[["data"]]["ysigma2",] ), max( ans.eq[["data"]]["ysigma2",] ), length.out=1024 )
				s2.lkl		<- chisqstretch.sulkl(s2, xn, sqrt(xsigma2), trafo=1 )
				y1			<- s2.lkl / sum(s2.lkl)				
				#				
				tmp			<- sapply(c(0.8,0.4,0.2,0.1,0.05), function(cu)
						{
							acc.naive				<- which( abs(tstat) <= cu )				
							tmp						<- density( ans.eq[["data"]]["ysigma2",acc.naive]	)
							y2						<- approx(tmp$x, tmp$y, s2, yleft=0, yright=0, rule=2)$y
							y2						<- y2 / sum(y2)										
							kl						<- log( y2/y1 )
							kl[ is.nan(kl) ]		<- 0
							kl[ is.infinite(kl) ]	<- NA
							kl[ y2==0 ]				<- 0
							kl						<- kl*y2	
							c(sum(kl), length(acc.naive) / ncol(ans.eq[["data"]]) )
						})				
				
				#get standard ABC
				tstat											<- ans.eq[["data"]]["sy2-sx2",]
				acc.naive										<- which( abs(tstat) <= quantile(abs(tstat), prob= length(acc.yneqxn) / ncol(ans.eq[["data"]])) )
				acc.naive										<- which( abs(tstat) <= 0.4 )
				acc.h.naive										<- project.nABC.movingavg.gethist(ans.eq[["data"]]["ysigma2",acc.naive], ans.eq[["xsigma2"]], nbreaks= 50, width= 0.5, plot=0, ylim=c(0,3.3))
				#
				#plot sigma2
				#
				cols		<- c(my.fade.col("black",0.6),my.fade.col("black",0.2),"black","black")
				ltys		<- c(1,1,3,4)				
				
				f.name	<- paste(dir.name,"/nABC.Chisq_mle_yyn_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,"_m",m,".pdf",sep='')
				print(f.name)
				pdf(f.name,version="1.4",width=4,height=5)
				par(mar=c(5,5,0.5,0.5))
				plot(1,1,type='n',bty='n',ylab=expression("numerical estimate of "*pi[abc]*'('*sigma^2*'|'*x*')'),xlab=expression(sigma^2),xlim=c(0,3),ylim=c(0,2.5))				
				lines(acc.h.ok$mids,acc.h.ok$density,type="S",col=cols[1], lty=ltys[1], lwd=1.5)				
				lines(acc.h.naive$mids,acc.h.naive$density,type="S",col=cols[2], lty=ltys[2], lwd=1.5)
				#lines(acc.h.yneqxn$mids,acc.h.yneqxn$density,type="S",col=cols[1], lty=ltys[3])
				#lines(acc.h.too$mids,acc.h.too$density,type="S",col=cols[1], lty=ltys[1])				
				sig2		<- seq(min(acc.h.ok$breaks),max(acc.h.ok$breaks),len=1000)
				su.lkl		<- chisqstretch.sulkl(sig2, length(x), sd(x), trafo=1, norm = 1, support= c(0,Inf), log=FALSE)
				lines(sig2,su.lkl,col=cols[3],lty=ltys[3], lwd=0.75)
				abline(v=ans.ok[["xsigma2"]],col=cols[4],lty=ltys[4], lwd=0.75)											
				legend("topright",fill=c("transparent","transparent",cols[1],"transparent","transparent","transparent","transparent",cols[2],"transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent"),lty=c(NA,NA,ltys[1],NA,NA,NA,NA,ltys[2],NA,NA,NA,NA,NA,ltys[3],NA,ltys[4],NA),border=NA,bty='n',legend=expression("n=60","","calibrated","ABC*",tau^'-'*"=0.589", tau^'+'*"=1.752","m=108","standard","ABC",c^'-'*"=-0.8",c^'+'*"=0.8","m=60","",pi*'('*sigma^2*'|'*x*')',"",argmax[sigma^2],pi*'('*sigma^2*'|'*x*')'))
				#legend("topright",fill=c("transparent","transparent",cols[1],"transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent"),lty=c(NA,NA,ltys[1],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,ltys[4],NA,ltys[3],NA),border=NA,bty='n',legend=expression("n=60","","calibrated","tolerances",tau^'-'*"=0.572", tau^'+'*"=1.808","m=97","","","","","","",pi*'('*sigma^2*'|'*x*')',"",argmax[sigma^2],pi*'('*sigma^2*'|'*x*')'))
				dev.off()
			}
		}
		
	}	
	else
	{
		print("HERE")
		stop()
		#load data
		cat(paste("\nnABC.Chisq",dir.name))
		f.name<- paste(dir.name,"/nABC.Chisq_mle_yn_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,".R",sep='')
		options(show.error.messages = FALSE, warn=1)		
		readAttempt<-try(suppressWarnings(load(f.name)))						
		options(show.error.messages = TRUE)		
		resume<- 0
		if(!resume || inherits(readAttempt, "try-error"))
		{				
			f.name<- list.files(dir.name, pattern=paste("^nABC.Chisq_mle_yn_",sep=''), full.names = TRUE)
			tmp<- sort(sapply(strsplit(f.name,'_',fixed=1),function(x)	as.numeric(substr(x[length(x)],2,nchar(x[length(x)])-2))		), index.return=1)
			f.name<- f.name[tmp$ix]				
			f.name2<- list.files(dir.name, pattern=paste("^nABC.Chisq_mle_yntoolarge_",sep=''), full.names = TRUE)
			tmp2<- sort(sapply(strsplit(f.name2,'_',fixed=1),function(x)	as.numeric(substr(x[length(x)],2,nchar(x[length(x)])-2))		), index.return=1)
			f.name2<- f.name2[tmp2$ix]
			f.name<- rbind( 	f.name[tmp$x%in%intersect(tmp$x,tmp2$x)], f.name2[tmp2$x%in%intersect(tmp$x,tmp2$x)]	)
			print(f.name)
			cat(paste("\nnABC.Chisq load data: ", ncol(f.name)))
			ans<- lapply(seq_len(ncol(f.name)),function(j)
					{
						out<- matrix(NA,2,4,dimnames=list(c("ok","large"),c("mean","hmode","dmode","xsigma2")))
						
						cat(paste("\nload",f.name[1,j]))
						readAttempt<-try(suppressWarnings(load( f.name[1,j] )))
						if(inherits(readAttempt, "try-error"))	stop("error at ok")
						#tmp fix bug (now resolved)
						#ans.ok[["data"]]["error",]<- ans.ok[["data"]]["error",]*59/60
						#accept if T in boundaries					
						acc.ok<- which( ans.ok[["data"]]["error",]<=ans.ok[["cir"]]  &  ans.ok[["data"]]["error",]>=ans.ok[["cil"]] )
						acc.h.ok<- project.nABC.movingavg.gethist(ans.ok[["data"]]["ysigma2",acc.ok], ans.ok[["xsigma2"]], nbreaks= 50, width= 0.5, plot=0)
						out["ok",]<- c(acc.h.ok[["mean"]],acc.h.ok[["hmode"]],acc.h.ok[["dmode"]],ans.ok[["xsigma2"]])
						
						cat(paste("\nload",f.name[2,j]))	
						#ans.naive<- ans.ok
						readAttempt<-try(suppressWarnings(load( f.name[2,j] )))
						if(inherits(readAttempt, "try-error"))	stop("error at toolarge")	
						#tmp fix bug (now resolved)
						#ans.naive[["cil"]]<- 0.5084666; ans.naive[["cir"]]<- 1.009202							
						acc.too<- which( ans.too[["data"]]["error",]<=ans.too[["cir"]]  &  ans.too[["data"]]["error",]>=ans.too[["cil"]] )
						acc.h.too<- project.nABC.movingavg.gethist(ans.too[["data"]]["ysigma2",acc.too], ans.too[["xsigma2"]], nbreaks= 50, width= 0.5, plot=0)
						out["large",]<- c(acc.h.too[["mean"]],acc.h.too[["hmode"]],acc.h.too[["dmode"]],ans.too[["xsigma2"]])
						#print(length(acc.too) / ncol(ans.too[["data"]]))											
						if(1 && j==1)
						{								
							cols<- c(my.fade.col("black",0.6),my.fade.col("black",0.2),"black")
							ltys<- c(1,1,4,3)
							require(pscl)
							#plot sigma2
							f.name<- paste(dir.name,"/nABC.Chisq_mle_yyn_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,"_m",j,".pdf",sep='')
							print(f.name)
							pdf(f.name,version="1.4",width=4,height=5)
							par(mar=c(5,5,0.5,0.5))
							plot(acc.h.too, col=cols[2],border=NA,main='',freq=0,ylab=expression("numerical estimate of "*pi[abc]*'('*sigma^2*'|'*x*')'),xlab=expression(sigma^2),xlim=c(0,3),ylim=c(0,3))
							plot(acc.h.ok, col=cols[1],border=NA,main='',add=1,freq=0)
							#plot(acc.h.ok, col=cols[1],border=NA,main='',freq=0,ylab=expression("numerical estimate of "*pi[abc]*'('*sigma^2*'|'*x*')'),xlab=expression(sigma^2),xlim=c(0,3),ylim=c(0,3))
							a		<- (xn-2)/2	 
							b		<- ans.ok[["xsigma2"]]*(xn)/2
							var.lkl	<- b*b/((a-1)*(a-1)*(a-2))		
							
							s.of.x	<- sqrt( ans.ok[["xsigma2"]]/(xn-1)*xn )								
							tmp		<- nabc.calibrate.m.and.tau.yesmxpw.yesKL("chisqstretch.kl", args = list(	n.of.x = xn, s.of.x = s.of.x, n.of.y = xn, 
											s.of.y = NA, mx.pw = 0.9, alpha = alpha, 
											calibrate.tau.u = T, for.mle=1, tau.u = 2), plot = F)																			
							print(tmp)
							
							yn		<- round(tmp[1]*3/100)*100
							tmp		<- chisqstretch.tau.lowup(0.9, 2, yn-1, alpha, for.mle=1)
							tmp		<- chisqstretch.kl(xn, s.of.x, yn, NA, 0.9, alpha, calibrate.tau.u = F, tau.u = tmp[2], plot = F)
							print(tmp)
							
							x<- seq(prior.l,prior.u,0.001)
							y<- densigamma(x,a,b) / diff(pigamma(c(prior.l,prior.u),a,b))							
							lines(x,y,col=cols[3],lty=ltys[4])
							abline(v=ans.ok[["xsigma2"]],col=cols[3],lty=ltys[3])								
							legend("topright",fill=c("transparent","transparent",cols[1],"transparent","transparent","transparent","transparent",cols[2],"transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent"),lty=c(NA,NA,ltys[1],NA,NA,NA,NA,ltys[2],NA,NA,NA,NA,NA,ltys[4],NA,ltys[3],NA),border=NA,bty='n',legend=expression("n=60","","calibrated","tolerances",tau^'-'*"=0.572", tau^'+'*"=1.808","m=97","calibrated","tolerances",tau^'-'*"=0.726",tau^'+'*"=1.392","m=300","",pi*'('*sigma^2*'|'*x*')',"",argmax[sigma^2],pi*'('*sigma^2*'|'*x*')'))
							#legend("topright",fill=c("transparent","transparent",cols[1],"transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent"),lty=c(NA,NA,ltys[1],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,ltys[4],NA,ltys[3],NA),border=NA,bty='n',legend=expression("n=60","","calibrated","tolerances",tau^'-'*"=0.572", tau^'+'*"=1.808","m=97","","","","","","",pi*'('*sigma^2*'|'*x*')',"",argmax[sigma^2],pi*'('*sigma^2*'|'*x*')'))
							dev.off()
							stop()
						}							
						out			
					})
			
			f.name<- paste(dir.name,"/nABC.Chisq_ynmean_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,".R",sep='')
			cat(paste("\nnABC.Chisq save 'ans' to ",f.name))				
			#save(ans,file=f.name)
			print(ans)
			#stop()
		}
	}
}
#------------------------------------------------------------------------------------------------------------------------
project.nABC.StretchedChi2<- function()
{	
	package.mkdir(DATA,"nABC.StretchedChisq")
	dir.name<- paste(DATA,"nABC.StretchedChisq",sep='/')
	subprog<- 7
	pdf.width<- 4
	pdf.height<-5
	
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,5),
									subp= return(as.numeric(substr(arg,6,nchar(arg)))),NA)	}))
		if(length(tmp)>0) subprog<- tmp[1]
	}
	
	#simulate N times from same ytau.u, simulate from same xsigma
	project.nABC.StretchedChi2.fix.x.uprior.ysig2<- function(N,tau.l,tau.u,prior.l,prior.u,alpha,x,yn,ymu, with.c=TRUE, for.mle=0)		
	{		
		if(tau.u<1)	stop("project.nABC.StretchedChi2.fix.x.uprior.ysig2: error at 1a")
		if(tau.l>1)	stop("project.nABC.StretchedChi2.fix.x.uprior.ysig2: error at 1b")
		if(prior.u<1)	stop("project.nABC.StretchedChi2.fix.x.uprior.ysig2: error at 1c")
		if(prior.l>1)	stop("project.nABC.StretchedChi2.fix.x.uprior.ysig2: error at 1d")
		
		
		args<- paste("chisqstretch",for.mle,tau.l,tau.l,tau.u,alpha,sep='/')
		#print(args)
		#perform one ABC - rejection run
		ans			<- vector("list",4)
		names(ans)	<- c("xsigma2","cil","cir","data")		
		ans[["xsigma2"]]	<- ifelse( for.mle, (length(x)-1)*var(x)/length(x), var(x) )		#either MLE or unbiased estimate of sig2
		if(with.c)
		{		
			tmp				<- chisqstretch(rnorm(yn,ymu,sd=sd(x)), var(x), args=args, verbose= 0)
			ans[["cil"]]	<- tmp["cil"]
			ans[["cir"]]	<- tmp["cir"]		
		}
		else
		{
			ans[["cil"]]	<- ans[["cir"]]<- NA
		}
		print(yn)
		ans[["data"]]		<- sapply(1:N,function(i)
				{					
					ysigma2	<- runif(1,prior.l,prior.u)
					y		<- rnorm(yn,ymu,sd=sqrt(ysigma2))
					#if(with.c)
					#{
					#	tmp	<- chisqstretch(y, var(x), args=args, verbose= 0)						
					#	tmp	<- c(  ysigma2,tmp[c("error","rho.mc")],var(y)-ans[["xsigma2"]]   )							
					#}
					#else 
					#{					
					tmp	<- c(ysigma2, var(y)/var(x), log( var(y)/var(x) ), var(y)-var(x) )
					#}
					names(tmp)<- c("ysigma2","error","rho.mc","sy2-sx2")
					tmp					
				})		
		ans
	}
	project.nABC.StretchedChi2.fix.x.stdprior.ysig2<- function(N,tau.l,tau.u,prior.l,prior.u,alpha,x,yn,ymu)		
	{		
		if(tau.u<1)	stop("project.nABC.StretchedChi2.fix.x.stdprior.ysig2: error at 1a")
		if(tau.l>1)	stop("project.nABC.StretchedChi2.fix.x.stdprior.ysig2: error at 1b")
		if(prior.u<1)	stop("project.nABC.StretchedChi2.fix.x.stdprior.ysig2: error at 1c")
		if(prior.l>1)	stop("project.nABC.StretchedChi2.fix.x.stdprior.ysig2: error at 1d")
		
		
		args<- paste("chisqstretch",tau.l,tau.u,alpha,sep='/')
		#perform one ABC - rejection run
		ans<- vector("list",4)
		names(ans)<- c("xsigma2","cil","cir","data")
		tmp<- chisqstretch(rnorm(yn,ymu,sd=sd(x)), var(x), args=args, verbose= 0)
		ans[["xsigma2"]]<- var(x)
		ans[["cil"]]<- tmp["cil"]
		ans[["cir"]]<- tmp["cir"]
		ans[["data"]]<- sapply(1:N,function(i)
				{					
					ysigma2<- exp(runif(1,log(prior.l),log(prior.u)))					
					y<- rnorm(yn,ymu,sd=sqrt(ysigma2))
					tmp<- chisqstretch(y, var(x), args=args, verbose= 0)					
					tmp<- c(ysigma2,tmp[c("error","rho.mc")])	
					names(tmp)<- c("ysigma2","error","rho.mc")
					tmp					
				})		
		ans
	}
	
	if(!is.na(subprog) && subprog==6)		#compare to summary likelihood
	{
		require(pscl)
		xn		<- yn<- 60
		dfx		<- xn-1
		dfy		<- yn-1
		alpha	<- 0.01		
		prior.u	<- 3
		prior.l	<- 1/3
		xsig2	<- 1
		
		xmu<- 0
		xsigma2<- 1
		for_mle=0
		
		tmp<- chisqstretch.tau.lowup(0.9, 2, yn-1, alpha, for.mle= for_mle)
		tau.l	<- tmp[1]
		tau.u	<- tmp[2]
		cil		<- tmp[5]
		cir		<- tmp[6]
		pw2		<- chisqstretch.pow(seq(prior.l,prior.u,by=0.001),yn,yn-1,cil,cir)
		
		x		<- rnorm(xn,xmu,sd=sqrt(xsigma2))	
		a		<- (xn-2)/2	 
		b		<- var(x)*(xn-1)/2
		
		var.lkl	<- b*b/((a-1)*(a-1)*(a-2))				
		tmp		<- chisqstretch.n.of.y(xn, sqrt(var.lkl), 0.9, alpha, tau.u.ub=2, for.mle= for_mle)
		yn		<- tmp[1]
		tau.l	<- tmp[2]
		tau.u	<- tmp[3]
		cil		<- tmp[4]
		cir		<- tmp[5]
		x2		<- seq(prior.l,prior.u,0.001)
		y		<- densigamma(x2,a,b) / diff(pigamma(c(prior.l,prior.u),a,b))							
		plot(x2,y,col="red",type='l')		
		
		lines(seq(prior.l,prior.u,by=0.001),pw2/(sum(pw2)*0.001),type='l',col="green")
		
		pw		<- chisqstretch.pow(seq(prior.l,prior.u,by=0.001),yn,yn-1,cil,cir)
		lines(seq(prior.l,prior.u,by=0.001),pw/(sum(pw)*0.001),type='l',col="blue")								
		abline(v=var(x)*(xn-1)/xn )	
		abline(v=1,col="green")
		stop()
		
		
		
		
		
		
		
		
		
		
		
		xsig2.mle<- yn / (yn-1) * xsig2
		#summary likelihood of sigma2 given sample mean and sum of squares
		rho		<- seq(prior.l,prior.u,length.out=1e3)
		shape	<- (xn-2)/2	 
		scale	<- xsig2*xn*xn/(xn-1)/2
		y		<- densigamma(rho, shape, scale)
		#mean	<- sum(rho*y)/sum(y)
		#print(c( sum((rho-mean)*(rho-mean)*y)/sum(y), beta*beta/((alpha-1)*(alpha-1)*(alpha-2))))
		var.Sx	<- scale*scale/((shape-1)*(shape-1)*(shape-2))
		mo.Sx	<- scale/(shape+1)
		print(c(mo.Sx,xsig2.mle,var.Sx))
		
		#abc summary likelihood of sigma2 that peaks at rho.star and tau.u sth max.pw is 0.9
		tmp		<- chisqstretch.tau.lowup(0.9, 2.5, dfy, alpha, for.mle=1)
		tau.l	<- tmp[1]
		tau.u	<- tmp[2]
		c.l		<- tmp[5]
		c.u		<- tmp[6]
		print(tmp)
		y2		<- chisqstretch.pow(rho, yn, dfy, c.l, c.u)
		
		#abc summary likelihood of sigma2 that peaks at rho.star and tau.u sth max.pw is 0.05
		tmp		<- chisqstretch.tau.lowup(0.05, 2.5, dfy, alpha, for.mle=1)
		tau.l	<- tmp[1]
		tau.u	<- tmp[2]
		c.l		<- tmp[5]
		c.u		<- tmp[6]
		print(tmp)
		y3		<- chisqstretch.pow(rho, yn, dfy, c.l, c.u)
		
		
		#abc summary likelihood of sigma2 that peaks at rho.star and tau.u sth max.pw is 0.9 and with matching variance
		tmp		<- chisqstretch.n.of.y(xn, sqrt(var.Sx), 0.9, alpha, tau.u.ub=tau.u, for.mle=1)
		print(tmp)
		yn		<- tmp[1]
		tau.l	<- tmp[2]
		tau.u	<- tmp[3]
		c.l		<- tmp[4]
		c.u		<- tmp[5]
		y4		<- chisqstretch.pow(rho, yn, yn-1, c.l, c.u)
		#print(y3); print(sum(y3))
		
		f.name<- paste(dir.name,"/nABC.Chisq_summarylkl.pdf",sep='')	
		cat(paste("\nwrite to",f.name))
		#pdf(f.name,version="1.4",width=4,height=6)				
		par(mar=c(4,4.25,.5,.5))
		th		<- rho*xsig2.mle
		plot(rho,y/mean(y),ylim=range(c(y/mean(y),y2/mean(y2),y3/mean(y3),y4/mean(y4))),type='l',xlab=expression(sigma^2),ylab="density")
		lines(th,y2/mean(y2),col="blue")
		#lines(th,y3/mean(y3),col="blue",lty=2)
		lines(th,y4/mean(y4),col="red")
		abline(v=mo.Sx,lty=3)
		abline(v=1,lty=3,col="red")
		#x2<- x[which(y2!=0)]
		#y2<- y2[which(y2!=0)]
		#y2<- y2/sum(diff(x2)*y2[-1])
		#dev.off()
	}
	if(!is.na(subprog) && subprog==2)		#check large n
	{
		xn<- yn<- NA
		if(exists("argv"))
		{
			tmp<- na.omit(sapply(argv,function(arg)
							{	switch(substr(arg,2,3),
										yn= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
			if(length(tmp)>0) xn<- yn<- tmp[1]
		}				
		df<- yn-1
		alpha<- 0.01		
		tau.u<- 2.2 						
		mx.pw<- 0.9
		
		ymu<- xmu<- 0
		xsigma2<- 1
		prior.u<- 4
		prior.l<- 0.2
		N<- 1e6
		
		resume<- 1		
		if(!is.na(yn))	#abc-run, but do not compute cl and cu because this will depend on tau.u, which we will choose later
		{		
			f.name<- paste(dir.name,"/nABC.Chisq_largensimu_",N,"_",xn,"_",prior.u,"_",prior.l,".R",sep='')
			cat(paste("\nnABC.Chisq: compute ",f.name))
			options(show.error.messages = FALSE, warn=1)		
			readAttempt<-try(suppressWarnings(load(f.name)))						
			options(show.error.messages = TRUE)						
			if(!resume || inherits(readAttempt, "try-error"))
			{
				x<- rnorm(xn,xmu,sd=sqrt(xsigma2))	
				x<- x/sd(x)		#enforce var(x)=1				
				ans.ok<- project.nABC.StretchedChi2.fix.x.uprior.ysig2(N,0.47,tau.u,prior.l,prior.u,alpha,x,yn,ymu, with.c=0 )
				cat(paste("\nnABC.Chisq: save ",f.name))
				save(ans.ok,file=f.name)				
			}
			else
				cat(paste("\nnABC.MA: resumed ",f.name))
		}	
		else
		{
			#load data					
			cat(paste("\nnABC.Chisq",dir.name))
			f.name<- paste(dir.name,"/nABC.Chisq_largen_",N,"_",prior.u,"_",prior.l,".R",sep='')			
			options(show.error.messages = FALSE, warn=1)		
			readAttempt<-try(suppressWarnings(load(f.name)))						
			options(show.error.messages = TRUE)		
			if(!resume || inherits(readAttempt, "try-error"))
			{				
				f.name<- list.files(dir.name, pattern=paste("^nABC.Chisq_largensimu_",sep=''), full.names = TRUE)
				f.name.yn<- sort(sapply(strsplit(f.name,'_',fixed=1),function(x)	as.numeric(x[length(x)-2])		), index.return=1)
				f.name<- f.name[f.name.yn$ix]				
				cat(paste("\nnABC.Chisq load data: ", length(f.name)))
				ans<- lapply(seq_along(f.name),function(j)
						{
							out<- matrix(NA,2,9,dimnames=list(c("fx.tau.u","fx.pw"),c("tau.l","tau.u","cl","cu","acc","mean","hmode","dmode","xsigma2")))							
							cat(paste("\nload",f.name[j]))
							readAttempt<-try(suppressWarnings(load( f.name[j] )))
							if(inherits(readAttempt, "try-error"))	stop("error at largen")														
							#determine tolerances for tau.u=2.2
							yn<- f.name.yn$x[j]							
							tau.l<- chisqstretch.tau.low(tau.u, yn-1, alpha)
							rej<- .Call("abcScaledChiSq",	c(yn-1,yn-1,tau.l,tau.u,alpha,1e-10,100,0.05)	)
							ans.ok[["cil"]]<- rej[1]
							ans.ok[["cir"]]<- rej[2]
							acc.ok<- which( ans.ok[["data"]]["error",]<=ans.ok[["cir"]]  &  ans.ok[["data"]]["error",]>=ans.ok[["cil"]] )
							acc.h.ok<- project.nABC.movingavg.gethist(ans.ok[["data"]]["ysigma2",acc.ok], ans.ok[["xsigma2"]], nbreaks= 50, width= 0.5, plot=0)
							out["fx.tau.u",]<- c(tau.l, tau.u, rej[1], rej[2], length(acc.ok)/ncol(ans.ok[["data"]]), acc.h.ok[["mean"]],acc.h.ok[["hmode"]],acc.h.ok[["dmode"]],ans.ok[["xsigma2"]])
							#determine tolerances sth TOST power is 0.95
							tmp<- chisqstretch.tau.lowup(mx.pw, 2.5, yn-1, alpha)
							#print(tmp)
							rej<- .Call("abcScaledChiSq",	c(yn-1,yn-1,tmp[1],tmp[2],alpha,1e-10,100,0.05)	)
							ans.pw<- ans.ok
							ans.pw[["cil"]]<- rej[1]
							ans.pw[["cir"]]<- rej[2]
							acc.pw<- which( ans.pw[["data"]]["error",]<=ans.pw[["cir"]]  &  ans.pw[["data"]]["error",]>=ans.pw[["cil"]] )
							acc.h.pw<- project.nABC.movingavg.gethist(ans.pw[["data"]]["ysigma2",acc.pw], ans.pw[["xsigma2"]], nbreaks= 50, width= 0.5, plot=0)
							out["fx.pw",]<- c(tmp[1], tmp[2], rej[1], rej[2], length(acc.pw)/ncol(ans.pw[["data"]]), acc.h.pw[["mean"]],acc.h.pw[["hmode"]],acc.h.pw[["dmode"]],ans.pw[["xsigma2"]])
							
							if(0 && j==1)
							{
								cols<- c(my.fade.col("black",0.2),my.fade.col("black",0.6),"black")
								ltys<- c(1,1,4)								
								#plot rho for fx.tau.u
								rho.h.ok<- project.nABC.movingavg.gethist(ans.ok[["data"]]["ysigma2",acc.ok]/ans.ok[["xsigma2"]], 1, nbreaks= 70, width= 0.5, plot=0)
								rho.h.pw<- project.nABC.movingavg.gethist(ans.pw[["data"]]["ysigma2",acc.pw]/ans.pw[["xsigma2"]], 1, nbreaks= 70, width= 0.5, plot=0)	
								f.name<- paste(dir.name,"/nABC.Chisq_",N,"_",yn,"_",prior.u,"_",prior.l,"_rho.pdf",sep='')
								pdf(f.name,version="1.4",width=4,height=5)
								par(mar=c(5,4.5,0.5,0.5))
								xlim<- c(0,2.5)	#range(c(rho.h.ok$breaks,rho.h.pw$breaks))
								plot(1,1,type='n',bty='n',xlim=xlim,ylim=range(c(rho.h.ok$density,rho.h.pw$density)),ylab=expression("n-ABC estimate of "*pi[tau]*'('*rho*'|'*x*')'),xlab=expression(rho))
								plot(rho.h.ok, col=cols[1], border=NA, main='',freq=0, add=1)
								plot(rho.h.pw, col=cols[2], border=NA, main='',freq=0, add=1)
								abline(v=1,col=cols[3],lty=ltys[3])
								legend("topright",fill=c("transparent","transparent",cols[1],"transparent","transparent","transparent",cols[2],"transparent","transparent","transparent","transparent","transparent"),lty=c(NA,NA,ltys[1],NA,NA,NA,ltys[2],NA,NA,NA,NA,ltys[3]),border=NA,bty='n',legend=expression("n=200","","calibrated","tolerances",tau^'-'*"=0.454", tau^'+'*"=2.2","calibrated","tolerances",tau^'-'*"=0.678",tau^'+'*"=1.5","",rho^symbol("\x2a")))
								dev.off()																
							}	
							#print(out)
							out			
						})
				names(ans)<- f.name.yn$x
				f.name<- paste(dir.name,"/nABC.Chisq_largen_",N,"_",prior.u,"_",prior.l,".R",sep='')
				cat(paste("\nnABC.Chisq save 'ans' to ",f.name))				
				#save(ans,file=f.name)
				#print(ans)
			}
			
			sample.size<- as.numeric(names(ans))
			cols<- c(my.fade.col("black",0.6),my.fade.col("black",0.2),"black")
			#cols<- c("black","gray50","black")
			ltys<- c(1,1,4)
			pdf.width<- 4
			pdf.height<-5
			#c("fx.tau.u","fx.pw"),c("tau.l","tau.u","cl","cu","acc","mean","hmode","dmode","xsigma2")))
			
			f.name<- paste(dir.name,"/nABC.Chisq_largen_",N,"_",prior.u,"_",prior.l,"_mode.pdf",sep='')			
			pdf(f.name,version="1.4",width=pdf.width,height=pdf.height)				
			par(mar=c(4,4.25,1,.5))
			ylim<- range(	sapply(ans, function(x) x[,"dmode"]) 	)
			plot(1,1,type='n',bty='n',log='x',xlim=range(sample.size),ylim=ylim, xlab="sample size n",ylab=expression("numerically estimated mode of "*pi[abc]*'('*sigma^2*'|'*x*')'))
			abline(h=0.1,lty=2, col="gray60")
			z<- sapply(ans, function(x) x["fx.pw","dmode"])	- 1 + (sample.size-1)/sample.size
			points(sample.size, z,		pch=18,cex=0.5,col=cols[1])
			z<- sapply(ans, function(x) x["fx.tau.u","dmode"]) - 1 + (sample.size-1)/sample.size
			points(sample.size, z,	pch=20,cex=0.75,col=cols[2])	
			legend("topleft",bty='n',border=NA,fill=c(cols[1],"transparent",cols[2],"transparent","transparent","transparent"),lty=c(ltys[1],NA,ltys[2],NA,NA,ltys[3]),legend=c("fixed power &",expression("decreasing "*tau^'+'),"increasing power &",expression("fixed "*tau^'+'),"",expression(scriptstyle(frac(n-1,n)*sigma[0]^2))))
			lines(sample.size, (sample.size-1)/sample.size, lty=4)
			dev.off()
			stop()
			f.name<- paste(dir.name,"/nABC.Chisq_largen_",N,"_",prior.u,"_",prior.l,"_acc.pdf",sep='')			
			pdf(f.name,version="1.4",width=pdf.width,height=pdf.height)				
			par(mar=c(4,4.25,1,.5))
			ylim<- range(	sapply(ans, function(x) x[,"acc"]) 	)
			plot(1,1,type='n',bty='n',log='x',xlim=range(sample.size),ylim=ylim, xlab="sample size n",ylab="acceptance prob")							
			points(sample.size,sapply(ans, function(x) x["fx.pw","acc"]),		pch=18,cex=1.5,col=cols[1])
			points(sample.size,sapply(ans, function(x) x["fx.tau.u","acc"]),	pch=20,cex=1.5,col=cols[2])	
			legend("topleft",bty='n',border=NA,fill=c(cols[1],"transparent",cols[2],"transparent","transparent","transparent"),lty=c(ltys[1],NA,ltys[2],NA,NA,ltys[3]),legend=c("fixed power &",expression("decreasing "*tau^'+'),"increasing power &",expression("fixed "*tau^'+'),"",expression(hat(nu)[x])))			
			dev.off()
			
			f.name<- paste(dir.name,"/nABC.Chisq_largen_",N,"_",prior.u,"_",prior.l,"_tauu.pdf",sep='')			
			pdf(f.name,version="1.4",width=pdf.width,height=pdf.height)				
			par(mar=c(4,4.25,1,.5))
			ylim<- range(	sapply(ans, function(x) x[,"tau.u"]) 	)
			plot(1,1,type='n',bty='n',log='x',xlim=range(sample.size),ylim=ylim, xlab="sample size n",ylab=expression(tau^'+'))							
			lines(sample.size,sapply(ans, function(x) x["fx.pw","tau.u"]),		lwd=2,col=cols[1])
			lines(sample.size,sapply(ans, function(x) x["fx.tau.u","tau.u"]),	lwd=2,col=cols[2])	
			legend(x=200,y=2,bty='n',border=NA,fill=c(cols[1],"transparent",cols[2],"transparent"),lty=c(ltys[1],NA,ltys[2],NA),legend=c("fixed power &",expression("decreasing "*tau^'+'),"increasing power &",expression("fixed "*tau^'+')))			
			dev.off()
			
			f.name<- paste(dir.name,"/nABC.Chisq_largen_",N,"_",prior.u,"_",prior.l,"_taul.pdf",sep='')			
			pdf(f.name,version="1.4",width=pdf.width,height=pdf.height)				
			par(mar=c(4,4.25,1,.5))
			ylim<- range(	sapply(ans, function(x) x[,"tau.l"]) 	)
			plot(1,1,type='n',bty='n',log='x',xlim=range(sample.size),ylim=ylim, xlab="sample size n",ylab=expression(tau^'-'))							
			points(sample.size,sapply(ans, function(x) x["fx.pw","tau.l"]),		pch=18,cex=1.5,col=cols[1])
			points(sample.size,sapply(ans, function(x) x["fx.tau.u","tau.l"]),	pch=20,cex=1.5,col=cols[2])	
			legend(x=500,y=2,bty='n',border=NA,fill=c(cols[1],"transparent",cols[2],"transparent","transparent","transparent"),lty=c(ltys[1],NA,ltys[2],NA,NA,ltys[3]),legend=c("fixed power &",expression("decreasing "*tau^'+'),"increasing power &",expression("fixed "*tau^'+'),"",expression(hat(nu)[x])))			
			dev.off()
		}
		#
	}
	if(!is.na(subprog) && subprog==5)		#compute naive ABC
	{
		m<- NA
		
		if(exists("argv"))
		{
			tmp<- na.omit(sapply(argv,function(arg)
							{	switch(substr(arg,2,2),
										m= return(as.numeric(substr(arg,3,nchar(arg)))),NA)	}))
			if(length(tmp)>0) m<- tmp[1]
		}
		
		xn<- yn<- 60
		df<- yn-1
		alpha<- 0.01		
		tau.u<- 2.2 		
		#tau.l<- chisqstretch.tau.low(tau.u, df, alpha)
		tau.h<- 0.65
		
		ymu<- xmu<- 0
		xsigma2<- 1
		prior.u<- 4
		prior.l<- 0.2
		N<- 1e6
		nbreaks<- 75
		
		resume<- 1
		if(!is.na(m))
		{		
			f.name<- paste(dir.name,"/nABC.sig2_stdabc_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,"_m",m,".R",sep='')
			cat(paste("\nnABC.Chisq: compute ",f.name))
			options(show.error.messages = FALSE, warn=1)		
			readAttempt<-try(suppressWarnings(load(f.name)))						
			options(show.error.messages = TRUE)						
			if(!resume || inherits(readAttempt, "try-error"))
			{
				x<- rnorm(xn,xmu,sd=sqrt(xsigma2))							
				f.name<- paste(dir.name,"/nABC.sig2_stdabc_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,"_m",m,".R",sep='')
				ans.ok<- project.nABC.StretchedChi2.fix.x.uprior.ysig2(N,tau.l,tau.u,prior.l,prior.u,alpha,x,yn,ymu,with.c=0)
				cat(paste("\nnABC.Chisq: save ",f.name))
				save(ans.ok,file=f.name)				
			}
			else
				cat(paste("\nnABC.sig2_stdabc: resumed ",f.name))
		}	
		else
		{
			#load data
			cat(paste("\nnABC.sig2_stdabc",dir.name))
			f.name<- paste(dir.name,paste("nABC.sig2_stdabc_",N,".R",sep=''),sep='/')
			options(show.error.messages = FALSE, warn=1)		
			readAttempt<-try(suppressWarnings(load(f.name)))						
			options(show.error.messages = TRUE)	
			if(!resume || inherits(readAttempt, "try-error"))
			{
				cat(paste("\nnABC.sig2_stdabc generate",paste(dir.name,paste("nABC.sig2_stdabc_",N,".R",sep=''),sep='/')))
				f.name<- list.files(dir.name, pattern=paste("^nABC.sig2_stdabc",'',sep=''), full.names = TRUE)				
				ans<- lapply(seq_along(f.name),function(i)
						{
							out<- matrix(NA,4,11,dimnames=list(c("mean","mode","acc","kl"),c("o0.8","o0.4","o0.2","o0.1","o0.05","l0.8","l0.4","l0.2","l0.1","l0.05","xsigma2")))
							cat(paste("\nload",f.name[i]))
							readAttempt<-try(suppressWarnings(load( f.name[i] )))
							if(inherits(readAttempt, "try-error"))	
								return(out)
							
							s2			<- seq( min( ans.ok[["data"]]["ysigma2",] ), max( ans.ok[["data"]]["ysigma2",] ), length.out=1024 )
							s2.lkl		<- chisqstretch.sulkl(s2, xn, sqrt(xsigma2), trafo=1 )
							y1			<- s2.lkl / sum(s2.lkl)
							
							cus<- c(0.8,0.4,0.2,0.1,0.05)							
							out[,1:5]<- sapply(cus,function(cu)
									{
										#print(cu)
										acc		<- which( ans.ok[["data"]]["sy2-sx2",]<=cu  &  ans.ok[["data"]]["sy2-sx2",]>=-cu )
										#get KL										
										tmp		<- density( ans.ok[["data"]]["ysigma2",acc]	)
										y2		<- approx(tmp$x, tmp$y, s2, yleft=0, yright=0, rule=2)$y
										y2		<- y2 / sum(y2)										
										kl		<- log( y2/y1 )
										kl[ is.nan(kl) ]		<- 0
										kl[ is.infinite(kl) ]	<- NA
										kl[ y2==0 ]				<- 0
										kl						<- kl*y2
										#get mode
										tmp<- project.nABC.movingavg.gethist( ans.ok[["data"]]["ysigma2",acc], xsigma2, nbreaks, width=.5, plot=0)
										#print(c(tmp[["mean"]],tmp[["dmode"]]))
										c(tmp[["mean"]],tmp[["dmode"]],length(acc)/ncol(ans.ok[["data"]]),sum(kl))
									})
							out[,6:10]<- sapply(cus,function(cu)
									{
										#print(cu)
										acc		<- which( ans.ok[["data"]]["rho.mc",]<=cu  &  ans.ok[["data"]]["rho.mc",]>=-cu )
										#get KL										
										tmp		<- density( ans.ok[["data"]]["ysigma2",acc]	)
										y2		<- approx(tmp$x, tmp$y, s2, yleft=0, yright=0, rule=2)$y
										y2		<- y2 / sum(y2)										
										kl		<- log( y2/y1 )
										kl[ is.nan(kl) ]		<- 0
										kl[ is.infinite(kl) ]	<- NA
										kl[ y2==0 ]				<- 0
										kl						<- kl*y2
										#get mode										
										tmp		<- project.nABC.movingavg.gethist( ans.ok[["data"]]["ysigma2",acc], xsigma2, nbreaks, width=.5, plot=0)
										#print(c(tmp[["mean"]],tmp[["dmode"]]))
										c(tmp[["mean"]],tmp[["dmode"]],length(acc)/ncol(ans.ok[["data"]]),sum(kl))
									})
#AAA							
							out[,11]<- ans.ok[["xsigma2"]]
							out
						})
				f.name<- paste(dir.name,paste("nABC.sig2_stdabc_",N,".R",sep=''),sep='/')
				cat(paste("\nnABC.StretchedF.sigmainference save 'ans' to ",f.name))				
				#save(ans,file=f.name)
				
			}
			else
				cat(paste("\nnABC.sig2_stdabc loaded",paste(dir.name,paste("nABC.sig2_stdabc_",N,".R",sep=''),sep='/')))
			
			cols<- c(my.fade.col("black",0.2),my.fade.col("black",0.6),my.fade.col("black",0.8),my.fade.col("black",1))
			ltys<- c(1,1,1,4)
			
			ok.idx	<- which(sapply(seq_along(ans),function(i) !any(is.na(ans[[i]]))	))
			ans		<- lapply(ok.idx,function(i) ans[[i]] )
			#extract xsigma2
			xsig2<- sapply(seq_along(ans),function(i) ans[[i]][1,"xsigma2"])
			cat(paste("\nmean xsig2",mean(xsig2)))
			
			#extract modes for S2y-S2x and log.S2y-log.S2x for cus
			#mean of densigamma is S^2(x) / (n-4)
			#mode of densigamma is S^2(x) / (n)			and we have S^2(x)/(n-1)
			cuidx	<- c(1,2,3,4)
			#print(ans[[1]])
			#print(ans[[1]][,cuidx])
			#print(ans[[1]][,cuidx+5])
			mo.s2	<- sapply(seq_along(ans),function(i) (xn-1)/(xn)*ans[[i]]["mode",cuidx])
			mo.ls2	<- sapply(seq_along(ans),function(i) (xn-1)/(xn)*ans[[i]]["mode",cuidx+5])
			me.s2	<- sapply(seq_along(ans),function(i) (xn-1)/(xn)*ans[[i]]["mean",cuidx])
			me.ls2	<- sapply(seq_along(ans),function(i) (xn-1)/(xn)*ans[[i]]["mean",cuidx+5])
			d.mo.s2	<- sapply(seq_along(ans),function(i) (xn-1)/(xn)*ans[[i]]["mode",cuidx]-(xn-1)/(xn)*ans[[i]][1,"xsigma2"])
			d.mo.ls2<- sapply(seq_along(ans),function(i) (xn-1)/(xn)*ans[[i]]["mode",cuidx+5]-(xn-1)/(xn)*ans[[i]][1,"xsigma2"])
			d.me.s2	<- sapply(seq_along(ans),function(i) (xn-1)/(xn)*ans[[i]]["mean",cuidx]-(xn-1)/(xn-4)*ans[[i]][1,"xsigma2"])
			d.me.ls2<- sapply(seq_along(ans),function(i) (xn-1)/(xn)*ans[[i]]["mean",cuidx+5]-(xn-1)/(xn-4)*ans[[i]][1,"xsigma2"])
			
			
			table	<- sapply(seq_along(cuidx),function(i)
					{
						
						means	<- sapply(list(mo.s2[i,],mo.ls2[i,],me.s2[i,],me.ls2[i,],d.mo.s2[i,],d.mo.ls2[i,],d.me.s2[i,],d.me.ls2[i,]),mean)
						names(means)<- c("mo.s2","mo.ls2","me.s2","me.ls2","d.mo.s2","d.mo.ls2","d.me.s2","d.me.ls2")
						print(colnames(ans[[1]])[i])
						print(means)
						round(means[c(1,5,3,7,2,6,4,8)],d=6)
					})
			print(table)
			require(xtable)
			print(xtable(table,digits=3))
			
			stop()
			h.xsig2	<- project.nABC.movingavg.gethist(xsig2, xsigma2, nbreaks= 50, width= 0.5, plot=0)
			h.mo.s2	<- lapply(seq_len(nrow(mo.s2)),function(i)	project.nABC.movingavg.gethist(mo.s2[i,], xsigma2, nbreaks= 50, width= 0.5, plot=0)	)	 
			h.mo.ls2<- lapply(seq_len(nrow(mo.ls2)),function(i)	project.nABC.movingavg.gethist(mo.ls2[i,], xsigma2, nbreaks= 50, width= 0.5, plot=0)	)
			
			if(0)
			{
				xlim	<- range(c(sapply(h.mo.ls2, function(x)	x$breaks), h.xsig2$breaks))
				ylim	<- range(c(sapply(h.mo.ls2, function(x)	x$counts), h.xsig2$counts))						
				f.name<- paste(dir.name,"/nABC.sig2_stdabc_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,"_modes.pdf",sep='')
				#pdf(f.name,version="1.4",width=4,height=5)			
				par(mar=c(5,4,0.5,0.5))						
				plot(1,1,bty='n',type='n',xlim=xlim,ylim=ylim,xlab=expression("mode of "*pi[tau]*'('*sigma^2*'|'*x*')'),ylab="n-ABC repetitions")
				sapply(seq_along(h.mo.ls2),function(i)	plot(h.mo.ls2[[i]], add=1, col=cols[i], border=NA)	)			
				h<- h.xsig2
				sapply(seq_along(h$breaks)[-1],function(i)
						{
							lines( c(h$breaks[i-1],h$breaks[i]), rep(h$counts[i-1],2), col=cols[4], lty=ltys[4] )
							lines( rep(h$breaks[i],2),c(h$counts[i-1],h$counts[i]), col=cols[4], lty=ltys[4] )						
						})
			}
			if(1)
			{
				xlim	<- range(c(sapply(h.mo.s2, function(x)	x$breaks), h.xsig2$breaks))
				ylim	<- range(c(sapply(h.mo.s2, function(x)	x$counts), h.xsig2$counts))						
				f.name<- paste(dir.name,"/nABC.sig2_stdabc_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,"_modes.pdf",sep='')
				#pdf(f.name,version="1.4",width=4,height=5)			
				par(mar=c(5,4,0.5,0.5))						
				plot(1,1,bty='n',type='n',xlim=xlim,ylim=ylim,xlab=expression("mode of "*pi[tau]*'('*sigma^2*'|'*x*')'),ylab="n-ABC repetitions")
				sapply(seq_along(h.mo.s2),function(i)	plot(h.mo.s2[[i]], add=1, col=cols[i], border=NA)	)
				
				h<- h.xsig2
				sapply(seq_along(h$breaks)[-1],function(i)
						{
							lines( c(h$breaks[i-1],h$breaks[i]), rep(h$counts[i-1],2), col=cols[4], lty=ltys[4] )
							lines( rep(h$breaks[i],2),c(h$counts[i-1],h$counts[i]), col=cols[4], lty=ltys[4] )						
						})
			}
			stop()
		}
	}
	if(!is.na(subprog) && subprog==8)		#check MLE, yn>xn
	{
		m<- NA
		
		if(exists("argv"))
		{
			tmp<- na.omit(sapply(argv,function(arg)
							{	switch(substr(arg,2,2),
										m= return(as.numeric(substr(arg,3,nchar(arg)))),NA)	}))
			if(length(tmp)>0) m<- tmp[1]
		}
		for.mle	<- 1
		xn		<- yn<- 60
		df		<- yn-1
		alpha	<- 0.01		
		tau.u	<- 2.2 		
		tau.l	<- chisqstretch.tau.low(tau.u, df, alpha, for.mle=for.mle)
		tau.h	<- 0.65
		
		ymu<- xmu<- 0
		xsigma2<- 1
		prior.u<- 4
		prior.l<- 0.2
		N<- 1e6
		
		resume<- 1
		if(!is.na(m))
		{		
			f.name<- paste(dir.name,"/nABC.Chisq_mle_yn_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,"_m",m,".R",sep='')
			cat(paste("\nnABC.Chisq: compute ",f.name))
			options(show.error.messages = FALSE, warn=1)		
			readAttempt<-try(suppressWarnings(load(f.name)))						
			options(show.error.messages = TRUE)						
			if(!resume || inherits(readAttempt, "try-error"))
			{
				x		<- rnorm(xn,xmu,sd=sqrt(xsigma2))
				#x 		<- (x - mean(x))/sd(x) * sqrt(xsigma2) + xmu
				#a		<- (xn-2)/2	 
				#b		<- var(x)*(xn-1)/2
				#var.lkl	<- b*b/((a-1)*(a-1)*(a-2))														
				tmp		<- nabc.calibrate.m.and.tau.yesmxpw.yesKL("chisqstretch.kl", args = list(	n.of.x = xn, s.of.x = sd(x), n.of.y = xn, 
								s.of.y = NA, mx.pw = 0.9, alpha = alpha, 
								calibrate.tau.u = T, for.mle=1, tau.u = tau.u), plot = F)
				print(tmp)
				yn		<- tmp[1]
				tau.l	<- tmp[3]
				tau.u	<- tmp[4]								
				#tmp		<- chisqstretch.n.of.y(xn, sqrt(var.lkl), 0.9, alpha, tau.u.ub=tau.u, for.mle=1)
				#print(tmp)
				#stop()
				#yn		<- tmp[1]
				#tau.l	<- tmp[2]
				#tau.u	<- tmp[3]
				f.name	<- paste(dir.name,"/nABC.Chisq_mle_yn_",N,"_",xn,"_",prior.u,"_",prior.l,"_m",m,".R",sep='')
				ans.ok	<- project.nABC.StretchedChi2.fix.x.uprior.ysig2(N,tau.l,tau.u,prior.l,prior.u,alpha,x,yn,ymu, for.mle=for.mle)
				cat(paste("\nnABC.Chisq: save ",f.name))
				save(ans.ok,file=f.name)				
				ans.ok	<- NULL				
				
				yn		<- round(yn*3/100)*100
				tmp		<- chisqstretch.tau.lowup(0.9, tau.u, yn-1, alpha, for.mle=for.mle)
				tau.l	<- tmp[1]
				tau.u	<- tmp[2]
				f.name	<- paste(dir.name,"/nABC.Chisq_mle_yntoolarge_",N,"_",xn,"_",prior.u,"_",prior.l,"_m",m,".R",sep='')				
				ans.too	<- project.nABC.StretchedChi2.fix.x.uprior.ysig2(N,tau.l,tau.u,prior.l,prior.u,alpha,x,yn,ymu, for.mle=for.mle)
				cat(paste("\nnABC.Chisq: save ",f.name))
				save(ans.too,file=f.name)				
				ans.too	<- NULL								
				#ans.naive<- project.nABC.StretchedChi2.fix.x.uprior.ysig2(N,xsigma2-tau.h,xsigma2+tau.h,prior.l,prior.u,alpha,x,yn,ymu, for.mle=for.mle)
				#f.name<- paste(dir.name,"/nABC.Chisq_mle_naive_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,"_m",m,".R",sep='')
				#cat(paste("\nnABC.Chisq: save ",f.name))
				#save(ans.naive,file=f.name)
				#ans.wprior<- project.nABC.StretchedChi2.fix.x.stdprior.ysig2(N,tau.l,tau.u,prior.l,prior.u,alpha,x,yn,ymu)
				#f.name<- paste(dir.name,"/nABC.Chisq_wprior_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,"_m",m,".R",sep='')
				#cat(paste("\nnABC.Chisq: save ",f.name))
				#save(ans.wprior,file=f.name)
			}
			else
				cat(paste("\nnABC.MA: resumed ",f.name))
		}	
		else
		{
			#load data
			cat(paste("\nnABC.Chisq",dir.name))
			f.name<- paste(dir.name,"/nABC.Chisq_mle_yn_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,".R",sep='')
			options(show.error.messages = FALSE, warn=1)		
			readAttempt<-try(suppressWarnings(load(f.name)))						
			options(show.error.messages = TRUE)		
			resume<- 0
			if(!resume || inherits(readAttempt, "try-error"))
			{				
				f.name<- list.files(dir.name, pattern=paste("^nABC.Chisq_mle_yn_",sep=''), full.names = TRUE)
				tmp<- sort(sapply(strsplit(f.name,'_',fixed=1),function(x)	as.numeric(substr(x[length(x)],2,nchar(x[length(x)])-2))		), index.return=1)
				f.name<- f.name[tmp$ix]				
				f.name2<- list.files(dir.name, pattern=paste("^nABC.Chisq_mle_yntoolarge_",sep=''), full.names = TRUE)
				tmp2<- sort(sapply(strsplit(f.name2,'_',fixed=1),function(x)	as.numeric(substr(x[length(x)],2,nchar(x[length(x)])-2))		), index.return=1)
				f.name2<- f.name2[tmp2$ix]
				f.name<- rbind( 	f.name[tmp$x%in%intersect(tmp$x,tmp2$x)], f.name2[tmp2$x%in%intersect(tmp$x,tmp2$x)]	)
				print(f.name)
				cat(paste("\nnABC.Chisq load data: ", ncol(f.name)))
				ans<- lapply(seq_len(ncol(f.name)),function(j)
						{
							out<- matrix(NA,2,4,dimnames=list(c("ok","large"),c("mean","hmode","dmode","xsigma2")))
							
							cat(paste("\nload",f.name[1,j]))
							readAttempt<-try(suppressWarnings(load( f.name[1,j] )))
							if(inherits(readAttempt, "try-error"))	stop("error at ok")
							#tmp fix bug (now resolved)
							#ans.ok[["data"]]["error",]<- ans.ok[["data"]]["error",]*59/60
							#accept if T in boundaries					
							acc.ok<- which( ans.ok[["data"]]["error",]<=ans.ok[["cir"]]  &  ans.ok[["data"]]["error",]>=ans.ok[["cil"]] )
							acc.h.ok<- project.nABC.movingavg.gethist(ans.ok[["data"]]["ysigma2",acc.ok], ans.ok[["xsigma2"]], nbreaks= 50, width= 0.5, plot=0)
							out["ok",]<- c(acc.h.ok[["mean"]],acc.h.ok[["hmode"]],acc.h.ok[["dmode"]],ans.ok[["xsigma2"]])
							
							cat(paste("\nload",f.name[2,j]))	
							#ans.naive<- ans.ok
							readAttempt<-try(suppressWarnings(load( f.name[2,j] )))
							if(inherits(readAttempt, "try-error"))	stop("error at toolarge")	
							#tmp fix bug (now resolved)
							#ans.naive[["cil"]]<- 0.5084666; ans.naive[["cir"]]<- 1.009202							
							acc.too<- which( ans.too[["data"]]["error",]<=ans.too[["cir"]]  &  ans.too[["data"]]["error",]>=ans.too[["cil"]] )
							acc.h.too<- project.nABC.movingavg.gethist(ans.too[["data"]]["ysigma2",acc.too], ans.too[["xsigma2"]], nbreaks= 50, width= 0.5, plot=0)
							out["large",]<- c(acc.h.too[["mean"]],acc.h.too[["hmode"]],acc.h.too[["dmode"]],ans.too[["xsigma2"]])
							#print(length(acc.too) / ncol(ans.too[["data"]]))											
							if(1 && j==1)
							{								
								cols<- c(my.fade.col("black",0.6),my.fade.col("black",0.2),"black")
								ltys<- c(1,1,4,3)
								require(pscl)
								#plot sigma2
								f.name<- paste(dir.name,"/nABC.Chisq_mle_yyn_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,"_m",j,".pdf",sep='')
								print(f.name)
								pdf(f.name,version="1.4",width=4,height=5)
								par(mar=c(5,5,0.5,0.5))
								plot(acc.h.too, col=cols[2],border=NA,main='',freq=0,ylab=expression("numerical estimate of "*pi[abc]*'('*sigma^2*'|'*x*')'),xlab=expression(sigma^2),xlim=c(0,3),ylim=c(0,3))
								plot(acc.h.ok, col=cols[1],border=NA,main='',add=1,freq=0)
								#plot(acc.h.ok, col=cols[1],border=NA,main='',freq=0,ylab=expression("numerical estimate of "*pi[abc]*'('*sigma^2*'|'*x*')'),xlab=expression(sigma^2),xlim=c(0,3),ylim=c(0,3))
								a		<- (xn-2)/2	 
								b		<- ans.ok[["xsigma2"]]*(xn)/2
								var.lkl	<- b*b/((a-1)*(a-1)*(a-2))		
								
								s.of.x	<- sqrt( ans.ok[["xsigma2"]]/(xn-1)*xn )								
								tmp		<- nabc.calibrate.m.and.tau.yesmxpw.yesKL("chisqstretch.kl", args = list(	n.of.x = xn, s.of.x = s.of.x, n.of.y = xn, 
												s.of.y = NA, mx.pw = 0.9, alpha = alpha, 
												calibrate.tau.u = T, for.mle=1, tau.u = 2), plot = F)																			
								print(tmp)
								
								yn		<- round(tmp[1]*3/100)*100
								tmp		<- chisqstretch.tau.lowup(0.9, 2, yn-1, alpha, for.mle=1)
								tmp		<- chisqstretch.kl(xn, s.of.x, yn, NA, 0.9, alpha, calibrate.tau.u = F, tau.u = tmp[2], plot = F)
								print(tmp)
								
								x<- seq(prior.l,prior.u,0.001)
								y<- densigamma(x,a,b) / diff(pigamma(c(prior.l,prior.u),a,b))							
								lines(x,y,col=cols[3],lty=ltys[4])
								abline(v=ans.ok[["xsigma2"]],col=cols[3],lty=ltys[3])								
								legend("topright",fill=c("transparent","transparent",cols[1],"transparent","transparent","transparent","transparent",cols[2],"transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent"),lty=c(NA,NA,ltys[1],NA,NA,NA,NA,ltys[2],NA,NA,NA,NA,NA,ltys[4],NA,ltys[3],NA),border=NA,bty='n',legend=expression("n=60","","calibrated","tolerances",tau^'-'*"=0.572", tau^'+'*"=1.808","m=97","calibrated","tolerances",tau^'-'*"=0.726",tau^'+'*"=1.392","m=300","",pi*'('*sigma^2*'|'*x*')',"",argmax[sigma^2],pi*'('*sigma^2*'|'*x*')'))
								#legend("topright",fill=c("transparent","transparent",cols[1],"transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent"),lty=c(NA,NA,ltys[1],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,ltys[4],NA,ltys[3],NA),border=NA,bty='n',legend=expression("n=60","","calibrated","tolerances",tau^'-'*"=0.572", tau^'+'*"=1.808","m=97","","","","","","",pi*'('*sigma^2*'|'*x*')',"",argmax[sigma^2],pi*'('*sigma^2*'|'*x*')'))
								dev.off()
								stop()
							}							
							out			
						})
				
				f.name<- paste(dir.name,"/nABC.Chisq_ynmean_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,".R",sep='')
				cat(paste("\nnABC.Chisq save 'ans' to ",f.name))				
				#save(ans,file=f.name)
				print(ans)
				#stop()
			}
		}
	}
	if(!is.na(subprog) && subprog==7)		#check MLE
	{
		m<- NA
		
		if(exists("argv"))
		{
			tmp<- na.omit(sapply(argv,function(arg)
							{	switch(substr(arg,2,2),
										m= return(as.numeric(substr(arg,3,nchar(arg)))),NA)	}))
			if(length(tmp)>0) m<- tmp[1]
		}
		for.mle	<- 1
		xn		<- yn<- 60
		df		<- yn-1
		alpha	<- 0.01		
		tau.u	<- 2.2 		
		tau.l	<- chisqstretch.tau.low(tau.u, df, alpha, for.mle=for.mle)
		tau.h	<- 0.65
		
		ymu<- xmu<- 0
		xsigma2<- 1
		prior.u<- 4
		prior.l<- 0.2
		N<- 1e6
		
		resume<- 1
		if(!is.na(m))
		{		
			f.name<- paste(dir.name,"/nABC.Chisq_mle_naive_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,"_m",m,".R",sep='')
			cat(paste("\nnABC.Chisq: compute ",f.name))
			options(show.error.messages = FALSE, warn=1)		
			readAttempt<-try(suppressWarnings(load(f.name)))						
			options(show.error.messages = TRUE)						
			if(!resume || inherits(readAttempt, "try-error"))
			{
				x<- rnorm(xn,xmu,sd=sqrt(xsigma2))				
				f.name<- paste(dir.name,"/nABC.Chisq_mle_ok_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,"_m",m,".R",sep='')
				ans.ok<- project.nABC.StretchedChi2.fix.x.uprior.ysig2(N,tau.l,tau.u,prior.l,prior.u,alpha,x,yn,ymu, for.mle=for.mle)
				cat(paste("\nnABC.Chisq: save ",f.name))
				save(ans.ok,file=f.name)
				ans.ok<- NULL				
				#ans.naive<- project.nABC.StretchedChi2.fix.x.uprior.ysig2(N,xsigma2-tau.h,xsigma2+tau.h,prior.l,prior.u,alpha,x,yn,ymu, for.mle=for.mle)
				#f.name<- paste(dir.name,"/nABC.Chisq_mle_naive_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,"_m",m,".R",sep='')
				#cat(paste("\nnABC.Chisq: save ",f.name))
				#save(ans.naive,file=f.name)
				#ans.wprior<- project.nABC.StretchedChi2.fix.x.stdprior.ysig2(N,tau.l,tau.u,prior.l,prior.u,alpha,x,yn,ymu)
				#f.name<- paste(dir.name,"/nABC.Chisq_wprior_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,"_m",m,".R",sep='')
				#cat(paste("\nnABC.Chisq: save ",f.name))
				#save(ans.wprior,file=f.name)
			}
			else
				cat(paste("\nnABC.MA: resumed ",f.name))
		}	
		else
		{
			#load data
			cat(paste("\nnABC.Chisq",dir.name))
			f.name<- paste(dir.name,"/nABC.Chisq_mle_modemean_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,".R",sep='')
			options(show.error.messages = FALSE, warn=1)		
			readAttempt<-try(suppressWarnings(load(f.name)))						
			options(show.error.messages = TRUE)
			resume<- 0
			if(!resume || inherits(readAttempt, "try-error"))
			{				
				f.name<- list.files(dir.name, pattern=paste("^nABC.Chisq_mle_ok_",sep=''), full.names = TRUE)
				tmp<- sort(sapply(strsplit(f.name,'_',fixed=1),function(x)	as.numeric(substr(x[length(x)],2,nchar(x[length(x)])-2))		), index.return=1)
				f.name<- f.name[tmp$ix]				
				#f.name2<- list.files(dir.name, pattern=paste("^nABC.Chisq_mle_naive_",sep=''), full.names = TRUE)
				#tmp2<- sort(sapply(strsplit(f.name2,'_',fixed=1),function(x)	as.numeric(substr(x[length(x)],2,nchar(x[length(x)])-2))		), index.return=1)
				#f.name2<- f.name2[tmp2$ix]
				#f.name3<- list.files(dir.name, pattern=paste("^nABC.Chisq_wprior_",sep=''), full.names = TRUE)
				#tmp3<- sort(sapply(strsplit(f.name3,'_',fixed=1),function(x)	as.numeric(substr(x[length(x)],2,nchar(x[length(x)])-2))		), index.return=1)
				#f.name3<- f.name3[tmp3$ix]											
				#f.name<- rbind( 	f.name[tmp$x%in%intersect(tmp$x,tmp2$x)], f.name2[tmp2$x%in%intersect(tmp$x,tmp2$x)]	)
				f.name<- t(as.matrix(f.name))
				#f.name<- rbind( 	f.name[tmp$x%in%intersect(intersect(tmp$x,tmp2$x),tmp3$x)], 
				#		f.name2[tmp2$x%in%intersect(intersect(tmp$x,tmp2$x),tmp3$x)],
				#		f.name3[tmp3$x%in%intersect(intersect(tmp$x,tmp2$x),tmp3$x)]	)	
				cat(paste("\nnABC.Chisq load data: ", ncol(f.name)))
				ans<- lapply(seq_len(ncol(f.name)),function(j)
						{
							out<- matrix(NA,2,4,dimnames=list(c("ok","naive"),c("mean","hmode","dmode","xsigma2")))
							
							cat(paste("\nload",f.name[1,j]))
							readAttempt<-try(suppressWarnings(load( f.name[1,j] )))
							if(inherits(readAttempt, "try-error"))	stop("error at ok")
							#tmp fix bug (now resolved)
							ans.ok[["data"]]["error",]<- ans.ok[["data"]]["error",]*59/60
							#accept if T in boundaries					
							acc.ok<- which( ans.ok[["data"]]["error",]<=ans.ok[["cir"]]  &  ans.ok[["data"]]["error",]>=ans.ok[["cil"]] )
							acc.h.ok<- project.nABC.movingavg.gethist(ans.ok[["data"]]["ysigma2",acc.ok], ans.ok[["xsigma2"]], nbreaks= 100, width= 0.5, plot=0)
							out["ok",]<- c(acc.h.ok[["mean"]],acc.h.ok[["hmode"]],acc.h.ok[["dmode"]],ans.ok[["xsigma2"]])
							
							
							#cat(paste("\nload",f.name[2,j]))	
							ans.naive<- ans.ok
							#readAttempt<-try(suppressWarnings(load( f.name[2,j] )))
							#if(inherits(readAttempt, "try-error"))	stop("error at naive")	
							#tmp fix bug (now resolved)
							ans.naive[["cil"]]<- 0.5084666
							ans.naive[["cir"]]<- 1.009202
							
							acc.naive<- which( ans.naive[["data"]]["error",]<=ans.naive[["cir"]]  &  ans.naive[["data"]]["error",]>=ans.naive[["cil"]] )
							acc.h.naive<- project.nABC.movingavg.gethist(ans.naive[["data"]]["ysigma2",acc.naive], ans.naive[["xsigma2"]], nbreaks= 100, width= 0.5, plot=0)
							out["naive",]<- c(acc.h.naive[["mean"]],acc.h.naive[["hmode"]],acc.h.naive[["dmode"]],ans.naive[["xsigma2"]])
							print(length(acc.naive) / ncol(ans.naive[["data"]]))			
#print(out)										
							if(1 && j==2)
							{
								require(pscl)
								cols<- c(my.fade.col("black",0.2),my.fade.col("black",0.6),"black")
								ltys<- c(1,1,4,3)
								
								#plot rho
								rho.h.ok	<- project.nABC.movingavg.gethist(ans.ok[["data"]]["ysigma2",acc.ok]/ans.ok[["xsigma2"]], 1, nbreaks= 100, width= 0.5, plot=1)
								
								pw<- chisqstretch.pow(seq(prior.l,prior.u,by=0.001),yn,df,ans.ok[["cil"]],ans.ok[["cir"]])
								print( c(seq(prior.l,prior.u,by=0.001)[ which.max(pw) ], sum(pw)*0.001 ) )
								lines(seq(prior.l,prior.u,by=0.001),pw/(sum(pw)*0.001),type='l',col="red")
								abline(v=1,col="red")
								
								#rho.h.ok	<- project.nABC.movingavg.gethist(ans.ok[["data"]]["error",], 1, nbreaks= 100, width= 0.5, plot=1)
								#print(rho.h.ok[["dmode"]])
								#stop()
								rho.h.naive	<- project.nABC.movingavg.gethist(ans.naive[["data"]]["ysigma2",acc.naive]/ans.naive[["xsigma2"]], 1, nbreaks= 50, width= 0.5, plot=0)								
								f.name<- paste(dir.name,"/nABC.Chisq_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,"_m",j,"_rho.pdf",sep='')
								#pdf(f.name,version="1.4",width=4,height=5)
								par(mar=c(5,5,0.5,0.5))
								plot(1,1,type='n',bty='n',ylab=expression("n-ABC estimate of "*pi[tau]*'('*rho*'|'*x*')'),xlab=expression(rho),xlim=c(0,3),ylim=range(c(rho.h.ok$density,rho.h.naive$density)))
								plot(rho.h.ok, col=cols[1],border=NA,main='',add=1,freq=0)
								plot(rho.h.naive, col=cols[2],border=NA,main='',add=1,freq=0)
								abline(v=1,col=cols[3],lty=ltys[3])
								legend("topright",fill=c("transparent","transparent",cols[1],"transparent","transparent","transparent",cols[2],"transparent","transparent","transparent","transparent","transparent"),lty=c(NA,NA,ltys[1],NA,NA,NA,ltys[2],NA,NA,NA,NA,ltys[3]),border=NA,bty='n',legend=expression("n=60","","calibrated","tolerances",tau^'-'*"=0.477", tau^'+'*"=2.2","naive","tolerances",tau^'-'*"=0.35",tau^'+'*"=1.65","",rho^symbol("\x2a")))
								#dev.off()
								#plot sigma2
								f.name<- paste(dir.name,"/nABC.Chisq_mle_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,"_m",j,".pdf",sep='')
								print(f.name)
								pdf(f.name,version="1.4",width=4,height=5)
								par(mar=c(5,5,0.5,0.5))
								plot(acc.h.ok, col=cols[1],border=NA,main='',freq=0,ylab=expression("numerical estimate of "*pi[abc]*'('*sigma^2*'|'*x*')'),xlab=expression(sigma^2),xlim=c(0,3),ylim=c(0,3))								
								#plot(acc.h.naive, col=cols[2],border=NA,main='',freq=0,ylab=expression("numerical estimate of "*pi[abc]*'('*sigma^2*'|'*x*')'),xlab=expression(sigma^2),xlim=c(0,2.75),ylim=c(0,3))
								plot(acc.h.naive, col=cols[2],border=NA,main='',freq=0,add=1)
								x<- seq(prior.l,prior.u,0.001)
								y<- densigamma(x,(xn-2)/2,ans.ok[["xsigma2"]]*xn/2) / diff(pigamma(c(prior.l,prior.u),(xn-2)/2,ans.ok[["xsigma2"]]*xn/2))							
								lines(x,y,col=cols[3],lty=ltys[4])
								abline(v=ans.ok[["xsigma2"]],col=cols[3],lty=ltys[3])								
								legend("topright",fill=c("transparent","transparent",cols[1],"transparent","transparent","transparent",cols[2],"transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent"),lty=c(NA,NA,ltys[1],NA,NA,NA,ltys[2],NA,NA,NA,NA,ltys[4],NA,ltys[3],NA),border=NA,bty='n',legend=expression("n=60","","calibrated","tolerances",tau^'-'*"=0.477", tau^'+'*"=2.2","naive","tolerances",tau^'-'*"=0.35",tau^'+'*"=1.65","",pi*'('*sigma^2*'|'*x*')',"",argmax[sigma^2],pi*'('*sigma^2*'|'*x*')'))
								#legend("topright",fill=c("transparent","transparent","transparent","transparent","transparent","transparent",cols[2],"transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent"),lty=c(NA,NA,NA,NA,NA,NA,ltys[2],NA,NA,NA,NA,ltys[4],NA,ltys[3],NA),border=NA,bty='n',legend=expression("n=60","","","","", "","naive","tolerances",tau^'-'*"=0.35",tau^'+'*"=1.65","",pi*'('*sigma^2*'|'*x*')',"",argmax[sigma^2],pi*'('*sigma^2*'|'*x*')'))
								dev.off()
								stop()
							}							
							out			
						})
				
				f.name<- paste(dir.name,"/nABC.Chisq_modemean_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,".R",sep='')
				cat(paste("\nnABC.Chisq save 'ans' to ",f.name))				
				#save(ans,file=f.name)
				print(ans)
				#stop()
			}			
			# 
			#compute means 
			est.theta0<- 	sapply(seq_along(ans),function(i)	c(ans[[i]]["ok","xsigma2"],ans[[i]]["ok","dmode"],ans[[i]]["ok","mean"],ans[[i]]["naive","dmode"],ans[[i]]["naive","mean"])	)
			rownames(est.theta0)<- c("xsigma2","ok.mo","ok.me","naive.mo","naive.me")
			cat("\n estimated means\n")
			print( apply(est.theta0,1,mean) )
			print( mean(est.theta0["ok.mo",]-est.theta0["xsigma2",]) )
			print( mean(est.theta0["ok.me",]-est.theta0["xsigma2",]*xn/(xn-4)) )
			print( mean(est.theta0["naive.mo",]-est.theta0["xsigma2",]) )
			print( mean(est.theta0["naive.me",]-est.theta0["xsigma2",]) )
			
		}
		stop()
	}	
	if(!is.na(subprog) && subprog==1)		#check unbiasedness
	{
		m<- NA
		
		if(exists("argv"))
		{
			tmp<- na.omit(sapply(argv,function(arg)
							{	switch(substr(arg,2,2),
										m= return(as.numeric(substr(arg,3,nchar(arg)))),NA)	}))
			if(length(tmp)>0) m<- tmp[1]
		}
		
		xn<- yn<- 60
		df<- yn-1
		alpha<- 0.01		
		tau.u<- 2.2 		
		tau.l<- chisqstretch.tau.low(tau.u, df, alpha)
		tau.h<- 0.65
		
		ymu<- xmu<- 0
		xsigma2<- 1
		prior.u<- 4
		prior.l<- 0.2
		N<- 1e6
		
		resume<- 1
		if(!is.na(m))
		{		
			f.name<- paste(dir.name,"/nABC.Chisq_wprior_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,"_m",m,".R",sep='')
			cat(paste("\nnABC.Chisq: compute ",f.name))
			options(show.error.messages = FALSE, warn=1)		
			readAttempt<-try(suppressWarnings(load(f.name)))						
			options(show.error.messages = TRUE)						
			if(!resume || inherits(readAttempt, "try-error"))
			{
				x<- rnorm(xn,xmu,sd=sqrt(xsigma2))							
				f.name<- paste(dir.name,"/nABC.Chisq_ok_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,"_m",m,".R",sep='')
				ans.ok<- project.nABC.StretchedChi2.fix.x.uprior.ysig2(N,tau.l,tau.u,prior.l,prior.u,alpha,x,yn,ymu)
				cat(paste("\nnABC.Chisq: save ",f.name))
				save(ans.ok,file=f.name)
				ans.ok<- NULL				
				ans.naive<- project.nABC.StretchedChi2.fix.x.uprior.ysig2(N,xsigma2-tau.h,xsigma2+tau.h,prior.l,prior.u,alpha,x,yn,ymu)
				f.name<- paste(dir.name,"/nABC.Chisq_naive_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,"_m",m,".R",sep='')
				cat(paste("\nnABC.Chisq: save ",f.name))
				save(ans.naive,file=f.name)
				ans.wprior<- project.nABC.StretchedChi2.fix.x.stdprior.ysig2(N,tau.l,tau.u,prior.l,prior.u,alpha,x,yn,ymu)
				f.name<- paste(dir.name,"/nABC.Chisq_wprior_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,"_m",m,".R",sep='')
				cat(paste("\nnABC.Chisq: save ",f.name))
				save(ans.wprior,file=f.name)
			}
			else
				cat(paste("\nnABC.MA: resumed ",f.name))
		}	
		else
		{
			#load data
			cat(paste("\nnABC.Chisq",dir.name))
			f.name<- paste(dir.name,"/nABC.Chisq_modemean_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,".R",sep='')
			options(show.error.messages = FALSE, warn=1)		
			readAttempt<-try(suppressWarnings(load(f.name)))						
			options(show.error.messages = TRUE)
			resume<- 0
			if(!resume || inherits(readAttempt, "try-error"))
			{				
				f.name<- list.files(dir.name, pattern=paste("^nABC.Chisq_ok_",sep=''), full.names = TRUE)
				tmp<- sort(sapply(strsplit(f.name,'_',fixed=1),function(x)	as.numeric(substr(x[length(x)],2,nchar(x[length(x)])-2))		), index.return=1)
				f.name<- f.name[tmp$ix]				
				f.name2<- list.files(dir.name, pattern=paste("^nABC.Chisq_naive_",sep=''), full.names = TRUE)
				tmp2<- sort(sapply(strsplit(f.name2,'_',fixed=1),function(x)	as.numeric(substr(x[length(x)],2,nchar(x[length(x)])-2))		), index.return=1)
				f.name2<- f.name2[tmp2$ix]
				f.name3<- list.files(dir.name, pattern=paste("^nABC.Chisq_wprior_",sep=''), full.names = TRUE)
				tmp3<- sort(sapply(strsplit(f.name3,'_',fixed=1),function(x)	as.numeric(substr(x[length(x)],2,nchar(x[length(x)])-2))		), index.return=1)
				f.name3<- f.name3[tmp3$ix]											
				f.name<- rbind( 	f.name[tmp$x%in%intersect(tmp$x,tmp2$x)], f.name2[tmp2$x%in%intersect(tmp$x,tmp2$x)]	)	
				#f.name<- rbind( 	f.name[tmp$x%in%intersect(intersect(tmp$x,tmp2$x),tmp3$x)], 
				#					f.name2[tmp2$x%in%intersect(intersect(tmp$x,tmp2$x),tmp3$x)],
				#					f.name3[tmp3$x%in%intersect(intersect(tmp$x,tmp2$x),tmp3$x)]	)		
				cat(paste("\nnABC.Chisq load data: ", ncol(f.name)))
				ans<- lapply(seq_len(ncol(f.name)),function(j)
						{
							out<- matrix(NA,3,4,dimnames=list(c("ok","naive","wprior"),c("mean","hmode","dmode","xsigma2")))
							
							cat(paste("\nload",f.name[1,j]))
							readAttempt<-try(suppressWarnings(load( f.name[1,j] )))
							if(inherits(readAttempt, "try-error"))	stop("error at ok")														
							
							acc.ok<- which( ans.ok[["data"]]["error",]<=ans.ok[["cir"]]  &  ans.ok[["data"]]["error",]>=ans.ok[["cil"]] )
							acc.h.ok<- project.nABC.movingavg.gethist(ans.ok[["data"]]["ysigma2",acc.ok], ans.ok[["xsigma2"]], nbreaks= 50, width= 0.5, plot=0)
							out["ok",]<- c(acc.h.ok[["mean"]],acc.h.ok[["hmode"]],acc.h.ok[["dmode"]],ans.ok[["xsigma2"]])
							
							cat(paste("\nload",f.name[2,j]))
							readAttempt<-try(suppressWarnings(load( f.name[2,j] )))
							if(inherits(readAttempt, "try-error"))	stop("error at naive")														
							acc.naive<- which( ans.naive[["data"]]["error",]<=ans.naive[["cir"]]  &  ans.naive[["data"]]["error",]>=ans.naive[["cil"]] )
							acc.h.naive<- project.nABC.movingavg.gethist(ans.naive[["data"]]["ysigma2",acc.naive], ans.naive[["xsigma2"]], nbreaks= 50, width= 0.5, plot=0)
							out["naive",]<- c(acc.h.naive[["mean"]],acc.h.naive[["hmode"]],acc.h.naive[["dmode"]],ans.naive[["xsigma2"]])
							
							#cat(paste("\nload",f.name[3,j]))
							#readAttempt<-try(suppressWarnings(load( f.name[3,j] )))
							#if(inherits(readAttempt, "try-error"))	stop("error at wprior")														
							#acc.wprior<- which( ans.wprior[["data"]]["error",]<=ans.wprior[["cir"]]  &  ans.wprior[["data"]]["error",]>=ans.wprior[["cil"]] )
							#acc.h.wprior<- project.nABC.movingavg.gethist(ans.wprior[["data"]]["ysigma2",acc.wprior], ans.wprior[["xsigma2"]], nbreaks= 50, width= 0.5, plot=0)
							#out["wprior",]<- c(acc.h.wprior[["mean"]],acc.h.wprior[["hmode"]],acc.h.wprior[["dmode"]],ans.wprior[["xsigma2"]])
							
							if(1 && j==10)
							{
								cols<- c(my.fade.col("black",0.2),my.fade.col("black",0.6),"black")
								ltys<- c(1,1,4)
								
								#plot rho
								rho.h.ok<- project.nABC.movingavg.gethist(ans.ok[["data"]]["ysigma2",acc.ok]/ans.ok[["xsigma2"]], 1, nbreaks= 50, width= 0.5, plot=0)
								rho.h.naive<- project.nABC.movingavg.gethist(ans.naive[["data"]]["ysigma2",acc.naive]/ans.naive[["xsigma2"]], 1, nbreaks= 50, width= 0.5, plot=0)	
								f.name<- paste(dir.name,"/nABC.Chisq_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,"_m",j,"_rho.pdf",sep='')
								pdf(f.name,version="1.4",width=4,height=5)
								par(mar=c(5,5,0.5,0.5))								
								plot(1,1,type='n',bty='n',ylab=expression("n-ABC estimate of "*pi[tau]*'('*rho*'|'*x*')'),xlab=expression(rho),xlim=c(0,3),ylim=range(c(rho.h.ok$density,rho.h.naive$density)))
								plot(rho.h.ok, col=cols[1],border=NA,main='',add=1,freq=0)
								plot(rho.h.naive, col=cols[2],border=NA,main='',add=1,freq=0)
								abline(v=1,col=cols[3],lty=ltys[3])
								legend("topright",fill=c("transparent","transparent",cols[1],"transparent","transparent","transparent",cols[2],"transparent","transparent","transparent","transparent","transparent"),lty=c(NA,NA,ltys[1],NA,NA,NA,ltys[2],NA,NA,NA,NA,ltys[3]),border=NA,bty='n',legend=expression("n=60","","calibrated","tolerances",tau^'-'*"=0.477", tau^'+'*"=2.2","naive","tolerances",tau^'-'*"=0.35",tau^'+'*"=1.65","",rho^symbol("\x2a")))
								dev.off()
								#plot sigma2
								f.name<- paste(dir.name,"/nABC.Chisq_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,"_m",j,".pdf",sep='')
								pdf(f.name,version="1.4",width=4,height=5)
								par(mar=c(5,5,0.5,0.5))
								plot(acc.h.ok, col=cols[1],border=NA,main='',freq=0,ylab=expression("n-ABC estimate of "*pi[tau]*'('*sigma^2*'|'*x*')'),xlab=expression(sigma^2),xlim=c(0,3),ylim=c(0,2))
								plot(acc.h.naive, col=cols[2],border=NA,main='',add=1,freq=0)
								abline(v=ans.ok[["xsigma2"]],col=cols[3],lty=ltys[3])								
								legend("topright",fill=c("transparent","transparent",cols[1],"transparent","transparent","transparent",cols[2],"transparent","transparent","transparent","transparent","transparent"),lty=c(NA,NA,ltys[1],NA,NA,NA,ltys[2],NA,NA,NA,NA,ltys[3]),border=NA,bty='n',legend=expression("n=60","","calibrated","tolerances",tau^'-'*"=0.477", tau^'+'*"=2.2","naive","tolerances",tau^'-'*"=0.35",tau^'+'*"=1.65","",scriptstyle(frac(1,n-1))*S^2(x)))								
								dev.off()
								stop()
							}							
							out			
						})
				
				f.name<- paste(dir.name,"/nABC.Chisq_modemean_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,".R",sep='')
				cat(paste("\nnABC.Chisq save 'ans' to ",f.name))				
				#save(ans,file=f.name)
				print(ans[,1:5])
				stop()
			}			
			#
			#compute means 
			est.theta0<- 	sapply(seq_along(ans),function(i)	c(ans[[i]]["ok","xsigma2"],ans[[i]]["wprior","xsigma2"],ans[[i]]["ok","dmode"],ans[[i]]["ok","mean"],ans[[i]]["naive","dmode"],ans[[i]]["naive","mean"],ans[[i]]["wprior","dmode"])	)
			rownames(est.theta0)<- c("xsigma2","wprior.xsigma2","ok.mo","ok.me","naive.mo","naive.me","wprior.mo")
			cat("\n estimated means\n")
			print( apply(est.theta0,1,mean) )
			print( mean(est.theta0["ok.mo",]-est.theta0["xsigma2",]*(xn-1)/xn) )
			print( mean(est.theta0["ok.me",]-est.theta0["xsigma2",]*(xn-1)/(xn-4)) )
			print( mean(est.theta0["naive.mo",]-est.theta0["xsigma2",]) )
			print( mean(est.theta0["naive.me",]-est.theta0["xsigma2",]*(xn-1)/(xn-4)) )
			print( mean(est.theta0["wprior.mo",]-est.theta0["wprior.xsigma2",]) )
#stop()
			xsigma2.h<- project.nABC.movingavg.gethist(est.theta0["xsigma2",], xsigma2, nbreaks= 50, width= 0.5, plot=0)
			ok.mo.h<- project.nABC.movingavg.gethist(est.theta0["ok.mo",], xsigma2, nbreaks= 50, width= 0.5, plot=0)
			ok.me.h<- project.nABC.movingavg.gethist(est.theta0["ok.me",], xsigma2, nbreaks= 50, width= 0.5, plot=0)
			naive.mo.h<- project.nABC.movingavg.gethist(est.theta0["naive.mo",], xsigma2, nbreaks= 50, width= 0.5, plot=0)
			xlim<- range(c(ok.mo.h$breaks,ok.me.h$breaks,xsigma2.h$breaks))
			ylim<- range(c(ok.mo.h$counts,ok.me.h$counts,xsigma2.h$counts))
			
			f.name<- paste(dir.name,"/nABC.Chisq_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,"_modes.pdf",sep='')
			pdf(f.name,version="1.4",width=4,height=5)			
			par(mar=c(5,4,0.5,0.5))
			cols<- c(my.fade.col("black",0.2),my.fade.col("black",0.6),my.fade.col("black",1))
			ltys<- c(1,1,4)
			plot(1,1,bty='n',type='n',xlim=xlim,ylim=ylim,xlab=expression("mode of "*pi[tau]*'('*sigma^2*'|'*x*')'),ylab="n-ABC repetitions")
			plot(ok.mo.h, add=1, col=cols[1], border=NA)			
			plot(naive.mo.h, add=1, col=cols[2], border=NA)
			
			h<- xsigma2.h
			sapply(seq_along(h$breaks)[-1],function(i)
					{
						lines( c(h$breaks[i-1],h$breaks[i]), rep(h$counts[i-1],2), col=cols[3], lty=ltys[3] )
						lines( rep(h$breaks[i],2),c(h$counts[i-1],h$counts[i]), col=cols[3], lty=ltys[3] )
						
					})
			legend("topright",bty='n',border=NA,fill=c("transparent","transparent",cols[1],"transparent",cols[2],"transparent","transparent","transparent"),lty=c(NA,NA,ltys[1],NA,ltys[2],NA,NA,ltys[3]),legend=expression("n=60","","calibrated","tolerances","naive","tolerances","",scriptstyle(frac(1,n-1))*S^2(x)))
			dev.off()						
		}	
	}
	stop()
}
#------------------------------------------------------------------------------------------------------------------------
#power in the rejection region
project.nABC.StretchedF.pow.sym<- function(rho, df, cu, cl= 1/cu) pf( cu / rho, df1= df, df2= df ) - pf( cl / rho, df1= df, df2= df )
#------------------------------------------------------------------------------------------------------------------------
ms.figure2A<- function()		#illustrate full calibrations of scaled ChiSquare
{
	library(devtools)
	require(roxygen2)
	code.dir	<- "/Users/Oliver/git/abc.star"
	roxygenize(code.dir)
	devtools::install(code.dir)
		
	require(abc.star)
	
	ftest.calibrate(n.of.x=n.of.x, t2.x=t2.x, p=p, what='KL', mx.pw=0.9, alpha=0.01, use.R= FALSE, plot=FALSE, verbose=FALSE)
	
	n.of.x	<- 60 
	n.of.y	<- 60
	outdir	<- '~/duke/2015_ABC_resubmission_figs'
	
	s.of.x	<- 60/(60-1)
	cali	<- vartest.calibrate(n.of.x=n.of.x, s.of.x=s.of.x, what='KL', alpha=0.01, mx.pw=0.9, df=0)
	outdir	<- '~/duke/2015_ABC_resubmission_figs'
	file	<- paste(outdir, '/Fig_varcalibrated.pdf', sep='') 
	ggsave(w=3, h=4, file=file)
	
}
#------------------------------------------------------------------------------------------------------------------------
ms.pipeline<- function()		#illustrate power of scaled ChiSquare
{
	cmd.hpcwrapper<- function(cmd, hpc.walltime=71, hpc.mem="600mb", hpc.nproc='1', hpc.q='pqeph')
	{
		wrap<- "#!/bin/sh"
		#hpcsys<- HPC.CX1.IMPERIAL
		tmp	<- paste("#PBS -l walltime=",hpc.walltime,":59:59,pcput=",hpc.walltime,":45:00",sep='')
		wrap<- paste(wrap, tmp, sep='\n')		
		tmp	<- paste("#PBS -l select=1:ncpus=",hpc.nproc,":mem=",hpc.mem,sep='')
		wrap<- paste(wrap, tmp, sep='\n')
		wrap<- paste(wrap, "#PBS -j oe", sep='\n')
		if(!is.na(hpc.q))
				wrap<- paste(wrap, paste("#PBS -q",hpc.q), sep='\n\n')
		wrap<- paste(wrap, "module load intel-suite mpi R/3.1.2 raxml examl/2013-05-09 beast/1.8.0 beagle-lib/2014-07-30", sep='\n')
		
		cmd<- lapply(seq_along(cmd),function(i){	paste(wrap,cmd[[i]],sep='\n')	})
		if(length(cmd)==1)
			cmd<- unlist(cmd)
		cmd	
	}
	cmd.hpccaller<- function(outdir, outfile, cmd)
	{
		if( nchar( Sys.which("qsub") ) )
		{
			file	<- paste(outdir,'/',outfile,'.qsub',sep='')
			cat(paste("\nwrite HPC script to",file,"\n"))
			cat(cmd,file=file)
			cmd		<- paste("qsub",file)
			cat( cmd )
			cat( system(cmd, intern=TRUE) )
			Sys.sleep(1)
		}
		else
		{
			file	<- paste(outdir,'/',outfile,'.sh',sep='')
			cat(paste("\nwrite Shell script to\n",file,"\nStart this shell file manually\n"))
			cat(cmd,file=file)
			Sys.chmod(file, mode = "777")	
			Sys.sleep(1)
		}		
	}
	
	if(0)
	{
		cmd		<- paste(CODE.HOME,'/misc/nabc.startme.R',' -exe=VARTESTPREC',sep='')
		cmd		<- cmd.hpcwrapper(cmd, hpc.walltime=71, hpc.mem="600mb", hpc.nproc='1', hpc.q='pqeph')
		outdir	<- paste(DATA,"nABC.vt",sep='/')		
		outfile	<- paste("vt",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
		cmd.hpccaller(outdir, outfile, cmd)
	}
	if(1)
	{
		cmd		<- paste(CODE.HOME,'/misc/nabc.startme.R',' -exe=VARTESTEVAL',sep='')
		cmd		<- cmd.hpcwrapper(cmd, hpc.walltime=71, hpc.mem="600mb", hpc.nproc='1', hpc.q='pqeph')
		outdir	<- paste(DATA,"nABC.vt",sep='/')		
		outfile	<- paste("vt",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
		cmd.hpccaller(outdir, outfile, cmd)
	}
	
}
#------------------------------------------------------------------------------------------------------------------------
ms.figure1<- function()		#illustrate power of scaled ChiSquare
{
	require(abc.star)
	
	n.of.x	<- 60 
	n.of.y	<- 60
	outdir	<- '~/duke/2015_ABC_resubmission_figs'
	cali	<- vartest.calibrate(n.of.x=n.of.x, n.of.y=n.of.y, tau.l=1/2, tau.u=2, what='CR', alpha=0.01)
	
	# problematic ABC tolerances for ABC inference: 
	# although power is not zero, does not plateau around 1, and still high around rho=1 (desirable),
	# the power is not maximised at the point of equality (rho=1).
	
	rho		<- seq(0.1, 3, len=1024)
	tmp		<- data.frame(rho=rho, power=vartest.pow(rho, n.of.x, n.of.y-1, cali['c.l'], cali['c.u']))
	ggplot(tmp,aes(x=rho,y=power)) + geom_line() + labs(y='Power\n(ABC acceptance probability)')
	
	#	tau.l decreasing
	tau.u	<- 2
	tmp	<- lapply(c(0.2,0.4,0.6,0.8),function(tau.l)
			{
				cali	<- vartest.calibrate(n.of.x=n.of.x, n.of.y=n.of.y, tau.l=tau.l, tau.u=tau.u, what='CR', alpha=0.01)
				data.table(rho=rho, power=vartest.pow(rho, n.of.x, n.of.y-1, cali['c.l'], cali['c.u']), tau.l=tau.l)			
			})
	tmp	<- do.call('rbind', tmp)
	set(tmp, NULL, 'tau.l', tmp[, factor(tau.l)])
	ggplot(tmp,aes(x=rho, y=power, linetype=tau.l, group=tau.l)) +
			geom_vline(xintercept=1, colour='grey30', lwd=1) +
			geom_line() + 
			labs(x=expression(rho), y='Power\n(ABC approximation to likelihood)', linetype=expression(tau^'-')) +
			scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.2), lim=c(-0.01, 1.01)) +
			theme_bw() +
			theme(legend.position=c(1,1), legend.justification=c(1,1)) 
	file	<- paste(outdir, '/Fig_vartauldecr.pdf', sep='') 
	ggsave(w=5, h=5, file=file)
	
	

	#	tau.u increasing
	tmp	<- lapply(seq(1.2,2.8,0.4),function(tau.u)
			{
				cali	<- vartest.calibrate(n.of.x=n.of.x, n.of.y=n.of.y, tau.l=1/tau.u, tau.u=tau.u, what='CR', alpha=0.01)
				data.table(rho=rho, power=vartest.pow(rho, n.of.x, n.of.y-1, cali['c.l'], cali['c.u']), tau.u=tau.u)			
			})
	tmp	<- do.call('rbind', tmp)
	set(tmp, NULL, 'tau.u', tmp[, factor(tau.u)])
	ggplot(tmp,aes(x=rho, y=power, linetype=tau.u, group=tau.u)) +
			geom_vline(xintercept=1, colour='grey30', lwd=1) +
			geom_line() + 
			labs(x=expression(rho), y='Power\n(ABC approximation to likelihood)', linetype=expression(tau^'+')) +
			scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.2), lim=c(-0.01, 1.01)) +
			theme_bw() +
			theme(legend.position=c(1,1), legend.justification=c(1,1)) 
	file	<- paste(outdir, '/Fig_vartauincr.pdf', sep='') 
	ggsave(w=5, h=5, file=file)
	
	
	#	m increasing
	tau.u	<- 1.6
	rho		<- seq(0.5, 2, len=1024)
	tmp	<- lapply(c(60, 70, 90, 120, 360),function(n.of.y)
			{
				cali	<- vartest.calibrate(n.of.x=n.of.x, n.of.y=n.of.y, tau.l=1/tau.u, tau.u=tau.u, what='CR', alpha=0.01)
				data.table(rho=rho, power=vartest.pow(rho, n.of.x, n.of.y-1, cali['c.l'], cali['c.u']), n.of.y=n.of.y)			
			})
	tmp	<- do.call('rbind', tmp)
	set(tmp, NULL, 'n.of.y', tmp[, factor(n.of.y)])
	ggplot(tmp,aes(x=rho, y=power, linetype=n.of.y, group=n.of.y)) +
			geom_vline(xintercept=1, colour='grey30', lwd=1) +
			geom_line() + 
			labs(x=expression(rho), y='Power\n(ABC approximation to likelihood)', linetype='m') +
			scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.2), lim=c(-0.01, 1.01)) +
			theme_bw() +
			theme(legend.position=c(1,1), legend.justification=c(1,1)) 
	file	<- paste(outdir, '/Fig_varmincr.pdf', sep='') 
	ggsave(w=5, h=5, file=file)			
}

#------------------------------------------------------------------------------------------------------------------------
project.nABC.StretchedF<- function()
{	
	#mean of pdf corresponding to power function
	project.nABC.StretchedF.pow.mean<- function(rho1,rho2, df, cl, cu)
	{
		rho.n<- 1e-3
		rho<- seq(rho1,rho2,by=rho.n)		
		pow<- project.nABC.StretchedF.pow.sym(rho, df, cu, cl)
		#plot(rho, pow, type='l'); print( sum( pow * rho.n ) );	print( sum( pow * rho.n * rho) )
		sum( pow * rho.n * rho)
	}
	#integral of pdf corresponding to power function
	project.nABC.StretchedF.pow.int<- function(rho1,rho2, df, cl, cu)
	{
		rho.n<- 1e-3
		rho<- seq(rho1,rho2,by=rho.n)		
		pow<- project.nABC.StretchedF.pow.sym(rho, df, cu, cl)		 
		marg.lkl<- sum( pow / rho * rho.n)
#plot(rho, pow / rho / marg.lkl, type='l');  lines(rho, pow / marg.lkl, col="red");		
#print(sum(pow / marg.lkl / rho * rho.n)); print(sum(pow / marg.lkl * rho.n)); print(log(rho2/rho1))
		sum( pow * rho.n ) / marg.lkl 
	}
	
	#simulate N times from same ysigma, simulate once from xsigma
	project.nABC.StretchedF.fix.x.fix.ysigma<- function(N,args,xn,xmu,xsigma,yn,ymu,ysigma)		
	{
		x<- rnorm(xn,xmu,xsigma)
		ans<- sapply(1:N,function(i)
				{
					y<- rnorm(yn,ymu,ysigma)
					ans<- get.dist.fstretch(y, x, args=args)
					ans[c("error","cil","cir","pval")]
				})		
		ans
	}
	#simulate N times from same ysigma, simulate once from xsigma and then resample
	project.nABC.StretchedF.resample.x.fix.ysigma<- function(N,args,xn,xmu,xsigma,yn,ymu,ysigma)		
	{
		x<- rnorm(xn,xmu,xsigma)
		ans<- sapply(1:N,function(i)
				{
					z<- sample(x,xn,replace=1)
					y<- rnorm(yn,ymu,ysigma)
					ans<- get.dist.fstretch(y, z, args=args)
					ans[c("error","cil","cir","pval")]
				})		
		ans
	}
	#simulate N times from same ysigma, simulate from same xsigma
	project.nABC.StretchedF.fix.xsigma.fix.ysigma<- function(N,args,xn,xmu,xsigma,yn,ymu,ysigma)		
	{		
		ans<- sapply(1:N,function(i)
				{
					x<- rnorm(xn,xmu,xsigma)
					y<- rnorm(yn,ymu,ysigma)
					ans<- get.dist.fstretch(y, x, args=args)
					ans[c("error","cil","cir","pval")]
				})		
		ans
	}
	#simulate N times from same ytau.u, simulate from same xsigma
	project.nABC.StretchedF.fix.xsigma.fix.ytau.u<- function(N,tau.l,tau.u,prior.l,prior.u,alpha,xn,xmu,xsigma,yn,ymu)		
	{		
		if(tau.u<1)	stop("project.nABC.StretchedF.fix.xsigma.fix.ytau.u: error at 1a")
		if(tau.l>1)	stop("project.nABC.StretchedF.fix.xsigma.fix.ytau.u: error at 1b")
		args<- paste("fstretch",tau.l,tau.u,alpha,sep='/')		
		ans<- sapply(1:N,function(i)
				{
					x<- rnorm(xn,xmu,xsigma)
					ysigma2<- runif(1,prior.l,prior.u)
					y<- rnorm(yn,ymu,sd=sqrt(ysigma2))
					ans<- get.dist.fstretch(x, y, args=args)					
					ans[11]<- ysigma2
					names(ans)[11]<- "ysigma2"
					ans[c("ysigma2","error","cil","cir","pval")]
				})		
		ans
	}
	#simulate N times from same ysigma sth log(ysigma)~U(log(prior.l),log(prior.u)), simulate from same xsigma
	project.nABC.StretchedF.fix.xsigma.fix.ynormtau.u<- function(N,tau.l,tau.u,prior.l,prior.u,alpha,xn,xmu,xsigma,yn,ymu)		
	{		
		if(tau.u<1)		stop("project.nABC.StretchedF.fix.xsigma.fix.ynormtau.u: error at 1a")
		if(tau.l>1)		stop("project.nABC.StretchedF.fix.xsigma.fix.ynormtau.u: error at 1b")
		if(prior.u<1)	stop("project.nABC.StretchedF.fix.xsigma.fix.ynormtau.u: error at 1c")
		if(prior.l>1)	stop("project.nABC.StretchedF.fix.xsigma.fix.ynormtau.u: error at 1d")
		
		args<- paste("fstretch",tau.l,tau.u,alpha,sep='/')		
		ans<- sapply(1:N,function(i)
				{
					x<- rnorm(xn,xmu,xsigma)
					ysigma2<- exp( rtnorm(1, mean=0, sd= diff(log(c(prior.l,prior.u))) / sqrt(12), lower=log(prior.l), upper=log(prior.u)) )
					y<- rnorm(yn,ymu,sd=sqrt(ysigma2))
					ans<- get.dist.fstretch(y, x, args=args)					
					ans[11]<- ysigma2
					names(ans)[11]<- "ysigma2"
					ans[c("ysigma2","error","cil","cir","pval")]
				})		
		ans
	}
	#simulate N times from same ysigma sth log(ysigma)~NORM(log(prior.l),log(prior.u)), simulate from same xsigma
	project.nABC.StretchedF.fix.xsigma.fix.ylogtau.u<- function(N,tau.l,tau.u,prior.l,prior.u,alpha,xn,xmu,xsigma,yn,ymu)		
	{		
		if(tau.u<1)		stop("project.nABC.StretchedF.fix.xsigma.fix.ylogtau.u: error at 1a")
		if(tau.l>1)		stop("project.nABC.StretchedF.fix.xsigma.fix.ylogtau.u: error at 1b")
		if(prior.u<1)	stop("project.nABC.StretchedF.fix.xsigma.fix.ylogtau.u: error at 1c")
		if(prior.l>1)	stop("project.nABC.StretchedF.fix.xsigma.fix.ylogtau.u: error at 1d")
		
		args<- paste("fstretch",tau.l,tau.u,alpha,sep='/')		
		ans<- sapply(1:N,function(i)
				{
					x<- rnorm(xn,xmu,xsigma)
					ysigma2<- exp( runif(1,log(prior.l),log(prior.u)) )
					y<- rnorm(yn,ymu,sd=sqrt(ysigma2))
					ans<- get.dist.fstretch(x, y, args=args)					
					ans[13]<- ysigma2
					names(ans)[13]<- "ysigma2"					
					ans[c("ysigma2","error","cil","cir","pval")]
				})		
		ans
	}
	#simulate N times from same ysigma sth log(ysigma)~NORM(log(prior.l),log(prior.u)), simulate from same xsigma
	project.nABC.naive.fix.xsigma.fix.ylogtau.u<- function(N,prior.l,prior.u,alpha,xn,xmu,xsigma,yn,ymu)		
	{		
		if(prior.u<1)	stop("project.nABC.naive.fix.xsigma.fix.ylogtau.u: error at 1c")
		if(prior.l>1)	stop("project.nABC.naive.fix.xsigma.fix.ylogtau.u: error at 1d")
		
		
		ans<- sapply(1:N,function(i)
				{
					x<- rnorm(xn,xmu,xsigma)
					ysigma2<- exp( runif(1,log(prior.l),log(prior.u)) )
					y<- rnorm(yn,ymu,sd=sqrt(ysigma2))
					
					tmp<- c(var(y),var(x),sd(y),sd(x))
					ans<- numeric(5)
					names(ans)<- c("ysigma2","sy2-sx2","log.sy2-log.sx2","sy-sx","log.sy-log.sx")
					ans<- c(ysigma2,diff(tmp[1:2]),diff(log(tmp[1:2])),diff(tmp[3:4]),diff(log(tmp[3:4])))					
					ans
				})		
		rownames(ans)<- c("ysigma2","sy2-sx2","log.sy2-log.sx2","sy-sx","log.sy-log.sx")
		ans
	}
	
	#simulate N times from same ytau.u, resample fixed x
	project.nABC.StretchedF.resample.x.fix.ytau.u<- function(N,tau,prior,alpha,x,yn,ymu)		
	{		
		if(tau<1)	stop("project.nABC.StretchedF.resample.x.fix.ytau.u: error at 1a")
		args<- paste("fstretch",tau,alpha,sep='/')		
		ans<- sapply(1:N,function(i)
				{					
					z<- sample(x,length(x),replace=1)
					ysigma2<- runif(1,1/prior,prior)
					y<- rnorm(yn,ymu,sd=sqrt(ysigma2))
					ans<- get.dist.fstretch(y, z, args=args)
					ans[11]<- ysigma2
					names(ans)[11]<- "ysigma2"
					ans[c("ysigma2","error","cil","cir","pval")]
				})		
		ans
	}
	#simulate N times from same ytau.u, take fixed x
	project.nABC.StretchedF.fix.x.fix.ytau.u<- function(N,tau,prior,alpha,x,yn,ymu)		
	{		
		if(tau<1)	stop("project.nABC.StretchedF.fix.x.fix.ytau.u: error at 1a")
		args<- paste("fstretch",tau,alpha,sep='/')		
		ans<- sapply(1:N,function(i)
				{					
					ysigma2<- runif(1,1/prior,prior)
					y<- rnorm(yn,ymu,sd=sqrt(ysigma2))
					ans<- get.dist.fstretch(y, x, args=args)
					ans[13]<- ysigma2
					names(ans)[13]<- "ysigma2"
					ans[c("ysigma","error","cil","cir","pval")]
				})		
		ans
	}
	
	project.nABC.StretchedF.pval.qq<- function(x, add, ...){
		x<- sort(x)
		print(c(x[1],x[length(x)]))
		e.cdf <- (seq_along(x)-0.5) / length(x)
		x<- c(0,x,1)
		e.cdf<- c(0,e.cdf,1)
		if(!add)
		{
			plot(1,1,type='n',xlim=c(0,1),ylim=c(0,1),xlab="expected quantiles under point H0",ylab="observed quantiles")
		}
		points(e.cdf,x,type='s',...)
	}
	
	project.nABC.StretchedF.gethist<- function(x, acc, nbreaks, plot=0)
	{
		#compute break points sth 1 is in the middle					
		breaks<- max( diff(c(1,max(x["ysigma2",acc])*1.1)) / nbreaks,	diff(c(min(x["ysigma2",acc])*0.9,1)) / nbreaks )						
		breaks<- c( rev(seq(from= 1-breaks/2, by=-breaks, length.out= nbreaks )), seq(from= 1+breaks/2, by=breaks, length.out= nbreaks ) )
		tmp<- hist(x["ysigma2",acc],breaks=breaks, plot= 0)
		if(plot)	plot(tmp)
		c( mean(x["ysigma2",acc]), tmp[["mids"]][which.max( tmp[["counts"]] )] )		
	}
	
	
	dir.name<- "/Users/Oliver/duke/2012_frequencyABC/sim.data"
	n<- m<- 100
	cu<- 1.5; 		cl<- 1/1.5				#typical ABC thresholds
	df<- n -1
	
	
	if(0)
	{
		file<- paste(CODE.HOME,paste("lib",paste("libphylodyr",.Platform$dynlib.ext,sep=''),sep='/'),sep='')		
		dyn.load(file)
		
		tau.u<- 2.77		#4 not OK, 3 OK.
		alpha<- 0.05
		err.cur<- err.tol<- 1e-10
		maxit<- 100
		incit<- 0.05
		xn<- yn<- 60
		rej<- numeric(2)
		verbose<- 1
		
		tau.up<- tau.u		
		#if pwi-crit<0 we want to increase tau.low
		#if pwi-crit>0 we want to decrease tau.low
		
		#check if calibration possible sth tau.l<1
		tau.low<- 1 
		tmp<- .Call("abcScaledF",	c(xn-1,yn-1,tau.low,tau.up,alpha,err.tol,maxit,0.05)	)
		if(tmp[4]>1e-10)	stop("error at 1a")
		err.cur<- project.nABC.StretchedF.pow.int(tau.low, tau.up, xn, tmp[1], tmp[2]) - 1
		if(verbose)		print(c(tau.low,err.cur))	
		if(err.cur<0)	stop("cannot calibrate sth tau.l<1")		
		
		#check if calibration possible sth tau.l>0 (this should always be OK)
		tau.low<- 1e-8
		tmp<- .Call("abcScaledF",	c(xn-1,yn-1,tau.low,tau.up,alpha,err.tol,maxit,0.05)	)
		if(tmp[4]>1e-10)	stop("error at 1b")
		err.cur<- project.nABC.StretchedF.pow.int(tau.low, tau.up, xn, tmp[1], tmp[2]) - 1
		if(verbose)		print(err.cur)
		if(err.cur>0)	stop("cannot calibrate sth tau.l>1e-8")
		
		
		#determine incit
		tau.low<- 1/tau.up
		tmp<- .Call("abcScaledF",	c(xn-1,yn-1,tau.low,tau.up,alpha,err.tol,maxit,0.05)	)
		if(tmp[4]>1e-10)	stop("error at 1c")
		rej<- tmp[1:2]			
		pwi<- project.nABC.StretchedF.pow.int(tau.low, tau.up, xn, rej[1], rej[2])
		err.cur<- pwi - 1
		incit<- - err.cur/maxit
		if(verbose)
		{
			cat(paste("tau.low",tau.low,"tau.up",tau.up,"cil",rej[1],"cir",rej[2],"pwi",pwi,"crit",1,'\n'))
			print(c(err.cur, incit ))
		}
		
		#if incit<0 we must increase tau.low 		until 'pwi' is smaller than required. 
		#if incit>0 we must decrease tau.low		until 'pwi' is larger than required.	
		while((incit<0 & err.cur>0) | (incit>0 & err.cur<0))
		{	
			if(incit>0 & tau.low + incit > 1)	
				incit<- incit/2
			if(incit<0 & tau.low + incit < 0)	
				incit<- incit/2			
			tau.low<- tau.low + incit
			tmp<- .Call("abcScaledF",	c(xn-1,yn-1,tau.low,tau.up,alpha,err.tol,maxit,0.05)	)
			if(tmp[4]>1e-10)	stop("xx: error at 1a")
			rej<- tmp[1:2]			
			pwi<- project.nABC.StretchedF.pow.int(tau.low, tau.up, xn, rej[1], rej[2])
			err.cur<- pwi - 1									
			if(verbose)	cat(paste("tau.low",tau.low,"tau.up",tau.up,"cil",rej[1],"cir",rej[2],"pwi",pwi,"crit",1,'\n'))	
		}
		if(incit<0)						#the correct tau.l must be between [ tau.low, tau.low-incit ]
			tau.up<- tau.low-incit		
		else							#the correct tau.l must be between [ tau.low-init, tau.low ]
		{	
			tau.up<- tau.low
			tau.low<- tau.low-incit
		}
		if(verbose)	print(c(tau.low,tau.up))		
		while(abs(err.cur)>=err.tol & maxit>0)
		{
			tau.l<- (tau.up + tau.low)/2
			tmp<- .Call("abcScaledF",	c(xn-1,yn-1,tau.l,tau.u,alpha,err.tol,maxit,0.05)	)
			if(tmp[4]>1e-10)	stop("error at 2a")
			rej<- tmp[1:2]			
			pwi<- project.nABC.StretchedF.pow.int(tau.l, tau.u, xn, rej[1], rej[2])
			err.cur<- pwi - 1
			if(verbose)	cat(paste("tau.l",tau.l,"cil",rej[1],"cir",rej[2],"pwi",pwi,"crit",1,'\n'))
			#if(round(tau.low,digits=4)==round(tau.up,digits=4))
			#	break
			if(err.cur>0)
				tau.up<- tau.l
			if(err.cur<0)
				tau.low<- tau.l
			maxit<- maxit-1
			if(verbose)	cat(paste("tau.l",tau.l,"err.cur",pwi - 1,"tau.low",tau.low,"tau.up",tau.up,"maxit",maxit,'\n'))			
		}	
		
		cat(paste("tau.l",tau.l,"err.cur",pwi - 1,"tau.low",tau.low,"tau.up",tau.up,"maxit",maxit,'\n'))		
		
		
		stop()
	}
	if(0)
	{
		#require(multicore, lib.loc= LIB.LOC)
		require(msm, lib.loc= LIB.LOC)
		file<- paste(CODE.HOME,paste("lib",paste("libphylodyr",.Platform$dynlib.ext,sep=''),sep='/'),sep='')		
		dyn.load(file)
		#ABC inference from within (-tau,tau) when observed data is resimulated, resampled or fixed
		
		N<- 1e6		#typically 5e4
		M<- 2e3		#typically 500
		m<- NA		
		if(exists("argv"))
		{
			tmp<- na.omit(sapply(argv,function(arg)
							{	switch(substr(arg,2,2),
										m= return(as.numeric(substr(arg,3,nchar(arg)))),NA)	}))
			if(length(tmp)>0) m<- tmp[1]
		}
		
		xn<- yn<- 61		
		xmu<- ymu<- 0
		xsigma<- 1 
		#tau.u<- 1.5;	tau.l<- 0.700134618892036
		prior.u<- 4
		prior.l<- 0.2
		tau.u<- 2.77
		tau.l<- 0.27 #calibration for post mean: 0.29085995330929 
		alpha<- 0.01
		fade<- 0.5
		nproc<- 8
		nbreaks<- 75
		package.mkdir(DATA,"nABC.StretchedF")
		dir.name<- paste(DATA,"nABC.StretchedF",sep='/')				
		resume<- 1
		
		if(!is.na(m))
		{		
			f.name<- paste(dir.name,"/nnABC.naive.sigmainference_",N,"_",m,".R",sep='')
			cat(paste("\nnABC.StretchedF: compute ",f.name))
			options(show.error.messages = FALSE, warn=1)		
			readAttempt<-try(suppressWarnings(load(f.name)))						
			options(show.error.messages = TRUE)						
			if(!resume || inherits(readAttempt, "try-error"))
			{
				ans<- project.nABC.naive.fix.xsigma.fix.ylogtau.u(N,prior.l,prior.u,alpha,xn,xmu,xsigma,yn,ymu)
				cat(paste("\nnABC.StretchedF: save ",f.name))
				save(ans,file=f.name)				
			}
			else
				cat(paste("\nnABC.MA: resumed ",f.name))
		}
		else
		{
			#load data
			cat(paste("\nnABC.StretchedF.sigmainference",dir.name))
			f.name<- paste(dir.name,paste("nABC.naive.sigmainference_",N,".R",sep=''),sep='/')
			options(show.error.messages = FALSE, warn=1)		
			readAttempt<-try(suppressWarnings(load(f.name)))						
			options(show.error.messages = TRUE)		
			if(!resume || inherits(readAttempt, "try-error"))
			{
				cat(paste("\nnABC.StretchedF.sigmainference generate",paste(dir.name,paste("nABC.naive.sigmainference_",N,".R",sep=''),sep='/')))
				f.name<- list.files(dir.name, pattern=paste("^nnABC.naive.sigmainference",'',sep=''), full.names = TRUE)
				ans<- sapply(seq_along(f.name),function(i)
						{
							cat(paste("\nload",f.name[i]))
							readAttempt<-try(suppressWarnings(load( f.name[i] )))
							if(inherits(readAttempt, "try-error"))	stop("error")
							out<- matrix(NA,2,10,dimnames=list(c("mean","mode"),c("o0.5","o0.4","o0.3","o0.2","o0.1","l0.5","l0.4","l0.3","l0.2","l0.1")))
							cu<- 0.5
							out[,1]<- project.nABC.StretchedF.gethist(	ans, 
									which( ans["sy2-sx2",]<=cu  &  ans["sy2-sx2",]>=-cu ), nbreaks)															
							cu<- 0.4
							out[,2]<- project.nABC.StretchedF.gethist(	ans, 
									which( ans["sy2-sx2",]<=cu  &  ans["sy2-sx2",]>=-cu ), nbreaks)
							cu<- 0.3
							out[,3]<- project.nABC.StretchedF.gethist(	ans, 
									which( ans["sy2-sx2",]<=cu  &  ans["sy2-sx2",]>=-cu ), nbreaks)
							cu<- 0.2
							out[,4]<- project.nABC.StretchedF.gethist(	ans, 
									which( ans["sy2-sx2",]<=cu  &  ans["sy2-sx2",]>=-cu ), nbreaks)															
							cu<- 0.1
							out[,5]<- project.nABC.StretchedF.gethist(	ans, 
									which( ans["sy2-sx2",]<=cu  &  ans["sy2-sx2",]>=-cu ), nbreaks)															
							
							cu<- 0.5
							#acc<- which( ans["log.sy2-log.sx2",]<=cu  &  ans["log.sy2-log.sx2",]>=-cu )
							#hist(ans["ysigma2",acc],breaks=100)
							out[,6]<- project.nABC.StretchedF.gethist(	ans, 
									which( ans["log.sy2-log.sx2",]<=cu  &  ans["log.sy2-log.sx2",]>=-cu ), nbreaks)
							cu<- 0.4
							out[,7]<- project.nABC.StretchedF.gethist(	ans, 
									which( ans["log.sy2-log.sx2",]<=cu  &  ans["log.sy2-log.sx2",]>=-cu ), nbreaks)
							cu<- 0.3							
							out[,8]<- project.nABC.StretchedF.gethist(	ans, 
									which( ans["log.sy2-log.sx2",]<=cu  &  ans["log.sy2-log.sx2",]>=-cu ), nbreaks)
							cu<- 0.2							
							out[,9]<- project.nABC.StretchedF.gethist(	ans, 
									which( ans["log.sy2-log.sx2",]<=cu  &  ans["log.sy2-log.sx2",]>=-cu ), nbreaks)
							cu<- 0.1							
							out[,10]<- project.nABC.StretchedF.gethist(	ans, 
									which( ans["log.sy2-log.sx2",]<=cu  &  ans["log.sy2-log.sx2",]>=-cu ), nbreaks)											
							out								
						})
				rownames(ans)<- c("me.o0.5","mo.o0.5","me.o0.4","mo.o0.4","me.o0.3","mo.o0.3","me.o0.2","mo.o0.2","me.o0.1","mo.o0.1","me.l0.5","mo.l0.5","me.l0.4","mo.l0.4","me.l0.3","mo.l0.3","me.l0.2","mo.l0.2","me.l0.1","mo.l0.1")
				
				f.name<- paste(dir.name,paste("nABC.naive.sigmainference_",N,".R",sep=''),sep='/')
				cat(paste("\nnABC.StretchedF.sigmainference save 'ans' to ",f.name))				
				save(ans,file=f.name)
			}
			
			ans.mean<- ans[seq(1,nrow(ans),2),]			
			ans.mode<- ans[seq(2,nrow(ans),2),]
			
			#print estimated mode for 'difference'
			ans<- ans.mode[c("mo.o0.5","mo.o0.3","mo.o0.1"), ]		
			nbreaks<- ceiling( nclass.Sturges(ans)*1 )+1
			breaks<- max( diff(c(1,max(c(ans,recursive=1))*1.1)) / nbreaks,	diff(c(min(c(ans,recursive=1))*0.9,1)) / nbreaks )		
			breaks<- c( rev(seq(from= 1-breaks/2, by=-breaks, length.out= nbreaks )), seq(from= 1+breaks/2, by=breaks, length.out= nbreaks ) )						
			m.h<- apply(ans,1,function(x){	hist(x, breaks=breaks, plot=0 )		})								
			xlim<- range(sapply(m.h,function(x){		range(x$breaks)		}))
			xlim<- c(0.6,1.45)
			ylim<- range(sapply(m.h,function(x){		range(c(0, x$counts))		}))		
			f.name<- paste(dir.name,"/nABC.StretchedF_naive_o_mo_",xn,".pdf",sep='')
			cat(paste("\nABC.StretchedF: write pdf to",f.name))
			pdf(f.name,version="1.4",width=5,height=6)
			par(mar=c(4,4,1,1))		
			plot(1,1,type='n',xlim=xlim,ylim=ylim,xlab=expression("estimated mode of "*sigma^2),ylab="n-ABC repetitions")
			cols<- c("red","cyan","blue")
			sapply(seq_along(m.h),function(i){		plot(m.h[[i]],col=my.fade.col(cols[i],0.4),border=my.fade.col(cols[i],0.7),add=TRUE, lty=1+i)		})		
			abline(v=1,lty=3)
			legend("topright",fill= c("transparent","transparent","red","transparent","cyan","transparent","blue"),legend=c(expression(S^2*'('*y*')'*-S^2*'('*x*')'),"",expression("["*c^'-'*","*c^'+'*"]=[-0.5,0.5]"),"",expression("["*c^'-'*","*c^'+'*"]=[-0.3,0.3]"),"",expression("["*c^'-'*","*c^'+'*"]=[-0.1,0.1]")),bty='n', border=NA)		
			dev.off()
			
			#print estimated mode for 'log.difference'
			ans<- ans.mode[c("mo.l0.5","mo.l0.3","mo.l0.1"), ]		
			nbreaks<- ceiling( nclass.Sturges(ans)*0.8 )+1
			breaks<- max( diff(c(1,max(c(ans,recursive=1))*1.1)) / nbreaks,	diff(c(min(c(ans,recursive=1))*0.9,1)) / nbreaks )		
			breaks<- c( rev(seq(from= 1-breaks/2, by=-breaks, length.out= nbreaks )), seq(from= 1+breaks/2, by=breaks, length.out= nbreaks ) )						
			m.h<- apply(ans,1,function(x){	hist(x, breaks=breaks, plot=0 )		})								
			xlim<- range(sapply(m.h,function(x){		range(x$breaks)		}))
			xlim<- c(0.6,1.45)
			ylim<- range(sapply(m.h,function(x){		range(c(0, x$counts))		}))		
			f.name<- paste(dir.name,"/nABC.StretchedF_naive_l_mo_",xn,".pdf",sep='')
			cat(paste("\nABC.StretchedF: write pdf to",f.name))
			pdf(f.name,version="1.4",width=5,height=6)
			par(mar=c(4,4,1,1))		
			plot(1,1,type='n',xlim=xlim,ylim=ylim,xlab=expression("estimated mode of "*sigma^2),ylab="n-ABC repetitions")
			cols<- c("red","cyan","blue")
			sapply(seq_along(m.h),function(i){		plot(m.h[[i]],col=my.fade.col(cols[i],0.4),border=my.fade.col(cols[i],0.7),add=TRUE, lty=1+i)		})		
			abline(v=1,lty=3)
			legend("topright",fill= c("transparent","transparent","red","transparent","cyan","transparent","blue"),legend=c(expression('log'*S^2*'('*y*')'*-'log'*S^2*'('*x*')'),"",expression("["*c^'-'*","*c^'+'*"]=[-0.5,0.5]"),"",expression("["*c^'-'*","*c^'+'*"]=[-0.3,0.3]"),"",expression("["*c^'-'*","*c^'+'*"]=[-0.1,0.1]")),bty='n', border=NA)		
			dev.off()
			
			#print estimated mean for 'log.difference'
			ans<- ans.mean[c("me.l0.5","me.l0.3","me.l0.1"), ]		
			nbreaks<- ceiling( nclass.Sturges(ans)*5 )+1
			breaks<- max( diff(c(1,max(c(ans,recursive=1))*1.1)) / nbreaks,	diff(c(min(c(ans,recursive=1))*0.9,1)) / nbreaks )		
			breaks<- c( rev(seq(from= 1-breaks/2, by=-breaks, length.out= nbreaks )), seq(from= 1+breaks/2, by=breaks, length.out= nbreaks ) )						
			m.h<- apply(ans,1,function(x){	hist(x, breaks=breaks, plot=0 )		})								
			xlim<- range(sapply(m.h,function(x){		range(x$breaks)		}))
			xlim<- c(0.88,1.1)
			ylim<- range(sapply(m.h,function(x){		range(c(0, x$counts))		}))		
			f.name<- paste(dir.name,"/nABC.StretchedF_naive_l_me_",xn,".pdf",sep='')
			cat(paste("\nABC.StretchedF: write pdf to",f.name))
			pdf(f.name,version="1.4",width=5,height=6)
			par(mar=c(4,4,1,1))		
			plot(1,1,type='n',xlim=xlim,ylim=ylim,xlab=expression("estimated mean of "*sigma^2),ylab="n-ABC repetitions")
			cols<- c("red","cyan","blue")
			sapply(seq_along(m.h),function(i){		plot(m.h[[i]],col=my.fade.col(cols[i],0.4),border=my.fade.col(cols[i],0.7),add=TRUE, lty=1+i)		})		
			abline(v=1,lty=3)
			legend("topleft",fill= c("transparent","transparent","red","transparent","cyan","transparent","blue"),legend=c(expression('log'*S^2*'('*y*')'*-'log'*S^2*'('*x*')'),"",expression("["*c^'-'*","*c^'+'*"]=[-0.5,0.5]"),"",expression("["*c^'-'*","*c^'+'*"]=[-0.3,0.3]"),"",expression("["*c^'-'*","*c^'+'*"]=[-0.1,0.1]")),bty='n', border=NA)		
			dev.off()
			
			#print estimated mean for 'difference'
			ans<- ans.mean[c("me.o0.5","me.o0.3","me.o0.1"), ]		
			nbreaks<- ceiling( nclass.Sturges(ans)*5 )+1
			breaks<- max( diff(c(1,max(c(ans,recursive=1))*1.1)) / nbreaks,	diff(c(min(c(ans,recursive=1))*0.9,1)) / nbreaks )		
			breaks<- c( rev(seq(from= 1-breaks/2, by=-breaks, length.out= nbreaks )), seq(from= 1+breaks/2, by=breaks, length.out= nbreaks ) )						
			m.h<- apply(ans,1,function(x){	hist(x, breaks=breaks, plot=0 )		})								
			xlim<- range(sapply(m.h,function(x){		range(x$breaks)		}))
			xlim<- c(0.88,1.1)
			ylim<- range(sapply(m.h,function(x){		range(c(0, x$counts))		}))
			f.name<- paste(dir.name,"/nABC.StretchedF_naive_o_me_",xn,".pdf",sep='')
			cat(paste("\nABC.StretchedF: write pdf to",f.name))
			pdf(f.name,version="1.4",width=5,height=6)
			par(mar=c(4,4,1,1))		
			plot(1,1,type='n',xlim=xlim,ylim=ylim,xlab=expression("estimated mean of "*sigma^2),ylab="n-ABC repetitions")
			cols<- c("red","cyan","blue")
			sapply(seq_along(m.h),function(i){		plot(m.h[[i]],col=my.fade.col(cols[i],0.4),border=my.fade.col(cols[i],0.7),add=TRUE, lty=1+i)		})		
			abline(v=1,lty=3)
			legend("topright",fill= c("transparent","transparent","red","transparent","cyan","transparent","blue"),legend=c(expression(S^2*'('*y*')'*-S^2*'('*x*')'),"",expression("["*c^'-'*","*c^'+'*"]=[-0.5,0.5]"),"",expression("["*c^'-'*","*c^'+'*"]=[-0.3,0.3]"),"",expression("["*c^'-'*","*c^'+'*"]=[-0.1,0.1]")),bty='n', border=NA)		
			dev.off()
			
			stop()
			
			ans<- ans.mode[c("u","u.naive"), ]		
			nbreaks<- ceiling( nclass.Sturges(ans)*0.6 )+1
			breaks<- max( diff(c(1,max(c(ans,recursive=1))*1.1)) / nbreaks,	diff(c(min(c(ans,recursive=1))*0.9,1)) / nbreaks )		
			breaks<- c( rev(seq(from= 1-breaks/2, by=-breaks, length.out= nbreaks )), seq(from= 1+breaks/2, by=breaks, length.out= nbreaks ) )						
			m.h<- apply(ans,1,function(x){	hist(x, breaks=breaks, plot=0 )		})								
			xlim<- range(sapply(m.h,function(x){		range(x$breaks)		}))
			xlim<- c(0.6,1.45)
			ylim<- range(sapply(m.h,function(x){		range(c(0, x$counts))		}))		
			cat(paste("\nABC.StretchedF: write pdf to",paste(dir.name,"/nABC.StretchedF_modethresh_",xn,".pdf",sep='')))
			pdf(paste(dir.name,"/nABC.StretchedF_modethresh_",xn,".pdf",sep=''),version="1.4",width=5,height=6)
			par(mar=c(4,4,1,1))		
			plot(1,1,type='n',xlim=xlim,ylim=ylim,xlab=expression("estimated mode of "*sigma^2),ylab="n-ABC repetitions")
			cols<- c("red","blue")
			sapply(seq_along(m.h),function(i){		plot(m.h[[i]],col=my.fade.col(cols[i],0.2),border=my.fade.col(cols[i],0.7),add=TRUE, lty=1+i)		})		
			abline(v=1,lty=3)
			legend("topright",fill= c("red","transparent","transparent","transparent","blue","transparent","transparent"),legend=c("symmetrized","tolerances",expression("["*tau^'-'*","*tau^'+'*"]=[1/2.77,2.77]"),"","rejection","region",expression("["*c^'-'*","*c^'+'*"]=[0.5,1.5]")),bty='n', border=NA)		
			dev.off()
			
		}
		stop()
	}
	if(1)
	{
		#require(multicore, lib.loc= LIB.LOC)
		require(msm, lib.loc= LIB.LOC)
		file<- paste(CODE.HOME,paste("lib",paste("libphylodyr",.Platform$dynlib.ext,sep=''),sep='/'),sep='')		
		dyn.load(file)
		#ABC inference from within (-tau,tau) when observed data is resimulated, resampled or fixed
		
		N<- 1e6		#typically 5e4
		M<- 2e3		#typically 500
		m<- NA		
		if(exists("argv"))
		{
			tmp<- na.omit(sapply(argv,function(arg)
							{	switch(substr(arg,2,2),
										m= return(as.numeric(substr(arg,3,nchar(arg)))),NA)	}))
			if(length(tmp)>0) m<- tmp[1]
		}
		
		xn<- yn<- 61		
		xmu<- ymu<- 0
		xsigma<- 1 
		#tau.u<- 1.5;	tau.l<- 0.700134618892036
		prior.u<- 4
		prior.l<- 0.2
		tau.u<- 2.77
		tau.l<- 0.27 #calibration for post mean: 0.29085995330929 
		alpha<- 0.01
		fade<- 0.5
		nproc<- 8
		nbreaks<- 75
		package.mkdir(DATA,"nABC.StretchedF")
		dir.name<- paste(DATA,"nABC.StretchedF",sep='/')				
		resume<- 1
		
		if(!is.na(m))
		{		
			f.name<- paste(dir.name,"/nnABC.StretchedF.sigmainference_",N,"_",m,".R",sep='')
			cat(paste("\nnABC.StretchedF: compute ",f.name))
			options(show.error.messages = FALSE, warn=1)		
			readAttempt<-try(suppressWarnings(load(f.name)))						
			options(show.error.messages = TRUE)						
			if(!resume || inherits(readAttempt, "try-error"))
			{
				ans.logu.naive<- 	project.nABC.StretchedF.fix.xsigma.fix.ylogtau.u(N,tau.l,tau.u,prior.l,prior.u,alpha,xn,xmu,xsigma,yn,ymu)
				ans.logu<- 			project.nABC.StretchedF.fix.xsigma.fix.ylogtau.u(N,1/tau.u,tau.u,prior.l,prior.u,alpha,xn,xmu,xsigma,yn,ymu)
				ans.u<- 		project.nABC.StretchedF.fix.xsigma.fix.ytau.u(N,1/tau.u,tau.u,prior.l,prior.u,alpha,xn,xmu,xsigma,yn,ymu)
				ans.u.naive<- 	project.nABC.StretchedF.fix.xsigma.fix.ytau.u(N,tau.l,tau.u,prior.l,prior.u,alpha,xn,xmu,xsigma,yn,ymu)
				
				cat(paste("\nnABC.StretchedF: save ",f.name))
				save(ans.logu,ans.logu.naive,ans.u,ans.u.naive,file=f.name)				
			}
			else
				cat(paste("\nnABC.MA: resumed ",f.name))
		}
		else
		{
			#load data
			cat(paste("\nnABC.StretchedF.sigmainference",dir.name))
			f.name<- paste(dir.name,paste("nABC.StretchedF.sigmainference_",N,".R",sep=''),sep='/')
			options(show.error.messages = FALSE, warn=1)		
			readAttempt<-try(suppressWarnings(load(f.name)))						
			options(show.error.messages = TRUE)		
			if(!resume || inherits(readAttempt, "try-error"))
			{
				cat(paste("\nnABC.StretchedF.sigmainference generate",paste(dir.name,paste("nABC.StretchedF.sigmainference_",N,".R",sep=''),sep='/')))
				f.name<- list.files(dir.name, pattern=paste("^nnABC.StretchedF.sigmainference",'',sep=''), full.names = TRUE)			
				ans<- sapply(seq_along(f.name),function(i)
						{
							cat(paste("\nload",f.name[i]))
							readAttempt<-try(suppressWarnings(load( f.name[i] )))
							if(inherits(readAttempt, "try-error"))	stop("error")
							out<- matrix(NA,2,7,dimnames=list(c("mean","mode"),c("logu","logu.ps","u","u.ps","u.h","u.naive","logu.naive")))
							out[,1]<- project.nABC.StretchedF.gethist(	ans.logu, 
									which( ans.logu["error",]<=ans.logu["cir",]  &  ans.logu["error",]>=ans.logu["cil",] ), nbreaks)
							out[,2]<- project.nABC.StretchedF.gethist(	ans.logu, 
									which( ans.u["ysigma2",]>=1/4 & ans.u["ysigma2",]<=4 & ans.logu["error",]<=ans.logu["cir",]  &  ans.logu["error",]>=ans.logu["cil",] ), nbreaks)															
							out[,3]<- project.nABC.StretchedF.gethist(	ans.u, 
									which( ans.u["error",]<=ans.u["cir",]  &  ans.u["error",]>=ans.u["cil",] ), nbreaks)
							out[,4]<- project.nABC.StretchedF.gethist(	ans.u, 
									which( ans.u["ysigma2",]>=1/4 & ans.u["ysigma2",]<=4 & ans.u["error",]<=ans.u["cir",]  &  ans.u["error",]>=ans.u["cil",] ), nbreaks)
							out[,5]<- project.nABC.StretchedF.gethist(	ans.u, 
									which( ans.u["ysigma2",]>=1/tau.u & ans.u["ysigma2",]<=4 & ans.u["error",]<=ans.u["cir",]  &  ans.u["error",]>=ans.u["cil",] ), nbreaks)																														
							out[,6]<- project.nABC.StretchedF.gethist(	ans.u.naive, 
									which( ans.u.naive["error",]<=ans.u.naive["cir",]  &  ans.u.naive["error",]>=ans.u.naive["cil",] ), nbreaks)
							out[,7]<- project.nABC.StretchedF.gethist(	ans.logu.naive, 
									which( ans.logu.naive["error",]<=ans.logu.naive["cir",]  &  ans.logu.naive["error",]>=ans.logu.naive["cil",] ), nbreaks)
							out									
						})
				cat(paste("\nnABC.StretchedF.sigmainference save 'ans' to ",paste(dir.name,paste("nABC.StretchedF.sigmainference_",N,".R",sep=''),sep='/')))
				save(ans,file=paste(dir.name,paste("nABC.StretchedF.sigmainference_",N,".R",sep=''),sep='/'))
			}
		}
		
		#print(ans)
		ans.mean<- ans[seq(1,nrow(ans),2),]
		rownames(ans.mean)<- c("logu","logu.ps","u","u.ps","u.h","u.naive","logu.naive")
		ans.mode<- ans[seq(2,nrow(ans),2),]
		rownames(ans.mode)<- c("logu","logu.ps","u","u.ps","u.h","u.naive","logu.naive")
		
		
		ans<- ans.mode[c("u","u.naive"), ]		
		nbreaks<- ceiling( nclass.Sturges(ans)*0.6 )+1
		breaks<- max( diff(c(1,max(c(ans,recursive=1))*1.1)) / nbreaks,	diff(c(min(c(ans,recursive=1))*0.9,1)) / nbreaks )		
		breaks<- c( rev(seq(from= 1-breaks/2, by=-breaks, length.out= nbreaks )), seq(from= 1+breaks/2, by=breaks, length.out= nbreaks ) )						
		m.h<- apply(ans,1,function(x){	hist(x, breaks=breaks, plot=0 )		})								
		xlim<- range(sapply(m.h,function(x){		range(x$breaks)		}))
		xlim<- c(0.6,1.45)
		ylim<- range(sapply(m.h,function(x){		range(c(0, x$counts))		}))		
		cat(paste("\nABC.StretchedF: write pdf to",paste(dir.name,"/nABC.StretchedF_modethresh_",xn,".pdf",sep='')))
		pdf(paste(dir.name,"/nABC.StretchedF_modethresh_",xn,".pdf",sep=''),version="1.4",width=5,height=6)
		par(mar=c(4,4,1,1))		
		plot(1,1,type='n',xlim=xlim,ylim=ylim,xlab=expression("estimated mode of "*sigma^2),ylab="n-ABC repetitions")
		cols<- c("black","gray60")
		sapply(seq_along(m.h),function(i){		plot(m.h[[i]],col=my.fade.col(cols[i],0.5),border=NA,add=TRUE)		})		
		abline(v=1,lty=2)
		legend("topright",fill= c("gray40","transparent","transparent","transparent","gray70","transparent","transparent"),legend=c("symmetrized","tolerances",expression("["*tau^'-'*","*tau^'+'*"]=[1/2.77,2.77]"),"","rejection","region",expression("["*c^'-'*","*c^'+'*"]=[0.5,1.5]")),bty='n', border=NA)		
		dev.off()
		
		
		ans<- ans.mode[c("u.ps","logu"), ]		
		nbreaks<- ceiling( nclass.Sturges(ans)*0.5 )+1
		breaks<- max( diff(c(1,max(c(ans,recursive=1))*1.1)) / nbreaks,	diff(c(min(c(ans,recursive=1))*0.9,1)) / nbreaks )		
		breaks<- c( rev(seq(from= 1-breaks/2, by=-breaks, length.out= nbreaks )), seq(from= 1+breaks/2, by=breaks, length.out= nbreaks ) )						
		m.h<- apply(ans,1,function(x){	hist(x, breaks=breaks, plot=0 )		})								
		xlim<- range(sapply(m.h,function(x){		range(x$breaks)		}))
		xlim<- c(0.6,1.45)
		ylim<- range(sapply(m.h,function(x){		range(c(0, x$counts))		}))		
		cat(paste("\nABC.StretchedF: write pdf to",paste(dir.name,"/nABC.StretchedF_modepriors_",xn,".pdf",sep='')))
		pdf(paste(dir.name,"/nABC.StretchedF_modepriors_",xn,".pdf",sep=''),version="1.4",width=5,height=6)
		par(mar=c(4,4,1,1))		
		plot(1,1,type='n',xlim=xlim,ylim=ylim,xlab=expression("estimated mode of "*sigma^2),ylab="n-ABC repetitions")
		sapply(seq_along(m.h),function(i){		plot(m.h[[i]],col=my.fade.col(cols[i],0.5), add=TRUE, border=NA)		})		
		abline(v=1,lty=2)			
		legend("topright",	fill= c("transparent","transparent","transparent","gray40","transparent","gray70"),
				c("symmetrized","tolerances and",'',expression(sigma^2*" ~ Unif[0.2,4]"),"",expression(1/sigma^2*" ~ Unif[0.2,4]")),
				bty='n', border=NA)
		dev.off()
		
		stop()		
	}	
	if(0)
	{
		require(multicore, lib.loc= LIB.LOC)
		require(msm, lib.loc= LIB.LOC)
		#ABC inference from within (-tau,tau) when observed data is resimulated, resampled or fixed
		
		N<- 5e4		#typically 5e4
		M<- 2e3	#typically 500
		xn<- yn<- 60		
		xmu<- ymu<- 0
		xsigma<- 1 
		tau.u<- 2.77
		tau.l<- 0.27
		alpha<- 0.05
		fade<- 0.5
		nproc<- 8
		f.name<- paste("/Users/Oliver/duke/2012_frequencyABC/sim.data/nABC.StretchedF.sigma.log.inference.ansh_",xn,".R",sep='')
		resume<- 1
		
		if(resume)
		{
			options(show.error.messages = FALSE)		
			readAttempt<-try(suppressWarnings(load(f.name)))
			if(!inherits(readAttempt, "try-error"))			cat(paste("\nnABC.StretchedF: resume ",f.name))
			options(show.error.messages = TRUE)
		}
		if(!resume || inherits(readAttempt, "try-error"))
		{
			package.mkdir(DATA,"nABC.StretchedF.sigmainference")
			dir.name<- paste(DATA,"nABC.StretchedF.sigmainference",sep='/')
			cat(paste("\nnABC.StretchedF:\ndir.name is ",dir.name,"\nN\t",N,"\nM\t",M,"\nxn\t",xn,"\ntau.l\t",tau.l,"\ntau.u\t",tau.u))
			
			#compute sigma histograms with symmetric tau and sample ysigma under log-uniform --> should give log.mean at 0
			ans.h.fix.xsigma.log<- mclapply(	seq_len(M), function(p)
					#lapply(	seq_len(M), function(p)
					{
						ans<- project.nABC.StretchedF.fix.xsigma.fix.ylogtau.u(N,1/tau.u,tau.u,1/tau.u,tau.u,alpha,xn,xmu,xsigma,yn,ymu)
						acc<- which( ans["error",]<=ans["cir",]  &  ans["error",]>=ans["cil",] )
						print(c(ncol(ans),length(acc)/ncol(ans),ans["cil",1],ans["cir",1]))
						#convert to log scale
						ans["ysigma",]<- log(ans["ysigma",])
						#compute break points sth 1 is in the middle					
						breaks<- max( diff(c(0,max(ans["ysigma",acc])*1.1)) / 20,	diff(c(min(ans["ysigma",acc])*1.1,0)) / 20 )						
						breaks<- c( rev(seq(from= -breaks/2, by=-breaks, length.out= 20 )), seq(from= breaks/2, by=breaks, length.out= 20 ) )
						ans.h<- hist(ans["ysigma",acc],breaks=breaks, plot= 0)
						ans.h[["log.mean"]]<- mean(ans["ysigma",acc])
						ans.h[["exp.mean"]]<- mean(exp(ans["ysigma",acc]))						
						ans.h
					},mc.preschedule= 1, mc.cores= nproc )
			
			#compute sigma histograms with symmetric tau and sample ysigma under lognormal --> should give log.mean at 0
			ans.h.fix.xsigma.lognorm<- mclapply(	seq_len(M), function(p)									
					{
						ans<- project.nABC.StretchedF.fix.xsigma.fix.ynormtau.u(N,1/tau.u,tau.u,1/tau.u,tau.u,alpha,xn,xmu,xsigma,yn,ymu)
						acc<- which( ans["error",]<=ans["cir",]  &  ans["error",]>=ans["cil",] )						
						print(c(ncol(ans),length(acc)/ncol(ans)))	
						ans["ysigma",]<- log(ans["ysigma",])
						#compute break points sth 1 is in the middle
						breaks<- max( diff(c(0,max(ans["ysigma",acc])*1.1)) / 20,	diff(c(min(ans["ysigma",acc])*1.1,0)) / 20 )						
						breaks<- c( rev(seq(from= -breaks/2, by=-breaks, length.out= 20 )), seq(from= breaks/2, by=breaks, length.out= 20 ) )
						ans.h<- hist(ans["ysigma",acc],breaks=breaks, plot= 0)
						ans.h[["log.mean"]]<- mean(ans["ysigma",acc])
						ans.h[["exp.mean"]]<- mean(exp(ans["ysigma",acc]))
						ans.h
					},
					mc.preschedule= 1, mc.cores= nproc )
			
			
			#compute sigma histograms with symmetric tau and sample ysigma under uniform --> should give log.mean larger than 0
			ans.h.fix.xsigma.unif<- mclapply(	seq_len(M), function(p)									
					{
						ans<- project.nABC.StretchedF.fix.xsigma.fix.ytau.u(N,1/tau.u,tau.u,1/tau.u,tau.u,alpha,xn,xmu,xsigma,yn,ymu)
						acc<- which( ans["error",]<=ans["cir",]  &  ans["error",]>=ans["cil",] )						
						print(c(ncol(ans),length(acc)/ncol(ans)))	
						ans["ysigma",]<- log(ans["ysigma",])
						#compute break points sth 1 is in the middle
						breaks<- max( diff(c(0,max(ans["ysigma",acc])*1.1)) / 20,	diff(c(min(ans["ysigma",acc])*1.1,0)) / 20 )						
						breaks<- c( rev(seq(from= -breaks/2, by=-breaks, length.out= 20 )), seq(from= breaks/2, by=breaks, length.out= 20 ) )						
						ans.h<- hist(ans["ysigma",acc],breaks=breaks, plot= 0)
						ans.h[["log.mean"]]<- mean(ans["ysigma",acc])
						ans.h[["exp.mean"]]<- mean(exp(ans["ysigma",acc]))
						ans.h
					},
					mc.preschedule= 1, mc.cores= nproc )
			
			
			
			ans.h<- list(3)
			ans.h[[1]]<- ans.h.fix.xsigma.unif
			ans.h[[2]]<- ans.h.fix.xsigma.log
			ans.h[[3]]<- ans.h.fix.xsigma.lognorm
			cat(paste("\nnABC.StretchedF: save file ",paste(dir.name,paste("nABC.StretchedF.sigma.log.inference.ansh_",xn,".R",sep=''),sep='/')))
			save(ans.h,file=paste(dir.name,paste("nABC.StretchedF.sigma.log.inference.ansh_",xn,".R",sep=''),sep='/'))									
		}
		#evaluate log.means
		means<- lapply(seq_along(ans.h),function(i)
				{
					sapply(ans.h[[i]],function(x){	x[["log.mean"]]		})					
				})
		nbreaks<- ceiling( nclass.Sturges(means[[1]])*2 )+1
		breaks<- max( diff(c(0,max(c(means,recursive=1))*1.1)) / nbreaks,	diff(c(min(c(means,recursive=1))*1.1,0)) / nbreaks )		
		breaks<- c( rev(seq(from= -breaks/2, by=-breaks, length.out= nbreaks )), seq(from= breaks/2, by=breaks, length.out= nbreaks ) )						
		m.h<- lapply(means,function(x){	hist(x, breaks=breaks, plot=0 )		})		
		
		xlim<- range(sapply(m.h,function(x){		range(x$breaks)		}))
		ylim<- range(sapply(m.h,function(x){		range(c(0, x$counts))		}))
		
		cat(paste("\nABC.StretchedF: write pdf to",paste(dir.name,"/nABC.StretchedF_meanslog_",xn,".pdf",sep='')))
		pdf(paste(dir.name,"/nABC.StretchedF_meanslog_",xn,".pdf",sep=''),version="1.4",width=5,height=6)
		par(mar=c(4,4,1,1))		
		plot(1,1,type='n',xlim=xlim,ylim=ylim,xlab=expression("estimated mean of log "*sigma^2),ylab="n-ABC repitions")
		plot(m.h[[1]],col=my.fade.col("blue",0.3),border=NA,add=TRUE)
		plot(m.h[[2]],col=my.fade.col("red",0.3),border=NA,add=TRUE)
		plot(m.h[[3]],col=my.fade.col("green",0.3),border=NA,add=TRUE)
		#legend("topleft",fill= c("blue","transparent","transparent","transparent","red","transparent","transparent"),legend=c("symmetrized","tolerances",expression("["*tau^'-'*","*tau^'+'*"]=[1/2.77,2.77]"),"","rejection","region",expression("["*c^'-'*","*c^'+'*"]=[0.5,1.5]")),bty='n', border=NA)
		legend("topleft",fill= c("blue","red","green"),legend=c(expression(rho*" ~ Unif"),expression("log "*rho*" ~ Unif"),expression("log "*rho*" ~ Norm")),bty='n')
		dev.off()
		
		stop()
	}		
	if(0)
	{
		require(multicore, lib.loc= LIB.LOC)
		#ABC inference from within (-tau,tau) when observed data is resimulated, resampled or fixed
		
		N<- 1e6		#typically 5e4
		M<- 2e3		#typically 500
		xn<- yn<- 61		
		xmu<- ymu<- 0
		xsigma<- 1 
		tau.u<- 2.77
		tau.l<- 0.27
		alpha<- 0.01
		fade<- 0.5
		nbreaks<- 75
		nproc<- 8
		dir.name<- paste(DATA,"nABC.StretchedF.sigmainference",sep='/')
		f.name<- paste("/Users/Oliver/duke/2012_frequencyABC/sim.data/nABC.StretchedF.sigma.asym.inference.ansh_",xn,".R",sep='')
		resume<- 1
		
		if(resume)
		{
			f.name<- paste(dir.name,paste("nABC.StretchedF.sigmainference.ansh_",xn,".R",sep=''),sep='/')
			options(show.error.messages = FALSE,warn=1)					
			readAttempt<-try(suppressWarnings(load(f.name)))		
			if(!inherits(readAttempt, "try-error"))			cat(paste("\nnABC.StretchedF: resume ",f.name))
			cat(paste("\nnABC.StretchedF: resume ",f.name))
			options(show.error.messages = TRUE)
		}
		if(!resume || inherits(readAttempt, "try-error"))
		{
			package.mkdir(DATA,"nABC.StretchedF.sigmainference")			
			cat(paste("\nnABC.StretchedF:\ndir.name is ",dir.name,"\nN\t",N,"\nM\t",M,"\nxn\t",xn,"\ntau.l\t",tau.l,"\ntau.u\t",tau.u))
			
			#compute sigma histograms with fix.sigma and asymmetric tau --> should give bias
			if(0)
			{
				ans.h.fix.xsigma.asy<- mclapply(	seq_len(M), function(p)
						#lapply(	seq_len(M), function(p)
						{
							tmp.fname<- paste("nABC.StretchedF.sigmainference.ansh_",xn,"_asy_",p,".R",sep='')
							options(show.error.messages = FALSE, warn=1)		
							readAttempt<-try(suppressWarnings(load(paste(dir.name,tmp.fname,sep='/'))))						
							options(show.error.messages = TRUE)						
							if(!inherits(readAttempt, "try-error"))
							{
								cat(paste("\nnABC.StretchedF: resumed ",tmp.fname))
							}			
							else
							{
								ans<- project.nABC.StretchedF.fix.xsigma.fix.ytau.u(N,tau.l,tau.u,tau.l,tau.u,alpha,xn,xmu,xsigma,yn,ymu)
								acc<- which( ans["error",]<=ans["cir",]  &  ans["error",]>=ans["cil",] )
								print(c(ncol(ans),length(acc)/ncol(ans)))						
								#compute break points sth 1 is in the middle
								breaks<- max( diff(c(1,max(ans["ysigma2",acc])*1.1)) / nbreaks,	diff(c(min(ans["ysigma2",acc])*0.9,1)) / nbreaks )		
								breaks<- c( rev(seq(from= 1-breaks/2, by=-breaks, length.out= nbreaks )), seq(from= 1+breaks/2, by=breaks, length.out= nbreaks ) )				
								ans.h<- hist(ans["ysigma2",acc],breaks=breaks, plot= 0)
								
								cat(paste("\nnABC.StretchedF: save file ",paste(dir.name,tmp.fname,sep='/')))
								save(ans.h,file=paste(dir.name,tmp.fname,sep='/'))	
							}										
							ans.h
						},mc.preschedule= 0, mc.cores= nproc )
			}
			#compute sigma histograms with fix.sigma and symmetrized tau
			if(1)
			{
				ans.h.fix.xsigma.sym<- mclapply(	seq_len(M), function(p)
						#lapply(	seq_len(M), function(p)
						{
							tmp.fname<- paste("nABC.StretchedF.sigmainference.ansh_",xn,"_sym_",p,".R",sep='')
							options(show.error.messages = FALSE)
							readAttempt<-try(suppressWarnings(load(paste(dir.name,tmp.fname,sep='/'))))
							options(show.error.messages = TRUE)						
							if(!inherits(readAttempt, "try-error"))
							{
								cat(paste("\nnABC.StretchedF: resume ",tmp.fname))
							}			
							else
							{							
								ans<- project.nABC.StretchedF.fix.xsigma.fix.ytau.u(N,1/tau.u,tau.u,1/tau.u,tau.u,alpha,xn,xmu,xsigma,yn,ymu)
								acc<- which( ans["error",]<=ans["cir",]  &  ans["error",]>=ans["cil",] )
								
								print(c(ncol(ans),length(acc)/ncol(ans)))	
								
								#compute break points sth 1 is in the middle
								breaks<- max( diff(c(1,max(ans["ysigma2",acc])*1.1)) / nbreaks,	diff(c(min(ans["ysigma2",acc])*0.9,1)) / nbreaks )		
								breaks<- c( rev(seq(from= 1-breaks/2, by=-breaks, length.out= nbreaks )), seq(from= 1+breaks/2, by=breaks, length.out= nbreaks ) )				
								ans.h<- hist(ans["ysigma2",acc],breaks=breaks, plot= 0)				
								cat(paste("\nnABC.StretchedF: save file ",paste(dir.name,tmp.fname,sep='/')))
								save(ans.h,file=paste(dir.name,tmp.fname,sep='/'))	
							}
							ans.h
						},mc.preschedule= 0, mc.cores= nproc )
			}	
			ans.h<- list(2)
			ans.h[[1]]<- ans.h.fix.xsigma.sym
			ans.h[[2]]<- ans.h.fix.xsigma.asy
			cat(paste("\nnABC.StretchedF: save file ",paste(dir.name,paste("nABC.StretchedF.sigmainference.ansh_",xn,".R",sep=''),sep='/')))
			save(ans.h,file=paste(dir.name,paste("nABC.StretchedF.sigmainference.ansh_",xn,".R",sep=''),sep='/'))									
		}		
		#evaluate modes
		modes<- lapply(seq_along(ans.h),function(i)
				{
					sapply(ans.h[[i]],function(x){	x[["mids"]][which.max( x[["counts"]] )]		})					
				})
		nbreaks<- ceiling( nclass.Sturges(modes[[1]])*0.7 )+1
		breaks<- max( diff(c(1,max(c(modes,recursive=1))*1.1)) / nbreaks,	diff(c(min(c(modes,recursive=1))*0.9,1)) / nbreaks )		
		breaks<- c( rev(seq(from= 1-breaks/2, by=-breaks, length.out= nbreaks )), seq(from= 1+breaks/2, by=breaks, length.out= nbreaks ) )						
		m.h<- lapply(modes,function(x){	hist(x, breaks=breaks, plot=0 )		})		
		
		xlim<- range(sapply(m.h,function(x){		range(x$breaks)		}))
		ylim<- range(sapply(m.h,function(x){		range(c(0, x$counts))		}))
		
		cat(paste("\nABC.StretchedF: write pdf to",paste(dir.name,"/nABC.StretchedF_modesasy_",xn,".pdf",sep='')))
		pdf(paste(dir.name,"/nABC.StretchedF_modesasy_",xn,".pdf",sep=''),version="1.4",width=5,height=6)
		par(mar=c(4,4,1,1))
		plot(1,1,type='n',xlim=xlim,ylim=ylim,xlab=expression("estimated mode of "*sigma^2),ylab="n-ABC repitions")
		abline(v=1,lty=3)
		plot(m.h[[1]],col=my.fade.col("blue",0.3),border=NA,add=TRUE)
		plot(m.h[[2]],col=my.fade.col("red",0.3),border=NA,add=TRUE)		
		legend("topright",fill= c("blue","transparent","transparent","transparent","red","transparent","transparent"),legend=c("symmetrized","tolerances",expression("["*tau^'-'*","*tau^'+'*"]=[1/2.77,2.77]"),"","rejection","region",expression("["*c^'-'*","*c^'+'*"]=[0.5,1.5]")),bty='n', border=NA)
		dev.off()
		
		stop()
	}	
	if(0)
	{
		require(multicore, lib.loc= LIB.LOC)
		#ABC inference from within (-tau,tau) when observed data is resimulated, resampled or fixed
		args<- "fstretch/1.75/0.01"
		N<- 1e3
		M<- 500
		xn<- yn<- 100		
		xmu<- ymu<- 0
		xsigma<- 1 
		tau<- 1.75
		alpha<- 0.01
		fade<- 0.5
		nproc<- 8
		f.name<- paste("/Users/Oliver/duke/2012_frequencyABC/sim.data/nABC.StretchedF.sigmainference.ansh_",xn,".R",sep='')
		resume<- 1
		
		if(resume)
		{
			options(show.error.messages = FALSE)		
			readAttempt<-try(suppressWarnings(load(f.name)))
			if(!inherits(readAttempt, "try-error"))			cat(paste("\nnABC.StretchedF: resume ",f.name))
			options(show.error.messages = TRUE)
		}
		if(!resume || inherits(readAttempt, "try-error"))
		{
			package.mkdir(DATA,"nABC.StretchedF.sigmainference")
			dir.name<- paste(DATA,"nABC.StretchedF.sigmainference",sep='/')
			cat(paste("\nnABC.StretchedF: dir.name is ",dir.name))
			
			x<- rnorm(xn,xmu,xsigma)
			
			#compute sigma histograms with fix.x			
			ans.h.fix.x<- mclapply(	seq_len(M), function(p)
					{
						ans<- project.nABC.StretchedF.fix.x.fix.ytau.u(N,tau,tau,alpha,x,yn,ymu)
						acc<- which( ans["error",]<=ans["cir",]  &  ans["error",]>=ans["cil",] )
						print(length(acc)/ncol(ans))	
						#compute break points sth 1 is in the middle
						breaks<- max( diff(c(1,max(ans["ysigma",acc])*1.1)) / 20,	diff(c(min(ans["ysigma",acc])*0.9,1)) / 20 )		
						breaks<- c( rev(seq(from= 1-breaks/2, by=-breaks, length.out= 20 )), seq(from= 1+breaks/2, by=breaks, length.out= 20 ) )				
						ans.h<- hist(ans["ysigma",acc],breaks=breaks, plot= 0)
						ans.h
					},
					mc.preschedule= 1, mc.cores= nproc )			
			
			#compute sigma histograms with resample.x
			ans.h.resample.x<- mclapply(	seq_len(M), function(p)
					{
						ans<- project.nABC.StretchedF.resample.x.fix.ytau.u(N,tau,tau,alpha,x,yn,ymu)
						acc<- which( ans["error",]<=ans["cir",]  &  ans["error",]>=ans["cil",] )
						print(length(acc)/ncol(ans))	
						
						#compute break points sth 1 is in the middle
						breaks<- max( diff(c(1,max(ans["ysigma",acc])*1.1)) / 20,	diff(c(min(ans["ysigma",acc])*0.9,1)) / 20 )		
						breaks<- c( rev(seq(from= 1-breaks/2, by=-breaks, length.out= 20 )), seq(from= 1+breaks/2, by=breaks, length.out= 20 ) )				
						ans.h<- hist(ans["ysigma",acc],breaks=breaks, plot= 0)
						ans.h
					},
					mc.preschedule= 1, mc.cores= nproc )
			
			#compute sigma histograms with fix.sigma
			ans.h.fix.xsigma<- mclapply(	seq_len(M), function(p)
					{
						ans<- project.nABC.StretchedF.fix.xsigma.fix.ytau.u(N,1/tau,tau,1/tau,tau,alpha,xn,xmu,xsigma,yn,ymu)
						acc<- which( ans["error",]<=ans["cir",]  &  ans["error",]>=ans["cil",] )
						print(length(acc)/ncol(ans))	
						
						#compute break points sth 1 is in the middle
						breaks<- max( diff(c(1,max(ans["ysigma",acc])*1.1)) / 20,	diff(c(min(ans["ysigma",acc])*0.9,1)) / 20 )		
						breaks<- c( rev(seq(from= 1-breaks/2, by=-breaks, length.out= 20 )), seq(from= 1+breaks/2, by=breaks, length.out= 20 ) )				
						ans.h<- hist(ans["ysigma",acc],breaks=breaks, plot= 0)
						ans.h
					},
					mc.preschedule= 1, mc.cores= nproc )
			
			ans.h<- list(3)
			ans.h[[1]]<- ans.h.fix.x
			ans.h[[2]]<- ans.h.fix.xsigma
			ans.h[[3]]<- ans.h.resample.x
			cat(paste("\nnABC.StretchedF: save file ",paste(dir.name,paste("nABC.StretchedF.sigmainference.ansh_",xn,".R",sep=''),sep='/')))
			save(ans.h,file=paste(dir.name,paste("nABC.StretchedF.sigmainference.ansh_",xn,".R",sep=''),sep='/'))									
		}
		
		#evaluate modes
		modes<- lapply(seq_along(ans.h),function(i)
				{
					sapply(ans.h[[i]],function(x){	x[["mids"]][which.max( x[["counts"]] )]		})					
				})
		nbreaks<- ceiling( nclass.Sturges(modes[[1]])*1.5 )+1
		breaks<- max( diff(c(1,max(c(modes,recursive=1))*1.1)) / nbreaks,	diff(c(min(c(modes,recursive=1))*0.9,1)) / nbreaks )		
		breaks<- c( rev(seq(from= 1-breaks/2, by=-breaks, length.out= nbreaks )), seq(from= 1+breaks/2, by=breaks, length.out= nbreaks ) )						
		m.h<- lapply(modes,function(x){	hist(x, breaks=breaks, plot=0 )		})		
		
		xlim<- range(sapply(m.h,function(x){		range(x$breaks)		}))
		ylim<- range(sapply(m.h,function(x){		range(c(0, x$counts))		}))
		
		pdf(paste(dir.name,"/nABC.StretchedF_modes_",xn,".pdf",sep=''),version="1.4",width=5,height=7)		
		plot(1,1,type='n',xlim=xlim,ylim=ylim,xlab="mode",ylab="counts")
		plot(m.h[[1]],col=my.fade.col("blue",0.3),border=NA,add=TRUE)
		plot(m.h[[2]],col=my.fade.col("red",0.3),border=NA,add=TRUE)
		plot(m.h[[3]],col=my.fade.col("green",0.3),border=NA,add=TRUE)
		legend("topleft",fill= c("blue","red","green"),legend=c("fixed obs","resimulated obs","resampled obs"),bty='n')
		dev.off()
	}
	if(0)
	{
		#compare qqplots when observed data is resimulated, resampled or fixed
		args<- "fstretch/1.75/0.01"
		N<- 1e3
		xn<- yn<- 100		
		xmu<- ymu<- 0
		xsigma<- 1 
		ysigma<- 1
		fade<- 0.5
		pdf(paste(dir.name,"/nABC.StretchedF_pvals.pdf",sep=''),version="1.4",width=5,height=5)
		for(i in 1:50)
		{
			ans<- project.nABC.StretchedF.fix.x.fix.ysigma(N,args,xn,xmu,xsigma,yn,ymu,ysigma)
			acc<- which( ans["error",]<=ans["cir",]  &  ans["error",]>=ans["cil",] )
			print(length(acc)/ncol(ans))		
			project.nABC.StretchedF.pval.qq(ans["pval",acc], ifelse(i==1,0,1), pch= 19, col=my.fade.col("blue", fade)) 
			
			ans<- project.nABC.StretchedF.resample.x.fix.ysigma(N,args,xn,xmu,xsigma,yn,ymu,ysigma)
			acc<- which( ans["error",]<=ans["cir",]  &  ans["error",]>=ans["cil",] )
			print(length(acc)/ncol(ans))		
			project.nABC.StretchedF.pval.qq(ans["pval",acc], 1, pch= 19, col=my.fade.col("green", fade)) 
			
			
			ans<- project.nABC.StretchedF.fix.xsigma.fix.ysigma(N,args,xn,xmu,xsigma,yn,ymu,ysigma)
			acc<- which( ans["error",]<=ans["cir",]  &  ans["error",]>=ans["cil",] )
			print(length(acc)/ncol(ans))		
			project.nABC.StretchedF.pval.qq(ans["pval",acc], 1, pch= 19, col=my.fade.col("red", fade))
			if(i==50)
			{
				abline(a=0, b=1, lty= 2)
				legend("topleft",fill= c("red","green","blue"),legend= c("with resimulation","with resampling","no resampling"), bty='n')
			}
		}
		dev.off()
	}
	if(0)		#plot mode for different choices of CL, CU
	{
		n<- m<- 60
		cu<- 1.5; 		cl<- 0.5				#typical ABC thresholds
		df<- n-1
		alpha<- 0.01
		#choose correct thresholds that are asymmetric on the normal scale (but symmetric on the log scale)
		rho.ok<- seq(1/cu,cu,by=0.001)
		pw.ok<- project.nABC.StretchedF.pow.sym(rho.ok,df, cu)
		#cat(paste("\npower in blue:",project.nABC.StretchedF.pow.to.one(rho.ok,df, cu)))
		
		rho.ok2<- seq(cl,1/cl,by=0.001)
		pw.ok2<- project.nABC.StretchedF.pow.sym(rho.ok2,df, 1/cl)
		#cat(paste("\npower in red:",project.nABC.StretchedF.pow.to.one(rho.ok2,df, cu)))
		
		#choose incorrect thresholds that are symmetric on the normal scale and max power off 1
		rho.sy<- seq(cl,cu,by=0.001)
		pw.sy<- project.nABC.StretchedF.pow.sym(rho.sy,df, cu, cl)
		
		#	pdf(paste(dir.name,"/nABC.StretchedF_poweroff_symmetrictol.pdf",sep=''),version="1.4",width=5,height=3)
		par(mar=c(5,4,0.5,0.5))
		plot(1,1,xlim= c(0,3),ylim= range(c(pw.sy,pw.ok,pw.ok2)), type='n', xlab= expression( sigma^2== sigma[0]^2 *" "* rho), ylab="power")
		#lines(rho.sy,pw.sy, col="blue")
		polygon( c(rho.sy,cu,cl), c(pw.sy,0,0), border= NA, col= "gray80")
		lines(rho.ok,pw.ok, col="blue", lty= 1.5)
		lines(rho.ok2,pw.ok2, col="red", lty= 1.5)
		#legend<- round(c(cl, 1/cu, cl),digits=2)
		#legend<- bquote(.(parse(text=paste(expression('['*c^'-'*", "*c^'+'*']'),"==~",legend,sep=""))))		
		#legend<- bquote(.(parse(text=paste(expression('['*c^'-'*", "*c^'+'*"]= [0.5,1.5]"),"==~",legend,sep=""))))
		
		legend<- expression('['*c^'-'*", "*c^'+'*"]= [0.5,1.5]", '['*c^'-'*", "*c^'+'*"]= [1/1.5,1.5]",'['*c^'-'*", "*c^'+'*"]= [0.5,1/0.5]")				
		legend("topright",fill=c("gray70","blue","red"),legend=legend, bty= 'n')
		abline(v=1, lty= 2)
		#	dev.off()
		stop()
	}
	if(0)	#illustrate scaled F test
	{
		xn<- yn<- 60		
		tau.low<- 1/2.77; tau.up<- 2.77
		df<- xn-1
		alpha<- 0.01
		verbose<- 1
		
		rej<- .Call("abcScaledF",	c(xn-1,yn-1,tau.low,tau.up,alpha,1e-10,100,0.05)	)
		xu<- seq(1/1.5,3,0.001)
		xl<- seq(0,1.5,0.001)
		yu<- df(xu/tau.up,xn-1,yn-1)
		yl<- df(xl/tau.low,xn-1,yn-1)
		
		pdf(paste(dir.name,"/nABC.StretchedF_example.pdf",sep=''),version="1.4",width=5,height=5)
		par(mar=c(5,4,0.5,0.5))
		plot(1,1,type='n',xlim=range(c(xu,xl)),ylim=range(c(yu,yl,-0.02)),ylab='',xlab="T")		
		cl<- seq(rej[1],1,0.001)
		cu<- seq(1,rej[2],0.001)
		polygon( c(cl,cu,rej[2:1]), c(df(cl/tau.low,xn-1,yn-1),df(cu/tau.up,xn-1,yn-1),0,0), col="gray30", border=NA)
		polygon( c(tau.low,tau.up,tau.up,tau.low), c(-0.03,-0.03,-1,-1), col= "gray80", border=NA )
		lines(xu,yu,lty=1, col="blue")		
		lines(xl,yl,lty=1, col= "red")
		abline(v=1,lty=2)
		
		legend<- expression( f[n-1*","*n-1](T/tau^'-'), f[n-1*","*n-1](T/tau^'+'),'['*c^'-'*", "*c^'+'*"]= ["*1.5^{-1}*",1.5]",'['*tau^'-'*", "*tau^'+'*"]= ["*2.77^{-1}*",2.77]","")						
		legend(x=1.55,y=0.4,fill=c("red","blue","gray30","gray80","transparent"),legend=legend, bty= 'n', border=NA)
		dev.off()
		stop()
		#my.fade.col("green", fade)
	}
	if(0)	#plot mean for different choices of CL, CU, assuming rho~ U(rho.l,rho.u) 	(which does not make much sense)
	{
		xn<- yn<- 100
		cu<- 1.5; 		cl<- 0.5				#typical ABC thresholds
		tau.low<- 0.271; tau.up<- 2.4
		df<- xn-1
		alpha<- 0.01
		verbose<- 1
		
		#tau.low<- 1/2.77; tau.up<- 2.77 	--> 6.663987e-01 1.500603e+00 
		#tau.low<- 0.271; tau.up<- 2.77		--> 5.002459e-01 1.500602e+00
		
		#choose correct thresholds that are asymmetric on the normal scale (but symmetric on the log scale)
		#BLUE
		tmp<- .Call("abcScaledF",	c(xn-1,yn-1,1/tau.up,tau.up,alpha,1e-10,100,0.05)	)
		rho.ok<- seq(1/tau.up,tau.up,by=0.001)
		pw.ok<- project.nABC.StretchedF.pow.sym(rho.ok, df, tmp[2], tmp[1])
		if(verbose)	cat( paste("tau.l",1/tau.up,"tau.up",tau.up,"cl",tmp[1],"cu",tmp[2],"\n") )
		plot(rho.ok,pw.ok,type='l')
		stop()
		
		
		#RED
		tmp<- .Call("abcScaledF",	c(xn-1,yn-1,tau.low,1/tau.low,alpha,1e-10,100,0.05)	)
		rho.ok2<- seq(tau.low,1/tau.low,by=0.001)
		pw.ok2<- project.nABC.StretchedF.pow.sym(rho.ok2, df, tmp[2], tmp[1])
		if(verbose)	cat( paste("tau.l",tau.low,"tau.up",1/tau.low,"cl",tmp[1],"cu",tmp[2],"\n") )
		
		#GREY
		#choose incorrect thresholds that are symmetric on the normal scale and max power off 1
		tmp<- .Call("abcScaledF",	c(xn-1,yn-1,tau.low,tau.up,alpha,1e-10,100,0.05)	)
		rho.sy<- seq(tau.low,tau.up,by=0.001)
		pw.sy<- project.nABC.StretchedF.pow.sym(rho.sy,df,tmp[2], tmp[1])
		if(verbose)	cat( paste("tau.l",tau.low,"tau.up",tau.up,"cl",tmp[1],"cu",tmp[2],"\n") )
		
		#BLUE
		xn<- yn<- 30
		tmp<- .Call("abcScaledF",	c(xn-1,yn-1,1/tau.up,tau.up,alpha,1e-10,100,0.05)	)
		rho.ok<- seq(1/tau.up,tau.up,by=0.001)
		pw.ok3<- project.nABC.StretchedF.pow.sym(rho.ok, df, tmp[2], tmp[1])
		if(verbose)	cat( paste("tau.l",1/tau.up,"tau.up",tau.up,"cl",tmp[1],"cu",tmp[2],"\n") )
		
		#BLUE
		xn<- yn<- 1500
		tmp<- .Call("abcScaledF",	c(xn-1,yn-1,1/tau.up,tau.up,alpha,1e-10,100,0.05)	)
		rho.ok<- seq(1/tau.up,tau.up,by=0.001)
		pw.ok4<- project.nABC.StretchedF.pow.sym(rho.ok, df, tmp[2], tmp[1])
		if(verbose)	cat( paste("tau.l",1/tau.up,"tau.up",tau.up,"cl",tmp[1],"cu",tmp[2],"\n") )
		
		
		#pdf(paste(dir.name,"/nABC.StretchedF_power.pdf",sep=''),version="1.4",width=5,height=5)
		par(mar=c(5,4,0.5,0.5))
		plot(1,1,xlim= c(0,3),ylim= range(c(pw.sy,pw.ok,pw.ok2)), type='n', xlab= expression( sigma^2== sigma[0]^2 *" "* rho), ylab="power")
		#lines(rho.sy,pw.sy, col="blue")
		polygon( c(rho.sy,tau.up,tau.low), c(pw.sy,0,0), border= NA, col= "gray80")
		lines(rho.ok,pw.ok, col="blue", lty= 1)
		lines(rho.ok,pw.ok3, col="blue", lty= 2 )
		lines(rho.ok,pw.ok4, col="blue", lty= 3 )
		#lines(rho.ok2,pw.ok2, col="red", lty= 1)
		legend<- expression('['*c^'-'*", "*c^'+'*"]= [0.5,1.5]",'['*tau^'-'*", "*tau^'+'*"]= [0.271,2.77]","", 
				'['*c^'-'*", "*c^'+'*"]= ["*1.5^{-1}*",1.5]",'['*tau^'-'*", "*tau^'+'*"]= ["*2.77^{-1}*",2.77]")#,"",
		# '['*c^'-'*", "*c^'+'*"]= [0.5,"*0.5^{-1}*"]",'['*tau^'-'*", "*tau^'+'*"]= [0.271,"*0.271^{-1}*"]")				
		legend("topright",fill=c("gray70","transparent","transparent","blue","transparent"
				#,"transparent","red","transparent"
				),legend=legend, bty= 'n', border=NA)
		abline(v=1, lty= 2)
		#dev.off()
		stop()
	}
	if(0)	#see if 'abcScaledF' works properly
	{					
		alpha <- 0.05
		tol <- 1.0e-10
		itmax <- 50
		ny1 <- 40
		ny2 <- 45
		rho1 <- 0.5625
		rho2 <- 1.7689
		itinc<- 0.05
		
		err <- -alpha
		c1 <- sqrt(rho1 * rho2)		#start C1 at midpoint and lower C1 by default by 0.05
		#if C1 is set to midpoint, then the level|rho1 must be smaller than alpha since C1,C2 are clearly too far on the right
		while (err < 0)
		{  c1 <- c1 - 0.05			
			h <- alpha + pf(c1/rho2,ny1,ny2)				
			c2 <- qf(h,ny1,ny2)*rho2									#C2 sth C2-C1 has area alpha |rho2		
			err <- pf(c2/rho1,ny1,ny2) - pf(c1/rho1,ny1,ny2) - alpha	
		}			
		c1L <- c1
		c1R <- c1 + 0.05
		it <- 0
		#now C1L,C2R are sth  rej|rho1 > alpha
		#thus the interval C1L,C2R contains the solution for C1 and the initial estimate C1
		
		
		while (abs(err) >= tol && it <= itmax)
		{  it <- it + 1
			c1 <- (c1L + c1R) / 2
			h <- alpha + pf(c1/rho2,ny1,ny2)
			c2 <- qf(h,ny1,ny2) * rho2
			err <- pf(c2/rho1,ny1,ny2) - pf(c1/rho1,ny1,ny2) - alpha			
			if (err <= 0) 
				c1R <- c1   #if interval too much on the right, then update C1R
			else
				c1L <- c1	#if interval too much on the left, then update C1L	
		}
		
		pow0 <- pf(c2,ny1,ny2) - pf(c1,ny1,ny2)
		
		cat(" ALPHA =",alpha,"   NY1 =",ny1,"   NY2 =",ny2,"   RHO1 =",rho1,"   RHO2 =",rho2,
				"   IT =",it,"\n","C1 =",c1,"   C2 =",c2,"   ERR =",err,"   POW0 =",pow0)
		
		if(!is.loaded("pdGetSummariesEpiPa"))
		{
			file<- paste(CODE.HOME,paste("lib",paste("libphylodyr",.Platform$dynlib.ext,sep=''),sep='/'),sep='')
			cat(paste("\nloading",file,'\n',sep=' '))
			dyn.load(file)
		}
		args<- c(ny1,ny2,rho1,rho2,alpha,tol,itmax,itinc)
		ans<- .Call("abcScaledF",args)
		print(ans)
	}
}