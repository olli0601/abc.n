ch11.mutost<- function()
{
	#
	require(roxygen2)	
	roxygenize('~/git/abc.star')
	#
	#
	outdir	<- '/Users/Oliver/duke/2014_ABCBookChapter/abcbook/chapters/chapter11/figures'
	#	example code for the ABC book chapter
	#
	#	TOST
	#
	xn<- 60; xmean<- 1; xsigma<- 1
	obs <- rnorm(xn, xmean, xsigma)
	obs <- (obs - mean(obs))/sd(obs) * xsigma + xmean	
	yn<- 60; ymean<- 1.85;	ysigma<- 1.2	
	sim <- rnorm(yn, ymean, ysigma)
	sim <- (sim - mean(sim))/sd(sim) * ysigma + ymean	
	#	type I error	
	mutost.calibrate(n.of.y=length(sim), s.of.y=sd(sim), tau.u=0.8, what='CR', alpha=0.01, plot=TRUE)		
	mutost.calibrate(n.of.y=length(sim), s.of.y=sd(sim), c.u=0.6, tau.u=0.8, what='ALPHA')
	#	power	
	rho 	<- seq(-0.7,0.7,0.01)	
	s.of.T 	<- sd(sim)/sqrt(length(sim))
	mutost.pow(rho, df=length(sim), s.of.T=s.of.T, c.u=0.2)
	mutost.pow(rho, df=length(sim), s.of.T=s.of.T, tau.u=0.7, alpha=0.01)
	
	#	power plot
	require(plyr)
	tmp <- ldply(c(0.3, 0.4, 0.5, 0.6, 0.9),function(tau.u)
			{
				return(data.frame(tau.u=as.factor(tau.u),rho=rho,power=mutost.pow(rho, df=length(sim), s.of.T=s.of.T, tau.u=tau.u, alpha=0.01)))
			})
	
	ggplot(tmp,aes(x=rho,y=power,colour=tau.u)) + geom_line() + labs(x=expression(rho), y='') +
			scale_colour_brewer(name=expression(tau^'+'), palette='Set1') +
			scale_y_continuous(breaks=seq(0,1,0.2)) +
			theme(legend.position=c(0,1), legend.justification=c(0,1))
	file	<- paste(outdir, '/mutost_power.pdf',sep='')
	ggsave(file, w=4, h=5)
	#
	mutost.calibrate(n.of.y=length(sim), s.of.y=sd(sim), tau.u=0.4, what='CR', alpha=0.01, plot=TRUE)
	mutost.calibrate(n.of.y=length(sim), s.of.y=sd(sim), tau.u=0.9, what='CR', alpha=0.01, plot=TRUE)
	mutost.calibrate(n.of.y=length(sim), s.of.y=sd(sim), tau.u=0.6, what='CR', alpha=0.01, plot=TRUE)
	#	power calibration	
	mutost.calibrate(n.of.y=length(sim), s.of.y=sd(sim), what='MXPW', plot=1, debug=1)
	mutost.calibrate(n.of.y=length(sim), s.of.y=sd(sim), what='MXPW', debug=0)
	file	<- paste(outdir, '/mutost_abc_cali_for_power.pdf',sep='')
	ggsave(file, w=5, h=5)
	#
	#	KL calibration
	#
	cali			<- mutost.calibrate(n.of.x=length(obs), s.of.x=sd(obs), s.of.y=sd(sim), what='KL', plot=0, debug=0)
	tau.u			<- cali['tau.u']
	n.of.y			<- cali['n.of.y']
	s.of.y			<- sd(sim)
	n.of.x			<- xn
	s.of.x			<- sd(obs)
	ssn 			<- s.of.x/sqrt(n.of.x)	
	df 				<- n.of.x - 1	
	#truncate pow and compute pow_norm
	pow_scale		<- 1.5
	pow_support 	<- c(-tau.u, tau.u) * pow_scale		
	pow_norm 		<- .Call("abc_mutost_integrate_pow", pow_support[1], pow_support[2],.Machine$double.eps^0.25,.Machine$double.eps^0.25,as.double(n.of.y-1),s.of.y/sqrt(n.of.y),tau.u,alpha,1,0)
	#compute the norm of lkl, given its support
	lkl_support 	<- pow_support
	lkl_norm 		<- diff(pt(lkl_support/ssn, df))	
	rho 			<- seq(lkl_support[1], lkl_support[2], length.out = 1000)
	df_lkl 			<- data.frame(x=rho, y=mutost.sulkl(rho, n.of.x, s.of.x, lkl_norm, lkl_support))
	df_lkl$dist 	<- "g"
	df_pow 			<- data.frame(x=rho, y=mutost.pow(rho, df=n.of.y-1, s.of.T=s.of.y/sqrt(n.of.y), tau.u=tau.u, alpha=alpha, norm=pow_norm, support= pow_support, log=FALSE))
	df_pow$dist 	<- "gabc"
	df 				<- rbind(df_pow, df_lkl)	
	ggplot(df, aes(x=x, y=y, colour=dist)) +
		geom_line() +
		scale_colour_manual(values=c("g"="black","gabc"="#E41A1C"), guide=FALSE) +
		labs(x=expression(rho), y="", colour='') 
	file			<- paste(outdir, '/mutost_abc_cali_for_KL.pdf',sep='')
	ggsave(file, w=4, h=5)
		
	#	time it
	system.time({ for(i in 1:1e4) mutost.calibrate(n.of.y=length(sim), s.of.y=sd(sim), tau.u.ub=2, what='MXPW', plot=0, debug=0) })
}


ch11.mutostabc<- function()
{
	outdir	<- '/Users/Oliver/duke/2014_ABCBookChapter/abcbook/chapters/chapter11/figures'
	abc.presim.uprior.mu<- function(abc.nit, xn, xmean, xsigma, prior.l, prior.u, ysigma, yn=NA )		
	{		
		ans			<- vector("list",5)
		names(ans)	<- c("x","xn","xmean","xsigma","sim")
		obs 		<- rnorm(xn, xmean, xsigma)
		obs 		<- (obs - mean(obs))/sd(obs) * xsigma + xmean
		ans[["x"]]			<- obs
		ans[["xmean"]]		<- xmean
		ans[["xsigma"]]		<- xsigma
		
		ans[["sim"]]		<- sapply(1:abc.nit, function(i)
				{					
					ymu		<- runif(1, prior.l, prior.u)
					y		<- rnorm(yn, ymu, sd=ysigma)
					tmp		<- c(yn, ymu, ysigma, mean(y), sd(y) )									
					tmp					
				})								
		rownames(ans[["sim"]])	<- c('m','ymu','ysigma','ysmean','yssd')
		ans
	}
	
	
	abc.presim14<- abc.presim.uprior.mu( 	abc.nit=1e7, xn=60, xmean=1.34, xsigma=1.4, 
											prior.l=1.34-5, prior.u=1.34+5, ysigma=1.4, yn=60 )
	abc.df14	<- as.data.table(t(abc.presim14$sim))								
	abc.df14[, it:=seq_len(nrow(abc.df14))]								
	tmp			<- abc.df14[,  as.list( mutost.calibrate(	n.of.y=m, s.of.y=yssd, tau.u.ub=1, what='MXPW', mx.pw=0.9, alpha=0.01)[1:4] ), by='it']
	abc.df14		<- merge(abc.df14, tmp, by='it')								
	save(file="/Users/Oliver/duke/2014_ABCBookChapter/data/abc.presim.mutost.sigma14.Rdata", abc.df14)
	#load("/Users/Oliver/duke/2014_ABCBookChapter/data/abc.presim.mutost.sigma14.Rdata")
	abc.df		<- copy(abc.df14)
	#
	# 	calibrated c.u, for predicted False Pos < 1% and mxpow=90%
	#
	abc.df[, T:= abc.df[, ysmean-1.34]]
	abc.df[, ABC_OK:= abc.df[, c.l<=T & T<=c.u]]
	abc.runs	<- data.table(RUN='ABC_OK', n.of.y=60, s.of.y=1.4, tau.u=mean(abc.df$tau.u), c.l=mean(abc.df$c.l), c.u=mean(abc.df$c.u) )
	#
	#	ABC fails to identify true mean: choose too large ABC tolerances that correspond to large tau.u and alpha=0.01
	#
	tmp			<- mutost.calibrate(	n.of.y=60, s.of.y=1.4, tau.u=1.25, what='CR', alpha=0.01, plot=1)
	abc.runs	<- rbind(abc.runs, as.data.table(c(RUN='ABC_TOOLARGE', n.of.y=60, s.of.y=1.4, tau.u=1.25, as.list(tmp))))
	#c.l        c.u 
	#-0.8178112  0.8178112 
	abc.df[, ABC_TOOLARGE:= abc.df[, tmp['c.l']<=T & T<=tmp['c.u']]]	
	#
	#	ABC with very small acceptance probability even at rho==0
	#
	tmp		<- mutost.calibrate(	n.of.y=60, s.of.y=1.4, tau.u=0.5, what='CR', alpha=0.01, plot=1)
	#c.l         c.u 
	#-0.06781116  0.06781116
	abc.runs<- rbind(abc.runs, as.data.table(c(RUN='ABC_TOOSMALL', n.of.y=60, s.of.y=1.4, tau.u=0.5, as.list(tmp))))
	abc.df[, ABC_TOOSMALL:= abc.df[, tmp['c.l']<=T & T<=tmp['c.u']]]
	#
	#	TABLES
	#	get acc prob, acc prob around rho*, alpha
	tmp		<- melt(abc.df, id.vars=c('it','ymu','tau.u'), measure.vars=c('ABC_TOOSMALL','ABC_TOOLARGE','ABC_OK'), variable.name='RUN')
	set(tmp, tmp[, which(RUN=='ABC_TOOSMALL')], 'tau.u', abc.runs[RUN=='ABC_TOOSMALL',c.u])
	set(tmp, tmp[, which(RUN=='ABC_TOOLARGE')], 'tau.u', abc.runs[RUN=='ABC_TOOLARGE',c.u])			
	abc.runs<- merge(abc.runs, tmp[, list( ACCPROB=mean(value)), by='RUN'], by='RUN')	
	abc.runs<- merge(abc.runs, subset(tmp, abs(ymu-1.34)<0.01)[, list( ACCPROB_RHOSTAR=mean(value)), by='RUN'], by='RUN')	
	abc.runs<- merge(abc.runs, subset(tmp, abs(ymu-1.34)>tau.u)[, list( FPR=mean(value)), by='RUN'], by='RUN')
	abc.runs<- merge(abc.runs, subset(tmp, abs(ymu-1.34)-tau.u>=0 & abs(ymu-1.34)-tau.u<0.005)[, list( ALPHA=mean(value)), by='RUN'], by='RUN')
	 
	#
	#	PLOTS
	#	illustrate how ABC depends on power properties
	#	show tolerances
	tmp		<- subset(abc.df, it<=500, select=c(it, c.l, c.u, tau.u))
	tmp[, RUN:='ABC_OK']
	tmp		<- rbind( tmp, abc.runs[, list(it=seq_len(500)), by=c('RUN','tau.u','c.l','c.u')], use.names=TRUE)
	tmp		<- merge(tmp, data.table(RUN=c('ABC_OK','ABC_TOOLARGE','ABC_TOOSMALL'), LEGEND=c('calibrated ABC with\npredicted max power ~ 0.9','standard ABC with\npredicted flat power ~ 1','standard ABC with\npredicted small power')), by='RUN')
	ggplot(tmp, aes(x=it, y=c.u, colour=LEGEND)) + geom_step(direction='vh',show_guide=FALSE) + theme_bw() +
			facet_wrap(~LEGEND, nrow=1) + labs(x='ABC iteration', y=expression(atop('Upper ABC tolerance',c^'+'))) +
			scale_colour_brewer(palette='Set2')	
	ggplot(tmp, aes(x=it, ymax=c.u, ymin=-c.u, fill=LEGEND)) + geom_ribbon(show_guide=FALSE) + theme_bw() +
			facet_wrap(~LEGEND, nrow=1) + labs(x='ABC iteration', y=expression(atop('ABC accept/reject region','[ '*h^'-'*', '*h^'+'*' ]'))) +
			scale_fill_brewer(palette='Set2')		
	ggsave(file=paste(outdir,'/mutostabc_tol.pdf',sep=''), w=9, h=3)
	#	show power
	abc.pw	<- abc.runs[, 	{
				rho			<- seq(-1.25,1.25,len=1024)
				pow_norm	<- mutost.pow.norm(n.of.y-1, s.of.y/sqrt(n.of.y), tau.u, 0.01, support=c(-tau.u*1.25, tau.u*1.25))
				list(rho=rho, pow=mutost.pow(rho, df=n.of.y-1, s.of.y/sqrt(n.of.y), tau.u=tau.u, alpha=0.01, support=c(-tau.u*1.25, tau.u*1.25)), pow_norm=pow_norm)
			}, by='RUN']					
	ggplot( abc.pw, aes(x=rho, ymax=pow, ymin=0, fill=RUN, group=RUN)) + geom_ribbon(show_guide=FALSE) +  
			scale_fill_brewer(palette='Set2') +
			scale_y_continuous(breaks=seq(0,1,0.2)) +
			theme_bw() + labs(x=expression(rho), y='power') + theme(strip.background = element_blank(), strip.text = element_blank()) +
			facet_wrap(~RUN, nrow=1)
	ggsave(file=paste(outdir,'/mutostabc_power.pdf',sep=''), w=9, h=3)
	#	show acceptances 
	abc.plot<- subset(melt(abc.df, id.vars=c('ymu','T'), measure.vars=c('ABC_OK','ABC_TOOLARGE','ABC_TOOSMALL'), variable.name='RUN'), value)
	ggplot( abc.plot, aes(fill=RUN)) +
			geom_histogram(aes(x=ymu-1.34, y= ..density..), breaks=seq(-1.25,1.25, length.out=40), show_guide=FALSE) +
			geom_line(data=abc.pw, aes(x=rho, y=pow/pow_norm), col='black') +			
			labs(x=expression(rho), y='ABC posterior density') + theme_bw() + theme(strip.background = element_blank(), strip.text = element_blank()) +
			scale_fill_brewer(palette='Set2') + 
			facet_wrap(~RUN, nrow=1)	
	ggsave(file=paste(outdir,'/mutostabc_accepted.pdf',sep=''), w=9, h=3)
	
}

ch11.tosz<- function()
{
	require(abc.star)
	outdir	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_ABCBookChapter/v160317'
	rho		<- seq(-2,2,0.01)
	tau.u	<- 1.5
	alpha	<- 0.01
	sigma	<- 0.5
	ztest.pow(rho, tau.u, alpha, sigma, norm = 1, support= c(-Inf,Inf), log=FALSE)
	
	sigma	<- 1.2
	n2s		<- function(n){ sigma/sqrt(n) }	#need formula to convert n.of.y into s.of.T, depends on application
	s2n		<- function(s){ round((sigma/s)^2) }	#need formula to convert s.of.T into n.of.y, depends on application	
	ztest.calibrate(n.of.x=80, n2s=n2s, s2n=s2n, what='KL', debug=FALSE, plot=TRUE, verbose=FALSE)
	
	#	power change tau-
	rho		<- seq(-1.2,1.2,0.01)
	tmp		<- lapply( c(0.35, 0.4, 0.5, 0.7, 0.9), function(tau.u)
			{
				data.table(rho=rho, tau.u=tau.u, power=ztest.pow(rho, tau.u, alpha=0.01, sigma=sigma/sqrt(80)))			
			})
	tmp	<- do.call('rbind', tmp)
	set(tmp, NULL, 'tau.u', tmp[, factor(tau.u)])
	ggplot(tmp,aes(x=rho,y=power,colour=tau.u, group=tau.u)) + geom_line() + labs(x=expression(rho), y='Power\n(ABC acceptance probability)') +
			scale_colour_brewer(name=expression(tau^'+'), palette='Set1') +
			scale_y_continuous(breaks=seq(0,1,0.2)) +
			theme(legend.position=c(1,1), legend.justification=c(1,1))
	file	<- paste(outdir, '/ztest_power.pdf',sep='')
	ggsave(file, w=4, h=5)
		
	#	power change m against likelihood
	rho			<- seq(-0.7,0.7,0.01)
	support		<- range(rho)
	s.of.x		<- sigma/sqrt(80)
	tmp			<- lapply(c(80, 120, 180, 240), function(ny)
			{
				cali		<- ztest.calibrate(n.of.y=ny, n2s=n2s, s2n=s2n, what='MXPW')
				s.of.T		<- sigma/sqrt(ny)
				pw.norm		<- ztest.pow.norm(cali['tau.u'], s.of.T, alpha=alpha, support=support)
				pw			<- ztest.pow(rho, cali['tau.u'], alpha, s.of.T, norm=pw.norm, support=support, log=FALSE)	
				lkl.norm	<- ztest.sulkl.norm(s.of.x, support=support)	
				lkl			<- ztest.sulkl(rho, s.of.x, norm=lkl.norm, support=support, log=FALSE)	
				df_lkl 		<- data.table(x=rho, y=lkl)
				df_lkl$dist <- "likelihood"
				df_pow 		<- data.table(x=rho, y=pw)
				df_pow$dist <- "ABC"
				gdf 		<- rbind(df_pow, df_lkl)
				gdf$ny		<- ny
				gdf
			})
	tmp	<- do.call('rbind', tmp)
	tmp	<- subset(tmp, dist=='ABC' | dist=='likelihood' & ny==80)
	set(tmp, NULL, 'ny', tmp[, factor(ny)])	
	ggplot(data=subset(tmp, dist=='ABC'), aes(x=x, y=y)) + geom_line(aes(colour=ny, group=ny)) +
			geom_line(data=subset(tmp, dist=='likelihood')) +
			labs(x=expression(rho), y='', colour='m') +
			theme(legend.position=c(1,1), legend.justification=c(1,1))
	file	<- paste(outdir, '/ztest_lkl.pdf',sep='')
	ggsave(file, w=4, h=5)
}

ch11.vartest<- function()
{
	n.of.x	<- 60
	n.of.y	<- 60
	outdir	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_ABCBookChapter/v160317'
	
	ans		<- vartest.calibrate(n.of.x=n.of.x, n.of.y=n.of.y, tau.u.ub=3, what='MXPW', mx.pw=0.9, alpha=0.01, plot=TRUE, verbose=FALSE)
	file	<- paste(outdir, '/vartest_abc_cali_for_power.pdf',sep='')
	ggsave(file, w=5, h=5)
	#
	#	POWER
	#
	
	rho		<- seq(0.1, 3, len=1024)
	tmp		<- lapply(c(1.4,1.7,2.0,2.4,3.0),function(tau.u)
			{
				cali	<- vartest.calibrate(n.of.x=n.of.x, n.of.y=n.of.y, tau.l=1/tau.u, tau.u=tau.u, what='CR', alpha=0.01)
				data.table(rho=rho, power=vartest.pow(rho, n.of.x, n.of.y-1, cali['c.l'], cali['c.u']), tau.u=tau.u)			
			})
	tmp		<- do.call('rbind', tmp)
	set(tmp, NULL, 'tau.u', tmp[, factor(tau.u)])
	tmp2	<- tmp[, list(rho= rho[which.max(power)], power=max(power)), by='tau.u']
	ggplot(tmp,aes(x=rho,y=power,colour=tau.u, group=tau.u)) + geom_line() + labs(x=expression(rho), y='') +
			geom_point(data=tmp2, colour='black') +
			scale_colour_brewer(name=expression(tau^'+'), palette='Set1') +
			scale_y_continuous(breaks=seq(0,1,0.2)) +
			theme(legend.position=c(1,1), legend.justification=c(1,1))
	file	<- paste(outdir, '/vartest_power.pdf',sep='')
	ggsave(file, w=4, h=5)
	#
	#	KL
	#
	s.of.x		<- 1.42	
	cali		<- vartest.calibrate(n.of.x=n.of.x, s.of.x=s.of.x, what='KL', mx.pw=0.9, alpha=0.01, plot=TRUE, verbose=FALSE)
	c.l			<- cali['c.l']
	c.u			<- cali['c.u']	
	tau.l		<- cali['tau.l']
	tau.u		<- cali['tau.u']
	n.of.y		<- cali['n.of.y']
	pow.scale	<- 1.5
	scale		<- n.of.x-1
	df			<- n.of.y-1
	#
	pow_support 	<- c(tau.l/pow_scale, tau.u*pow_scale) 	
	pow_norm 		<- vartest.pow.norm(scale, df, c.l, c.u, trafo=1, support=pow_support) 
	lkl_support		<- pow_support	
	lkl_norm		<- vartest.su.lkl.norm(n.of.x, support=lkl_support)
	
	rho_lkl 		<- seq(lkl_support[1], lkl_support[2], length.out = 1000)
	df_lkl 			<- data.frame(x=rho_lkl, y=vartest.sulkl(rho_lkl, n.of.x, norm=lkl_norm, support=lkl_support))
	df_lkl$dist 	<- "g"
	rho_pow	 		<- seq(pow_support[1], pow_support[2], length.out = 1000) 		
	df_pow 			<- data.frame(x=rho_pow, y=vartest.pow(rho_pow, scale, df, c.l, c.u, trafo= 1, norm=pow_norm))
	df_pow$dist 	<- "gabc"	
	gdf 			<- rbind(df_pow, df_lkl)
	
	ggplot(data = gdf, aes(x=x, y=y, colour=dist)) +
			geom_line() +
			scale_colour_manual(values=c("g"="black","gabc"="#E41A1C"), guide=FALSE) +
			theme(legend.position='bottom') +
			labs(x=expression(rho))	
	file			<- paste(outdir, '/vartest_KL.pdf',sep='')
	ggsave(file, w=4, h=5)
}

ch11.ratetestabc<- function()
{
	abc.presim.uprior.beta<- function(abc.nit, n.of.x, beta.x, prior.l, prior.u, n.of.y=NA ){
		
		ans  		<- vector("list",4)
		names(ans)	<- c("x","n.of.x","beta.x","sim")
		
		# xn is the number of generated samples
		obs 		<- rexp(n.of.x, rate=1/beta.x)
		#obs 		<- (obs - mean(obs))/sd(obs) * xsigma + xmean  # normalization of the generated normal samples
		
		ans[["x"]]			<- obs
		ans[['n.of.x']] <- n.of.x
		ans[["beta.x"]]		<- beta.x
		
		# 'sim' is for simulating simulated samples for given observed x samples
		# 'abc.nit' is number of ABC iterations
		ans[["sim"]]		<- sapply(1:abc.nit,  function(i){
					
					ybeta		<- runif(1, prior.l, prior.u) # uniform prior
					y		<- rexp(n.of.y, 1/ybeta)  # generate simulated samples based on uniform prior
					tmp		<- c( n.of.y,  ybeta,  mean(y) )									
					tmp					
				}   )
		
		rownames(ans[["sim"]])	<- c('m', 'ybeta', 'ysmean')
		
		ans
	}
	#==========================
	ratetest.plot<- function(n.of.y, c.l, c.u, tau.l, tau.u, pow_scale=1.5){
		
		pow_support <- c(tau.l/pow_scale, tau.u*pow_scale)   
		pow_norm 	<- ratetest.pow.norm(c.l=c.l, c.u=c.u, m=n.of.y, trafo=1, support=pow_support)	
		
		tmp			<- data.frame(rho=seq(pow_support[1], pow_support[2], length.out = 1024))	
		tmp$power	<- ratetest.pow(tmp$rho, c.l=c.l, c.u=c.u, m=n.of.y, norm=pow_norm, trafo= 1)*pow_norm	
		
		p	<- ggplot(tmp, aes(x=rho, y=power)) + geom_line() + labs(x=expression(rho), y='Power\n(ABC acceptance probability)') +
				scale_y_continuous(breaks=seq(0,1,0.2), limits=c(0,1)) +
				scale_x_continuous(limits=c(0,4)) +
				geom_vline(xintercept = c(tau.l, tau.u), linetype = "dotted") +
				geom_vline(xintercept = c(c.l, c.u), linetype = "dashed") +
				ggtitle(paste("n.of.y=", n.of.y, "\ntau.l=", round(tau.l,d=5), " tau.u=", round(tau.u,d=5), "\nc.l=", round(c.l,d=5), " c.u=", round(c.u,d=5)))
		print(p)
	}
	
	start.time <- Sys.time()
	abc.presim <- abc.presim.uprior.beta(abc.nit=100000, n.of.x=40, beta.x=2, prior.l=0.1, prior.u=4, n.of.y=40)
	end.time <- Sys.time()
	time.taken <- end.time - start.time
	time.taken
	
	
	abc.df  <- as.data.table(t(abc.presim$sim))	
	abc.df[ , it:=seq_len(nrow(abc.df))]	
	print(abc.df)
	
	start.time <- Sys.time()
	tmp  <- abc.df[ , as.list(ratetest.calibrate(n.of.y = m, mx.pw=0.9, tau.u.ub=5, what='MXPW')[1:4]), by='it']
	end.time <- Sys.time()
	time.taken <- end.time - start.time
	time.taken
	
	print(tmp[1:10,])
	
	abc.df14  	<- merge(abc.df, tmp, by='it')		
	
	abc.df  	<- copy(abc.df14)
#abc.df <- abc.df14
	print(abc.df[1:10,])
	is.data.table(abc.df)
	
	abc.df <- as.data.table(abc.df)
	
	abc.df[,T:=abc.df[,ysmean/2]]
	print(abc.df[1:10,])
	
	abc.df[, ABC_OK := abc.df[, c.l <=T & T <= c.u]]
	print(abc.df[1:100,])
	
	abc.ok<- data.table(RUN='ABC_OK', n.of.y=40, tau.l=mean(abc.df$tau.l), tau.u=mean(abc.df$tau.u), c.l=mean(abc.df$c.l), c.u=mean(abc.df$c.u) )
	print(abc.ok)
	
	ratetest.plot(n.of.y=40, c.l=0.7601569, c.u=1.28276, tau.l=0.5413799, tau.u=1.916706)
	
	sel.seq <- abc.df[which(abc.df[['ABC_OK']]=='TRUE'), ybeta]
	
	hist(sel.seq, freq=F, breaks=20)
	
	sell <- sel.seq/2
	d <- density(sel.seq)
	plot(d,ylim=c(0,1), main='')
	
	d <- density(sell)
	plot(d, ylim=c(0,2), xlim=c(0,3), main='')
	abline(v=1, lty='dotted')
}