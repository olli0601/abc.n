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
	yn<- 60; ymean<- 1.85;	ysigma<- 0.7	
	sim <- rnorm(yn, ymean, ysigma)
	sim <- (sim - mean(sim))/sd(sim) * ysigma + ymean	
	#	type I error	
	ans <- mutost.calibrate(n.of.y=length(sim), s.of.y=sd(sim), tau.u=0.8, what='CR', alpha=alpha, plot=TRUE)		
	ans	<- mutost.calibrate(n.of.y=length(sim), s.of.y=sd(sim), c.u=0.6, tau.u=0.7, what='ALPHA')
	#	power	
	rho 	<- seq(-0.7,0.7,0.01)	
	s.of.T 	<- sd(sim)/sqrt(length(sim))
	accprob	<- mutost.pow(rho, df=length(sim), s.of.T=s.of.T, tau.u=0.3, alpha=0.01)
	#	power plot
	require(plyr)
	tmp <- ldply(c(0.23,0.25, 0.3, 0.4, 0.6),function(tau.u)
			{
				return(data.frame(tau.u=as.factor(tau.u),rho=rho,power=mutost.pow(rho, df=length(sim), s.of.T=s.of.T, tau.u=tau.u, alpha=0.01)))
			})
	
	ggplot(tmp,aes(x=rho,y=power,colour=tau.u)) + geom_line() + labs(x=expression(rho), y='Power\n(ABC acceptance probability)') +
			scale_colour_brewer(name=expression(tau^'+'), palette='Set1') +
			scale_y_continuous(breaks=seq(0,1,0.2)) +
			theme(legend.position=c(0,1), legend.justification=c(0,1))
	file	<- paste(outdir, '/mutost_power.pdf',sep='')
	ggsave(file, w=5, h=5)
	#
	mutost.calibrate(n.of.y=length(sim), s.of.y=sd(sim), tau.u=0.23, what='CR', alpha=alpha, plot=TRUE)
	mutost.calibrate(n.of.y=length(sim), s.of.y=sd(sim), tau.u=0.6, what='CR', alpha=alpha, plot=TRUE)
	mutost.calibrate(n.of.y=length(sim), s.of.y=sd(sim), tau.u=0.4, what='CR', alpha=alpha, plot=TRUE)
	#	power calibration	
	mutost.calibrate(n.of.y=length(sim), s.of.y=sd(sim), tau.u.ub=2, what='MXPW', plot=1, debug=1)
	file	<- paste(outdir, '/mutost_abc_cali_for_power.pdf',sep='')
	ggsave(file, w=5, h=5)
	#	time it
	system.time({ for(i in 1:1e4) mutost.calibrate(n.of.y=length(sim), s.of.y=sd(sim), tau.u.ub=2, what='MXPW', plot=0, debug=0) })	
}

ch11.vartest<- function()
{
	n.of.x	<- 60
	n.of.y	<- 60
	outdir	<- '/Users/Oliver/duke/2014_ABCBookChapter/abcbook/chapters/chapter11/figures'

	ans	<- vartest.calibrate(n.of.x=n.of.x, n.of.y=n.of.y, tau.u.ub=3, what='MXPW', mx.pw=0.9, alpha=0.01, plot=TRUE, verbose=FALSE)
	file	<- paste(outdir, '/vartest_abc_cali_for_power.pdf',sep='')
	ggsave(file, w=5, h=5)
	
	rho		<- seq(0.1, 3, len=1024)
	tmp		<- lapply(c(1.4,1.7,2.0,2.4,3.0),function(tau.u)
			{
				cali	<- vartest.calibrate(n.of.x=n.of.x, n.of.y=n.of.y, tau.l=1/tau.u, tau.u=tau.u, what='CR', alpha=0.01)
				data.table(rho=rho, power=vartest.pow(rho, n.of.x, n.of.y-1, cali['c.l'], cali['c.u']), tau.u=tau.u)			
			})
	tmp	<- do.call('rbind', tmp)
	set(tmp, NULL, 'tau.u', tmp[, factor(tau.u)])
	ggplot(tmp,aes(x=rho,y=power,colour=tau.u, group=tau.u)) + geom_line() + labs(x=expression(rho), y='Power\n(ABC acceptance probability)') +
			scale_colour_brewer(name=expression(tau^'+'), palette='Set1') +
			scale_y_continuous(breaks=seq(0,1,0.2)) +
			theme(legend.position=c(1,1), legend.justification=c(1,1))
	file	<- paste(outdir, '/vartest_power.pdf',sep='')
	ggsave(file, w=5, h=5)
	
	
	
}