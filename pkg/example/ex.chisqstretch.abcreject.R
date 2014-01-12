#
#	ABC* rejection sampler for the normal variance example 
#	
xn			<- yn	<- 60
df			<- yn-1
ymu			<- 0 
xmu			<- 0
xsigma2		<- 1
prior.u		<- 4
prior.l		<- 0.2
N			<- 1e6
#	function to precompute ABC* simulations
chi2stretch.simu<- function(N, prior.l, prior.u, x, yn, ymu)		
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
	rownames(ans[["data"]])	<- c("ysigma2","T","log.sy2/sx2","sy2-sx2")
	ans
}
#	function to produce histogram
chi2stretch.hist<- function(x, theta, nbreaks= 20, breaks= NULL, width= 0.5, plot=0, rtn.dens=0,...)
{
	#compute break points sth theta is in the middle
	if(is.null(breaks))
	{
		breaks<- max(abs( theta - x ))*1.1 / nbreaks								
		breaks<- c( rev(seq(from= theta-breaks/2, by=-breaks, length.out= nbreaks )), seq(from= theta+breaks/2, by=breaks, length.out= nbreaks ) )		
	}
	ans.h<- hist(x, breaks=breaks, plot= 0)
	ans.h[["mean"]]<- mean(x)		
	ans.h[["hmode"]]<- mean(ans.h[["breaks"]][seq(which.max(ans.h[["counts"]]),length.out=2)])	
	tmp<- density(x, kernel="biweight",from=breaks[1],to=breaks[length(breaks)],width = max(EPS,width*diff(summary(x)[c(2,5)])))
	ans.h[["dmode"]]<- tmp[["x"]][which.max( tmp[["y"]])]
	if(rtn.dens)
		ans.h[["dens"]]<- tmp
	if(plot)
	{
		plot(ans.h, freq=0,...)
		lines(tmp)
	}
	ans.h
}
#	pseudo data
x				<- rnorm(xn,xmu,sd=sqrt(xsigma2))
x 				<- (x - mean(x))/sd(x) * sqrt(xsigma2) + xmu
#	calibrated ABC* parameters
abc.param		<- chisqstretch.calibrate(length(x), sd(x), mx.pw=0.9, alpha=0.01, max.it=100, debug=FALSE, plot=FALSE)
#	get ABC* simulations
simu			<- chi2stretch.simu(N, prior.l, prior.u, x, abc.param['n.of.y'], ymu)
#	evaluate ABC* simulations
tstat			<- simu[["data"]]["T",] 
acc				<- which( tstat>=abc.param["cl"]  &  tstat<=abc.param["cu"] )
chi2stretch.hist(simu[["data"]]["ysigma2",acc], simu[["xsigma2"]], nbreaks= 50, width= 0.5, plot=0, ylim=c(0,2.25))

