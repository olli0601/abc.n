#------------------------------------------------------------------------------------------------------------------------
# nABC PROJECTS
#------------------------------------------------------------------------------------------------------------------------
project.nABC.pval<- function()
{
	#suppose location shift
	#-- try reproduce 'good fit': shift equals c+
	#-- increase c+ and shift and explore qq plot
	#--
	dir.name<- "/Users/olli0601/duke/2012_frequencyABC/sim.data"
	fade<- 0.5
	N<- 1e2
	alpha<- 0.01
	x.sigma<- 2
	x.mu<- 0
	x.n<- 100
	tau<- 0.75
	shift<- 0
	x<- rnorm(x.n, x.mu, x.sigma)
	x<- (( x - mean(x) ) / sd(x) + x.mu )* x.sigma

	project.nABC.pval.simtost<- function(N,x.n,x.mu,x.sigma,shift,tau,alpha=0.01,resample=0)
	{
		sapply(seq_len(N),function(i){
				y<- rnorm(x.n, x.mu+shift, x.sigma)
				if(resample)
					z<-	sample(x,x.n,replace=1)
				else
					z<- x
				ans<- nabc.mutost(y, z, args= NA, verbose= 0, alpha= alpha, tau.u= tau )
				ans["pval"]<- ans["pval"]*( 1 - ans["ar"] ) + ans["ar"]/2		#my code automatically rescales the pvals
				ans["sim.mean"]<- mean(y)
				ans
			})
	}
	project.nABC.pval.simtostu<- function(N,x.n,x.mu,x.sigma,shift,tau,thresh,alpha=0.01,resample=0)
	{
		shift<- runif(N,-thresh,thresh)+shift
		sapply(seq_len(N),function(i){
					y<- rnorm(x.n, x.mu+shift[i], x.sigma)
					if(resample)
						z<-	sample(x,x.n,replace=1)
					else
						z<- x
					ans<- nabc.mutost(y, z, args= NA, verbose= 0, alpha= alpha, tau.u= tau )
					ans["pval"]<- ans["pval"]*( 1 - ans["ar"] ) + ans["ar"]/2		#my code automatically rescales the pvals
					ans["sim.mean"]<- mean(y)
					ans
				})
	}
	project.nABC.pval.qq<- function(x, add, ...){
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


	if(1)
	{
		#compare pvals for different tau when simus are from uniform prior
		pdf(paste(dir.name,"/TOST_pval4b.pdf",sep=''),version="1.4",width=4,height=4)
		par(mar=c(4.5,4.5,0.5,0.5))

		N<- 1e3
		x.sigma<- 0.3
		shift<- c(0.5,2, 3)
		tau<- c(2, 2, 2)
		alpha<- 0.01
		for(i in 1:40)
		{
			ans<- project.nABC.pval.simtostu(N,x.n,x.mu,x.sigma,shift[1],tau[1],tau[1],alpha=alpha,resample=1)
			acc<- which(ans["error",]<=ans["cir",])
			print(length(acc) / ncol(ans))
			#'al' prints the point null t statistic of accepted samples against the t density. clearly, these T are more broadly distributed, which explains the U shape
			#hist( ans["al",acc], freq=0 ); t<- seq(-6,6,0.01); lines(t, dt(t,103))
			project.nABC.pval.qq((ans["pval",acc]-ans["ar",acc]/2) / (1-ans["ar",acc]), ifelse(i==1,0,1), col=my.fade.col("red", fade))

			ans<- project.nABC.pval.simtostu(N,x.n,x.mu,x.sigma,shift[2],tau[2],tau[2],alpha=alpha,resample=1)
			acc<- which(ans["error",]<=ans["cir",])
			print(length(acc) / ncol(ans))
			project.nABC.pval.qq((ans["pval",acc]-ans["ar",acc]/2) / (1-ans["ar",acc]), 1, col=my.fade.col("green", fade))

			ans<- project.nABC.pval.simtostu(N,x.n,x.mu,x.sigma,shift[3],tau[3],tau[3],alpha=alpha,resample=1)
			acc<- which(ans["error",]<=ans["cir",])
			print(length(acc) / ncol(ans))
			project.nABC.pval.qq((ans["pval",acc]-ans["ar",acc]/2) / (1-ans["ar",acc]), 1, col=my.fade.col("blue", fade))

			if(i==40)
			{
				abline(a=0,b=1, lty= 2)
				legend("topleft",fill= c("red","green","blue"),legend= bquote(.(parse(text=paste("sh==~",shift,sep="")))), bty='n')
			}
		}
		dev.off()
	}
	if(0)
	{
		#compare pvals for different tau when simus are from uniform prior
		pdf(paste(dir.name,"/TOST_pval4.pdf",sep=''),version="1.4",width=4,height=4)
		par(mar=c(4.5,4.5,0.5,0.5))

		N<- 1e3
		x.sigma<- 0.3
		tau<- c(0.5, 1, 2)
		alpha<- 0.01
		for(i in 1:40)
		{
			ans<- project.nABC.pval.simtostu(N,x.n,x.mu,x.sigma,0,tau[1],tau[1],alpha=alpha,resample=1)
			acc<- which(ans["error",]<=ans["cir",])
			print(length(acc) / ncol(ans))
			#'al' prints the point null t statistic of accepted samples against the t density. clearly, these T are more broadly distributed, which explains the U shape
			#hist( ans["al",acc], freq=0 ); t<- seq(-6,6,0.01); lines(t, dt(t,103))
			project.nABC.pval.qq((ans["pval",acc]-ans["ar",acc]/2) / (1-ans["ar",acc]), ifelse(i==1,0,1), col=my.fade.col("red", fade))

			ans<- project.nABC.pval.simtostu(N,x.n,x.mu,x.sigma,0,tau[2],tau[2],alpha=alpha,resample=1)
			acc<- which(ans["error",]<=ans["cir",])
			print(length(acc) / ncol(ans))
			project.nABC.pval.qq((ans["pval",acc]-ans["ar",acc]/2) / (1-ans["ar",acc]), 1, col=my.fade.col("green", fade))

			ans<- project.nABC.pval.simtostu(N,x.n,x.mu,x.sigma,0,tau[3],tau[3],alpha=alpha,resample=1)
			acc<- which(ans["error",]<=ans["cir",])
			print(length(acc) / ncol(ans))
			project.nABC.pval.qq((ans["pval",acc]-ans["ar",acc]/2) / (1-ans["ar",acc]), 1, col=my.fade.col("blue", fade))

			if(i==40)
			{
				abline(a=0,b=1, lty= 2)
				legend("topleft",fill= c("red","green","blue"),legend= bquote(.(parse(text=paste(expression(tau),"==~",tau,sep="")))), bty='n')
			}
		}
		dev.off()
	}
	if(0)
	{
		#compare pvals across different shifts and different tau
		pdf(paste(dir.name,"/TOST_pval3.pdf",sep=''),version="1.4",width=4,height=4)
		par(mar=c(4.5,4.5,0.5,0.5))
		for(i in 1:40)
		{
			shift<- 0.2
			ans<- project.nABC.pval.simtost(N,x.n,x.mu,x.sigma,shift,shift+0.6,resample=1)
			acc<- which(ans["error",]<=ans["cir",])
			print(length(acc))
			project.nABC.pval.qq((ans["pval",acc]-ans["ar",acc]/2) / (1-ans["ar",acc]), ifelse(i==1,0,1), col=my.fade.col("red", fade))

			shift<- 0.5
			ans<- project.nABC.pval.simtost(N,x.n,x.mu,x.sigma,shift,shift+0.6,resample=1)
			acc<- which(ans["error",]<=ans["cir",])
			print(length(acc))
			project.nABC.pval.qq((ans["pval",acc]-ans["ar",acc]/2) / (1-ans["ar",acc]), 1, col=my.fade.col("green", fade))


			shift<- 1
			ans<- project.nABC.pval.simtost(N,x.n,x.mu,x.sigma,shift,shift+0.6,resample=1)
			acc<- which(ans["error",]<=ans["cir",])
			print(length(acc))
			project.nABC.pval.qq((ans["pval",acc]-ans["ar",acc]/2) / (1-ans["ar",acc]), 1, col=my.fade.col("blue", fade))

			if(i==40)
			{
				abline(a=0,b=1, lty= 2)
				legend("topleft",fill= c("blue","green","red"),legend= bquote(.(parse(text=paste("sh==~",c(0.2,0.5,1),sep="")))), bty='n')
			}
		}
		dev.off()
	}
	if(0)
	{
		#compare raw and adjusted pvalues when the model is true, case tau<<sigma
		pdf(paste(dir.name,"/TOST_pval1.pdf",sep=''),version="1.4",width=4,height=4)
		par(mar=c(4.5,4.5,0.5,0.5))
		for(i in 1:40)
		{
			ans<- project.nABC.pval.simtost(N,x.n,x.mu,x.sigma,shift,tau, resample=1)
			acc<- which(ans["error",]<=ans["cir",])
			print(length(acc))
			#qq plot on all and accepted
			#project.nABC.pval.qq(ans["pval",], ifelse(i==1,0,1))
			y<- ans["pval",acc]
			project.nABC.pval.qq(y, ifelse(i==1,0,1), col=my.fade.col("blue", fade))
			project.nABC.pval.qq((y-ans["ar",acc]/2) / (1-ans["ar",acc]), 1, col=my.fade.col("red", fade))
			if(i==40)
			{
				abline(a=0,b=1, lty= 2)
				legend("topleft",fill= c("blue","red"),legend= c("pval", "rescaled p-val"), bty='n')
			}
		}
		dev.off()
	}
	if(0)
	{
		resample<- 1
		#compare raw and adjusted pvalues when the model is true, case tau>>sigma
		#if resample==0, then this does not give a uniform distribution !!
		pdf(paste(dir.name,"/TOST_pval1b.pdf",sep=''),version="1.4",width=4,height=4)
		tau<- 2
		x.sigma<- 0.75
		par(mar=c(4.5,4.5,0.5,0.5))
		for(i in 1:40)
		{
			ans<- project.nABC.pval.simtost(N,x.n,x.mu,x.sigma,shift,tau, resample=resample)
			acc<- which(ans["error",]<=ans["cir",])
			print(length(acc))
			#qq plot on all and accepted
			#project.nABC.pval.qq(ans["pval",], ifelse(i==1,0,1))
			y<- ans["pval",acc]
			project.nABC.pval.qq(y, ifelse(i==1,0,1), col=my.fade.col("blue", fade))
			project.nABC.pval.qq((y-ans["ar",acc]/2) / (1-ans["ar",acc]), 1, col=my.fade.col("red", fade))
			if(i==40)
			{
				abline(a=0,b=1, lty= 2)
				legend("topleft",fill= c("blue","red"),legend= c("pval", "rescaled p-val"), bty='n')
			}
		}
		dev.off()
	}
	if(0)
	{
		#compare pvals when averaged over obs data and when this is not done
		pdf(paste(dir.name,"/TOST_pval2.pdf",sep=''),version="1.4",width=4,height=4)
		tau<- 2
		x.sigma<- 0.75
		par(mar=c(4.5,4.5,0.5,0.5))
		print(shift)
		for(i in 1:40)
		{
			p<- sapply(seq_len(N),function(i){
						y<- rnorm(x.n, x.mu+shift, x.sigma)
						z<- rnorm(x.n, x.mu, x.sigma)
						t.test(z,y)$p.value
					})
			project.nABC.pval.qq(p,  ifelse(i==1,0,1), col=my.fade.col("red", fade))
			p<- sapply(seq_len(N),function(i){
						y<- rnorm(x.n, x.mu+shift, x.sigma)
						t.test(x,y)$p.value
					})
			project.nABC.pval.qq(p,  1, col=my.fade.col("blue", fade))
			p<- sapply(seq_len(N),function(i){
						y<- rnorm(x.n, x.mu+shift, x.sigma)
						z<- sample(x,x.n,replace=1)
						t.test(z,y)$p.value
					})
			project.nABC.pval.qq(p,  1, col=my.fade.col("green", fade))

			if(i==40)
			{
				abline(a=0,b=1, lty= 2)
				legend("topleft",fill= c("blue","red","green"),legend= c("fixed obs", "random obs","resampled obs"), bty='n')
			}
		}
		dev.off()
	}

	#print(ans)
	stop()
}
#------------------------------------------------------------------------------------------------------------------------
project.nABC.plotGAUSSresults<- function()
{
	dir.name<- "/Users/olli0601/duke/2012_frequencyABC/sim.data"
	if(0)
	{
		#compare estimates of sigma^2 when SAME confidence threshold
		h<- 0.04
		res.stabc<- c(1.29, 0.56, 0.48, 2.47,		1.25,0.48,0.52,2.22,		1.16,0.42,0.54,2.03,		1.13,0.34,0.6,1.86,		1.09,0.28,0.65,1.69,		1.06,0.21,0.7,1.47,		1.05,0.17,0.78,1.42)
		res.fabc<- c(1.27,0.53,0.48,2.36,		1.23,0.46,0.52,2.13,		1.18,0.4,0.55,1.98,		1.14,0.34,0.6,1.86,		1.09,0.27,0.66,1.67,		1.06,0.21,0.71,1.52,		1.04,0.16,0.77,1.38)
		res.stabc<- matrix( res.stabc, nrow=4,ncol=7, dimnames=list(c("mean","sd","95l","95u"),seq(3,1.5,-0.25))  )
		res.fabc<- matrix( res.fabc, nrow=4,ncol=7, dimnames=list(c("mean","sd","95l","95u"),seq(3,1.5,-0.25))  )
		print(res.stabc)

		tau<- as.numeric(colnames(res.stabc))
		xlim<- range(tau)+c(-0.25,0.25)
		ylim<- c( min(c(res.stabc["95l",],res.fabc["95l",])),  max(c(res.stabc["95u",],res.fabc["95u",])) )
		pdf(paste(dir.name,"/GAUSS_sametolerance.pdf",sep=''),version="1.4",width=6,height=6)
		plot(1,1,type='n',xlim=xlim,ylim=ylim,xlab="threshold",ylab=expression(paste("posterior ",sigma^2," : mean, 95%")))
		abline(h=1,col="blue",lty=1)
		m<- res.stabc
		h2<- -0.015
		sapply(seq_len(ncol(m)),function(i)
				{
					#polygon(tau[i]+c(-h,h,h,-h), c(m["95u",i],m["95u",i],m["95l",i],m["95l",i]),border="red")
					lines(tau[i]+h2+c(-h,h) ,  c(m["95u",i],m["95u",i]))
					lines(tau[i]+h2+c(-h,h) ,  c(m["95l",i],m["95l",i]))
					lines(h2+c(tau[i],tau[i]),	c(m["95l",i],m["95u",i]), lty=2	)
					points(h2+tau[i],m["mean",i],pch=19)
				})
		m<- res.fabc
		h2<- 0.015
		sapply(seq_len(ncol(m)),function(i)
				{
					#polygon(tau[i]+c(-h,h,h,-h), c(m["95u",i],m["95u",i],m["95l",i],m["95l",i]),border="red")
					lines(tau[i]+h2+c(-h,h) ,  c(m["95u",i],m["95u",i]),col="red")
					lines(tau[i]+h2+c(-h,h) ,  c(m["95l",i],m["95l",i]),col="red")
					lines(h2+c(tau[i],tau[i]),	c(m["95l",i],m["95u",i]), lty=2,col="red"	)
					points(h2+tau[i],m["mean",i],pch=19,col="red")
				})
		legend("topleft",bty='n',fill=c("black","red"),legend=c("Std-ABC","n-ABC"))
		dev.off()
	}
	if(1)
	{
		require(fields)
		res.fabc<- c(0.11,0.31,0.36,0.4,0.41,0.42,0.43,0.44,0.44,0.44,		0.05,0.2,0.27,0.3,0.32,0.33,0.34,0.35,0.36,0.36,			0.01,0.06,0.12,0.15,0.16,0.18,0.2,0.21,0.21,0.21,			0.005,0.01,0.04,0.06,0.08,0.09,0.1,0.11,0.11,0.12)
		res.fabc<- matrix(res.fabc, nrow=10,ncol=4, dimnames=list(c(10,20,30,40,50,60,70,80,90,100),c(1.75,1.5,1.25,1.15)))
		res.stabc<- c(0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,			0.015,0.0125,0.0175,0.015,0.0125,0.0125,0.01,0.0125,0.015,0.015,			0.0425,0.0525,0.0525,0.06,0.0575,0.055,0.055,0.06,0.06,0.055,				0.1175,0.1325,0.1375,0.1325,0.1325,0.175,0.175,0.1375,0.1375,0.135	)
		res.stabc<- matrix(res.stabc, nrow=10,ncol=4, dimnames=list(c(10,20,30,40,50,60,70,80,90,100),c(1.15,1.25,1.5,1.75)))
		zlim= range(c(res.stabc,res.fabc))

		pdf(paste(dir.name,"/GAUSS_acc_stdabc.pdf",sep=''),version="1.4",width=6,height=6)
		par(mar=c(5,4,4,4))
		m<- res.stabc
		grid<- cbind( rep( as.numeric(rownames(m)), ncol(m) ), rep( as.numeric(colnames(m)), each=nrow(m) ) )
		tmp<- as.image(as.vector(m),x=grid,nrow=nrow(m),ncol=2*ncol(m))
		tmp<- image.smooth(tmp, theta=0.2)
		image.plot(tmp, xlab=expression(sqrt(n)),ylab="threshold",legend.lab="acceptance",zlim=zlim, legend.mar=4.1)
		contour(tmp$x,tmp$y,tmp$z,add=TRUE, nlevels= 5, lwd=2,labcex=1.5,col="white", levels=c(0.1))
		dev.off()

		m<- res.fabc
		pdf(paste(dir.name,"/GAUSS_acc_fabc.pdf",sep=''),version="1.4",width=6,height=6)
		par(mar=c(5,4,4,4))
		grid<- cbind( rep( as.numeric(rownames(m)), ncol(m) ), rep( as.numeric(colnames(m)), each=nrow(m) ) )
		tmp<- as.image(as.vector(m),x=grid,nrow=nrow(m),ncol=2*ncol(m))
		tmp<- image.smooth(tmp, theta=0.2)
		image.plot(tmp, xlab=expression(sqrt(n)),ylab="threshold",legend.lab="acceptance",zlim=zlim, legend.mar=4.1)
		contour(tmp$x,tmp$y,tmp$z,add=TRUE, nlevels= 5, lwd=2,labcex=1.5,col="white", levels=c(0.2,0.3))
		dev.off()
	}
	stop()
}
#------------------------------------------------------------------------------------------------------------------------
project.nABC.movingavg.get.2D.mode<- function(x,y,xlim=NA,ylim=NA,n.hists=5,nbin=2, nlevels=5, width.infl=0.25, gridsize=c(100,100), method="kde",plot=0)
{
	if(!method%in%c("kde","ash"))	stop("project.nABC.movingavg.get.2D.mode: error at 1a")
	if(any(is.na(xlim)))	xlim<- range(x)*1.05
	if(any(is.na(ylim)))	ylim<- range(y)*1.05
	
	if(method=="kde")
	{
		require(ash)
		bins<- bin2(cbind(x, y), ab=rbind(xlim,ylim),nbin=nbin*c(nclass.Sturges(x),nclass.Sturges(y)))
		f <- ash2(bins,rep(n.hists,2))
		mxidx<- c( (which.max(f$z)-1)%%nrow(f$z)+1, (which.max(f$z)-1)%/%ncol(f$z)+1 ) #row, col
		mx<- c(		mean( f$x[ c(mxidx[1],ifelse(mxidx[1]<length(f$x),mxidx[1]+1,mxidx[1])) ] ),
					mean( f$y[ c(mxidx[2],ifelse(mxidx[2]<length(f$y),mxidx[2]+1,mxidx[2])) ] )		)
		if(plot)
		{
			image(f$x,f$y,f$z, col=head( rev(gray(seq(0,.95,len=trunc(50*1.4)))), 50))				
			contour(f$x,f$y,f$z,add=TRUE, nlevels= nlevels)
			points(mx, col="red", pch=19)
		}
	}
	else if(method=="ash")
	{
		require(KernSmooth)
		require(fields)
		x.bw<- width.infl*diff(summary(x)[c(2,5)])
		y.bw<- width.infl*diff(summary(y)[c(2,5)])
		f <- bkde2D(cbind(x, y), range.x=list(xlim,ylim), bandwidth=c(x.bw,y.bw), gridsize=gridsize)
		mxidx<- c( (which.max(f$fhat)-1)%%nrow(f$fhat)+1, (which.max(f$fhat)-1)%/%ncol(f$fhat)+1 ) #row, col	
		#print(mxidx)
		#print(f$x2[ c(mxidx[2],ifelse(mxidx[2]<length(f$x2),mxidx[2]+1,mxidx[2])) ] )
		mx<- c(		mean( f$x1[ c(mxidx[1],ifelse(mxidx[1]<length(f$x1),mxidx[1]+1,mxidx[1])) ] ),
					mean( f$x2[ c(mxidx[2],ifelse(mxidx[2]<length(f$x2),mxidx[2]+1,mxidx[2])) ] )		)
		if(plot)
		{
			image(f$x1,f$x2,f$fhat, col=head( rev(gray(seq(0,.95,len=trunc(50*1.4)))), 50))
			contour(f$x1, f$x2, f$fhat, nlevels= nlevels, add=1)			
		}
	}	
	mx
}
#------------------------------------------------------------------------------------------------------------------------
project.nABC.movingavg.gethist<- function(x, theta, nbreaks= 20, width= 0.5, plot=0, rtn.dens=0)
{
	#compute break points sth theta is in the middle
	breaks<- max(abs( theta - x ))*1.1 / nbreaks								
	breaks<- c( rev(seq(from= theta-breaks/2, by=-breaks, length.out= nbreaks )), seq(from= theta+breaks/2, by=breaks, length.out= nbreaks ) )
	ans.h<- hist(x,breaks=breaks, plot= 0)
	ans.h[["mean"]]<- mean(x)		
	ans.h[["hmode"]]<- mean(ans.h[["breaks"]][seq(which.max(ans.h[["counts"]]),length.out=2)])	
	tmp<- density(x, kernel="biweight",from=breaks[1],to=breaks[length(breaks)],width = max(EPS,width*diff(summary(x)[c(2,5)])))
	ans.h[["dmode"]]<- tmp[["x"]][which.max( tmp[["y"]])]
	if(rtn.dens)
		ans.h[["dens"]]<- tmp
	if(plot)
	{
		plot(ans.h, freq=0)
		lines(tmp)
	}
	ans.h
}
#------------------------------------------------------------------------------------------------------------------------
project.nABC.movingavg.a2rho<- function(x) 
{
	atanh( x/(1+x*x) )
}
#------------------------------------------------------------------------------------------------------------------------
project.nABC.movingavg.rho2a<- function(x)
{
	x<- tanh(x)
	x[x>.5]<- .5
	x[x<-.5]<- -.5
	(1-sqrt(1-4*x*x))/(2*x)
} 
#------------------------------------------------------------------------------------------------------------------------
project.nABC.movingavg.estimateTheta0<- function(m, theta.names,links.names )
{
	verbose<- 1
	lsets<- lapply(links.names,function(x)
			{
				#nABC.getlevelset.2d(m, x, theta.names, rho.eq=0, rho.eq.sep=35, rho.eq.q=0.05, theta.sep=250, plot=0, method="quantile", verbose=verbose)
				nABC.getlevelset.2d(m, x, theta.names, rho.eq=0, rho.eq.sep=15, rho.eq.q=0.005, theta.sep=250, plot=0, method="fixed", verbose=verbose)
			})
	names(lsets)<- links.names
	#print(lsets)
	theta.intersection<- nABC.getlevelsetintersection.2d(lsets,theta.names, 50, plot=0, verbose=verbose)				
	cat(paste("\nfinal number of theta in intersection",ncol(theta.intersection),"\n"))
	apply(theta.intersection,1,mean)
}
#------------------------------------------------------------------------------------------------------------------------
project.nABC.movingavg.detJac<- function(a,sig2,ax,vx)
{
	( 1-a^4 ) / (( 1 + a*a + a^4 )*vx )
}
#------------------------------------------------------------------------------------------------------------------------
project.nABC.movingavg.avgdetJac<- function(a.tau.l, a.tau.u, ax, vx, s, alpha, empirical.rho= NULL, empirical.links=NULL, plot=0)
{
	rhox<- 	project.nABC.movingavg.a2rho(ax)
	#rhox<- ax	#the ax stored is actually hat(rho)_x
	#ax	<- project.nABC.movingavg.rho2a(rhox)
#print(c(rhox,ax,vx))

	sigma2<- 1
	rho<- seq(rhox+a.tau.l, rhox+a.tau.u, by=0.001)	
	a<- seq(project.nABC.movingavg.rho2a(rhox+a.tau.l),project.nABC.movingavg.rho2a(rhox+a.tau.u), by=0.001)
#print(range(a))
	detjac<- abs(project.nABC.movingavg.detJac( a, sigma2, ax, vx))	
	ax.idx	<- which.min(abs(a-ax))
	#print(a[ax.idx])	
	#print(mean(detjac[1:ax.idx]-detjac[ax.idx]))
	#print(mean(detjac[ax.idx:length(detjac)]-detjac[ax.idx]))	
	ans		<- c(ax, mean(detjac[ax.idx:length(detjac)]-detjac[ax.idx])-mean(detjac[1:ax.idx]-detjac[ax.idx]))	
	#pw		<- nabc.acf.equivalence.pow(rho, a.tau.u, alpha, s)
	dens.a	<- nabc.acf.equivalence.pow(project.nABC.movingavg.a2rho(a)-rhox, a.tau.u, alpha, s)
	dens.a	<- dens.a * detjac
	dens.a	<- dens.a / sum(dens.a)		
	ans		<- c(ans, a[which.max(dens.a)]-ax )
	names(ans)<- c("ax","mean","pow")
	if(!is.null(empirical.links))
	{
		detjac.j<- nabc.estimate.jac( empirical.links, cbind(th1=a, th2=rep(sigma2,length(a))), c(2e-3, 2e-2), ax, vx )
		#print(detjac.j[150:160]); print(detjac[150:160]); print(detjac.j-detjac)

		dens.a.j<- nabc.acf.equivalence.pow(project.nABC.movingavg.a2rho(a)-rhox, a.tau.u, alpha, s)
		dens.a.j<- dens.a.j * detjac.j
		dens.a.j<- dens.a.j / sum(dens.a.j, na.rm=1)		
		tmp			<- names(ans)
		ans			<- c(ans, a[which.max(dens.a.j)]-ax )
		names(ans)<- c(tmp,"lhatpow")
	}
	if(!is.null(empirical.rho))
	{
		#empirical.rho.d	<- density(empirical.rho, kernel="biweight",width = max(EPS,2*diff(summary(empirical.rho)[c(2,5)])))
		dens.a.e	<- approx(empirical.rho$x,empirical.rho$y,xout=project.nABC.movingavg.a2rho(a)-rhox,rule=2)$y
		dens.a.e	<- dens.a.e * detjac
		dens.a.e	<- dens.a.e / sum(dens.a.e)		
		tmp			<- names(ans)
		ans			<- c(ans, a[which.max(dens.a.e)]-ax )
		names(ans)<- c(tmp,"empirical")
	}		
	#plot<- 1
	if(plot)
	{	
		#plot(empirical.links[[2]])
		#stop()
		plot(1,1,type='n',xlim=range(a),ylim=range(c(detjac.j,detjac), na.rm=1),xlab="blar")
		lines(a,detjac.j,lty=4)
		lines(a,detjac,lty=2)
		stop()
		plot(a,dens.a,type='l', ylim=range(dens.a,dens.a.e,dens.a.j, na.rm=1), lty=2)
		lines(a,dens.a.e,lty=3)
		lines(a,dens.a.j,lty=4)
		abline(v=ax)
		abline(v=a[which.max(dens.a)],lty=2)
		abline(v=a[which.max(dens.a.e)],lty=3)
		abline(v=a[which.max(dens.a.j)],lty=4)
		stop()
	}
	ans		
}
#------------------------------------------------------------------------------------------------------------------------
project.nABC.changeofvariables<- function()
{
	dir.name<- paste(DATA,"nABC.cov",sep='/')
	my.mkdir(DATA,"nABC.cov")
	
	cv.link<- function(x,a=0.2){  tmp<- x; tmp[x<=a]<- 2*x[x<=a]; tmp[x>a]<- 0.5*x[x>a]+2*a-0.5*a; tmp	}	
	cov.gethist<- function(a, r.n, r.sigma)
	{
		rho<- rnorm(r.n,0, r.sigma)	
		theta<- cv.link(rho,a=a)
		tmp<- hist(theta, breaks=100, plot=0)
		tmp
	}
	
	if(0)
	{
		r.sigma<- c(0.05, 0.2, 0.5)
		theta.h<- lapply(r.sigma, function(x)		cov.gethist(0.2,5e6,x)		)
		plot(theta.h[[1]], freq=0)
		plot(theta.h[[2]], freq=0)
	}
	if(1)
	{
		a0<- 0.1
		sigma20<- 1		
		rho0<- 	project.nABC.movingavg.a2rho(a0)
		
		#plot det of Jacobian for MA(1) process
		a<- seq(-0.1,0.3,0.001)
		sigma2<- 1
		cov<- project.nABC.movingavg.detJac(a, sigma2, a0, sigma20)
		#plot(a,cov,type='l')

		#plot difference in det(J) for increasing tau, assuming that most of the mass in 
		#the power is withint -tau/2,tau/2
		tau.u<- seq(0.01,0.2,0.005)
		
		cov<- sapply(tau.u, function(x)
				{					
					project.nABC.movingavg.detJac(c(project.nABC.movingavg.rho2a(rho0-x/2), a0, project.nABC.movingavg.rho2a(rho0+x/2)), sigma2, a0, sigma20)	
				})				
		cov<- rbind(cov[1,]-cov[2,], cov[2,]-cov[3,])		
		#get the max difference with increasing tau
		cov<- sapply(seq_len(ncol(cov)),function(i)
				{
					apply(cov[,1:i,drop=0],1,max)
				})
		colnames(cov)<- tau.u
		print(cov)
		print(cov[1,]+cov[2,])
	}
	if(0)
	{
		a0<- 0.2
		r.n<- 5e6
		r.sigma<- 0.3
		r.sigmas<- seq(0.05,0.2,0.025)
		tau.u<- 0.2
		verbose<- 0
		resume<- 0
		dir.name<- paste(DATA,"nABC.cov",sep='/')
		
		#collect differences in differentials and modes 
		cat(paste("\nnABC.cov dirname is",dir.name))
		f.name<- paste(dir.name,paste("nABC.cov_MA1_",a0,"_",tau.u,"_",r.n,".R",sep=''),sep='/')
		options(show.error.messages = FALSE, warn=1)		
		readAttempt<-try(suppressWarnings(load(f.name)))						
		options(show.error.messages = TRUE)		
		if(!resume || inherits(readAttempt, "try-error"))
		{
			cov<- sapply(r.sigmas,function(r.sigma)
			{
				#r.sigma<- 0.2
				#get link function
				rho0<- 	project.nABC.movingavg.a2rho(a0)		#rho(x)
				tau.u<- 1*r.sigma
				cat(paste("\nprocess r.sigma",r.sigma,"tau.u",tau.u))
				if(verbose) cat(paste("\ntrue rho0",rho0))		
				if(verbose) plot( seq(-0.423,0.423,0.001), project.nABC.movingavg.rho2a(seq(-0.423,0.423,0.001)), type='l')
				
				#evaluate difference quotient at boundary of interval hypothesis
				a.differential<- project.nABC.movingavg.rho2a(c(rho0-tau.u,rho0,rho0+tau.u))
				a.differential<- c(-tau.u/(a.differential[2]-a.differential[1]), -tau.u/(a.differential[3]-a.differential[2]))
				a.differential<- c(a.differential,diff(a.differential))
		
				#estimate mode of 'rho' for power centering on rho0 
				rho<-  	rnorm(r.n,rho0, r.sigma)					#rho		
				rho<-	rho[-0.423<=rho & rho<=0.423 & (rho0-tau.u)<=rho & rho<=(rho0+tau.u)]
				rho.h<- project.nABC.movingavg.gethist(rho, rho0, nbreaks= 100)
				if(verbose) plot(rho.h)
				if(verbose) cat(paste("\ndmode of rho",rho.h[["dmode"]]))
				if(verbose) cat(paste("\nlength of rho",length(rho)))
				
				#estimate mode of 'a' for power centering on rho0
				a<- project.nABC.movingavg.rho2a( rho )
				a.h<- project.nABC.movingavg.gethist(a, a0, nbreaks= 100)
				if(verbose) plot(a.h, col=my.fade.col("black",0.3),border=NA)
				if(verbose) abline(v=a0,lty=2)
				if(verbose) cat(paste("\ndmode of a",a.h[["dmode"]]))				
				
				ans<- c(a.differential,rho0, rho.h[["dmode"]],a0,a.h[["dmode"]])
				names(ans)<- c("diff.l","diff.u","diff.d","rho0","rho.dmode","a0","a.dmode")
				ans
			})
			colnames(cov)<- r.sigmas
			f.name<- paste(dir.name,paste("nABC.cov_MA1_",a0,"_",tau.u,"_",r.n,".R",sep=''),sep='/')
			cat(paste("\nsave 'cov' to",f.name))
			save(cov,file=f.name)	
		}
		print(cov)
	}
	stop()
}
#------------------------------------------------------------------------------------------------------------------------
project.nABC.movingavg.get.fixed.ts<- function(n,mu,a,sd, verbose=0, tol=1e-3)
{		
	verbose<- 1
	m		<- ceiling( max(1, n / 5000) )
	m		<- 1e4 / m
	rho0	<- project.nABC.movingavg.a2rho(a)
	error	<- 1
	while(error>tol)
	{
		x<-	rnorm(m*n+1,mu,sd=sd)
		x<- x[-1] + x[-(m*n+1)]*a
		x<- matrix(x,ncol=m)		
		x<- sapply(seq_len(ncol(x)),function(i)		x[,i]/sd(x[,i])*sqrt((1+a*a)*sd*sd)	)		#get correct variance
		#rhox	<- apply(x,2,function(col)	atanh(acf(ts(col), plot=0,lag=1)[["acf"]][2,1,1]) )
		rhox	<- apply(x,2,function(col)	nabc.acf.equivalence.cor(col, leave.out=2)["z"] )				
		error	<- abs(rhox - rho0 )
		ans		<- which.min(error)
		error	<- error[ans]
		if(verbose)	cat(paste("\ncurrent error is",error,"and should be <",tol))						
	}	
	ans<- x[,ans]
	if(verbose)	cat(paste("\nrhox of time series is",nabc.acf.equivalence.cor(ans, leave.out=2)["z"],"and should be",rho0))	
	if(verbose)	cat(paste("\nvar of time series is",var(ans),"and should be",(1+a*a)*sd*sd))
	ans
}
#------------------------------------------------------------------------------------------------------------------------
project.nABC.movingavg<- function()
{
	my.mkdir(DATA,"nABC.movingavg_mode")
	dir.name<- paste(DATA,"nABC.movingavg_mode",sep='/')
	subprog<- 1
	resume<- 1
	pdf.width<- 4
	pdf.height<-5
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,5),
									subp= return(as.numeric(substr(arg,6,nchar(arg)))),NA)	}))
		if(length(tmp)>0) subprog<- tmp[1]
	}
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,4),
									res= return(as.numeric(substr(arg,5,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
	}
	
	project.nABC.movingavg.a2nu<- function(x) x/(1+x*x) 
	project.nABC.movingavg.nu2a<- function(x)
	{
		x[x>.5]<- .5
		x[x<-.5]<- -.5		
		(1-sqrt(1-4*x*x))/(2*x)
	} 
	project.nABC.movingavg.getaccr<- function(c,x,acc,method="sym")
	{
		if(!method%in%c("sym","inv"))	stop("project.nABC.movingavg.getaccr: error at 1a")
		if(method=="sym")
			return(which(x<=c & x>=-c)/length(x)-acc)
		else
			return(which(x<=c & x>=1/c)/length(x)-acc)		
	}
	fstretch.solvefor.tau<- function(tau.u, cir, n, alpha) pf( cir / tau.u,  n, n) - pf( 1 / (cir*tau.u),  n, n) - alpha
	
	
	project.nABC.movingavg.fixx.unifsigma.unifma<- function(	N,alpha,ma.n,x,
																xmapa.tau.l,xmapa.tau.u,xmapa.prior.l,xmapa.prior.u,
																xsig2.tau.l,xsig2.tau.u,xsig2.prior.l,xsig2.prior.u, 
																xmapa.leave.out=0, xsig2.leave.out=0	)
	{
		if(xmapa.tau.u<0 || xmapa.tau.l>0)	stop("project.nABC.movingavg.fixsigma.unifma: error at 1a")
		if(xsig2.tau.u<1 || xsig2.tau.l>1)	stop("project.nABC.movingavg.fixsigma.unifma: error at 1b")									
		ans<- vector("list",7)
		names(ans)<- c("xv","xa","v.cil","v.cir","a.cil","a.cir","data")		
		ans[["xv"]]<- 	var(x[seq.int(1,length(x),by=1+xsig2.leave.out)])
		ans[["xa"]]<-	nabc.acf.equivalence.cor(x, leave.out=xmapa.leave.out)["z"]
		tmp<- .Call("abcScaledChiSq",	c(	floor(length(x) / (1+xsig2.leave.out)) - 1, floor(length(x) / (1+xsig2.leave.out)) - 1,	xsig2.tau.l,	xsig2.tau.u,	alpha,1e-10,100,0.05)	)			
		ans[["v.cil"]]<- tmp[1]
		ans[["v.cir"]]<- tmp[2]
		tmp<- nabc.acf.equivalence(x, x, args=paste("acfequiv",xmapa.leave.out,xmapa.tau.l,xmapa.tau.u,alpha,sep='/')	)
		ans[["a.cil"]]<- tmp["al"]
		ans[["a.cir"]]<- tmp["ar"]	
		ans[["data"]]<- sapply(1:N,function(i)
				{
					#cat(paste("\nproject.nABC.movingavg.unifsigma.unifma iteration",i))					
					ymapa<- project.nABC.movingavg.rho2a( runif(1, project.nABC.movingavg.a2rho( xmapa.prior.l ), project.nABC.movingavg.a2rho( xmapa.prior.u )) )					
					ysigma2<- runif(1,xsig2.prior.l, xsig2.prior.u)
					y<-	rnorm(ma.n+1,0,sd=sqrt(ysigma2))
					y<- y[-1] + y[-(ma.n+1)]*ymapa
					tmp<- nabc.acf.equivalence.cor(y, leave.out=xmapa.leave.out)				
					out.a<- c(ymapa, project.nABC.movingavg.a2rho(ymapa), tmp["z"],	( tmp["z"] - ans[["xa"]] )*sqrt(tmp["n"]-3))
					names(out.a)<- c("a.theta","a.link","a.link.mc","a.error")
#print(tmp["n"]); print(ans[["xa"]]); print(tmp); print(( tmp["z"] - ans[["xa"]] )*sqrt(tmp["n"]-3)); print(out.a); stop()					
					out.v<- c(ysigma2,	(1+ymapa*ymapa)*ysigma2,	var(y[seq.int(1,length(y),by=1+xsig2.leave.out)])/ans[["xv"]])
					names(out.v)<- c("v.theta","v.link","v.error")
					c(out.a,out.v)					
				})							
		ans
	}	
	
			
	
	if(0)	#other potential summaries
	{
		N<- 1e4
		alpha<- 0.01
		nbreaks<- 10	
		ma.n<- 1e4
		xma.pa<- 0.4		
		xsig2<- 1
		yma.pa<- 0		
		ysig2<- 1
		zma.pa<- -0.4		
		zsig2<- 1
		
		other.sus<- function(n,a,s){
			x<-	rnorm(n+1,0,s)
			x<- x[-1] + x[-(n+1)]*a
			y<- diff(x)
			y<- y[seq.int(1,n-1,by=2)]			
			z<- x[seq.int(1,n,by=2)]
			c(sd(x),sd(y),sd(z))
		}
		 
		ans1<- replicate(N,other.sus(ma.n,xma.pa,xsig2))		
		ans2<- replicate(N,other.sus(ma.n,yma.pa,ysig2))
		ans3<- replicate(N,other.sus(ma.n,zma.pa,zsig2))
		
		print(apply(ans1,1,mean))
		print(apply(ans2,1,mean))
		print(apply(ans3,1,mean))
		stop()
	}
	if(!is.na(subprog) && subprog==1)		#perform nABC repeatedly for small enough tau.u(a)
	{
		N<- 2e6
		M<- 2e3
		m<- NA
		if(exists("argv"))
		{
			tmp<- na.omit(sapply(argv,function(arg)
							{	switch(substr(arg,2,2),
										m= return(as.numeric(substr(arg,3,nchar(arg)))),NA)	}))
			if(length(tmp)>0) m<- tmp[1]
		}
		xa<- 0.1 
		r.xa<- project.nABC.movingavg.a2nu(xa)		#r for xa
		z.xa<- project.nABC.movingavg.a2rho(xa)		#r for xa
		xsigma2<- 1	#sqrt(2)
		alpha<- 0.01
		nbreaks<- 20		
		xn<- 5e3		
		
		tau.u<- 0.1
		tau.l<- -tau.u
		xsig2.tau.u<- 1.1
		xsig2.tau.l<- 1/xsig2.tau.u
		xsig2.prior.u<- 1.15 
		xsig2.prior.l<- 0.8
		
		verbose<- 1
		if(verbose)	cat(paste("true xmapa, correlation scale",r.xa,"true xmapa, test scale",z.xa,"\n"))						
		prior.u<- project.nABC.movingavg.rho2a( .423 )	#project.nABC.movingavg.rho2a( z.xa+tau.u )		
		prior.l<- project.nABC.movingavg.rho2a( -.423 )	#project.nABC.movingavg.rho2a( z.xa+tau.l )
		if(verbose)	cat(paste("prior mapa thresholds from test scale",prior.l,prior.u,"\n"))
		if(verbose)	cat(paste("sym prior mapa thresholds from test scale",project.nABC.movingavg.rho2a( z.xa+tau.l ),project.nABC.movingavg.rho2a( z.xa+tau.u ),"\n"))

		if(!is.na(m))
		{		
			f.name<- paste(dir.name,"/nABC.MA1_ok_",N,"_",xn,"_",round(prior.l,d=2),"_",round(prior.u,d=2),"_",round(tau.u,d=2),"_",round(xsig2.prior.l,d=2),"_",round(xsig2.prior.u,d=2),"_",round(xsig2.tau.u,d=2),"_m",m,".R",sep='')
			cat(paste("\nnABC.MA: compute ",f.name))
			options(show.error.messages = FALSE, warn=1)		
			readAttempt<-try(suppressWarnings(load(f.name)))						
			options(show.error.messages = TRUE)						
			if(!resume || inherits(readAttempt, "try-error"))
			{
				x<-	rnorm(xn+1,0,sd=sqrt(xsigma2))
				x<- x[-1] + x[-(xn+1)]*xa				
				ans<- project.nABC.movingavg.fixx.unifsigma.unifma(	N,alpha,xn,x,
																	tau.l,tau.u,prior.l,prior.u,
																	xsig2.tau.l,xsig2.tau.u,xsig2.prior.l,xsig2.prior.u,		
																	xmapa.leave.out=2, xsig2.leave.out=1)					
				cat(paste("\nnABC.MA: save ",f.name))
				save(ans,file=f.name)				
			}
			else
				cat(paste("\nnABC.MA: resumed ",f.name))
		}	
		else
		{
			#collect all abc runs, estimate 2D mode 
			cat(paste("\nnABC.MA",dir.name))
			save.f.name<- paste(dir.name,"/nABC.MA1_modemean_",N,"_",xn,"_",round(prior.l,d=2),"_",round(prior.u,d=2),"_",round(tau.u,d=2),"_",round(xsig2.prior.l,d=2),"_",round(xsig2.prior.u,d=2),"_",round(xsig2.tau.u,d=2),".R",sep='')
			options(show.error.messages = FALSE, warn=1)		
			readAttempt<-try(suppressWarnings(load(save.f.name)))						
			options(show.error.messages = TRUE)
			
			#tau.u<- 0.097
			#tau.l<- -tau.u 
			#xsig2.tau.u<- 1.047
			#xsig2.tau.l<- nabc.chisqstretch.tau.low(xsig2.tau.u, floor(5e3 / 2) - 1, alpha)
			#print(xsig2.tau.l)
			#xsig2.tau.l<- 1/xsig2.tau.u																
			#v.rej<- .Call("abcScaledChiSq",	c(	floor(xn / 2) - 1,	floor(xn / 2) - 1, xsig2.tau.l,	xsig2.tau.u,	alpha,1e-10,100,0.05)	)			
			#a.rej<- nabc.acf.equivalence(rnorm(xn), rnorm(xn), args=paste("acfequiv",3,tau.l,tau.u,alpha,sep='/')	)
			
			if(!resume || inherits(readAttempt, "try-error"))
			{
				#accept if both SD and ACF ok
				cat(paste("\nnABC.MA generate",save.f.name))
				f.name<- list.files(dir.name, pattern=paste("^nABC.MA1_ok_",'',sep=''), full.names = TRUE)
				xa.symu<- project.nABC.movingavg.rho2a( z.xa+tau.u )
				xa.syml<- project.nABC.movingavg.rho2a( z.xa+tau.l )
				modes<- sapply(seq_along(f.name),function(i)
						{
							cat(paste("\nnABC.MA load",f.name[i]))
							readAttempt<-try(suppressWarnings(load( f.name[i] )))
							if(inherits(readAttempt, "try-error"))	return(rep(NA,6))
							out<- c(ans[["xa"]],ans[["xv"]])
							if(0)
							{
								ans[["v.cil"]]<- v.rej[1]
								ans[["v.cir"]]<- v.rej[2]								
								ans[["a.cil"]]<- a.rej["al"]
								ans[["a.cir"]]<- a.rej["ar"]								
							}														
							#ans.sd
							acc<- which( ans[["data"]]["v.error",]<=ans[["v.cir"]]  &  ans[["data"]]["v.error",]>=ans[["v.cil"]] )
							tmp<- project.nABC.movingavg.get.2D.mode(ans[["data"]]["a.theta",acc],ans[["data"]]["v.theta",acc], xlim= c(-0.4,0.4),ylim=c(0.8,1.2),plot=0, nbin=10, nlevels=20)
							#abline(h=ans[["xv"]],col="red"); abline(v=ans[["xa"]],col="red")
							out<- c(out,tmp)
							#ans.sdacf
							acc<- which( 	ans[["data"]]["a.error",]<=ans[["a.cir"]]  &  ans[["data"]]["a.error",]>=ans[["a.cil"]]	& 
											ans[["data"]]["v.error",]<=ans[["v.cir"]]  &  ans[["data"]]["v.error",]>=ans[["v.cil"]] )
							tmp<- project.nABC.movingavg.get.2D.mode(ans[["data"]]["a.theta",acc],ans[["data"]]["v.theta",acc], xlim= c(-0.2,0.4),ylim=c(0.8,1.2), plot=0, nbin=10, nlevels=20)
							#tmp<- project.nABC.movingavg.get.2D.mode(ans[["data"]]["a.link",acc],ans[["data"]]["v.link",acc], xlim= c(-0.2,0.4),ylim=c(0.8,1.2), plot=1, nbin=10, nlevels=20)
							#abline(h=ans[["xv"]],col="red"); abline(v=ans[["xa"]],col="red")
							out<- c(out,tmp)
							names(out)<- c("xa","xv","ya.dmode.sd","yv.dmode.sd","ya.dmode.sdacf","yv.dmode.sdacf")
							#stop()
							if(0)
							{
								print(ans[["xa"]])
								print(ans[["xv"]])
								print(ans[["data"]][,1:5])
								print(length(acc)/ncol(ans[["data"]]))
								print(ncol(ans[["data"]]))
								tmp<- project.nABC.movingavg.gethist(ans[["data"]]["a.link",acc], ans[["xa"]], nbreaks= 50, width= 0.5, plot=1)
								print(tmp[["dmode"]])
								tmp<- project.nABC.movingavg.gethist(ans[["data"]]["v.link",acc], ans[["xv"]], nbreaks= 50, width= 0.5, plot=1)
								print(tmp[["dmode"]])
								stop()
								print(range(ans[["data"]]["a.theta",acc]))
								print(range(ans[["data"]]["v.theta",acc]))								
							}
							out
						})					
				cat(paste("\nnABC.MA: save file ",save.f.name))
				save(modes,file=save.f.name)
				stop()
			}
			
			print(modes[,1:5])			
			modes<- modes[, !is.na(modes[1,]) ]
			cat(paste("\nlength of ABC repetitions is",ncol(modes)))
			ax		<- project.nABC.movingavg.rho2a(modes["xa",])
			sig2x	<- modes["xv",]/(1+project.nABC.movingavg.rho2a(modes["xa",])^2)
			err		<- 	sqrt(	(modes["ya.dmode.sdacf",]-ax)^2	+	(modes["yv.dmode.sdacf",]-sig2x)^2	)
			cat(paste("\n mean v1",mean(modes["xv",])))
			cat(paste("\n mean v2",mean(modes["xa",])))
			cat(paste("\n mean ax",mean(ax)))			
			cat(paste("\n mean sig2x",mean(sig2x)))
			cat(paste("\n mean mode of a",mean(modes["ya.dmode.sdacf",])))
			cat(paste("\n mean mode of sig2",mean(modes["yv.dmode.sdacf",])))					
			cat(paste("\n mean mode-(ax,sig2x)",mean(err)))
#bookmark_tableMA			
			stop()
			require(ash)
			xlim<- c(-0.1,0.2)
			ylim<- c(0.8,1.2)
			nbin<- c(100,100)			
			bins<- 		bin2(cbind(modes["ya.dmode.sdacf",],modes["yv.dmode.sdacf",]),nbin=c(10,10))
			f.sdacf<- 	ash2(bins,rep(5,2))
			bins<- 		bin2(cbind(modes["ya.dmode.sd",],modes["yv.dmode.sd",]), ab=rbind(xlim,ylim),nbin=nbin)
			f.sd<- 		ash2(bins,rep(5,2))
			bins<- 		bin2(cbind(modes["xa",],modes["xv",]), ab=rbind(xlim,ylim),nbin=nbin)
			f.x<- 		ash2(bins,rep(5,2))
			
			f.name<- paste(dir.name,"/nABC.MA1_modes_sdacf_vs_acf_",N,"_",xn,"_",m,".pdf",sep='')
			#cat(paste("\nnABC.MA: write pdf to",f.name))
			pdf(f.name,version="1.4",width=pdf.width,height=pdf.height)			
			par(mar=c(4.25,4.25,1,1))
			cols<- c(my.fade.col("black",1),my.fade.col("black",0.2),my.fade.col("black",0.6))
			ltys<- c(2,3,4)
			#plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim,xlab=expression("estimated mode of "*a),ylab=expression("estimated mode of "*sigma^2))
			#points(modes["xa",],modes["xv",],col=cols[1],pch=22)
			image(f.x$x,f.x$y,f.x$z, col=head( rev(gray(seq(0,1,len=trunc(50*1.4)))), 50), ,xlab=expression("estimated mode of "*a),ylab=expression("estimated mode of "*sigma^2))
			#contour(f.x$x,f.x$y,f.x$z, add=1, nlevels= 5,col=cols[1],lty=ltys[1],drawlabels=0 )
			
			#image(f.sd$x,f.sd$y,f.sd$z, col=cols,xlab=expression("estimated mode of "*a),ylab=expression("estimated mode of "*sigma^2))
			#points(modes["ya.dmode.sd",],modes["yv.dmode.sd",],col=cols[2],pch=19)
			contour(f.sd$x,f.sd$y,f.sd$z, add=1, nlevels= 4,col=cols[2],lty=ltys[1],drawlabels=0, lwd=1.5 )		
			#image(f.sdacf$x,f.sdacf$y,f.sdacf$z, col=cols,add=1)
			#points(modes["ya.dmode.sdacf",],modes["yv.dmode.sdacf",],col=cols[3],pch=19)
			contour(f.sdacf$x,f.sdacf$y,f.sdacf$z, add=1, nlevels= 5,lty=ltys[3],col=cols[3],drawlabels=1, lwd=1.5 )
			legend("topleft", bty='n',box.col="white",border=NA, lty=c(ltys[1],ltys[3]),fill=c(cols[2],cols[3]),legend=expression("only "*nu[1],"both "*nu[1]*" and "*nu[2]))
			abline(h=1,lty=3)
			abline(v=0.1,lty=3)
			dev.off()			
		}			
		stop()
	}
	if(!is.na(subprog) && subprog==2)			#large n
	{
		N<- 2e6
		xn<- NA
		if(exists("argv"))
		{
			tmp<- na.omit(sapply(argv,function(arg)
							{	switch(substr(arg,2,3),
										xn= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
			if(length(tmp)>0) xn<- tmp[1]
		}
		if(exists("argv"))
		{
			tmp<- na.omit(sapply(argv,function(arg)
							{	switch(substr(arg,2,2),
										N= return(as.numeric(substr(arg,3,nchar(arg)))),NA)	}))
			if(length(tmp)>0) N<- tmp[1]
		}
		xa<- 0.1 
		r.xa<- project.nABC.movingavg.a2nu(xa)		#r for xa
		z.xa<- project.nABC.movingavg.a2rho(xa)		#r for xa
		xsigma2<- 1	#sqrt(2)
		alpha<- 0.01
		nbreaks<- 20		
		xmapa.leave.out<- 2
		xsig2.leave.out<- 1
		mx.pw<- 0.9
		
		tau.u<- 0.1
		tau.l<- -tau.u
		prior.u<- project.nABC.movingavg.rho2a( z.xa+.199 )		
		prior.l<- project.nABC.movingavg.rho2a( z.xa-.199 )
		
		xsig2.tau.u<- 1.1
		xsig2.tau.l<- 1/xsig2.tau.u
		xsig2.prior.u<- 1.7 
		xsig2.prior.l<- 1/xsig2.prior.u
		
		verbose<- 1
		if(verbose)	cat(paste("true xmapa, correlation scale",r.xa,"true xmapa, test scale",z.xa,"\n"))						
		if(verbose)	cat(paste("prior mapa thresholds from test scale",prior.l,prior.u,"\n"))
		if(verbose)	cat(paste("sym prior mapa thresholds from test scale",project.nABC.movingavg.rho2a( z.xa+tau.l ),project.nABC.movingavg.rho2a( z.xa+tau.u ),"\n"))
		
		if(!is.na(xn))
		{		
			f.name<- paste(dir.name,"/nABC.MA1_largensimu_",N,"_",xn,"_",round(prior.l,d=2),"_",round(prior.u,d=2),"_",round(tau.u,d=2),"_",round(xsig2.prior.l,d=2),"_",round(xsig2.prior.u,d=2),"_",round(xsig2.tau.u,d=2),".R",sep='')
			cat(paste("\nnABC.MA: compute ",f.name))
			options(show.error.messages = FALSE, warn=1)		
			readAttempt<-try(suppressWarnings(load(f.name)))						
			options(show.error.messages = TRUE)						
			if(!resume || inherits(readAttempt, "try-error"))
			{
				x<- project.nABC.movingavg.get.fixed.ts(xn,0,xa,xsigma2, verbose=0)
				ans<- project.nABC.movingavg.fixx.unifsigma.unifma(	N,alpha,xn,x,
																	tau.l,tau.u,prior.l,prior.u,
																	xsig2.tau.l,xsig2.tau.u,xsig2.prior.l,xsig2.prior.u,		
																	xmapa.leave.out=xmapa.leave.out, xsig2.leave.out=xsig2.leave.out)					
				cat(paste("\nnABC.MA: save ",f.name))
				save(ans,file=f.name)				
			}
			else
				cat(paste("\nnABC.MA: resumed ",f.name))
		}	
		else
		{
			#collect all abc runs, estimate 2D mode 
			cat(paste("\nnABC.MA",dir.name))
			save.f.name<- paste(dir.name,"/nABC.MA1_largen_",N,"_",round(prior.l,d=2),"_",round(prior.u,d=2),"_",round(tau.u,d=2),"_",round(xsig2.prior.l,d=2),"_",round(xsig2.prior.u,d=2),"_",round(xsig2.tau.u,d=2),".R",sep='')
			options(show.error.messages = FALSE, warn=1)		
			readAttempt<-try(suppressWarnings(load(save.f.name)))						
			options(show.error.messages = TRUE)
			
			#set up fixed tau			
			xn<- 500
			fixed.tau<- matrix(NA,2,4,dimnames=list(c("a","v"),c("tau.l","tau.u","c.l","c.u")))						
			fixed.tau["a",1:2]<- 	nabc.acf.equivalence.tau.lowup(mx.pw, 0.35, floor(xn / (1+xmapa.leave.out)), alpha)[1:2]		
			fixed.tau["a",3:4]<-	nabc.acf.equivalence.abctol(fixed.tau["a","tau.l"], fixed.tau["a","tau.u"], floor(xn / (1+xmapa.leave.out)), alpha)			
			fixed.tau["v",1:2]<-	nabc.chisqstretch.tau.lowup(mx.pw, 2, floor(xn / (1+xsig2.leave.out))-1, alpha)[1:2]
			fixed.tau["v",3:4]<-	.Call("abcScaledChiSq",	c(	floor(xn / (1+xsig2.leave.out))-1, floor(xn / (1+xsig2.leave.out))-1,	fixed.tau["v","tau.l"], fixed.tau["v","tau.u"],	alpha,1e-10,100,0.05)	)[1:2]
			if(!resume || inherits(readAttempt, "try-error"))
			{
				#accept if both SD and ACF ok
				cat(paste("\nnABC.MA generate",save.f.name))
				f.name<- list.files(dir.name, pattern=paste("^nABC.MA1_largensimu_",'',sep=''), full.names = TRUE)
				f.name.n<- sort(as.numeric(sapply(strsplit(f.name,'_'), function(x)	x[length(x)-6])), index.return=1)
				f.name<- f.name[f.name.n$ix]
				ans<- lapply(seq_along(f.name),function(i)
						{
							out<- matrix(NA,2,16,dimnames=list(c("fx.tau.u","fx.pw"),c("acf.tau.l","acf.tau.u","acf.cl","acf.cu","sd.tau.l","sd.tau.u","sd.cl","sd.cu","acc","a.hmode","a.hmode.uniest","v.hmode","acf.hmode","sd.hmode","xa","xv")))
							
							cat(paste("\nnABC.MA load",f.name[i]))
							readAttempt<-try(suppressWarnings(load( f.name[i] )))							
							if(inherits(readAttempt, "try-error"))	return(out)
							xn<- f.name.n$x[i]				
							
							out[2,c("xa","xv")]<- out[1,c("xa","xv")]<- c(project.nABC.movingavg.rho2a( ans[["xa"]] ),	ans[["xv"]])
							#consider fixed tau							
							out[1,c("acf.cl","acf.cu")]		<- nabc.acf.equivalence.abctol(fixed.tau["a","tau.l"],fixed.tau["a","tau.u"], floor(xn / (1+xmapa.leave.out)), alpha)
							out[1,c("sd.cl","sd.cu")]		<- .Call("abcScaledChiSq",	c(	floor(xn / (1+xsig2.leave.out))-1, floor(xn / (1+xsig2.leave.out))-1,	fixed.tau["v","tau.l"], fixed.tau["v","tau.u"],	alpha,1e-10,100,0.05)	)[1:2]
							acc<- which( 	ans[["data"]]["a.error",]<=out[1,"acf.cu"]  &  ans[["data"]]["a.error",]>=out[1,"acf.cl"]	& 
											ans[["data"]]["v.error",]<=out[1,"sd.cu"]  &  ans[["data"]]["v.error",]>=out[1,"sd.cl"] )																									
							out[1,c("acf.tau.l","acf.tau.u","sd.tau.l","sd.tau.u","acc")]<-	c(fixed.tau["a",c("tau.l","tau.u")],fixed.tau["v",c("tau.l","tau.u")],length(acc)/ncol(ans[["data"]]))
							out[1,c("a.hmode","v.hmode")]	<- project.nABC.movingavg.get.2D.mode(ans[["data"]]["a.theta",acc],ans[["data"]]["v.theta",acc], xlim= c(-0.2,0.4),ylim=c(1/1.7,1.7), plot=0, nbin=10, nlevels=20)
							#abline(h=ans[["xv"]],col="red"); abline(v=ans[["xa"]],col="red")
							out[1,c("acf.hmode","sd.hmode")]<- project.nABC.movingavg.get.2D.mode(ans[["data"]]["a.link",acc],ans[["data"]]["v.link",acc], xlim= c(-0.2,0.4),ylim=c(1/1.7,1.7), plot=0, nbin=10, nlevels=20)
							out[1,"a.hmode.uniest"]			<- project.nABC.movingavg.gethist(ans[["data"]]["a.theta",acc], out[1,"xa"], nbreaks= 50, width= 0.5, plot=1)[["dmode"]]
							alink.fxtau						<- ans[["data"]]["a.link",acc]-ans[["xa"]]							
							link.fxtau						<- nabc.get.locfit.links(2, as.data.frame(cbind( a=ans[["data"]]["a.theta",1:50000], sigma2=ans[["data"]]["v.theta",1:50000], VAR=ans[["data"]]["v.error",1:50000], ACF=ans[["data"]]["a.error",1:50000]/sqrt(floor(xn / (1+xmapa.leave.out)) - 3))), th.thin=1, th.sep=40	)
							#now fixed power							
							out[2,c("acf.tau.l","acf.tau.u")]<-	nabc.acf.equivalence.tau.lowup(mx.pw, 0.35, floor(xn / (1+xmapa.leave.out)), alpha)[1:2]
							out[2,c("acf.cl","acf.cu")]<-	nabc.acf.equivalence.abctol(out[2,"acf.tau.l"],out[2,"acf.tau.u"], floor(xn / (1+xmapa.leave.out)), alpha)							
							out[2,c("sd.tau.l","sd.tau.u")]<-	nabc.chisqstretch.tau.lowup(mx.pw, 2, floor(xn / (1+xsig2.leave.out))-1, alpha)[1:2]
							out[2,c("sd.cl","sd.cu")]<-	.Call("abcScaledChiSq",	c(	floor(xn / (1+xsig2.leave.out))-1, floor(xn / (1+xsig2.leave.out))-1,	out[2,"sd.tau.l"],out[2,"sd.tau.u"],	alpha,1e-10,100,0.05)	)[1:2]			
							acc<- which( 	ans[["data"]]["a.error",]<=out[2,"acf.cu"]  &  ans[["data"]]["a.error",]>=out[2,"acf.cl"]	& 
											ans[["data"]]["v.error",]<=out[2,"sd.cu"]  &  ans[["data"]]["v.error",]>=out[2,"sd.cl"] )							
							out[2,"acc"]<-	length(acc)/ncol(ans[["data"]])														
							out[2,c("a.hmode","v.hmode")]<- project.nABC.movingavg.get.2D.mode(ans[["data"]]["a.theta",acc],ans[["data"]]["v.theta",acc], xlim= c(-0.2,0.4),ylim=c(1/1.7,1.7), plot=0, nbin=10, nlevels=20)
							#abline(h=ans[["xv"]],col="red"); abline(v=ans[["xa"]],col="red")
							out[2,c("acf.hmode","sd.hmode")]<- project.nABC.movingavg.get.2D.mode(ans[["data"]]["a.link",acc],ans[["data"]]["v.link",acc], xlim= c(-0.2,0.4),ylim=c(1/1.7,1.7), plot=0, nbin=10, nlevels=20)
							out[2,"a.hmode.uniest"]	<- project.nABC.movingavg.gethist(ans[["data"]]["a.theta",acc], out[1,"xa"], nbreaks= 50, width= 0.5, plot=1)[["dmode"]]
							alink.fxpw				<- ans[["data"]]["a.link",acc]-ans[["xa"]]
							link.fxpw				<- nabc.get.locfit.links(2, as.data.frame(cbind( a=ans[["data"]]["a.theta",1:50000], sigma2=ans[["data"]]["v.theta",1:50000], VAR=ans[["data"]]["v.error",1:50000], ACF=ans[["data"]]["a.error",1:50000]/sqrt(floor(xn / (1+xmapa.leave.out)) - 3))), th.thin=1, th.sep=40	)														
							
							#jac.fxpw				<- nabc.get.locfit.jac(2, as.data.frame(cbind( a=ans[["data"]]["a.theta",1:50000], sigma2=ans[["data"]]["v.theta",1:50000], VAR=ans[["data"]]["v.error",1:50000], ACF=ans[["data"]]["a.error",1:50000]/sqrt(floor(xn / (1+xmapa.leave.out)) - 3))), th.thin=1, th.sep=40	)
							
							#tmp<- as.data.frame(cbind( a=ans[["data"]]["a.theta",1:50000], sigma2=ans[["data"]]["v.theta",1:50000], VAR=ans[["data"]]["v.error",1:50000], ACF=ans[["data"]]["a.error",1:50000]/sqrt(floor(xn / (1+xmapa.leave.out)) - 3)))
							
							#tmp<- matrix(c(seq(0,0.2,by=0.01), rep(1,21)),ncol=2,byrow=0)
							#print( tmp )
							#print( jac.fxpw[[1]][[2]] )
							#print( jac.fxpw[[2]][[1]] )
							#tmp2<- apply(tmp,1,function(x) checkjac(x[1],x[2],out[1,"xa"], out[1,"xv"]) )
							#print( tmp2[3,] )
							#print( predict(jac.fxpw[[2]][[1]], tmp ) )
							#stop()
							#print( tmp2[1,]-predict(jac.fxpw[[1]][[1]], tmp ) )
							#print( tmp2[2,]-predict(jac.fxpw[[1]][[2]], tmp ) )
							#print( tmp2[3,]-predict(jac.fxpw[[2]][[1]], tmp ) )
							#print( tmp2[4,]-predict(jac.fxpw[[2]][[2]], tmp ) )
							
							#tmp2<- apply(tmp,1,function(x) checklink(x[1],x[2],out[1,"xa"], out[1,"xv"]) )
							#print( tmp2[1,]-predict(link.fxpw[[1]], tmp ) )
							#print( tmp2[2,]-predict(link.fxpw[[2]], tmp ) )

							#print( link.fxpw[[1]] )
							#print( link.fxpw[[2]] )
							#tmp2<- apply(tmp,1,function(x) checklink(x[1],x[2],out[1,"xa"], out[1,"xv"]) )
							#print( tmp2[1,]-predict(link.fxpw[[1]], tmp ) )
							#print( tmp2[2,]-predict(link.fxpw[[2]], tmp ) )
							stop()
							if(0)
							{
								print(ans[["xa"]])
								print(ans[["xv"]])
								print(ans[["data"]][,1:5])
								print(length(acc)/ncol(ans[["data"]]))
								print(ncol(ans[["data"]]))
								tmp<- project.nABC.movingavg.gethist(ans[["data"]]["a.link",acc]-ans[["xa"]], 0, nbreaks= 50, width= 0.5, plot=1)
								abline(v=0)
								#plot(	seq(out[1,"acf.tau.l"],out[1,"acf.tau.u"],by=0.001), 
								#		nabc.acf.equivalence.pow(seq(out[1,"acf.tau.l"],out[1,"acf.tau.u"],by=0.001), out[1,"acf.tau.u"], alpha, 1/sqrt(floor(xn / (1+xmapa.leave.out))-3)), type='l')
								tmp<- project.nABC.movingavg.gethist(ans[["data"]]["a.theta",acc], ans[["xa"]], nbreaks= 50, width= 0.5, plot=1)
								print(tmp[["dmode"]])
								tmp<- project.nABC.movingavg.gethist(ans[["data"]]["v.link",acc], ans[["xv"]], nbreaks= 50, width= 0.5, plot=1)
								print(tmp[["dmode"]])
								stop()
								print(range(ans[["data"]]["a.theta",acc]))
								print(range(ans[["data"]]["v.theta",acc]))								
							}									
							#f.name<- paste(dir.name,"/nABC.MA1_largen_",N,"_",round(prior.l,d=2),"_",round(prior.u,d=2),"_",round(tau.u,d=2),"_",round(xsig2.prior.l,d=2),"_",round(xsig2.prior.u,d=2),"_",round(xsig2.tau.u,d=2),"_",xn,".R",sep='')
							#cat(paste("\nnABC.MA: save file ",f.name))
							#save(link.fxtau, link.fxpw, file=f.name, compress="xz")
							#stop()
							list(	out			= out, 
									alink.fxtau	= density(alink.fxtau, kernel="biweight",from=min(alink.fxtau),to=max(alink.fxtau),width = max(EPS,0.5*diff(summary(alink.fxtau)[c(2,5)])),n=512), 
									alink.fxpw 	= density(alink.fxpw, kernel="biweight",from=min(alink.fxpw),to=max(alink.fxpw),width = max(EPS,0.5*diff(summary(alink.fxpw)[c(2,5)])),n=512),
									link.fxpw	= link.fxpw,
									link.fxtau	= link.fxtau
									)
						})	
				names(ans)<- f.name.n$x
				cat(paste("\nnABC.MA: save file ",save.f.name))
				#save(ans,file=save.f.name)		
stop()
			}
		}
		sample.size	<- as.numeric(names(ans))
		ok.idx		<- which(sapply(seq_along(ans),function(i)		!is.matrix( ans[[i]] )	))
		#ok.idx		<- which(sapply(seq_along(ans),function(i)		all(!is.na(ans[[i]][["out"]]))	))
		alink.fxpw	<- lapply(ok.idx,function(i)		ans[[i]][["alink.fxpw"]]	)
		link.fxpw	<- lapply(ok.idx,function(i)		ans[[i]][["link.fxpw"]]	)
		alink.fxtau	<- lapply(ok.idx,function(i)		ans[[i]][["alink.fxtau"]]	)
		link.fxtau	<- lapply(ok.idx,function(i)		ans[[i]][["link.fxtau"]]	)
		ans			<- lapply(ok.idx,function(i)		ans[[i]][["out"]]	)
		sample.size	<- sample.size[ ok.idx ]
		names(link.fxpw)<- names(link.fxtau)<- names(alink.fxpw)<- names(alink.fxtau)<- sample.size
		#ans			<- lapply(ok.idx,function(i)		ans[[i]]	)		
		f.name<- paste(dir.name,"/nABC.MA1_largen_",N,"_",xn,"_",round(prior.l,d=2),"_",round(prior.u,d=2),"_",round(tau.u,d=2),"_",round(xsig2.prior.l,d=2),"_",round(xsig2.prior.u,d=2),"_",round(xsig2.tau.u,d=2),"_detjacs.R",sep='')
		options(show.error.messages = FALSE, warn=1)		
		readAttempt<-try(suppressWarnings(load(f.name)))						
		options(show.error.messages = TRUE)
		resume<- 1
		if(!resume || inherits(readAttempt, "try-error"))
		{
			cat(paste("\ncreate",f.name))
			detjac.fxpow<- sapply(seq_along(ans), function(i)
					{
						#print(sample.size[i])
						c(ans[[i]][2,"acf.tau.u"], sample.size[i], project.nABC.movingavg.avgdetJac(	ans[[i]][2,"acf.tau.l"],	ans[[i]][2,"acf.tau.u"],	ans[[i]][2,"xa"],	ans[[i]][2,"xv"], 1/sqrt(floor(sample.size[i] / (1+xmapa.leave.out))-3), alpha, empirical.rho= alink.fxpw[[i]], empirical.links=link.fxpw[[i]]))					
					})
			
			detjac.fxtau<- sapply(seq_along(ans), function(i)
					{
						#print(sample.size[i])
						c(ans[[i]][1,"acf.tau.u"], sample.size[i], project.nABC.movingavg.avgdetJac(	ans[[i]][1,"acf.tau.l"],	ans[[i]][1,"acf.tau.u"],	ans[[i]][1,"xa"],	ans[[i]][1,"xv"], 1/sqrt(floor(sample.size[i] / (1+xmapa.leave.out))-3), alpha=alpha, empirical.rho= alink.fxtau[[i]], empirical.links=link.fxtau[[i]]))					
					})
			f.name<- paste(dir.name,"/nABC.MA1_largen_",N,"_",xn,"_",round(prior.l,d=2),"_",round(prior.u,d=2),"_",round(tau.u,d=2),"_",round(xsig2.prior.l,d=2),"_",round(xsig2.prior.u,d=2),"_",round(xsig2.tau.u,d=2),"_detjacs.R",sep='')
			cat(paste("\nsave",f.name))
			save(detjac.fxpow,detjac.fxtau, file=f.name)
		}
#print(detjac.fxpow)

		detjac.fxtau.e	<- detjac.fxtau["empirical",]		
		detjac.fxtau.p	<- detjac.fxtau["pow",]
		detjac.fxtau.l	<- detjac.fxtau["lhatpow",]

		detjac.fxpow.e	<- detjac.fxpow["empirical",]
		detjac.fxpow.p	<- detjac.fxpow["pow",]
		detjac.fxpow.l	<- detjac.fxpow["lhatpow",]
		
		cols		<- c(my.fade.col("black",0.6),my.fade.col("black",0.2))
		ltys		<- c(1,1,4)
		
		if(0)
		{
		f.name<- paste(dir.name,"/nABC.MA1_largen_",N,"_",xn,"_",round(prior.l,d=2),"_",round(prior.u,d=2),"_",round(tau.u,d=2),"_",round(xsig2.prior.l,d=2),"_",round(xsig2.prior.u,d=2),"_",round(xsig2.tau.u,d=2),"_detjac.pdf",sep='')
		#pdf(f.name,version="1.4",width=pdf.width,height=pdf.height)				
		par(mar=c(4,4.5,1,.5))
		ylim<- range(	c(detjac.fxpow.e,detjac.fxtau.e,detjac.fxpow.p,detjac.fxtau.p,detjac.fxpow.l,detjac.fxtau.l)	)
		ylim<- c(-0.15,0.02)
		plot(1,1,type='n',bty='n',log='x',xlim=range(sample.size),ylim=ylim, xlab="sample size n",ylab=expression("predicted bias in a"))
		abline(h=0,lty=ltys[3])
		points(sample.size,detjac.fxtau.e,pch=20,cex=1.5,col=cols[2])
		lines(sample.size,detjac.fxtau.p,col=cols[2],lwd=2)
		lines(sample.size,detjac.fxtau.l,col=cols[2],lwd=2, lty=2)
		points(sample.size,detjac.fxpow.e,pch=18,cex=1.5,col=cols[1])
		lines(sample.size,detjac.fxpow.p,col=cols[1],lwd=2)
		lines(sample.size,detjac.fxpow.l,col=cols[1],lwd=2, lty=2)
		legend("bottomleft",bty='n',border=NA,fill=c(cols[1],"transparent",cols[2],"transparent","transparent","transparent"),lty=c(ltys[1],NA,ltys[2],NA,NA,NA),legend=c("fixed power &",expression("decreasing "*tau^'+'),"increasing power &",expression("fixed "*tau^'+'),"",""))		
		#dev.off()
		}
	

		
		f.name<- paste(dir.name,"/nABC.MA1_largen_",N,"_",xn,"_",round(prior.l,d=2),"_",round(prior.u,d=2),"_",round(tau.u,d=2),"_",round(xsig2.prior.l,d=2),"_",round(xsig2.prior.u,d=2),"_",round(xsig2.tau.u,d=2),"_amode.pdf",sep='')		
		pdf(f.name,version="1.4",width=pdf.width,height=pdf.height)				
		par(mar=c(4,4.5,1,.75))
		ylim<- range(	sapply(ans, function(x) x[,"a.hmode"]-x[1,"xa"]) 	)
		#ylim<- c(-0.04,0.015)		
		plot(1,1,type='n',bty='n',log='x',xlim=range(sample.size),ylim=ylim, xlab="sample size n",ylab=expression("bias in a: mode of "*hat(pi)[tau]*'('*a*'|'*x*') - '*a[x]))
		abline(h=0,lty=ltys[3])
		points(sample.size,sapply(ans, function(x) x[1,"a.hmode"]-x[1,"xa"]),pch=20,cex=1.5,col=cols[2])		
		points(sample.size,sapply(ans, function(x) x[2,"a.hmode"]-x[2,"xa"]),pch=23,cex=1.25,col=cols[1])
		lines(sample.size,detjac.fxtau.p,col=cols[2],lwd=2)
		lines(sample.size,detjac.fxpow.p,col=cols[1],lwd=2)
		legend("bottomleft",bty='n',border=NA,fill=c(cols[1],"transparent",cols[2],"transparent","transparent","transparent","transparent"),legend=c("fixed power &",expression("decreasing "*tau^'+'),"increasing power &",expression("fixed "*tau^'+'),"","dots: simulated",expression("lines: expected")))		
		dev.off()
stop()

		f.name<- paste(dir.name,"/nABC.MA1_largen_",N,"_",xn,"_",round(prior.l,d=2),"_",round(prior.u,d=2),"_",round(tau.u,d=2),"_",round(xsig2.prior.l,d=2),"_",round(xsig2.prior.u,d=2),"_",round(xsig2.tau.u,d=2),"_vmode.pdf",sep='')		
		pdf(f.name,version="1.4",width=pdf.width,height=pdf.height)				
		par(mar=c(4,4.5,1,.5))
		ylim<- range(	sapply(ans, function(x) x[,"v.hmode"]-x[1,"xv"]) 	)
		plot(1,1,type='n',bty='n',xlim=range(sample.size),ylim=ylim, xlab="sample size n",ylab=expression("mode of "*sigma^2*" - "*hat(sigma)[x]^2))
		abline(h=0,lty=ltys[3])
		points(sample.size,sapply(ans, function(x) x[1,"v.hmode"]-x[1,"xv"]),pch=20,cex=1.5,col=cols[2])		
		points(sample.size,sapply(ans, function(x) x[2,"v.hmode"]-x[2,"xv"]),pch=18,cex=1.5,,col=cols[1])		
		legend("topleft",bty='n',border=NA,fill=c(cols[1],"transparent",cols[2],"transparent","transparent","transparent"),lty=c(ltys[1],NA,ltys[2],NA,NA,NA),legend=c("fixed power &",expression("decreasing "*tau^'+'),"increasing power &",expression("fixed "*tau^'+'),"",""))		
		dev.off()
		
		stop()
		
		
		f.name<- paste(dir.name,"/nABC.MA.test_SDACF_pw_",xma.tau.u,"_",xsig2.tau.u,"_als.pdf",sep='')
		pdf(f.name,version="1.4",width=5,height=4.5)				
		par(mar=c(4,4.25,.5,.5))
		ylim<- range(	sapply(ans, function(x) x[,"a.star"]) 	)
		ylim<- c(0.08,0.12)
		plot(1,1,type='n',xlim=range(sample.size),ylim=ylim, xlab="sample size n",ylab="level set estimate of a")
		abline(h=0.1,lty=2, col="gray60")
		points(sample.size,sapply(ans, function(x) x[1,"a.star"]),pch=20,cex=1.5,col="gray50")		
		points(sample.size,sapply(ans, function(x) x[2,"a.star"]),pch=18,cex=1.5)		
		legend("bottomleft",bty='n',legend=c(expression("decreasing "*tau),expression("constant "*tau)),pch=c(18,20),col=c("black","gray50"))		
		dev.off()				
		
		
		f.name<- paste(dir.name,"/nABC.MA.test_SDACF_pw_",xma.tau.u,"_",xsig2.tau.u,"_vls.pdf",sep='')
		pdf(f.name,version="1.4",width=5,height=4.5)						
		par(mar=c(4,4.25,.5,.5))
		ylim<- range(	sapply(ans, function(x) x[,"v.star"]) 	)
		ylim<- c(0.7,1.3)
		plot(1,1,type='n',xlim=range(sample.size),ylim=ylim, xlab="sample size n",ylab=expression("level set estimate of "*sigma^2))
		abline(h=1,lty=2, col="gray60")
		points(sample.size,sapply(ans, function(x) x[1,"v.star"]),pch=20,cex=1.5,col="gray50")		
		points(sample.size,sapply(ans, function(x) x[2,"v.star"]),pch=18,cex=1.5)		
		legend("bottomleft",bty='n',legend=c(expression("decreasing "*tau),expression("constant "*tau)),pch=c(18,20),col=c("black","gray50"))		
		dev.off()
		
		
		
		f.name<- paste(dir.name,"/nABC.MA.test_SDACF_pw_",xma.tau.u,"_",xsig2.tau.u,"_vmode.pdf",sep='')
		pdf(f.name,version="1.4",width=5,height=4.5)				
		par(mar=c(4,4.25,.5,.5))
		ylim<- range(	sapply(ans, function(x) x[,"v.mode.d"]) 	)
		plot(1,1,type='n',xlim=range(sample.size),ylim=ylim, xlab="sample size n",ylab=expression("mode of "*sigma^2))
		abline(h=1,lty=2, col="gray60")
		points(sample.size,sapply(ans, function(x) x[1,"v.mode.d"]),pch=20,cex=1.5,col="gray50")		
		points(sample.size,sapply(ans, function(x) x[2,"v.mode.d"]),pch=18,cex=1.5)		
		legend("topleft",bty='n',legend=c(expression("decreasing "*tau),expression("constant "*tau)),pch=c(18,20),col=c("black","gray50"))		
		dev.off()
		
		f.name<- paste(dir.name,"/nABC.MA.test_SDACF_pw_",xma.tau.u,"_",xsig2.tau.u,"_accrate.pdf",sep='')
		pdf(f.name,version="1.4",width=5,height=4.5)		
		par(mar=c(4,4.25,.5,.5))
		ylim<- range(	c(0,sapply(ans, function(x) x[,"accr"])) 	)
		plot(1,1,type='n',xlim=range(sample.size),ylim=ylim, xlab="sample size n",ylab="acceptance prob")
		points(sample.size,sapply(ans, function(x) x[1,"accr"]),pch=20,cex=1.5,col="gray50")		
		points(sample.size,sapply(ans, function(x) x[2,"accr"]),pch=18,cex=1.5)		
		legend("topleft",bty='n',legend=c(expression("decreasing "*tau),expression("constant "*tau)),pch=c(18,20),col=c("black","gray50"))		
		dev.off()
		
		
		#xma.tau.u<- 0.2		
		#ac<- xma.tau.u*(sample.size-4)+qnorm(alpha)		
		a.tau<- sapply(seq_along(ans), function(i) (ans[[i]][2,"ac"]-qnorm(alpha))/(sample.size[i]-4))			
		#print(a.tau)				
		v.tau<- sapply(seq_along(ans), function(i)
				{
					uniroot(	fstretch.solvefor.tau,	c(1.01,2.5), tol=.Machine$double.eps^0.5, cir=exp(ans[[i]][2,"vc"]), n=sample.size[i]-1, alpha=alpha)$root			
				})
		#print(v.tau)
		
		f.name<- paste(dir.name,"/nABC.MA.test_SDACF_pw_",xma.tau.u,"_",xsig2.tau.u,"_tau.pdf",sep='')
		pdf(f.name,version="1.4",width=5,height=4.5)
		par(mar=c(4,4.25,.5,.5))
		plot(1,1,type='n',xlim=range(sample.size),ylim=c(0,1), xlab="sample size n",ylab=expression('relative decrease in '*tau))
		points(sample.size,a.tau/xma.tau.u,pch=18,cex=1.5,col="black")
		points(sample.size,v.tau/xsig2.tau.u,pch=5,cex=1,col="black")		
		points(sample.size,rep(1,length(sample.size)),pch=20,cex=1.5,col="gray50")
		
		legend("bottomleft",bty='n',legend=c(	expression("constant "*tau[1]*"=2.4, constant "*tau[2]*"=0.2"),
						expression("decreasing "*tau[1]*" for "*sigma^2),
						expression("decreasing "*tau[2]*" for "*a)),
				pch=c(20,5,18),col=c("gray50","black","black"))
		dev.off()
		stop()
		
		
		#NOTE: the CIL and CIR may increase for constant TAU as sample size increases
		f.name<- paste(dir.name,"/nABC.MA.test_SDACF_pw_",xma.tau.u,"_",xsig2.tau.u,"_c.pdf",sep='')
		#pdf(f.name,version="1.4",width=5,height=4.5)		
		par(mar=c(4,4.25,.5,4))
		ylim<- apply(	sapply(ans, function(x) x[1,c("ac","vc")]), 1, max 	)		
		plot(1,1,type='n',xlim=range(sample.size),ylim=c(0,1), xlab="sample size n",ylab='',yaxt='n')
		axis(2,at=c(50,100,150,200)/ylim["ac"],labels=c(50,100,150,200))		
		mtext(expression("tolerance "*c[2]^'+'),2,line=2.5)
		axis(4,at=log(c(1.2,1.4,1.6,1.8))/ylim["vc"],labels=c(1.2,1.4,1.6,1.8))
		mtext(expression("tolerance "*c[1]^'+'),4,line=2.5)		
		points(sample.size,sapply(ans, function(x) x[1,"ac"])/ylim["ac"],pch=20,cex=1.5,col="gray50")
		points(sample.size,sapply(ans, function(x) x[1,"vc"])/ylim["vc"],pch=1,cex=1.1,col="gray50")		
		points(sample.size,sapply(ans, function(x) x[2,"ac"])/ylim["ac"],pch=18,cex=1.5)
		points(sample.size,sapply(ans, function(x) x[2,"vc"])/ylim["vc"],pch=5,cex=1.1)
		legend("topleft",bty='n',legend=c(expression("decreasing "*tau),expression("constant "*tau)),pch=c(18,20),col=c("black","gray50"))
		stop()
		abline(h=0.1,lty=3)		
		lines(sample.size,ans["v.mode.h",],lty=2,col="blue")
		lines(sample.size,ans["v.mode.d",],lty=3,col="blue")
		
		stop()
		

		if(1)
		{
			cat(paste("\nload",f.names[100]))
			load(f.names[100])
			print(range(ans["v.ytheta",]))
			#enforce symmetric prior around true value in a
			ans<- ans[,which(ans["a.ytheta",]<=xma.prior.u & ans["a.ytheta",]>=xma.prior.l)]				
			acc<- which( 	ans["a.error",]<=ans["a.cir",]  &  ans["a.error",]>=ans["a.cil",] &
							ans["v.error",]<=ans["v.cir",]  &  ans["v.error",]>=ans["v.cil",]		)		
			a.ftau.h<- project.nABC.movingavg.gethist(ans["a.ytheta",acc], 0.1, nbreaks= 50)				
			tmp<- density(ans["a.ytheta",acc], kernel="biweight",width = max(EPS,0.5*diff(summary(ans["a.ytheta",acc])[c(2,5)])))
			a.ftau.m<- tmp[["x"]][which.max( tmp[["y"]])]
			tmp<- project.nABC.movingavg.a2rho(ans["a.ytheta",acc]) - project.nABC.movingavg.a2rho(0.1)
			r.ftau.h<- project.nABC.movingavg.gethist(tmp, 0, nbreaks= 50)
			tmp<- density(tmp, kernel="biweight",width = max(EPS,0.5*diff(summary(tmp)[c(2,5)])))
			r.ftau.m<- tmp[["x"]][which.max( tmp[["y"]])]																		
			
			a.c<- quantile( abs(ans["a.error",]), probs=c(.1) )						
			v.c<- quantile( abs(log(ans["v.error",which( abs(ans["a.error",])<=a.c )])), probs=c(.1) )						
			acc<- which( 	ans["a.error",]<=a.c  &  ans["a.error",]>=-a.c &
							ans["v.error",]<=exp(v.c)  &  ans["v.error",]>=exp(-v.c)		)														
			a.facc.h<- project.nABC.movingavg.gethist(ans["a.ytheta",acc], 0.1, nbreaks= 50)								
			tmp<- density(ans["a.ytheta",acc], kernel="biweight",width = max(EPS,0.5*diff(summary(ans["a.ytheta",acc])[c(2,5)])))
			a.facc.m<- tmp[["x"]][which.max( tmp[["y"]])]
			tmp<- project.nABC.movingavg.a2rho(ans["a.ytheta",acc]) - project.nABC.movingavg.a2rho(0.1)
			r.facc.h<- project.nABC.movingavg.gethist(tmp, 0, nbreaks= 50)
			tmp<- density(tmp, kernel="biweight",width = max(EPS,0.5*diff(summary(tmp)[c(2,5)])))
			r.facc.m<- tmp[["x"]][which.max( tmp[["y"]])]																		
			
			f.name<- paste(dir.name,"/nABC.MA.test_SDACF_pw_",xma.tau.u,"_",xsig2.tau.u,"_example1.pdf",sep='')
			pdf(f.name,version="1.4",width=5,height=4.5)						
			par(mar=c(4,4.25,.5,.5))
			plot(1,1,type='n',xlim=range(c(r.ftau.h[["breaks"]],r.facc.h[["breaks"]])),ylim=range(c(r.ftau.h[["density"]],r.facc.h[["density"]])),xlab=expression(rho[2]),ylab="density")
			cols<- c("gray60","gray30")
			plot(r.ftau.h,freq=0,add=1,col=my.fade.col(cols[1],0.8),border=my.fade.col(cols[1],0.8))
			plot(r.facc.h,freq=0,add=1,col=my.fade.col(cols[2],0.8),border=my.fade.col(cols[2],0.8))
			abline(v=r.ftau.m,lty=2,col=cols[1])
			abline(v=r.facc.m,lty=2)
			legend("topright",bty='n',legend=c(expression("decreasing "*tau),expression("constant "*tau)),fill=c(cols[2],cols[1]),border=NA)
			dev.off()
			
			f.name<- paste(dir.name,"/nABC.MA.test_SDACF_pw_",xma.tau.u,"_",xsig2.tau.u,"_example2.pdf",sep='')
			pdf(f.name,version="1.4",width=3.6,height=5)								
			par(mar=c(4,3.85,0.2,0.2))
			plot(1,1,type='n',xlim=range(c(a.ftau.h[["breaks"]],a.facc.h[["breaks"]])),ylim=range(c(a.ftau.h[["density"]],a.facc.h[["density"]])),xlab="a",ylab="density")
			cols<- c("gray60","gray30")
			plot(a.ftau.h,freq=0,add=1,col=my.fade.col(cols[1],0.8),border=my.fade.col(cols[1],0.8))
			plot(a.facc.h,freq=0,add=1,col=my.fade.col(cols[2],0.8),border=my.fade.col(cols[2],0.8))
			abline(v=a.ftau.m,lty=2,col=cols[1])
			abline(v=a.facc.m,lty=2)
			legend("topright",bty='n',legend=c(expression("decreasing "*tau),expression("constant "*tau)),fill=c(cols[2],cols[1]),border=NA)
			dev.off()
		}
		
		
	}
	if(0)	#check mode estimation & link function for estimating both sigma2 and mapa	USE ACF & SD(lag2)
	{
		require(locfit)
		my.mkdir(DATA,"nABC.movingavg_mode")
		dir.name<- paste(DATA,"nABC.movingavg_mode",sep='/')
		resume<- 1
		verbose<- 1
		N<- 1.2e3
		alpha<- 0.01
		nbreaks<- 20		
		xma.pa<- 0.1
		#xma.tau.u<- 0.025
		xma.tau.u<- 0.2
		xsig2<- 1
		#xsig2.tau.u<- 1.13
		xsig2.tau.u<- 2.4
		#the round ma.n's ie 500 1000 5000 used a unif prior on mapa space
		#the odd ma.n's ie 501 1001 5001 used a unif prior on rho space		
		ma.n<- 10001
				
		if(exists("argv"))
		{
			tmp<- na.omit(sapply(argv,function(arg)
							{	switch(substr(arg,2,2),
										m= return(as.numeric(substr(arg,3,nchar(arg)))),NA)	}))
			if(length(tmp)>0) ma.n<- tmp[1]
			tmp<- na.omit(sapply(argv,function(arg)
							{	switch(substr(arg,2,2),
										M= return(as.numeric(substr(arg,3,nchar(arg)))),NA)	}))
			if(length(tmp)>0) N<- tmp[1]
			tmp<- na.omit(sapply(argv,function(arg)
							{	switch(substr(arg,2,2),
										t= return(as.numeric(substr(arg,3,nchar(arg)))),NA)	}))
			if(length(tmp)>0) xma.tau.u<- tmp[1]
			tmp<- na.omit(sapply(argv,function(arg)
							{	switch(substr(arg,2,2),
										s= return(as.numeric(substr(arg,3,nchar(arg)))),NA)	}))
			if(length(tmp)>0) xsig2.tau.u<- tmp[1]
		}
		
		
		r.xma.pa<- project.nABC.movingavg.a2nu(xma.pa)		#r for xma.pa
		z.xma.pa<- project.nABC.movingavg.a2rho(xma.pa)		#r for xma.pa						
		xma.tau.l<- -xma.tau.u
		xsig2.tau.l<- 1/xsig2.tau.u
		if(verbose)	cat(paste("true xmapa, correlation scale",r.xma.pa,"true xmapa, test scale",z.xma.pa,"\n"))
		
		if(0)
		{
			xsig2.prior.u<- xsig2.tau.u 
			xsig2.prior.l<- xsig2.tau.l			
			xma.prior.u<- project.nABC.movingavg.rho2a( z.xma.pa+xma.tau.u )		
			xma.prior.l<- project.nABC.movingavg.rho2a( z.xma.pa+xma.tau.l )			
		}
		else if(0)
		{
			xsig2.prior.u<- xsig2.tau.u 
			xsig2.prior.l<- xsig2.tau.l			
			xma.prior.u<- project.nABC.movingavg.rho2a( z.xma.pa+xma.tau.u*15 )		
			xma.prior.l<- project.nABC.movingavg.rho2a( z.xma.pa+xma.tau.l*15 )			
		}
		else{
			xsig2.prior.u<- 4#1.2 
			xsig2.prior.l<- 0.2#0.8			
			xma.prior.u<- project.nABC.movingavg.rho2a( .423 )		
			xma.prior.l<- project.nABC.movingavg.rho2a( -.423 )
		}
		if(verbose)	cat(paste("prior mapa thresholds from test scale",xma.prior.l,xma.prior.u,"\n"))	
			
		f.name<- paste(dir.name,"/nABC.MA.test_SDACF_dpw_",N,"_",round(xma.tau.u,d=6),"_",xsig2.tau.u,"_",ma.n,".R",sep='')
		#f.name<- paste(dir.name,"/nABC.MA.test_SDACF_",N,"_",round(xma.tau.u,d=6),"_",xsig2.tau.u,"_",ma.n,".R",sep='')
		cat(paste("\nnABC.StretchedF: attempt to read ",f.name))
	
		options(show.error.messages = FALSE, warn=1)		
		readAttempt<-try(suppressWarnings(load(f.name)))						
		options(show.error.messages = TRUE)						
		if(!resume || inherits(readAttempt, "try-error"))
		{			
			sparse<- substr(f.name,1,nchar(f.name)-2)	
			ans<- project.nABC.movingavg.unifsigma.unifma(	N,alpha,ma.n,xma.pa,xsig2,
															xma.tau.l,xma.tau.u,xma.prior.l,xma.prior.u,
															xsig2.tau.l,xsig2.tau.u,xsig2.prior.l,xsig2.prior.u, sparse=sparse		)																										
			if(is.na(sparse))
			{
				cat(paste("\nnABC.StretchedF: save ",f.name))
				save(ans,file=f.name)
			}
		}
		else
			cat(paste("\nnABC.StretchedF: resumed ",f.name))
stop()			
		acc<- which( 	ans["a.error",]<=ans["a.cir",]  &  ans["a.error",]>=ans["a.cil",] &
						ans["v.error",]<=ans["v.cir",]  &  ans["v.error",]>=ans["v.cil",]		)	
		
		print(c(ncol(ans),length(acc)/ncol(ans),ans["a.cil",1],ans["a.cir",1],ans["v.cil",1],ans["v.cir",1]))
		print(ans[,1:10])
		print(length(acc)/ncol(ans))			
		stop()
		#m<- cbind( ans["a.ytheta",acc], ans["v.ytheta",acc], ans["a.link.mc",acc]-z.xma.pa, log(ans["v.link.mc",acc])-0 )
		m<- cbind( ans["a.ytheta",acc], ans["v.ytheta",acc], ans["a.link.mc",acc]-z.xma.pa, log(ans["v.link.mc",acc])-log(1.01) )
		colnames(m)<- c("a","sigma2","ACF","VAR")
		m<- as.data.frame(m)
		theta.names<- c("a","sigma2")
		links.names<- c("ACF","VAR")
		theta0<- project.nABC.movingavg.estimateTheta0(m, c("a","sigma2"), c("ACF","VAR"))
		print(theta0)
		
		
		stop()
		acc<- which( 	ans["a.error",]<=ans["a.cir",]  &  ans["a.error",]>=ans["a.cil",] &
						ans["v.error",]<=ans["v.cir",]  &  ans["v.error",]>=ans["v.cil",]		)	
		#plot histogram on theta-space.
		f.name<- paste(dir.name,"/nABC.MA.test_SDACF_",N,"_",round(xma.tau.u,d=6),"_",xsig2.tau.u,"_",ma.n,"_abcfit.pdf",sep='')
		pdf(f.name,version="1.4",width=4,height=4)
		par(mar=c(4,4.5,0.5,0.75))
		plot.2D.dens(ans["a.ytheta",acc],ans["v.ytheta",acc],xlab="a",ylab=expression(sigma^2),xlim= range(ans["a.ytheta",]), ylim= range(ans["v.ytheta",]),method="ash",zero.abline=0, palette="gray")
		abline(v=xma.pa,lty=3)
		abline(h=xsig2,lty=3)
		points(xma.pa,xsig2,col="red", pch=16)
		dev.off()
		
		#plot histogram on rho-space. should have max at z.xma.pa
		ans.ylink<- project.nABC.movingavg.gethist(ans["a.ylink",acc], z.xma.pa, nbreaks= nbreaks)		
		print(ans.ylink)
		plot(ans.ylink, main=N)
		abline(v=z.xma.pa, col="red")
		#stop()		
		ans.ylink<- project.nABC.movingavg.gethist((1+xma.pa*xma.pa)*ans["v.ylink",acc], (1+xma.pa*xma.pa)*xsig2, nbreaks= nbreaks)
		plot(ans.ylink, main=N)
		abline(v=(1+xma.pa*xma.pa)*xsig2, col="red")	
		#stop()
#bookmark_linkMA		

		#reconstruct the link function with local linear regression
		f.name<- paste(dir.name,"/nABC.MA.test_SDACF_",N,"_",round(xma.tau.u,d=6),"_",xsig2.tau.u,"_",ma.n,"_locfit.pdf",sep='')
		pdf(f.name,version="1.4",width=4,height=4)
		thin<- 2000
		lnk.df<- data.frame(mapa=ans["a.ytheta",], sig2=ans["v.ytheta",], rhomc=ans["a.link.mc",] )		
		lnk.df<- lnk.df[seq.int(1,nrow(lnk.df),by=thin),]
		x<- locfit(rhomc~mapa:sig2, data=lnk.df)
		out<- plot.persplocfit(x, pv= c("mapa","sig2"), xlab= "a", ylab= expression(sigma^2), zlab= expression(rho2), theta=-50, phi=10	)					
		lines (trans3d(out$x, min(out$y), z= project.nABC.movingavg.a2rho(out$x), pmat = out$pmat), col = "black", lty=4)
		lines (trans3d(min(out$x), out$y, z= project.nABC.movingavg.a2rho(min(out$x)), pmat = out$pmat), col = "black", lty=4)
		lines (trans3d(project.nABC.movingavg.rho2a(z.xma.pa), out$y, z= z.xma.pa, pmat = out$pmat), col = "red", lty=1, lwd=1.5)		
		dev.off()
stop()		
		tmp<- project.nABC.movingavg.get.2D.mode(ans["a.ytheta",acc],ans["v.ytheta",acc], xlim= c(-0.2,0.4),ylim=c(0.8,1.2), plot=0)
		print(tmp)
stop()
					
	}
	if(0)	#check mode estimation & link function for estimating both sigma2 and mapa	USE SD(lag2)
	{
		require(locfit)
		my.mkdir(DATA,"nABC.movingavg_mode")
		dir.name<- paste(DATA,"nABC.movingavg_mode",sep='/')
		resume<- 1
		verbose<- 1
		N<- 1e5
		alpha<- 0.01
		nbreaks<- 25		
		xma.pa<- 0.1
		xma.tau.u<- 0.025
		xsig2<- 1
		xsig2.tau.u<- 1.13
		#the round ma.n's ie 500 1000 5000 used a unif prior on mapa space
		#the odd ma.n's ie 501 1001 5001 used a unif prior on rho space		
		ma.n<- 10006
		
		if(exists("argv"))
		{
			tmp<- na.omit(sapply(argv,function(arg)
							{	switch(substr(arg,2,2),
										m= return(as.numeric(substr(arg,3,nchar(arg)))),NA)	}))
			if(length(tmp)>0) ma.n<- tmp[1]
			tmp<- na.omit(sapply(argv,function(arg)
							{	switch(substr(arg,2,2),
										M= return(as.numeric(substr(arg,3,nchar(arg)))),NA)	}))
			if(length(tmp)>0) N<- 1e5*tmp[1]
			tmp<- na.omit(sapply(argv,function(arg)
							{	switch(substr(arg,2,2),
										t= return(as.numeric(substr(arg,3,nchar(arg)))),NA)	}))
			if(length(tmp)>0) xma.tau.u<- tmp[1]
			tmp<- na.omit(sapply(argv,function(arg)
							{	switch(substr(arg,2,2),
										s= return(as.numeric(substr(arg,3,nchar(arg)))),NA)	}))
			if(length(tmp)>0) xsig2.tau.u<- tmp[1]
		}
		
		
		r.xma.pa<- project.nABC.movingavg.a2nu(xma.pa)		#r for xma.pa
		z.xma.pa<- project.nABC.movingavg.a2rho(xma.pa)		#r for xma.pa
		xsig2.tau.l<- 1/xsig2.tau.u
		xma.tau.l<- -xma.tau.u		
		if(verbose)	cat(paste("true xmapa, correlation scale",r.xma.pa,"true xmapa, test scale",z.xma.pa,"true sigma",(1+xma.pa*xma.pa)*xsig2,"\nN is",N,"\nma.n is",ma.n))
		
		if(1)
		{
			xma.prior.u<- project.nABC.movingavg.rho2a( .423 )			
			xma.prior.l<- project.nABC.movingavg.rho2a( -.423 )
			
			xsig2.prior.u<- 1.2 
			xsig2.prior.l<- 0.8
		}	
		else
		{
			xma.prior.u<- project.nABC.movingavg.rho2a( z.xma.pa+xma.tau.u*15 )		
			xma.prior.l<- project.nABC.movingavg.rho2a( z.xma.pa+xma.tau.l*15 )
			xsig2.prior.u<- xsig2.tau.u 
			xsig2.prior.l<- xsig2.tau.l
		}
		if(verbose)	cat(paste("\nnABC.movingavg_mode prior mapa thresholds from test scale",xma.prior.l,xma.prior.u,"\n"))	
		if(verbose)	cat(paste("\nnABC.movingavg_modeprior sig2 thresholds from test scale",xsig2.prior.l,xsig2.prior.u,"\n"))
		
		f.name<- paste(dir.name,"/nABC.MA.test_SD_",N,"_",round(xma.tau.u,d=6),"_",xsig2.tau.u,"_",ma.n,".R",sep='')
		cat(paste("\nnABC.movingavg_mode: attempt to read ",f.name))	
		options(show.error.messages = FALSE, warn=1)		
		readAttempt<-try(suppressWarnings(load(f.name)))						
		options(show.error.messages = TRUE)						
		if(!resume || inherits(readAttempt, "try-error"))
		{
			ans.lag20<- project.nABC.movingavg.unifsigma.unifma(		N,alpha,ma.n,xma.pa,xsig2,
																xma.tau.l,xma.tau.u,xma.prior.l,xma.prior.u,
																xsig2.tau.l,xsig2.tau.u,xsig2.prior.l,xsig2.prior.u, lag2=0		)
			ans<- project.nABC.movingavg.unifsigma.unifma(		N,alpha,ma.n,xma.pa,xsig2,
																xma.tau.l,xma.tau.u,xma.prior.l,xma.prior.u,
																xsig2.tau.l,xsig2.tau.u,xsig2.prior.l,xsig2.prior.u		)											

			cat(paste("\nnABC.movingavg_mode: save ",f.name))
			save(ans, ans.lag20,file=f.name)
		}
		else
			cat(paste("\nnABC.movingavg_mode: resumed ",f.name))
		if(1)
		{
			#plot histogram of p-values
			nbreaks<- 40
			cols<- c("gray60","gray30")		
			acc<- which( 	ans.lag20["v.error",]<=ans.lag20["v.cir",]  &  ans.lag20["v.error",]>=ans.lag20["v.cil",] 	)
			h.ts.lag20<- hist( ans.lag20["v.ts.pval",acc], breaks=nbreaks, plot=0 )
			h.pf.lag20<- hist( ans.lag20["v.pfam.pval",acc], breaks=nbreaks, plot=0 )
			acc<- which( 	ans["v.error",]<=ans["v.cir",]  &  ans["v.error",]>=ans["v.cil",] 	)		
			h.ts.lag21<- hist( ans["v.ts.pval",acc], breaks=nbreaks, plot=0)
			h.pf.lag21<- hist( ans["v.pfam.pval",acc], breaks=nbreaks, plot=0)
			#plot histogram of p-values: ts-pval		
			f.name<- paste(dir.name,"/nABC.MA.test_SD_tspval_",N,"_",round(xma.tau.u,d=6),"_",xsig2.tau.u,"_",ma.n,"_abcfit.pdf",sep='')
			pdf(f.name,version="1.4",width=4,height=4)
			par(mar=c(4.25,4.25,0.5,0.5))
			#ylim<- range(c(h.ts.lag20$density,h.ts.lag21$density))
			ylim<- range(c(h.ts.lag20$density,h.ts.lag21$density)); ylim[1]<- 0; ylim[2]<- ylim[2]
			plot(1,1,type='n',xlim=c(0,1),ylim=ylim,xlab="p-value",ylab="density",main='')		
			plot(h.ts.lag21,col=my.fade.col(cols[1],0.8),border=my.fade.col(cols[1],0.8),add=TRUE,freq=FALSE)
			plot(h.ts.lag20,col=my.fade.col(cols[2],0.8),border=my.fade.col(cols[2],0.8),add=TRUE,freq=FALSE)
			legend(x=0.1,y=ylim[2]*1,fill=c(cols[1],cols[2]),legend=c("thinned time series","unthinned time series"),bty='n',border=NA)
			dev.off()
			
			f.name<- paste(dir.name,"/nABC.MA.test_SD_pfampval_",N,"_",round(xma.tau.u,d=6),"_",xsig2.tau.u,"_",ma.n,"_abcfit.pdf",sep='')
			pdf(f.name,version="1.4",width=4,height=4)
			par(mar=c(4.25,4.25,0.5,0.5))		
			ylim<- range(c(h.pf.lag20$density,h.pf.lag21$density)); ylim[1]<- 0; ylim[2]<- ylim[2]*1.3
			#ylim<- range(c(h.pf.lag20$counts,h.pf.lag21$counts))
			plot(1,1,type='n',xlim=c(0,1),ylim=ylim,xlab="p-value",ylab="density",main='')		
			plot(h.pf.lag21,col=my.fade.col(cols[1],0.8),border=my.fade.col(cols[1],0.8),add=TRUE,freq=0)
			plot(h.pf.lag20,col=my.fade.col(cols[2],0.8),border=my.fade.col(cols[2],0.8),add=TRUE,freq=0)
			legend(x=0.1,y=ylim[2]*1,fill=c(cols[1],cols[2]),legend=c("thinned time series","unthinned time series"),bty='n',border=NA)
			dev.off()
			stop()
		}
		acc<- which( 	ans["v.error",]<=ans["v.cir",]  &  ans["v.error",]>=ans["v.cil",] 	)		
		print(c(ncol(ans),length(acc)/ncol(ans),ans["v.cil",1],ans["v.cir",1]))
		nbreaks<- 25
		#plot histogram on theta-space. should have max at z.xma.pa
		f.name<- paste(dir.name,"/nABC.MA.test_SD_",N,"_",round(xma.tau.u,d=6),"_",xsig2.tau.u,"_",ma.n,"_abcfit.pdf",sep='')
		pdf(f.name,version="1.4",width=4,height=4)
		par(mar=c(4,4.5,0.5,0.75))
		plot.2D.dens(ans["a.ytheta",acc],ans["v.ytheta",acc],xlab="a",ylab=expression(sigma^2),xlim= range(ans["a.ytheta",]), ylim= range(ans["v.ytheta",]),method="ash",zero.abline=0, palette="gray")
		abline(v=xma.pa,lty=3)
		abline(h=xsig2,lty=3)
		tmp<- seq(min(ans["a.ytheta",]),max(ans["a.ytheta",]),0.001)
		lines(tmp,(1+xma.pa*xma.pa)*xsig2/(1+tmp*tmp),type='l',col="red", lwd=1.5, lty=1)
		dev.off()
#bookmark_linkMA
		#plot histogram on rho-space. should have max at z.xma.pa
		ans.ylink<- project.nABC.movingavg.gethist(ans["v.ylink",acc], (1+xma.pa*xma.pa)*xsig2, nbreaks= nbreaks)		
		print(ans.ylink)
		plot(ans.ylink, main=N)
		abline(v=(1+xma.pa*xma.pa)*xsig2, col="red")
		
		#reconstruct the link function with local linear regression
		f.name<- paste(dir.name,"/nABC.MA.test_SD_",N,"_",round(xma.tau.u,d=6),"_",xsig2.tau.u,"_",ma.n,"_locfit.pdf",sep='')
		pdf(f.name,version="1.4",width=4,height=4)
		thin<- 2000
		lnk.df<- data.frame(mapa=ans["a.ytheta",], sig2=ans["v.ytheta",], rhomc=ans["v.link.mc",] )		
		lnk.df<- lnk.df[seq.int(1,nrow(lnk.df),by=thin),]
		x<- locfit(rhomc~mapa:sig2, data=lnk.df)
		out<- plot.persplocfit(x, pv= c("mapa","sig2"), xlab= "a", ylab= expression(sigma^2), zlab= expression(rho1), palette="gray"	)					
		lines(trans3d(out$x, min(out$y), z= (1+out$x*out$x)*min(out$y), pmat = out$pmat), col = "black",lty=4)
		lines(trans3d(max(out$x), out$y, z= (1+max(out$x)^2)*out$y, pmat = out$pmat), col = "black",lty=4)
		tmp<- seq(min(out$x),sqrt((1+xma.pa*xma.pa)*xsig2 / min(out$y) - 1)*0.95,0.001)
		lines(trans3d(x=tmp, y=(1+xma.pa*xma.pa)*xsig2/(1+tmp*tmp), z= (1+xma.pa*xma.pa)*xsig2, pmat = out$pmat), col = "red", lwd=1.5, lty=1)
		dev.off()
		#points( trans3d(lnk.df[["mapa"]],lnk.df[["sig2"]],lnk.df[["rhomc"]],pmat=out$pmat), col=my.fade.col("red",0.5), pch=16)		
		stop()																	
	}
	if(0)	#check mode estimation & link function for tau=0.025 and various ma.n
	{
		my.mkdir(DATA,"nABC.MA_test")
		dir.name<- paste(DATA,"nABC.MA_test",sep='/')
		resume<- 1
		N<- 5e5
		tau.u<- 0.025		
		xma.pa<- 0.4
		#the round ma.n's ie 500 1000 5000 used a unif prior on mapa space
		#the odd ma.n's ie 501 1001 5001 used a unif prior on rho space		 		
		ma.n<- 10001
		alpha<- 0.01
		nbreaks<- 10
		xma.sigma<- sqrt(2)
		verbose<- 1
		if(exists("argv"))
		{
			tmp<- na.omit(sapply(argv,function(arg)
							{	switch(substr(arg,2,2),
										m= return(as.numeric(substr(arg,3,nchar(arg)))),NA)	}))
			if(length(tmp)>0) ma.n<- tmp[1]
			tmp<- na.omit(sapply(argv,function(arg)
							{	switch(substr(arg,2,2),
										M= return(as.numeric(substr(arg,3,nchar(arg)))),NA)	}))
			if(length(tmp)>0) N<- 1e5*tmp[1]
			tmp<- na.omit(sapply(argv,function(arg)
							{	switch(substr(arg,2,2),
										t= return(as.numeric(substr(arg,3,nchar(arg)))),NA)	}))
			if(length(tmp)>0) tau.u<- tmp[1]
		}
		
		
		r.xma.pa<- project.nABC.movingavg.a2nu(xma.pa)		#r for xma.pa
		z.xma.pa<- project.nABC.movingavg.a2rho(xma.pa)		#r for xma.pa						
		tau.l<- -tau.u		
		if(verbose)	cat(paste("true xmapa, correlation scale",r.xma.pa,"true xmapa, test scale",z.xma.pa,"\n"))		
		
		#prior.u<- project.nABC.movingavg.nu2a( r.xma.pa+tau.u )		
		#prior.l<- project.nABC.movingavg.nu2a( r.xma.pa+tau.l )
		#if(verbose)	cat(paste("prior mapa thresholds from correlation scale",prior.l,prior.u,"\n"))
		
		prior.u<- project.nABC.movingavg.rho2a( z.xma.pa+tau.u )		
		prior.l<- project.nABC.movingavg.rho2a( z.xma.pa+tau.l )
		if(verbose)	cat(paste("prior mapa thresholds from test scale",prior.l,prior.u,"\n"))	
		
		prior.u<- 0.4 + 0.035
		prior.l<- 0.4 - 0.035
		tau.u<- project.nABC.movingavg.a2rho(prior.u)-z.xma.pa
		tau.l<- project.nABC.movingavg.a2rho(prior.l)-z.xma.pa
		print(c(project.nABC.movingavg.a2rho( prior.l )-z.xma.pa, project.nABC.movingavg.a2rho( prior.u )-z.xma.pa))
		
		f.name<- paste(dir.name,"/nABC.MA.test_",N,"_",round(tau.u,d=6),"_",ma.n,".R",sep='')
		cat(paste("\nnABC.StretchedF: attempt to read ",f.name))
		options(show.error.messages = FALSE, warn=1)		
		readAttempt<-try(suppressWarnings(load(f.name)))						
		options(show.error.messages = TRUE)						
		if(!resume || inherits(readAttempt, "try-error"))
		{
			ans3<- project.nABC.movingavg.fixsigma.unifma(N,tau.l,tau.u,prior.l,prior.u,alpha,ma.n,xma.sigma,xma.pa)
			cat(paste("\nnABC.StretchedF: save ",f.name))
			save(ans3,file=f.name)
		}
		else
			cat(paste("\nnABC.StretchedF: resumed ",f.name))

		ans<- ans3
		print(ans[,1:10])
		acc<- which( ans["error",]<=ans["cir",]  &  ans["error",]>=ans["cil",] )		
		print(c(ncol(ans),length(acc)/ncol(ans),ans["cil",1],ans["cir",1]))
		
		#plot histogram on rho-space. should have max at z.xma.pa
		ans.ylink<- project.nABC.movingavg.gethist(ans["ylink",acc], z.xma.pa, nbreaks= nbreaks)
		ans.d<- density(ans["ylink",acc], kernel="biweight",width = max(EPS,diff(summary(ans["ylink",acc])[c(2,5)])))
		ans.d.mode<- ans.d[["x"]][which.max( ans.d[["y"]] )]
		print(ans.ylink)
		print(ans.d.mode)
		plot(ans.ylink, main=N)
		abline(v=z.xma.pa, col="red")
		abline(v=ans.d.mode, col="blue")
		abline(v=ans.ylink[["mean"]], col="green")
stop()
		
		#reconstruct the link function with local linear regression
		plot(ans["ytheta",acc], ans["link.mc",acc], pch=19, col=my.fade.col("gray80",0.25))
		abline(v=xma.pa,lty=3)
		abline(h=z.xma.pa,lty=3)		
		tmp<- seq(min(ans["ytheta",acc]),max(ans["ytheta",acc]),0.001)
		lines(tmp, project.nABC.movingavg.a2rho( tmp ),col="red")
		library(KernSmooth)
		h <- dpill(ans["ytheta",acc], ans["link.mc",acc])
		lines( locpoly(ans["ytheta",acc], ans["link.mc",acc], bandwidth = h, degree=1), col="blue" )
		stop()				
	}
	if(!is.na(subprog) && subprog==3)		#analytical power
	{
		#TODOO
		xa<- 0.1
		xsig2<- 1
		xv<- (1+xa*xa)*xsig2 
		
		ltys<- c(1,2,4)
		cols<- c(my.fade.col("black",1),my.fade.col("black",1),my.fade.col("black",0.3))
		a<- seq(-0.5,0.5,0.001)
		nu2<- project.nABC.movingavg.a2nu(a)
		rho2<- project.nABC.movingavg.a2rho(a)
		jac<- checkjac(a,xsig2,xa,xv)
		jac<- apply(matrix(jac,ncol=length(a),byrow=1),2,function(c)	det(matrix(c,ncol=2,byrow=1))	)
		jac<- -jac
		jac.max<- max(jac)
		print(range(jac))
		jac<- jac-jac.max + max(rho2)*0.8
		jac.ax<- c(0,0.2,0.4)-max(rho2)*0.8+jac.max
		print(jac.ax)
		print(range(jac))
		print(length(jac))
		print(length(a))
		
		ylim<- range(c(nu2,rho2))
		f.name<- paste(dir.name,"/nABC.MA1_link.pdf",sep='')
		cat(paste("\nsave ",f.name))
		pdf(f.name,version="1.4",width=pdf.width,height=pdf.height)		
		par(mar=c(4,4,.5,4))
		plot(1,1,xlim=range(a),ylim=ylim,type='n',bty='n',xlab="a",ylab=expression(nu[2]*" and "*rho[2]))
		axis(4,at=c(0,0.2,0.4),labels=round(jac.ax,d=2))
		mtext(expression("| "*partialdiff*"L |"),side=4,at=0.2,line=3)
		lines(a,nu2, lty=ltys[1], col=cols[1])
		lines(a,rho2, lty=ltys[2], col=cols[2])
		lines(a,jac, lty=ltys[3], col=cols[3])
		abline(v=0, col=my.fade.col("black",0.2))
		legend("bottomright",bty='n',border=NA,lty=c(ltys[1],ltys[2],ltys[3]),fill=c(cols[1],cols[2],cols[3]),legend= expression(nu[2],rho[2],"|"*partialdiff*"L|"))
		dev.off()
	}
	if(!is.na(subprog) && subprog==4)		#analytical power
	{
		xma.pa<- 0.1
		xma.sigma<- 1
		ma.n<- 5e3
		tau.u<- 0.09			#take tau= 0.1 on the test statistic scale
		tau.l<- -tau.u
		prior.u<- xma.pa+tau.u		
		prior.l<- xma.pa+tau.l				
		N<- 1e3
		alpha<- 0.01
	
		
		#power of test stat
		#print(project.nABC.movingavg.a2nu(c(prior.u,xma.pa,prior.l)))
		#print(project.nABC.movingavg.a2rho(c(prior.u,xma.pa,prior.l)))
		print(project.nABC.movingavg.rho2a(c(prior.l,prior.u)))
		tmp<- nabc.acf.equivalence.cor(rnorm(ma.n,0,1), leave.out=2)['n']
		print(tmp)
		rho<- seq(tau.l,tau.u,0.001)
		y<- nabc.acf.equivalence.pow(rho, tau.u, alpha, 1/sqrt(tmp-3))
		
		plot(rho,y,type='l',ylim=c(0,1))
		
stop()
		
		#why is there such a bad correspondence betw link and link.mc ?
		#simulations suggest we need a lot of ma.n to get a reliable monte carlo estimate of the z score (rho transformed autocorrelation)
		link.mc<- sapply(1:N,function(i)
				{
					x<-	rnorm(ma.n+1,0,xma.sigma)
					x<- x[-1] + x[-(ma.n+1)]*xma.pa 
					cor.sim<- cor(x[-1],x[-length(x)])
					z.sim<- .5 * log( (1+cor.sim)/(1-cor.sim) )					
				})
		hist(link.mc)
		abline(v=project.nABC.movingavg.a2rho(xma.pa))
		
stop()		
		
		args<- paste("acfequiv",tau,alpha,sep='/')
		
		x<-	rnorm(ma.n+1,0,xma.sigma)
		x<- x[-1] + x[-(ma.n+1)]*ma.pa[1] 
		y<-	rnorm(ma.n+1,0,xma.sigma)
		y<- y[-1] + y[-(ma.n+1)]*ma.pa[1] 
						
		plot(x, type='l')
		print(c( cor(x[-1],x[-length(x)]), var(x) ))				
		nabc.acf.equivalence(x,y,args)
		
		x<- seq(-1,1,0.001)
		y<- .5 * log( (1+x)/(1-x) )
		plot(x,y,type='l')
		
	}

	#TODO 	run nABC with cors and mallows, and show that it performs better because of the link function?
	#TODO 	estimate link function and show how well the estimate captures the truth
	stop()
}
#------------------------------------------------------------------------------------------------------------------------
project.nABC.compareSEIRS<- function()
{	
	my.mkdir(DATA,"nABC.SEIIRScompare")
	d.name<- paste(DATA,"nABC.SEIIRScompare",sep='/')

	if(0)	#plot stuff for paper
	{		
		phyloType<- "GENDIST_TIER2"
		dataType<- "NBH3N2_NL_EU101027"				
		suType<- 8
		cex<- 1.5
		os<- ABC.obsinit(DATAT[dataType])
		
		f.name<- "simuobs_SBRI_STIRS_FDWET_1,6e+07_3,5_1,8_50_0,81_10_0,41_0,482_0_4000_1,5_-6_0,08_0_8_1_1_2_0,9_2_1,7_0,006_4_0,01_1_1_0_0_15_42_0,5_1_1_10_0,01_-1_0_0,5_-1_0.R"		
		cat(paste("load",paste(d.name,f.name,sep='/')))
		load(paste(d.name,f.name,sep='/'))
		
		ws<- epi.whichSummaries(strsplit(names(SUMODTYPE[suType]), '_')[[1]])
		print(names(simuobs[["data"]][[1]]))
		#print(simuobs[["data"]][[1]][["TOTAL.ILI"]])
		
		COLS<- c(my.fade.col("black",0.5),my.fade.col("black",0.8))
		
		
		#plot H3N2 ILI
		oss<- list( os[["TOTAL.ILI"]], simuobs[["data"]][[1]][["TOTAL.ILI"]] )
		xlab<- "year"
		ylab<- "ILI / week / 100 000"
		itupyr<- 1/52
		ref.yr<- 1965
		sim.log<- 1
		
		sim.t.ranges<- sapply(oss,function(os){ os[[2]] })		#row 1: emergenceT, row 2: deathT
		z<- rep(0, diff(range(sim.t.ranges)) )
		sim.t.ranges.pyr<- sim.t.ranges*itupyr	+	ref.yr
		sim.log.pyr<- sim.log*itupyr
		ylim<- range( c(450,sapply(oss,function(os){ range( os[[3]] ) })),na.rm=TRUE )
		cat(paste("\nplot",paste(d.name,"ts_ILI.pdf",sep='/')))
		pdf(paste(d.name,"ts_ILI.pdf",sep='/'),version="1.4",width=10,height=5)
		par(mar=c(5,4.5,0.5,1.5))
		plot(1,1,	xlim=range( as.vector(sim.t.ranges.pyr) ),	ylim=ylim,	type='n',bty='n',	xlab= xlab,	ylab=ylab,		las=1, cex.axis= cex, cex.lab= cex)
		#plot first obss as polygon
		x<- seq(sim.t.ranges.pyr[1,1],	by=sim.log.pyr,	length=		length(oss[[1]][[3]]))
		y<- oss[[1]][[3]]
		polygon( c(x,rev(x)), c(y,rep(0, length(y))), col=COLS[1],border=NA )		
		#plot second obss as polygon
		x<- seq(sim.t.ranges.pyr[1,2],	by=sim.log.pyr,	length=		length(oss[[2]][[3]]))
		y<- oss[[2]][[3]]
		polygon( c(x,rev(x)), c(y,rep(0, length(y))), col="transparent",border=COLS[2] )
		legend("topleft",legend= c("H3N2 regression estimate, Netherlands","H3N2 ILI under compartment model"), fill= c(COLS[1],"white"), bty='n', cex= cex)		
		dev.off()
		
		
		#plot case report attack rate
		oss<- list( os[["ANN.ATT.R"]], simuobs[["data"]][[1]][["ANN.ATT.R"]] )		
		ylab<- "ILI attack rate"
		h<- 0.8
		h2<- c(0, 0, 0.15)
		ref.yr<- 1965
		solstice<- 172
		fade<- 0.8
		
		obs.t.pyr<- lapply(oss, function(x){		as.numeric(names(x))			}	)
		xlim<- range( c( obs.t.pyr, recursive= TRUE)  )
		xlim[2]<- xlim[2]+1
		xlim<- ref.yr + xlim
		ylim<- c(0, 1.05*max(  sapply(oss, max )  )  )
		cat(paste("\nplot",paste(d.name,"ts_ATTR.pdf",sep='/')))
		pdf(paste(d.name,"ts_ATTR.pdf",sep='/'),version="1.4",width=10,height=3.5)		
		par(mar=c(5,5.5,0.5,1.5), mgp= c(3,1,0))
		plot(1,1,xlim=xlim,ylim=ylim,type='n',xlab="season",xaxt= 'n', bty='n',ylab= '', las= 1, cex.axis= cex, cex.lab= cex)
		mtext(side = 2, text = ylab, line = 4, cex= cex)
		xlim<- seq(xlim[1],xlim[2],5)
		axis(1, at=xlim+solstice/DAYSPYEAR, labels= sapply(seq_along(xlim),function(i){	paste(substr(xlim[i]-1,3,nchar(xlim[i])),'/',substr(xlim[i],3,nchar(xlim[i])),sep='')	}), cex.axis= cex, cex.lab= cex		)
		#plot polygon for first data set
		sapply(seq_along(oss[[1]]),function(j)
				{
					polygon(	ref.yr+c(obs.t.pyr[[1]][j]+1-h,  obs.t.pyr[[1]][j]+h,  obs.t.pyr[[1]][j]+h,  obs.t.pyr[[1]][j]+1-h),
								c(oss[[1]][j],oss[[1]][j],0,0)	,
								col=COLS[1], 		border=NA			)
				})
		#plot empty polygon for all other data sets		
		sapply(seq_along(oss[[2]]),function(j)
				{
					polygon(	ref.yr+c(obs.t.pyr[[2]][j]+h2[2]+1-h,  obs.t.pyr[[2]][j]+h2[2]+h,  obs.t.pyr[[2]][j]+h2[2]+h,  obs.t.pyr[[2]][j]+h2[2]+1-h),		
								c(oss[[2]][j],oss[[2]][j],0,0)	,			col="transparent", 		border=COLS[2]			)
				})
		dev.off()
		
		sus<- list(oss[[1]][seq.int(1,length(oss[[1]]),2)], oss[[2]][seq.int(1,length(oss[[2]]),2)] )
		
		x<- seq(min(sus[[1]])*1,max(sus[[1]])*1.5,by=diff(range(sus[[1]]))/500)
		y<- dnorm(x,mean(sus[[1]]), sd=sd(sus[[1]]))
		cat(paste("\nplot",paste(d.name,"ts_ATTR_hist_Neth.pdf",sep='/')))
		pdf(paste(d.name,"ts_ATTR_hist_Neth.pdf",sep='/'),version="1.4",width=4,height=6)				
		par(mar=c(5,4.5,0.5,1.5))
		plot(1,1,xlim=range(x),ylim=range(c(0,y)),type='n',bty='n',xlab='',ylab="density", cex.axis= cex, cex.lab= cex)
		mtext(side = 1, text = ylab, line = 4, cex= cex)
		lines(x,y,lty=1)
		points(sus[[1]],rep(0,length(sus[[1]])),pch=19,col=COLS[1], cex=cex)
		dev.off()
		
		x2<- seq(min(sus[[2]])/1.03,max(sus[[2]])*1.03,by=diff(range(sus[[2]]))/500)
		y2<- dnorm(x2,mean(sus[[2]]), sd=sd(sus[[2]]))
		cat(paste("\nplot",paste(d.name,"ts_ATTR_hist_sim.pdf",sep='/')))
		pdf(paste(d.name,"ts_ATTR_hist_sim.pdf",sep='/'),version="1.4",width=4,height=6)						
		par(mar=c(5,4.5,0.5,1.5))
		plot(1,1,xlim=range(x2),ylim=range(c(0,y2)),type='n',bty='n',xlab='',ylab="density", cex.axis= cex, cex.lab= cex)
		mtext(side = 1, text = ylab, line = 4, cex= cex)
		lines(x2,y2,lty=1)
		points(sus[[2]],rep(0,length(sus[[2]])),pch=21,col=COLS[2], cex=cex)
		dev.off()
		
		
		#plot logs of consecutive case report attack rates
		oss<- list( os[["ANN.ATT.R"]], simuobs[["data"]][[1]][["ANN.ATT.R"]] )
		oss<- lapply(oss, function(x)
				{
					x[x==0]<- min(x[x!=0])
					x<- log(	x[-1]  /  x[ -length(x) ]		)
					x
				})
		ylab<- "log 1st order diff'ces in\n ILI attack rate"
		
		obs.t.pyr<- lapply(oss, function(x){		as.numeric(names(x))			}	)
		xlim<- range( c( obs.t.pyr, recursive= TRUE)  )
		xlim[2]<- xlim[2]+1
		xlim<- ref.yr + xlim
		ylim<- 1.1*range(  sapply(oss, range )  )		
		cat(paste("\nplot",paste(d.name,"ts_FDATTR.pdf",sep='/')))
		pdf(paste(d.name,"ts_FDATTR.pdf",sep='/'),version="1.4",width=10,height=3.5)		
		par(mar=c(5,7,0.5,1.5))
		plot(1,1,xlim=xlim,ylim=ylim,type='n',xlab="between this and following season",xaxt= 'n',bty='n', ylab= '', las= 1, cex.axis= cex, cex.lab= cex)
		mtext(side = 2, text = ylab, line = 4, cex= cex)
		xlim<- seq(xlim[1],xlim[2],5)
		axis(1, at=xlim+solstice/DAYSPYEAR, labels= sapply(seq_along(xlim),function(i){	paste(substr(xlim[i]-2,3,nchar(xlim[i])),'/',substr(xlim[i]-1,3,nchar(xlim[i])),sep='')	}), cex.axis= cex, cex.lab= cex		)
		#plot polygon for first data set
		sapply(seq_along(oss[[1]]),function(j)
				{
					polygon(	ref.yr+c(obs.t.pyr[[1]][j]+1-h,  obs.t.pyr[[1]][j]+h,  obs.t.pyr[[1]][j]+h,  obs.t.pyr[[1]][j]+1-h),
								c(oss[[1]][j],oss[[1]][j],0,0)	,
								col=COLS[1], 		border=NA			)
				})		
		#plot empty polygon for all other data sets		
		sapply(seq_along(oss[[2]]),function(j)
			{
				polygon(	ref.yr+c(obs.t.pyr[[2]][j]+h2[2]+1-h,  obs.t.pyr[[2]][j]+h2[2]+h,  obs.t.pyr[[2]][j]+h2[2]+h,  obs.t.pyr[[2]][j]+h2[2]+1-h),		
							c(oss[[2]][j],oss[[2]][j],0,0)	,			col="transparent", 		border=COLS[2]			)
			})
		abline(h=0)
		dev.off()
		
		sus<- list(oss[[1]][seq.int(2,length(oss[[1]]),2)], oss[[2]][seq.int(2,length(oss[[2]]),2)] )
		
		x<- seq(min(sus[[1]])*1.5,max(sus[[1]])*1.5,by=diff(range(sus[[1]]))/500)
		y<- dnorm(x,mean(sus[[1]]), sd=sd(sus[[1]]))
		cat(paste("\nplot",paste(d.name,"ts_FDATTR_hist_Neth.pdf",sep='/')))
		pdf(paste(d.name,"ts_FDATTR_hist_Neth.pdf",sep='/'),version="1.4",width=4,height=6)				
		par(mar=c(5,4.5,0.5,1.5))
		plot(1,1,xlim=range(x),ylim=range(c(0,y)),type='n',bty='n',xlab='',ylab="density", cex.axis= cex, cex.lab= cex)
		mtext(side = 1, text = ylab, line = 4, cex= cex)
		lines(x,y,lty=1)
		points(sus[[1]],rep(0,length(sus[[1]])),pch=19,col=COLS[1], cex=cex)
		dev.off()
		
		
		x2<- seq(min(sus[[2]])/1.2,max(sus[[2]])*1.2,by=diff(range(sus[[2]]))/500)
		y2<- dnorm(x2,mean(sus[[2]]), sd=sd(sus[[2]]))
		cat(paste("\nplot",paste(d.name,"ts_FDATTR_hist_sim.pdf",sep='/')))
		pdf(paste(d.name,"ts_FDATTR_hist_sim.pdf",sep='/'),version="1.4",width=4,height=6)						
		par(mar=c(5,4.5,0.5,1.5))
		plot(1,1,xlim=range(x2),ylim=range(c(0,y2)),type='n',bty='n',xlab='',ylab="density", cex.axis= cex, cex.lab= cex)
		mtext(side = 1, text = ylab, line = 4, cex= cex)
		lines(x2,y2,lty=1)
		points(sus[[2]],rep(0,length(sus[[2]])),pch=21,col=COLS[2], cex=cex)
		dev.off()

		
		#plot XATT
		oss<- list( NULL, simuobs[["data"]][[1]][["TOTAL.INC.ATT.R.0"]] )
		ylab<- "population-level attack rate"
		sus<- list(NULL, oss[[2]][seq.int(1,length(oss[[2]]),2)] )
		x2<- seq(min(sus[[2]])/1.04,max(sus[[2]])*1.04,by=diff(range(sus[[2]]))/500)
		y2<- dnorm(x2,mean(sus[[2]]), sd=sd(sus[[2]]))
		cat(paste("\nplot",paste(d.name,"ts_XATTR_hist_sim.pdf",sep='/')))
		pdf(paste(d.name,"ts_XATTR_hist_sim.pdf",sep='/'),version="1.4",width=4,height=6)						
		par(mar=c(5,4.5,0.5,1.5))
		plot(1,1,xlim=range(x2),ylim=range(c(0,y2)),type='n',bty='n',xlab='',ylab="density", cex.axis= cex, cex.lab= cex)
		mtext(side = 1, text = ylab, line = 4, cex= cex)
		lines(x2,y2,lty=1)
		points(sus[[2]],rep(0,length(sus[[2]])),pch=21,col=COLS[2], cex=cex)		
		dev.off()
		
		stop()						
	}
	if(0)		#compare power of ABC procedure on simulated data
	{
		alpha<- 0.01
		f.name<- paste(d.name,"power_schuir_simu_SEIRS.R",sep='/')
		cat(paste("\nload file",f.name))
		load(f.name)
		pw<- POWER
		pw<- list(pw[[1]],pw[[2]],pw[[3]],pw[[6]])
		names(pw)<- c(expression(mu*"-attack"),"pop-attack",expression(sigma*"-attack"),"mx-peaks")
		
		
		x<- pw[[2]]
		xname<- "XATT"
		tau<- c( 0.05, 0.02, 0.01 )
		cols<- c("black","gray20","gray40")
		pw.tost<- sapply(seq_along(tau),function(i)
				{
					se<- x[3]
					nm<- x[4]
					df<- x[5]
					PowerTOST:::.power.TOST(alpha=alpha, -tau[i], tau[i], seq(-tau[i], tau[i], length.out= 1e3), se, nm, df, bk = 4)	
				})
		f.name<- paste(d.name,"/nABC.compareSEIRS_simu_pw_",xname,".pdf",sep='')
		cat(paste("\nABC.StretchedF: write pdf to",f.name))
		pdf(f.name,version="1.4",width=4,height=5)
		par(mar=c(4,4,.5,.5))				
		plot(1,1,type='n',xlim=c(-max(tau),max(tau)),ylim=c(0,1),ylab="power",xlab=expression(rho["pop-attack"]))
		sapply(seq_along(tau),function(i)
				{
					x<- seq(-tau[i], tau[i], length.out= 1e3)
					lines(x,pw.tost[,i],col=cols[i],lty=length(tau)-i+1)
				})
		#legend("topleft",bty='o',border=NA,bg="white",box.col="white",lty=seq_along(tau),legend=c(expression(tau["pop-attack"]^'+'==0.01),expression(tau["pop-attack"]^'+'==0.02),expression(tau["pop-attack"]^'+'==0.05)))
		legend("bottomleft",bty='o',border=NA,bg="white",box.col="white",lty=seq_along(tau),legend=c(expression(tau^'+'==0.01),expression(tau^'+'==0.02),expression(tau^'+'==0.05)))
		dev.off()		
	}
	if(0)		#compare posterior histograms of real data
	{
		grace.after.annealing<- 1
		resume<- 1
		if(resume)
		{
			f.name<- paste(d.name,"real_allruns.R",sep='/')
			options(show.error.messages = FALSE, warn=1)		
			readAttempt<-try(suppressWarnings(load(f.name)))						
			options(show.error.messages = TRUE)
		}									
		if(!resume || inherits(readAttempt, "try-error"))
		{
			f.name<- c("abc.ci.mcmc.anneal.SEIIRS_PTPR_NL.B1_fit_Tier1_Ferguson_9pa_v1.01")
			tmp<- paste("nABC.SEIIRScompare",f.name,sep='/')
			
			print(tmp)
			
			post<- lapply(tmp,function(x)
					{
						cat(paste("\nprocess ABC run",x))
						abc.core<- ABC.load( x )
						mabc<- ABCMU.MMCMC.init( x )
						acc<- ABC.MMCMC.get.acceptance(mabc, grace.after.annealing= grace.after.annealing)
						samples<- ABC.MMCMC.getsamples(mabc, grace.after.annealing= grace.after.annealing)
						list(abc.core=abc.core,acc=acc,samples=samples)
					})
			f.name<- paste(d.name,"real_allruns.R",sep='/')
			cat("\nsave runs to ",f.name)
			save(post,file=f.name)
		}
		else
			cat(paste("\nproject.nABC.compareSEIRS: resumed ",f.name))
		
		xname<- "R0"
		yname<- "durImm"
		cols<- c("black","gray20","gray40")
		post<- list(post[[1]])
		xs<- lapply(seq_along(post),function(i)
				{
					c(sapply(seq_along(post[[i]][["samples"]]),function(r)	post[[i]][["samples"]][[r]][,xname]), recursive=1)	
				})	
		ys<- lapply(seq_along(post),function(i)
				{
					c(sapply(seq_along(post[[i]][["samples"]]),function(r)	post[[i]][["samples"]][[r]][,yname]), recursive=1)
				})	
		hxs<- lapply(seq_along(xs),function(i)
				{
					project.nABC.movingavg.gethist(xs[[i]],theta=3.5,nbreaks=100)
				})	
		hys<- lapply(seq_along(ys),function(i)
				{
					project.nABC.movingavg.gethist(ys[[i]],theta=10,nbreaks=37)
				})
		#compare R0 histograms
		f.name<- paste(d.name,"/nABC.compareSEIRS_simu_",xname,".pdf",sep='')
		cat(paste("\nABC.StretchedF: write pdf to",f.name))
		#pdf(f.name,version="1.4",width=4,height=5)
		par(mar=c(4,4,.5,.5))		
		plot(1,1,type='n',xlim=range(sapply(xs,range)),ylim=range(sapply(seq_along(hxs),function(i) hxs[[i]][["density"]] )),ylab="density",xlab=expression(R[0]))	
		sapply(seq_along(post),function(i)
				{			
					plot(hxs[[i]],add=1,freq=0, border=NA, col=my.fade.col(cols[i],0.5))
					abline(v=hxs[[i]][["dmode"]],lty=2,col=my.fade.col(cols[i],0.75))				
				})
		legend("topright",bty='n',border=NA,fill=c("gray30","gray50","gray70"),legend=c(expression(tau^'+'==0.01),expression(tau^'+'==0.02),expression(tau^'+'==0.05)))
		#dev.off()
		
	}
	if(1)		#compare posterior histograms of simulated data
	{
		require(ash)
		grace.after.annealing<- 1
		resume<- 0
		if(resume)
		{
			f.name<- paste(d.name,"simu_allruns.R",sep='/')
			options(show.error.messages = FALSE, warn=1)		
			readAttempt<-try(suppressWarnings(load(f.name)))						
			options(show.error.messages = TRUE)
		}									
		if(!resume || inherits(readAttempt, "try-error"))
		{
			f.name<- c("abc.ci.mcmc.anneal.SEIIRS_PTPR_NL_1_fit_Tier1_Ferguson_9pa_v7","abc.ci.mcmc.anneal.SEIIRS_PTPR_NL_1_fit_Tier1_Ferguson_9pa_v6","abc.ci.mcmc.anneal.SEIIRS_PTPR_NL_1_fit_Tier1_Ferguson_9pa_v0.03","abc.ci.mcmc.anneal.SEIIRS_PTPR_NL_1_fit_Tier1_Ferguson_9pa_v0.025","abc.ci.mcmc.anneal.SEIIRS_PTPR_NL_1_fit_Tier1_Ferguson_9pa_v0.02","abc.ci.mcmc.anneal.SEIIRS_PTPR_NL_1_fit_Tier1_Ferguson_9pa_v0.01")
			f.name<- c("abc.ci.mcmc.anneal.SEIIRS_PTPR_NL_1_fit_Tier1_Ferguson_9pa_v0.015","abc.ci.mcmc.anneal.SEIIRS_PTPR_NL_1_fit_Tier1_Ferguson_9pa_v0.011","abc.ci.mcmc.anneal.SEIIRS_PTPR_NL_1_fit_Tier1_Ferguson_9pa_v0.009")
			f.name<- c("abc.ci.mcmc.anneal.SEIIRS_PTPR_NL_1_fit_Tier1_Ferguson_9pa_v0.009","abc.ci.mcmc.anneal.SEIIRS_PTPR_NL_1_fit_Tier1_Ferguson_9pa_v1")			
			tmp<- paste("nABC.SEIIRScompare",f.name,sep='/')
			
			print(tmp)
			
			post<- lapply(tmp,function(x)
					{
						cat(paste("\nprocess ABC run",x))
						abc.core<- ABC.load( x )
						mabc	<- ABCMU.MMCMC.init( x )
						acc		<- ABC.MMCMC.get.acceptance(mabc, grace.after.annealing= grace.after.annealing)
						samples	<- ABC.MMCMC.getsamples(mabc, grace.after.annealing= grace.after.annealing)
						links	<- ABC.MMCMC.getsamples(mabc, only.nonconst=FALSE, grace.after.annealing= grace.after.annealing, what= "rho.mc")
						list(abc.core=abc.core,acc=acc,samples=samples, links=links)
					})
			f.name<- paste(d.name,"simu_allruns.R",sep='/')
			cat("\nsave runs to ",f.name)
			save(post,file=f.name)
		}
		else
			cat(paste("\nproject.nABC.compareSEIRS: resumed ",f.name))
		
		post	<- list(post[[1]],post[[2]])
		cols	<- c("black","gray20","gray40")
		xnames	<- c("R0","durImm","repProb")
		xlab	<- expression(R[0],1/nu,rho)
		xtrue	<- c(3.5,10,0.08)
		lapply(seq_along(xnames),function(j)
						{
							xname<- xnames[j]
							cat(paste("\nprocess",xname))
							xs<- lapply(seq_along(post),function(i)
									{
										c(sapply(seq_along(post[[i]][["samples"]]),function(r)	post[[i]][["samples"]][[r]][,xname]), recursive=1)	
									})	
							print(range(xs))
							hxs<- lapply(seq_along(xs),function(i)
									{
										project.nABC.movingavg.gethist(xs[[i]],theta=xtrue[j],nbreaks=37)
									})	
							
							f.name<- paste(d.name,"/nABC.compareSEIRS_simu_",xname,".pdf",sep='')
							cat(paste("\nABC.StretchedF: write pdf to",f.name))
							pdf(f.name,version="1.4",width=4,height=5)
							par(mar=c(4,4,.5,.5))		
							plot(1,1,type='n',xlim=range(sapply(xs,range)),ylim=range(sapply(seq_along(hxs),function(i) hxs[[i]][["density"]] )),ylab="density",xlab=xlab[j])
							abline(v=xtrue[j],lty=1,col="red")	
							sapply(seq_along(post),function(i)
									{			
										plot(hxs[[i]],add=1,freq=0, border=NA, col=my.fade.col(cols[i],0.5))
										abline(v=hxs[[i]][["dmode"]],lty=2,col=my.fade.col(cols[i],0.75))
										abline(v=hxs[[i]][["mean"]],lty=3,col=my.fade.col(cols[i],0.75))
									})
							legend("topright",bty='n',border=NA,fill=c("gray30","gray50","gray70"),legend=c(expression(tau^'+'==0.01),expression(tau^'+'==0.02),expression(tau^'+'==0.05)))		
							dev.off()														
						})		
		xnames	<- c("MED.ANN.ATT.R","AMED.FD.ATT.R","MED.INC.ATT.R")
		xlab	<- expression("MED.ANN.ATT.R","AMED.FD.ATT.R","MED.INC.ATT.R")
		xtrue	<- c(0,0,0)				
		lapply(seq_along(xnames),function(j)
						{
							xname<- xnames[j]
							cat(paste("\nprocess",xname))
							xs<- lapply(seq_along(post),function(i)
									{
										c(sapply(seq_along(post[[i]][["links"]]),function(r)	post[[i]][["links"]][[r]][,xname]), recursive=1)	
									})	
							print(range(xs))
							hxs<- lapply(seq_along(xs),function(i)
									{
										project.nABC.movingavg.gethist(xs[[i]],theta=xtrue[j],nbreaks=37)
									})	
							
							f.name<- paste(d.name,"/nABC.compareSEIRS_simu_",xname,".pdf",sep='')
							cat(paste("\nABC.StretchedF: write pdf to",f.name))
							pdf(f.name,version="1.4",width=4,height=5)
							par(mar=c(4,4,.5,.5))		
							plot(1,1,type='n',xlim=range(sapply(xs,range)),ylim=range(sapply(seq_along(hxs),function(i) hxs[[i]][["density"]] )),ylab="density",xlab=xlab[j])
							abline(v=xtrue[j],lty=1,col="red")	
							sapply(seq_along(post),function(i)
									{			
										plot(hxs[[i]],add=1,freq=0, border=NA, col=my.fade.col(cols[i],0.5))
										abline(v=hxs[[i]][["dmode"]],lty=2,col=my.fade.col(cols[i],0.75))
										abline(v=hxs[[i]][["mean"]],lty=3,col=my.fade.col(cols[i],0.75))
									})
							legend("topright",bty='n',border=NA,fill=c("gray30","gray50","gray70"),legend=c(expression(tau^'+'==0.01),expression(tau^'+'==0.02),expression(tau^'+'==0.05)))		
							dev.off()														
						})				
		stop()
	}
}
#------------------------------------------------------------------------------------------------------------------------
project.nABC.StretchedChi2<- function()
{	
	my.mkdir(DATA,"nABC.StretchedChisq")
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
			tmp				<- nabc.chisqstretch(rnorm(yn,ymu,sd=sqrt(ans[["xsigma2"]])), ans[["xsigma2"]], args=args, verbose= 0)
			ans[["cil"]]	<- tmp["cil"]
			ans[["cir"]]	<- tmp["cir"]		
		}
		else
		{
			ans[["cil"]]	<- ans[["cir"]]<- NA
		}
		ans[["data"]]		<- sapply(1:N,function(i)
								{					
									ysigma2	<- runif(1,prior.l,prior.u)
									y		<- rnorm(yn,ymu,sd=sqrt(ysigma2))
									if(with.c)
									{
										tmp	<- nabc.chisqstretch(y, ans[["xsigma2"]], args=args, verbose= 0)						
										tmp	<- c(  ysigma2,tmp[c("error","rho.mc")],var(y)-ans[["xsigma2"]]   )							
									}
									else 
									{					
										tmp	<- c(ysigma2, var(y)/ans[["xsigma2"]], log( var(y)/ans[["xsigma2"]] ), var(y)-ans[["xsigma2"]])
									}
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
		tmp<- nabc.chisqstretch(rnorm(yn,ymu,sd=sd(x)), var(x), args=args, verbose= 0)
		ans[["xsigma2"]]<- var(x)
		ans[["cil"]]<- tmp["cil"]
		ans[["cir"]]<- tmp["cir"]
		ans[["data"]]<- sapply(1:N,function(i)
				{					
					ysigma2<- exp(runif(1,log(prior.l),log(prior.u)))					
					y<- rnorm(yn,ymu,sd=sqrt(ysigma2))
					tmp<- nabc.chisqstretch(y, var(x), args=args, verbose= 0)					
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
		tmp		<- nabc.chisqstretch.tau.lowup(0.9, 2.5, dfy, alpha, for.mle=1)
		tau.l	<- tmp[1]
		tau.u	<- tmp[2]
		c.l		<- tmp[5]
		c.u		<- tmp[6]
		print(tmp)
		y2		<- nabc.chisqstretch.pow(rho, yn, dfy, c.l, c.u)
		
		#abc summary likelihood of sigma2 that peaks at rho.star and tau.u sth max.pw is 0.05
		tmp		<- nabc.chisqstretch.tau.lowup(0.05, 2.5, dfy, alpha, for.mle=1)
		tau.l	<- tmp[1]
		tau.u	<- tmp[2]
		c.l		<- tmp[5]
		c.u		<- tmp[6]
		print(tmp)
		y3		<- nabc.chisqstretch.pow(rho, yn, dfy, c.l, c.u)
		
		
		#abc summary likelihood of sigma2 that peaks at rho.star and tau.u sth max.pw is 0.9 and with matching variance
		tmp		<- nabc.chisqstretch.n.of.y(xn, sqrt(var.Sx), 0.9, alpha, tau.u.ub=tau.u, for.mle=1)
		print(tmp)
		yn		<- tmp[1]
		tau.l	<- tmp[2]
		tau.u	<- tmp[3]
		c.l		<- tmp[4]
		c.u		<- tmp[5]
		y4		<- nabc.chisqstretch.pow(rho, yn, yn-1, c.l, c.u)
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
							tau.l<- nabc.chisqstretch.tau.low(tau.u, yn-1, alpha)
							rej<- .Call("abcScaledChiSq",	c(yn-1,yn-1,tau.l,tau.u,alpha,1e-10,100,0.05)	)
							ans.ok[["cil"]]<- rej[1]
							ans.ok[["cir"]]<- rej[2]
							acc.ok<- which( ans.ok[["data"]]["error",]<=ans.ok[["cir"]]  &  ans.ok[["data"]]["error",]>=ans.ok[["cil"]] )
							acc.h.ok<- project.nABC.movingavg.gethist(ans.ok[["data"]]["ysigma2",acc.ok], ans.ok[["xsigma2"]], nbreaks= 50, width= 0.5, plot=0)
							out["fx.tau.u",]<- c(tau.l, tau.u, rej[1], rej[2], length(acc.ok)/ncol(ans.ok[["data"]]), acc.h.ok[["mean"]],acc.h.ok[["hmode"]],acc.h.ok[["dmode"]],ans.ok[["xsigma2"]])
							#determine tolerances sth TOST power is 0.95
							tmp<- nabc.chisqstretch.tau.lowup(mx.pw, 2.5, yn-1, alpha)
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
			plot(1,1,type='n',bty='n',log='x',xlim=range(sample.size),ylim=ylim, xlab="sample size n",ylab=expression("mode of "*pi[tau]*'('*sigma^2*'|'*x*')'))
			abline(h=0.1,lty=2, col="gray60")				
			points(sample.size,sapply(ans, function(x) x["fx.pw","dmode"]),		pch=18,cex=0.5,col=cols[1])
			points(sample.size,sapply(ans, function(x) x["fx.tau.u","dmode"]),	pch=20,cex=0.75,col=cols[2])	
			legend("topleft",bty='n',border=NA,fill=c(cols[1],"transparent",cols[2],"transparent","transparent","transparent"),lty=c(ltys[1],NA,ltys[2],NA,NA,ltys[3]),legend=c("fixed power &",expression("decreasing "*tau^'+'),"increasing power &",expression("fixed "*tau^'+'),"",expression(hat(nu)[x])))
			abline(h=1,lty=4)
			dev.off()

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
		tau.l<- nabc.chisqstretch.tau.low(tau.u, df, alpha)
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
							out<- matrix(NA,2,11,dimnames=list(c("mean","mode"),c("o0.8","o0.4","o0.2","o0.1","o0.05","l0.8","l0.4","l0.2","l0.1","l0.05","xsigma2")))
							cat(paste("\nload",f.name[i]))
							readAttempt<-try(suppressWarnings(load( f.name[i] )))
							if(inherits(readAttempt, "try-error"))	
								return(out)
							
							cus<- c(0.8,0.4,0.2,0.1,0.05)							
							out[,1:5]<- sapply(cus,function(cu)
									{
										#print(cu)
										acc<- which( ans.ok[["data"]]["sy2-sx2",]<=cu  &  ans.ok[["data"]]["sy2-sx2",]>=-cu )
										#print(length(acc))
										tmp<- project.nABC.movingavg.gethist( ans.ok[["data"]]["ysigma2",acc], xsigma2, nbreaks, .5, plot=0)
										#print(c(tmp[["mean"]],tmp[["dmode"]]))
										c(tmp[["mean"]],tmp[["dmode"]])
									})
							out[,6:10]<- sapply(cus,function(cu)
									{
										#print(cu)
										acc<- which( ans.ok[["data"]]["rho.mc",]<=cu  &  ans.ok[["data"]]["rho.mc",]>=-cu )
										#print(length(acc))
										tmp<- project.nABC.movingavg.gethist( ans.ok[["data"]]["ysigma2",acc], xsigma2, nbreaks, .5, plot=0)
										#print(c(tmp[["mean"]],tmp[["dmode"]]))
										c(tmp[["mean"]],tmp[["dmode"]])
									})
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
		tau.l	<- nabc.chisqstretch.tau.low(tau.u, df, alpha, for.mle=for.mle)
		tau.h	<- 0.65

		ymu<- xmu<- 0
		xsigma2<- 1
		prior.u<- 4
		prior.l<- 0.2
		N<- 2e6
		
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
#bookmark7
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
			if(!resume || inherits(readAttempt, "try-error"))
			{				
				f.name<- list.files(dir.name, pattern=paste("^nABC.Chisq_mle_ok_",sep=''), full.names = TRUE)
				tmp<- sort(sapply(strsplit(f.name,'_',fixed=1),function(x)	as.numeric(substr(x[length(x)],2,nchar(x[length(x)])-2))		), index.return=1)
				f.name<- f.name[tmp$ix]				
				f.name2<- list.files(dir.name, pattern=paste("^nABC.Chisq_mle_naive_",sep=''), full.names = TRUE)
				tmp2<- sort(sapply(strsplit(f.name2,'_',fixed=1),function(x)	as.numeric(substr(x[length(x)],2,nchar(x[length(x)])-2))		), index.return=1)
				f.name2<- f.name2[tmp2$ix]
				#f.name3<- list.files(dir.name, pattern=paste("^nABC.Chisq_wprior_",sep=''), full.names = TRUE)
				#tmp3<- sort(sapply(strsplit(f.name3,'_',fixed=1),function(x)	as.numeric(substr(x[length(x)],2,nchar(x[length(x)])-2))		), index.return=1)
				#f.name3<- f.name3[tmp3$ix]											
				f.name<- rbind( 	f.name[tmp$x%in%intersect(tmp$x,tmp2$x)], f.name2[tmp2$x%in%intersect(tmp$x,tmp2$x)]	)
				#f.name<- t(as.matrix(f.name))
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
							#accept if T in boundaries					
							acc.ok<- which( ans.ok[["data"]]["error",]<=ans.ok[["cir"]]  &  ans.ok[["data"]]["error",]>=ans.ok[["cil"]] )
							acc.h.ok<- project.nABC.movingavg.gethist(ans.ok[["data"]]["ysigma2",acc.ok], ans.ok[["xsigma2"]], nbreaks= 100, width= 0.5, plot=0)
							out["ok",]<- c(acc.h.ok[["mean"]],acc.h.ok[["hmode"]],acc.h.ok[["dmode"]],ans.ok[["xsigma2"]])
							
							
							cat(paste("\nload",f.name[2,j]))	
							#ans.naive<- ans.ok
							readAttempt<-try(suppressWarnings(load( f.name[2,j] )))
							if(inherits(readAttempt, "try-error"))	stop("error at naive")	
							#tmp fix bug (now resolved)
				ans.naive[["cil"]]<- 0.5084666
				ans.naive[["cir"]]<- 1.009202
				
							acc.naive<- which( ans.naive[["data"]]["error",]<=ans.naive[["cir"]]  &  ans.naive[["data"]]["error",]>=ans.naive[["cil"]] )
							acc.h.naive<- project.nABC.movingavg.gethist(ans.naive[["data"]]["ysigma2",acc.naive], ans.naive[["xsigma2"]], nbreaks= 100, width= 0.5, plot=0)
							out["naive",]<- c(acc.h.naive[["mean"]],acc.h.naive[["hmode"]],acc.h.naive[["dmode"]],ans.naive[["xsigma2"]])
				print(length(acc.naive) / ncol(ans.naive[["data"]]))			
#print(out)										
							if(0 && j==2)
							{
								require(pscl)
								cols<- c(my.fade.col("black",0.2),my.fade.col("black",0.6),"black")
								ltys<- c(1,1,4,3)
								
								#plot rho
								rho.h.ok	<- project.nABC.movingavg.gethist(ans.ok[["data"]]["ysigma2",acc.ok]/ans.ok[["xsigma2"]], 1, nbreaks= 100, width= 0.5, plot=1)
								
								pw<- nabc.chisqstretch.pow(seq(prior.l,prior.u,by=0.001),yn,df,ans.ok[["cil"]],ans.ok[["cir"]])
								print( c(seq(prior.l,prior.u,by=0.001)[ which.max(pw) ], sum(pw)*0.001 ) )
								lines(seq(prior.l,prior.u,by=0.001),pw/(sum(pw)*0.001),type='l',col="red")
								abline(v=1,col="red")
								
								#rho.h.ok	<- project.nABC.movingavg.gethist(ans.ok[["data"]]["error",], 1, nbreaks= 100, width= 0.5, plot=1)
								print(rho.h.ok[["dmode"]])
								stop()
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
								#pdf(f.name,version="1.4",width=4,height=5)
								par(mar=c(5,5,0.5,0.5))
								plot(acc.h.ok, col=cols[1],border=NA,main='',freq=0,ylab=expression("n-ABC estimate of "*pi[tau]*'('*sigma^2*'|'*x*')'),xlab=expression(sigma^2),xlim=c(0,3),ylim=c(0,2))
								plot(acc.h.naive, col=cols[2],border=NA,main='',add=1,freq=0)
								x<- seq(prior.l,prior.u,0.001)
								y<- densigamma(x,(xn-2)/2,ans.ok[["xsigma2"]]*xn/2) / diff(pigamma(c(prior.l,prior.u),(xn-2)/2,ans.ok[["xsigma2"]]*xn/2))							
								lines(x,y,col=cols[3],lty=ltys[4])
								abline(v=ans.ok[["xsigma2"]],col=cols[3],lty=ltys[3])								
								legend("topright",fill=c("transparent","transparent",cols[1],"transparent","transparent","transparent",cols[2],"transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent"),lty=c(NA,NA,ltys[1],NA,NA,NA,ltys[2],NA,NA,NA,NA,ltys[4],NA,ltys[3],NA),border=NA,bty='n',legend=expression("n=60","","calibrated","tolerances",tau^'-'*"=0.477", tau^'+'*"=2.2","naive","tolerances",tau^'-'*"=0.35",tau^'+'*"=1.65","",pi*'('*sigma^2*'|'*x*')',"",argmax[sigma^2],pi*'('*sigma^2*'|'*x*')'))								
								#dev.off()
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
			print( mean(est.theta0["ok.me",]-est.theta0["xsigma2",]) )
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
		tau.l<- nabc.chisqstretch.tau.low(tau.u, df, alpha)
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
				#f.name<- rbind( 	f.name[tmp$x%in%intersect(tmp$x,tmp2$x)], f.name2[tmp2$x%in%intersect(tmp$x,tmp2$x)]	)	
				f.name<- rbind( 	f.name[tmp$x%in%intersect(intersect(tmp$x,tmp2$x),tmp3$x)], 
									f.name2[tmp2$x%in%intersect(intersect(tmp$x,tmp2$x),tmp3$x)],
									f.name3[tmp3$x%in%intersect(intersect(tmp$x,tmp2$x),tmp3$x)]	)				
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
							
							cat(paste("\nload",f.name[3,j]))
							readAttempt<-try(suppressWarnings(load( f.name[3,j] )))
							if(inherits(readAttempt, "try-error"))	stop("error at wprior")														
							acc.wprior<- which( ans.wprior[["data"]]["error",]<=ans.wprior[["cir"]]  &  ans.wprior[["data"]]["error",]>=ans.wprior[["cil"]] )
							acc.h.wprior<- project.nABC.movingavg.gethist(ans.wprior[["data"]]["ysigma2",acc.wprior], ans.wprior[["xsigma2"]], nbreaks= 50, width= 0.5, plot=0)
							out["wprior",]<- c(acc.h.wprior[["mean"]],acc.h.wprior[["hmode"]],acc.h.wprior[["dmode"]],ans.wprior[["xsigma2"]])
							
							if(0 && j==10)
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
								plot(acc.h.ok, col=cols[1],border=NA,main='',freq=1,ylab="n-ABC accepted samples",xlab=expression(sigma^2),xlim=c(0,3))
								plot(acc.h.naive, col=cols[2],border=NA,main='',add=1,freq=1)
								abline(v=ans.ok[["xsigma2"]],col=cols[3],lty=ltys[3])								
								legend("topright",fill=c("transparent","transparent",cols[1],"transparent","transparent","transparent",cols[2],"transparent","transparent","transparent","transparent","transparent"),lty=c(NA,NA,ltys[1],NA,NA,NA,ltys[2],NA,NA,NA,NA,ltys[3]),border=NA,bty='n',legend=expression("n=60","","calibrated","tolerances",tau^'-'*"=0.477", tau^'+'*"=2.2","naive","tolerances",tau^'-'*"=0.35",tau^'+'*"=1.65","",hat(nu)[x]))								
								dev.off()
								stop()
							}							
							out			
						})
								
				f.name<- paste(dir.name,"/nABC.Chisq_modemean_",N,"_",xn,"_",prior.u,"_",prior.l,"_",tau.u,".R",sep='')
				cat(paste("\nnABC.Chisq save 'ans' to ",f.name))				
				save(ans,file=f.name)
				print(ans[,1:5])
				stop()
			}			
			#
			#compute means 
			est.theta0<- 	sapply(seq_along(ans),function(i)	c(ans[[i]]["ok","xsigma2"],ans[[i]]["wprior","xsigma2"],ans[[i]]["ok","dmode"],ans[[i]]["ok","mean"],ans[[i]]["naive","dmode"],ans[[i]]["naive","mean"],ans[[i]]["wprior","dmode"])	)
			rownames(est.theta0)<- c("xsigma2","wprior.xsigma2","ok.mo","ok.me","naive.mo","naive.me","wprior.mo")
			cat("\n estimated means\n")
			print( apply(est.theta0,1,mean) )
			print( mean(est.theta0["ok.mo",]-est.theta0["xsigma2",]) )
			print( mean(est.theta0["ok.me",]-est.theta0["xsigma2",]) )
			print( mean(est.theta0["naive.mo",]-est.theta0["xsigma2",]) )
			print( mean(est.theta0["naive.me",]-est.theta0["xsigma2",]) )
			print( mean(est.theta0["wprior.mo",]-est.theta0["wprior.xsigma2",]) )

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
			legend("topright",bty='n',border=NA,fill=c("transparent","transparent",cols[1],"transparent",cols[2],"transparent","transparent","transparent"),lty=c(NA,NA,ltys[1],NA,ltys[2],NA,NA,ltys[3]),legend=expression("n=60","","calibrated","tolerances","naive","tolerances","",hat(nu)[x]))
			dev.off()						
		}	
	}
	if(!is.na(subprog) && subprog==3)		#illustrate power of scaled ChiSquare
	{
		print("illustrate power")		
		alpha<- 0.01		
		tau.up<- 2.2 #2.95				
		plot.pdf<- 1 
		verbose<- 1
		
		tau.up<- 2.2 #1.09
		yn<- 60 #5e3
		scale<- yn
		df<- yn-1
		tau.low<- nabc.chisqstretch.tau.low(tau.up, df, alpha, for.mle=1)
		#tau.low<- 0.35
		rej<- .Call("abcScaledChiSq",	c(scale,df,tau.low,tau.up,alpha,1e-10,100,0.05)	)
		pw<- nabc.chisqstretch.pow(seq(tau.low,tau.up,by=0.001),scale,df,rej[1],rej[2])
		print( seq(tau.low,tau.up,by=0.001)[ which.max(pw) ] )
		plot(seq(tau.low,tau.up,by=0.001),pw,type='l')
		print(c(tau.low,rej))
		stop()
		
		#plot height of power
		yns<- c(10,30,60,100)
		tau.low<- sapply(yns, function(x)	nabc.chisqstretch.tau.low(tau.up, x-1, alpha)	)
		pw<- lapply(seq_along(yns),function(i)
				{
					rej<- .Call("abcScaledChiSq",	c(yns[i]-1,yns[i]-1,tau.low[i],tau.up,alpha,1e-10,100,0.05)	)
					nabc.chisqstretch.pow(seq(tau.low[i],tau.up,by=0.001),yns[i]-1,yns[i]-1,rej[1],rej[2])					
				})	
		yn<- 60
		tau.ups<- c(1.5,1.7,2,3)
		tau.low2<- sapply(tau.ups, function(x)	nabc.chisqstretch.tau.low(x, yn-1, alpha)	)
		pw2<- lapply(seq_along(tau.ups),function(i)
				{
					rej<- .Call("abcScaledChiSq",	c(yn-1,yn-1,tau.low2[i],tau.ups[i],alpha,1e-10,100,0.05)	)
					nabc.chisqstretch.pow(seq(tau.low2[i],tau.ups[i],by=0.001),yn-1,yn-1,rej[1],rej[2])					
				})	
				
		#xlim<- range(c(tau.low,tau.up,tau.ups))
		xlim<- range(c(tau.low,tau.up))
		ylim<- range(c(pw,pw2))
		ltys<- c(2,3,1,4)
		if(plot.pdf)
		{
			f.name<- paste(dir.name,"/nABC.Chisq_power2.pdf",sep='')
			pdf(f.name,version="1.4",width=pdf.width,height=pdf.height)
		}
		par(mar=c(5,4,0.5,0.5))
		plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim,xlab=expression(rho),ylab="power")
		sapply(seq_along(yns),function(i)
				{
					rho<- seq(tau.low[i],tau.up,by=0.001)
					lines(rho,pw[[i]],lty=ltys[i])
				})
		
		legend("topright",bty='n',lty=c(NA,NA,NA,ltys),legend=expression("calibrated",chi^2*"-test","","n=10","n=30","n=60","n=100"))
		#sapply(seq_along(tau.ups),function(i)
		#		{
		#			rho<- seq(tau.low2[i],tau.ups[i],by=0.001)
		#			lines(rho,pw2[[i]], lty=4)
		#		})
		if(plot.pdf)	
			dev.off()
		
		
		yn<- 60
		df<- yn-1
		tau.low<- nabc.chisqstretch.tau.low(tau.up, df, alpha)
		#plot location of power
		rej<- .Call("abcScaledChiSq",	c(df,df,tau.low,tau.up,alpha,1e-10,100,0.05)	)						
		rho<- seq(tau.low,tau.up,by=0.001)
		pw<- nabc.chisqstretch.pow(rho,df,df,rej[1],rej[2])
		if(verbose)	cat(paste("\ncase: rho.max at 1\ntau.low is",tau.low,"\ttau.up is",tau.up,"\tcl is",rej[1],"\tcu is",rej[2]))
		if(verbose)	cat(paste("\nrho.max is",rho[ which.max(pw) ],"\tpow.max is",pw[ which.max(pw) ]))
				
		pw.f<- project.nABC.StretchedF.pow.sym(rho, df, rej[2], cl= rej[1])
		
		rej<- .Call("abcScaledChiSq",	c(df,df,1/tau.up,tau.up,alpha,1e-10,100,0.05)	)
		rho.sym<- seq(1/tau.up,tau.up,by=0.001)		
		pw.sym<- nabc.chisqstretch.pow(rho.sym,df,df,rej[1],rej[2])
		if(verbose)	cat(paste("\n\n\ncase: tau.low=1/tau.up\ntau.low is",1/tau.up,"\ttau.up is",tau.up,"\tcl is",rej[1],"\tcu is",rej[2]))
		if(verbose)	cat(paste("\nrho.max is",rho.sym[ which.max(pw.sym) ],"\tpow.max is",pw.sym[ which.max(pw.sym) ]))
		
		h<- 0.65	#diff(c(tau.low,tau.up))/2
		rej<- .Call("abcScaledChiSq",	c(df,df,1-h,1+h,alpha,1e-10,100,0.05)	)
		rho.diff<- seq(1-h,1+h,by=0.001)		
		pw.diff<- nabc.chisqstretch.pow(rho.diff,df,df,rej[1],rej[2])
		if(verbose)	cat(paste("\n\n\ncase: +-h\ntau.low is",1-h,"\ttau.up is",1+h,"\tcl is",rej[1],"\tcu is",rej[2]))
		if(verbose)	cat(paste("\nrho.max is",rho.diff[ which.max(pw.diff) ],"\tpow.max is",pw.diff[ which.max(pw.diff) ]))
		
		cols<- c(my.fade.col("black",0.2),my.fade.col("black",0.6),"black")		
		ltys<- c(1,1,4)
		xlim<- range(c(rho.diff,rho.sym,rho))
		ylim<- range(c(pw,pw.sym,pw.diff))
		if(plot.pdf)
		{
			f.name<- paste(dir.name,"/nABC.Chisq_power.pdf",sep='')
			pdf(f.name,version="1.4",width=pdf.width,height=pdf.height)
		}
		par(mar=c(5,4,0.5,0.5))		
		plot(1,1,xlim=xlim,ylim=ylim,type='n',bty='n',xlab=expression(rho),ylab="power")
		#lines(rho.diff,pw.diff,col=cols[2],lty=ltys[2])
		polygon(c(rho.diff,rho.diff[length(rho.diff)],rho.diff[1]), c(pw.diff,0,0), border=NA,col=cols[2] )
		polygon(c(rho,rho[length(rho)],rho[1]), c(pw,0,0), border=NA,col=cols[1] )
		#lines(rho,pw,col=cols[1],lty=ltys[1])
		#lines(rho.sym,pw.sym,col="gray20")
		
		lines(rho,pw.f,col=cols[3],lty=ltys[3])		
		lines(rep(rho[which.max(pw.f)],2), c(-1,pw.f[which.max(pw.f)]), col=cols[3],lty=ltys[3])
		
		lines(rep(rho[which.max(pw)],2), c(-1,pw[which.max(pw)]), col=cols[1],lty=ltys[1])
		lines(rep(rho.diff[which.max(pw.diff)],2), c(-1,pw.diff[which.max(pw.diff)]), col=cols[2],lty=ltys[2])
		
		#legend<- expression('['*c^'-'*", "*c^'+'*"]= [0.5,1.5]", '['*c^'-'*", "*c^'+'*"]= [1/1.5,1.5]",'['*c^'-'*", "*c^'+'*"]= [0.5,1/0.5]")				
		legend("topright",fill=c(cols[1],"transparent",cols[2],"transparent",cols[3],"transparent"),lty=c(ltys[1],NA,ltys[2],NA,ltys[3],NA),legend=expression("calibrated",chi^2*"-test","naive",chi^2*"-test","F-test",""), bty= 'n', border=NA)
		if(plot.pdf)	
			dev.off()				
	}
	if(!is.na(subprog) && subprog==4)		#illustrate scaled ChiSquare
	{
		yn<- 60		
		tau.up<- 2.2
		df<- yn-1
		alpha<- 0.01
		tau.low<- nabc.chisqstretch.tau.low(tau.up, df, alpha) 		
		verbose<- 1
		
		rej<- .Call("abcScaledChiSq",	c(df,df,tau.low,tau.up,alpha,1e-10,100,0.05)	)
		#print(rej)
		
		cols<- c("gray60","gray80","black","black","black")
		xu<- seq(1/1.5,3,0.001)
		xl<- seq(0,1.5,0.001)		
		yu<- dchisq(xu/tau.up*df,df)
		#yu<- yu / mean(yu)
		yl<- dchisq(xl/tau.low*df,df)
		#yl<- yl / mean(yl)
	
		f.name<- paste(dir.name,"/nABC.Chisq_example.pdf",sep='')
		pdf(f.name,version="1.4",width=pdf.width,height=pdf.height)
		par(mar=c(5,4,0.5,0.5))
		plot(1,1,type='n',bty='n',xlim=range(c(xu,xl)),ylim=range(c(yu,yl,-0.001)),ylab="scaled density",xlab="T")
		cl<- seq(rej[1],1,0.001)
		cu<- seq(1,rej[2],0.001)
		polygon( c(cl,cu,rej[2:1]), c(dchisq(cl/tau.low*df,df),dchisq(cu/tau.up*df,df),0,0), col=cols[1], border=NA)
		polygon( c(tau.low,tau.up,tau.up,tau.low), c(-0.0005,-0.0005,-1,-1), col=cols[2], border=NA )		
		lines(xu,yu, col=cols[3], lty=4)		
		lines(xl,yl, col=cols[4], lty=1)
		#abline(v=1,lty=2)
		
		legend<- expression( f[n-1](T/tau^'-'), f[n-1](T/tau^'+'),'['*c^'-'*", "*c^'+'*"]",'['*tau^'-'*", "*tau^'+'*"]","")						
		legend(x=1.38,y=0.01,lty=c(1,4,NA,NA,NA),fill=c(cols[3],cols[4],cols[1],cols[2],"transparent"),legend=legend, bty= 'n', border=NA)
		dev.off()
		
		if(verbose)	cat(paste("\ncase: rho.max at 1\ntau.low is",tau.low,"\ttau.up is",tau.up,"\tcl is",rej[1],"\tcu is",rej[2]))				
	}
	
	
	stop()
}
#------------------------------------------------------------------------------------------------------------------------
#power in the rejection region
project.nABC.StretchedF.pow.sym<- function(rho, df, cu, cl= 1/cu) pf( cu / rho, df1= df, df2= df ) - pf( cl / rho, df1= df, df2= df )
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
		my.mkdir(DATA,"nABC.StretchedF")
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
		my.mkdir(DATA,"nABC.StretchedF")
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
			my.mkdir(DATA,"nABC.StretchedF.sigmainference")
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
			my.mkdir(DATA,"nABC.StretchedF.sigmainference")			
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
			my.mkdir(DATA,"nABC.StretchedF.sigmainference")
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
#------------------------------------------------------------------------------------------------------------------------
project.nABC.TOST<- function()
{
	require(PowerTOST)	
	my.mkdir(DATA,"nABC.mutost")
	dir.name<- paste(DATA,"nABC.mutost",sep='/')
	subprog<- 5
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
	project.nABC.mutost.fix.x.uprior.ysig2<- function(N,tau,prior,alpha,x,yn,stdize)		
	{		
		if(!is.matrix(tau)	|| nrow(tau)!=1	|| ncol(tau)!=2)	
			stop("project.nABC.mutost.fix.x.uprior.ysig2: error at 1a")		
		if(!is.matrix(prior)	|| nrow(prior)!=2	|| ncol(prior)!=2)	
			stop("project.nABC.mutost.fix.x.uprior.ysig2: error at 1b")
		
		args<- paste("mutost",stdize,tau["mu","l"],tau["mu","u"],alpha,sep='/')
		#perform one ABC - rejection run
		ans				<- vector("list",5)
		names(ans)		<- c("xmu","xsigma2","cil","cir","data")		
		ans[["xmu"]]	<- mean(x)
		ans[["xsigma2"]]<- var(x)		
		tmp				<- nabc.mutost.onesample(rnorm(yn,mean(x),sd=sd(x)), mean(x), std.sd=sd(x), args= args, verbose= 0)
		ans[["cil"]]	<- tmp[["cil"]]
		ans[["cir"]]	<- tmp[["cir"]]		
		ans[["data"]]	<- sapply(1:N,function(i)
				{					
					ymu			<- runif(1,prior["mu","l"],prior["mu","u"])
					ysigma2		<- runif(1,prior["sig2","l"],prior["sig2","u"])
					tmp			<- nabc.mutost.onesample(rnorm(yn,ymu,sd=sqrt(ysigma2)), mean(x), std.sd=sd(x), args= args, verbose= 0)
					tmp			<- c(ymu,ysigma2,tmp[c("error","rho.mc")])
					names(tmp)	<- c("ymu","ysigma2","error","rho.mc")					
					tmp					
				})		
		ans
	}	
	project.nABC.mutost.fix.x.fix.ysig2<- function(N,tau,prior,alpha,x,yn,stdize)		
	{		
		if(!is.matrix(tau)	|| nrow(tau)!=1	|| ncol(tau)!=2)	
			stop("project.nABC.mutost.fix.x.uprior.ysig2: error at 1a")		
		if(!is.matrix(prior)	|| nrow(prior)!=2	|| ncol(prior)!=2)	
			stop("project.nABC.mutost.fix.x.uprior.ysig2: error at 1b")
		
		if(stdize!=2)
			args<- paste("mutost",stdize,tau["mu","l"],tau["mu","u"],alpha,sep='/')
		else
			args<- paste("mutost",stdize,1,0.9,alpha,sep='/')
		#perform one ABC - rejection run
		ans				<- vector("list",5)
		names(ans)		<- c("xmu","xsigma2","cil","cir","data")		
		ans[["xmu"]]	<- mean(x)
		ans[["xsigma2"]]<- var(x)
		tmp				<- nabc.mutost.onesample(x, mean(x), std.sd=sd(x), args= args, verbose= 0)
		ans[["cil"]]	<- tmp[["cil"]]
		ans[["cir"]]	<- tmp[["cir"]]		
		ans[["data"]]	<- sapply(1:N,function(i)
				{					
					ymu			<- runif(1,prior["mu","l"],prior["mu","u"])
					tmp			<- rnorm(yn,0,1)
					tmp			<- tmp/sd(tmp)*sqrt(ans[["xsigma2"]])+ymu
					tmp			<- nabc.mutost.onesample(tmp, mean(x), std.sd=sd(x), args= args, verbose= 0)
					tmp			<- c(ymu,ans[["xsigma2"]],tmp[c("error","rho.mc")])
					names(tmp)	<- c("ymu","ysigma2","error","rho.mc")
					tmp					
				})		
		ans
	}	
	
	if(!is.na(subprog) && subprog==2)
	{
		xn<- yn	<- 60
		alpha	<- 0.01	
		tau		<- 0.5
		tau		<- matrix(c(-tau,tau),ncol=2,dimnames=list(c("mu"),c("l","u"))) 		
		prior	<- matrix(c(-5,5,0.5,4),ncol=2,byrow=1,dimnames=list(c("mu","sig2"),c("l","u")))
		
		xmu		<- 0.5
		xsigma2	<- 2
		N		<- 2e6
		stdize	<- 2
		m		<- NA
		resume	<- 0
		if(!is.na(m))
		{		
			f.name<- paste(dir.name,"/nABC.mutost_unbiasedrepeat_",N,"_",xn,"_",stdize,"_",prior[1,1],"_",prior[1,2],"_",prior[2,1],"_",prior[2,2],"_",tau[1,2],"_m",m,".R",sep='')
			cat(paste("\nnABC.mutost: compute ",f.name))
			options(show.error.messages = FALSE, warn=1)		
			readAttempt<-try(suppressWarnings(load(f.name)))						
			options(show.error.messages = TRUE)						
			if(!resume || inherits(readAttempt, "try-error"))
			{
				x	<- rnorm(xn,0,1)
				x	<- x/sd(x)*sqrt(xsigma2)+xmu						
				ans	<- project.nABC.mutost.fix.x.uprior.ysig2(N,tau,prior,alpha,x,yn,stdize)				
				cat(paste("\nnABC.mutost: save ",f.name))
				save(ans,file=f.name)				
			}
			else
				cat(paste("\nnABC.mutost: resumed ",f.name))
		}
		else
		{
			#load data					
			cat(paste("\nnABC.mutost",dir.name))
			f.name<- paste(dir.name,"/nABC.mutost_unbiased_",N,"_",xn,"_",prior[1,1],"_",prior[1,2],"_",prior[2,1],"_",prior[2,2],"_",tau[1,2],".R",sep='')			
			options(show.error.messages = FALSE, warn=1)		
			readAttempt<-try(suppressWarnings(load(f.name)))						
			options(show.error.messages = TRUE)		
			if(!resume || inherits(readAttempt, "try-error"))
			{				
				f.name<- list.files(dir.name, pattern=paste("^nABC.mutost_unbiasedrepeat_",sep=''), full.names = TRUE)
				f.name.yn<- sort(sapply(strsplit(f.name,'_',fixed=1),function(x)	as.numeric(x[length(x)-2])		), index.return=1)
				f.name<- f.name[f.name.yn$ix]				
				cat(paste("\nnABC.mutost load data: ", length(f.name)))
				ans<- lapply(seq_along(f.name),function(j)
						{														
							cat(paste("\nload",f.name[j]))
							readAttempt<-try(suppressWarnings(load( f.name[j] )))
							if(inherits(readAttempt, "try-error"))	stop("error at unbiased")
							#print(ans)
							#print(any(is.na(ans[["data"]])))
							#print(which(apply(ans[["data"]],2,function(x) any(is.na(x)))))
							#print(c(ans[["xmu"]],ans[["cil"]],ans[["cir"]]))
							#determine tolerances for tau.u=2.2
							acc<- which( ans[["data"]]["error",]<=ans[["cir"]]  &  ans[["data"]]["error",]>=ans[["cil"]] )
							print(length(acc)/ncol(ans[["data"]]))
							print(length(acc))
							hist<- project.nABC.movingavg.gethist(ans[["data"]]["ymu",acc], ans[["xmu"]], nbreaks= 50, width= 0.5, plot=1)
							abline(v=ans[["xmu"]],col="red")
stop()
							out["fx.tau.u",]<- c(tau.l, tau.u, rej[1], rej[2], length(acc.ok)/ncol(ans.ok[["data"]]), acc.h.ok[["mean"]],acc.h.ok[["hmode"]],acc.h.ok[["dmode"]],ans.ok[["xsigma2"]])
							#determine tolerances sth TOST power is 0.95
							tmp<- nabc.chisqstretch.tau.lowup(mx.pw, 2.5, yn-1, alpha)
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
		}
		stop()
	}
	if(!is.na(subprog) && subprog==5)		#check effect of n != m
	{
		xn		<- 60
		yn		<- 60
		alpha	<- 0.01	
		tau		<- 0.75
		tau		<- matrix(c(-tau,tau),ncol=2,dimnames=list(c("mu"),c("l","u"))) 		
		prior	<- matrix(c(-5,5,0.5,4),ncol=2,byrow=1,dimnames=list(c("mu","sig2"),c("l","u")))
		
		xmu		<- 0.5
		xsigma2	<- 2
		N		<- 1e6
		stdize	<- 0
		m		<- 1
		resume	<- 1
		if(1)
		{		
			x	<- rnorm(xn,0,1)
			x	<- x/sd(x)*sqrt(xsigma2)+xmu
			
			f.name<- paste(dir.name,"/nABC.mutost_samplesrepeat_",N,"_",xn,"_",yn,"_",stdize,"_",prior[1,1],"_",prior[1,2],"_",prior[2,1],"_",prior[2,2],"_",tau[1,2],"_m",m,".R",sep='')
			cat(paste("\nnABC.mutost: compute ",f.name))
			options(show.error.messages = FALSE, warn=1)		
			readAttempt<-try(suppressWarnings(load(f.name)))						
			options(show.error.messages = TRUE)						
			if(!resume || inherits(readAttempt, "try-error"))
			{										   
				ans.equal	<- project.nABC.mutost.fix.x.fix.ysig2(N,tau,prior,alpha,x,yn,stdize)				
				cat(paste("\nnABC.mutost: save ",f.name))
				save(ans.equal,file=f.name)												
			}
			else
				cat(paste("\nnABC.mutost: resumed ",f.name))
			
			yn<- 8*xn
			f.name<- paste(dir.name,"/nABC.mutost_samplesrepeat_",N,"_",xn,"_",yn,"_",stdize,"_",prior[1,1],"_",prior[1,2],"_",prior[2,1],"_",prior[2,2],"_",tau[1,2],"_m",m,".R",sep='')
			cat(paste("\nnABC.mutost: compute ",f.name))
			options(show.error.messages = FALSE, warn=1)		
			readAttempt<-try(suppressWarnings(load(f.name)))						
			options(show.error.messages = TRUE)						
			if(!resume || inherits(readAttempt, "try-error"))
			{										   
				ans.more	<- project.nABC.mutost.fix.x.fix.ysig2(N,tau,prior,alpha,x,yn,stdize)				
				cat(paste("\nnABC.mutost: save ",f.name))
				save(ans.more,file=f.name)												
			}
			else
				cat(paste("\nnABC.mutost: resumed ",f.name))
			
			h<- list()
			ans<- ans.equal
			acc<- which( ans[["data"]]["error",]<=ans[["cir"]]  &  ans[["data"]]["error",]>=ans[["cil"]] )
			print(length(acc)/ncol(ans[["data"]]))	
			h[[1]]<- project.nABC.movingavg.gethist(ans[["data"]]["ymu",acc], ans[["xmu"]], nbreaks= 50, width= 0.5, plot=0, rtn.dens=1)
			ans<- ans.more
			acc<- which( ans[["data"]]["error",]<=ans[["cir"]]  &  ans[["data"]]["error",]>=ans[["cil"]] )
			print(length(acc)/ncol(ans[["data"]]))			
			h[[2]]<- project.nABC.movingavg.gethist(ans[["data"]]["ymu",acc], ans[["xmu"]], nbreaks= 50, width= 0.5, plot=0, rtn.dens=1)
			
			ltys<- c(2,3)
			x<- seq(prior[1,1],prior[1,2],length.out=1e3)
			y<-	dnorm(x,xmu,sqrt(xsigma2/xn))
			y<- y / diff(pnorm(prior[1,],xmu,sqrt(xsigma2/xn)))
			y2<- nabc.mutost.pow(x, xn-1, tau, sqrt(xsigma2/xn), alpha)
			x2<- x[which(y2!=0)]
			y2<- y2[which(y2!=0)]
			y2<- y2/sum(diff(x2)*y2[-1])			
							  
			tmp	<- nabc.mutost.onesample.n.of.y(xn, sqrt(xsigma2/xn), 0.9, sqrt(xsigma2), alpha, tau.u.ub=2*tau[1,2] )
			yn	<- tmp[1]
			tau	<- tmp[3]						
			y3<- nabc.mutost.pow(x, yn-1, tau, sqrt(xsigma2/yn), alpha)
			x3<- x[which(y3!=0)]
			y3<- y3[which(y3!=0)]
			y3<- y3/sum(diff(x3)*y3[-1])
			
			y4<- nabc.mutost.pow(x, xn-1, tau, sqrt(xsigma2/xn), alpha)
			x4<- x[which(y4!=0)]
			y4<- y4[which(y4!=0)]
			y4<- y4/sum(diff(x4)*y4[-1])
			
			
			xlim<- range(sapply(h,function(x) range(x$dens$x)))
			ylim<- range(c(y,sapply(h,function(x) range(x$dens$y))))
			plot(1,1,type='n',xlim=xlim,ylim=ylim,xlab=expression(mu))
			lines(x,y,col="red")
			lines(x2+xmu,y2,col="blue")			
			lines(x3+xmu,y3,col="green")
			lines(x4+xmu,y4,col="pink")
			sapply(seq_along(h),function(i)
					{
						lines(h[[i]][["dens"]][["x"]], h[[i]][["dens"]][["y"]],lty=ltys[i])
					})
			
			abline(v=xmu,col="red")
			#print(hist)
			stop()
		}		
	}
	if(!is.na(subprog) && subprog==3)	#TOST power
	{
		cols	<- c("black","blue","DarkViolet","red")				
		xn		<- 7	#60
		alpha	<- 0.01
		s		<- 0.1163	#sqrt(2)
		
		s.of.T	<- s/sqrt(xn)
		df		<- xn-1				
		tau.us	<- c(0.1,0.16,0.18,0.2,0.225)#c(0.5) 
		ltys	<- seq_along(tau.us)
		rho		<- lapply(tau.us,function(tau.u)	seq(-tau.u,tau.u,length.out=1e3))
	#print(rho[[1]])	
		pw		<- lapply(seq_along(tau.us),function(i)
				{
					nabc.mutost.pow(rho[[i]], df, tau.us[i], s.of.T, alpha)								
				})
	#print(pw[[1]])	
		plot(1,1,type='n',xlim=range(unlist(rho)),ylim=range(unlist(pw)))
		sapply(seq_along(tau.us),function(i)
				{
					lines(rho[[i]],pw[[i]],lty=ltys[i])					
				})
		#stop()
		#main issue with TOST: tau.u cannot be made arbitrarily small
		tau.u<- seq(0.1,0.225,length.out=1e3)
		y<- tau.u/s.of.T+qt( alpha, df )
		y[y<=0]<- 0
		plot(tau.u,y,type='l')
		#legend("topleft", legend= bquote(.(parse(text=paste("n==~",n,sep=""))) ), bty= 'n', text.col=cols)
		#legend("topright",legend= expression(paste(tau,"=0.005")), bty= 'n')
		#abline(v=0, col="gray30",lty=2)
		#dev.off()
		stop()
	}
	
	if(0)	#TOST power for obs SEIRS model
	{
		alpha<- 0.01
		#load("~/workspace_sandbox/phylody/data/nABC.SEIIRScompare/power_schuir_obs_SEIRS.R")		
		dir.name	<- paste(DATA,"nABC.SEIIRScompare",sep='/')
		f.name		<- paste(dir.name,"power_schuir_obs_SEIRS_MATTlag2.1.R",sep='/')
		cat(paste("\nload",f.name))
		load(f.name)
		pw			<- POWER
		pw			<- list(pw[[1]],pw[[3]],pw[[2]],pw[[5]])		
		names.pw	<- expression(mu*"-attack",mu*"-fd-attack",mu*"-pop-attack",mu*"-t-pk")		
		tau.us		<- c(0.015,2.5,0.015,0.2)
		tau.us		<- c(0.015,4.5,0.017,0.2)		
		ltys		<- c(1,3,4,5)
		cols		<- c(my.fade.col("black",1),my.fade.col("black",0.8),my.fade.col("black",0.6),my.fade.col("black",0.4))
		f.name		<- paste(dir.name,"power_schuir_obs_SEIRS.pdf",sep='/')
		cat(paste("\nplot",f.name))
		pdf(f.name,version="1.4",width=3,height=6)
		par(mar=c(4,4,.5,.5))
		plot(1,1,type='n',bty='n',xlim=c(-1.1,1.1),ylim=c(0,1.2),xaxt='n',ylab="power",xlab=expression(rho))
		axis(1, at=c(-1,-0.5,0,0.5,1),labels=expression(tau^'-',tau^'-'/2,0,tau^'+'/2,tau^'+'))
		sapply(seq_along(pw),function(k)
				{				
					tau.u	<- tau.us[k]
					tau.l	<- -tau.u				
					tau		<- seq(tau.l, tau.u, length.out=5e2)
					se		<- pw[[k]][3]
					nm		<- pw[[k]][4]
					df		<- pw[[k]][5]
					y		<- PowerTOST:::.power.TOST(alpha=alpha, tau.l, tau.u, tau, se, nm, df, bk = 4)
					lines(seq(-1,1,length.out=5e2),y, lty=ltys[k], col=cols[k], lwd=1.5)								
				})
		legend("topright",bty='n',legend=names.pw,fill=cols,lty=ltys,border=NA)
		dev.off()
		stop()
		pw<- POWER
		pw<- list(pw[[1]],pw[[2]],pw[[3]],pw[[5]])		
		names(pw)<- c("MATT","XATT","AMEDFDA","SDTPK")
		print(pw)

		sapply(seq_along(pw),function(k)
				{				
					x<- pw[[k]]
					tau.l<- x[1]
					tau.u<- x[2]
					if(k==2)	tau.l<- -tau.u
					tau<- seq(tau.l, tau.u, length.out= 1e3)
					se<- x[3]
					nm<- x[4]
					df<- x[5]
					y.cur<- PowerTOST:::.power.TOST(alpha=alpha, tau.l, tau.u, tau, se, nm, df, bk = 4)
					y1a<- PowerTOST:::.power.TOST(alpha=alpha, tau.l/1.5, tau.u/1.5, tau, se, nm, df, bk = 4)
					y2<- PowerTOST:::.power.TOST(alpha=alpha, tau.l/2.5, tau.u/2.5, tau, se, nm, df, bk = 4)
					y4<- PowerTOST:::.power.TOST(alpha=alpha, tau.l/6, tau.u/6, tau, se, nm, df, bk = 4)
					y6<- PowerTOST:::.power.TOST(alpha=alpha, tau.l/8, tau.u/8, tau, se, nm, df, bk = 4)
					y8<- PowerTOST:::.power.TOST(alpha=alpha, tau.l/16, tau.u/16, tau, se, nm, df, bk = 4)
					y24<- PowerTOST:::.power.TOST(alpha=alpha, tau.l/36, tau.u/36, tau, se, nm, df, bk = 4)
					plot(tau,y.cur,type='l',xlab="rho",ylab="power",ylim=c(0,1.2))
					lines(tau,y1a,lty=2)
					lines(tau,y2,lty=3)
					lines(tau,y4,lty=4)
					lines(tau,y6,lty=5)
					lines(tau,y8,lty=6)
					lines(tau,y24,lty=7)
					legend("topright",bty='n',names(pw)[k])
					legend("topleft",lty=1:7,legend=c("/1","/1.5","/2.5","/6","/8","/16","/36"))
					if(k==3)	stop()
				})
	}
	if(!is.na(subprog) && subprog==4)	#TOST power for simu SEIRS model
	{
		alpha		<- 0.01
		dir.name	<- paste(DATA,"nABC.SEIIRScompare",sep='/')
		f.name		<- paste(dir.name,"power_schuir_simu_SEIRS_MATTlag2.1.R",sep='/')
		cat(paste("\nload",f.name))
		load(f.name)
		pw			<- POWER
		pw			<- list(pw[[1]],pw[[3]],pw[[2]])	
		print(pw)
		names.pw	<- expression(mu*"-attack",mu*"-fd-attack",mu*"-pop-attack")		
		tau.us		<- c(0.0005,0.5,0.02)
		tau.us		<- c(0.000275,0.225,0.01)
		ltys		<- c(1,3,4)
		cols		<- c(my.fade.col("black",0.8),my.fade.col("black",0.6),my.fade.col("black",0.4))
		f.name		<- paste(dir.name,"power_schuir_simu_SEIRS.pdf",sep='/')
		cat(paste("\nplot",f.name))
		pdf(f.name,version="1.4",width=3,height=4)
		par(mar=c(4,4,.5,.5))
		plot(1,1,type='n',bty='n',xlim=c(-1.1,1.1),ylim=c(0,1.2),xaxt='n',ylab="power",xlab=expression(rho))
		axis(1, at=c(-1,-0.5,0,0.5,1),labels=expression(tau^'-',tau^'-'/2,0,tau^'+'/2,tau^'+'))
		sapply(seq_along(pw),function(k)
			{				
				tau.u	<- tau.us[k]
				tau.l	<- -tau.u				
				rho		<- seq(tau.l, tau.u, length.out=5e2)
				s.of.T	<- pw[[k]][6]
				df		<- pw[[k]][5]
				y		<- nabc.mutost.pow(rho, df, tau.u, s.of.T, alpha)
				lines(seq(-1,1,length.out=5e2),y, lty=ltys[k], col=cols[k], lwd=1.5)								
			})
		legend("topright",bty='n',legend=names.pw,fill=cols,lty=ltys,border=NA)
		dev.off()
		
		#plot also tau.u for any of these values
		f.name		<- paste(dir.name,"power_schuir_simu_SEIRS_cus.pdf",sep='/')
		cat(paste("\nplot",f.name))
		pdf(f.name,version="1.4",width=3,height=4)
		par(mar=c(4,4,.5,.5))		
		plot(1,1,type='n',xlim=c(0,1.5),ylim=c(0,2.5),xaxt='n')
		axis(1, at=c(0,0.5,1,1.5),labels=expression(0,0.5*tau^'+',tau^'+',1.5*tau^'+'))
		sapply(seq_along(pw),function(k)
			{
				tau.u	<- seq(0,1.5*tau.us[k],length.out=5e2)
				df		<- pw[[k]][5]
				s.of.T	<- pw[[k]][6]
				c.u		<- tau.u/s.of.T+qt( alpha, df )
				c.u[c.u<=0]<- 0
				lines(seq(0,1.5,length.out=5e2), c.u, lty=ltys[k], col=cols[k], lwd=1.5)
			})
		legend("topright",bty='n',legend=names.pw,fill=cols,lty=ltys,border=NA)
		dev.off()
		stop()
	}
	if(0)	#TOST power as n increases for different ABC thresholds
	{
		taur<- c(0.05, 0.045, 0.042)
		tau0<- 0.04
		se<- 0.01
		n<-  seq(10, 500, 1)
		cols<- c("black","blue","red")
		pow<- sapply(seq_along(taur),function(i){		PowerTOST:::.power.TOST(alpha = 0.05, -taur[i], taur[i], tau0, se, 2*n, 2*n-2, bk = 4)			})		#assume variance the same

		pdf(paste(dir.name,"/loc_tost_power.pdf",sep=''),version="1.4",width=4,height=4)
		par(mar=c(4,4,1,1))
		plot(1,1,type='n', xlim= range(n), ylim= c(0,1),xlab="n",ylab=expression(paste("power  P( R | ",H[1]," )")), bty='n')
		sapply(1:ncol(pow),function(j)
				{
					lines( n, pow[,j], col=cols[j] )
				})
		legend(x=5,y=0.3,expression(paste(tau,"= 0.042")),bty='n', text.col=cols[3])
		legend(x=-30,y=0.5,expression(paste(tau,"= 0.045")),bty='n', text.col=cols[2])
		legend(x=-50,y=0.7,expression(paste(tau,"= 0.05")),bty='n', text.col=cols[1])
		legend("bottomright", legend= expression(paste(nu[k](x)-nu[k](theta),"=0.04")), bty= 'n')
		abline(h=0.8, col="black",lty=2)
		dev.off()
	}
	if(0)
	{
		n<- c(100)
		se<- 0.01
		cols<- c("blue")
		#tight taur appropriate for large n but not small n
		taur<- 0.005
		tau0<- seq(-taur,taur,by=taur/50)
		pow<- sapply(seq_along(n),function(i){		PowerTOST:::.power.TOST(alpha = 0.05, -taur, taur, tau0, se, 2*n[i], 2*n[i]-2, bk = 4)			})		#assume variance the same
		pdf(paste(dir.name,"/loc_tost_change_n100.pdf",sep=''),version="1.4",width=6,height=6)
		plot(1,1,type='n', xlim= range(tau0), ylim= c(0,1),xlab=expression(nu(theta)-nu(x)),ylab="power")
		sapply(1:ncol(pow),function(j)
				{
					lines( tau0, pow[,j], lty= 1, col=cols[j] )
				})

		legend("topleft", legend= bquote(.(parse(text=paste("n==~",n,sep=""))) ), bty= 'n', text.col=cols)
		legend("topright",legend= expression(paste(tau,"=0.005")), bty= 'n')
		abline(v=0, col="gray30",lty=2)
		dev.off()
	}
	if(0)	#TOST power as n changes
	{
		n<- c(10, 30, 60, 200)
		se<- 0.01
		cols<- c("black","blue","DarkViolet","red")
		#tight taur appropriate for large n but not small n
		taur<- 0.005
		tau0<- seq(-taur,taur,by=taur/50)
		pow<- sapply(seq_along(n),function(i){		PowerTOST:::.power.TOST(alpha = 0.05, -taur, taur, tau0, se, 2*n[i], 2*n[i]-2, bk = 4)			})		#assume variance the same

		pdf(paste(dir.name,"/loc_tost_change_n.pdf",sep=''),version="1.4",width=6,height=6)
		plot(1,1,type='n', xlim= range(tau0), ylim= c(0,1),xlab=expression(nu(theta)-nu(x)),ylab="power")
		sapply(1:ncol(pow),function(j)
				{
					lines( tau0, pow[,j], lty= 1, col=cols[j] )
				})

		legend("topleft", legend= bquote(.(parse(text=paste("n==~",n,sep=""))) ), bty= 'n', text.col=cols)
		legend("topright",legend= expression(paste(tau,"=0.005")), bty= 'n')
		abline(v=0, col="gray30",lty=2)
		dev.off()
	}
	if(0)	#TOST power as tau changes
	{
		n<- 30
		se<- 0.01
		taur<- c(0.001, 0.005, 0.01, 0.015)
		pow<- sapply(seq_along(taur),function(i)
					{
						tau0<- seq(-taur[i],taur[i],by=taur[i]/50)
						PowerTOST:::.power.TOST(alpha = 0.05, -taur[i], taur[i], tau0, se, 2*n, 2*n-2, bk = 4)
					})		#assume variance the same
		pdf(paste(dir.name,"/loc_tost_change_tau.pdf",sep=''),version="1.4",width=6,height=6)
		plot(1,1,type='n', xlim= range(c(-max(taur),max(taur))), ylim= c(0,1),xlab=expression(nu(theta)-nu(x)),ylab="power")
		sapply(1:ncol(pow),function(i)
				{
					tau0<- seq(-taur[i],taur[i],by=taur[i]/50)
					lines( tau0, pow[,i], lty= i )
				})
		legend("topleft",legend= "n=30", bty= 'n')
		legend("topright", legend= taur, lty= seq_along(taur), bty= 'n')
		abline(v=0, col="blue")
		dev.off()
		#abline(h=0.05, col="red")
	}
	if(1)	#TOST power as sd changes
	{
		n<-  30
		taur<- 0.01
		#power for rho=0 and increasing variance
		tau0<- 0
		se<- seq(1e-3, 0.05, by= 1e-4)					
		pow<- sapply(seq_along(taur),function(i){		PowerTOST:::.power.TOST(alpha = 0.05, -taur, taur, tau0, se, 2*n, 2*n-2, bk = 4)			})		#assume variance the same
		
		#power for rho=taur/2 and increasing variance
		tau0<- taur/3
		se<- seq(1e-3, 0.05, by= 1e-4)
		pow2<- sapply(seq_along(taur),function(i){		PowerTOST:::.power.TOST(alpha = 0.05, -taur, taur, tau0, se, 2*n, 2*n-2, bk = 4)			})		#assume variance the same
		
		plot(1,1, type='n', xlim= range(se), ylim=c(0,1))
		lines(se, pow )
		lines(se, pow2, col="blue")
		#plot(1,1,type='n', xlim= range(n), ylim= c(0,1),xlab="n",ylab=expression(paste("power  P( R | ",H[1]," )")), bty='n')
		
	}
	if(0)#one sided t - test 	H0: mu1-mu2<=tau
	{
		n<- c(10, 30, 200)
		alpha<- 0.05
		se<- 0.01
		taur<- 0.001
		tau0<- seq(taur, 0.015, by= 0.015/50)
		pow<- sapply(seq_along(n), function(i){ 1 - pt( qt(1-alpha, 2*n[i]-2) - tau0/se*sqrt(n[i]), 2*n[i]-2)	})

		pdf(paste(dir.name,"/loc_ost_change_n.pdf",sep=''),version="1.4",width=6,height=6)
		plot(1,1,type='n', xlim= range(c(-0.005,tau0)), ylim= c(0,1),xlab=expression(nu(theta)-nu(x)),ylab="power")
		sapply(1:ncol(pow),function(j)
				{
					lines( tau0, pow[,j], lty= j )
				})
		legend("topleft", legend= n, lty= seq_along(n), bty= 'n')
		legend("bottomright",legend= expression(paste(tau,"=0.001")), bty= 'n')
		abline(v=0, col="blue")
		dev.off()
	}
	stop()
}
#------------------------------------------------------------------------------------------------------------------------
project.nABC.ISBA2012talk<- function()
{
	require(sfsmisc)
	dir.name<- "/Users/olli0601/duke/2012_frequencyABC/sim.data"
	n<- 5000; sigmax<- 3; sigmay<-1.5; mux<- 5; muy<- 0
	x<- rnorm(n,mux,sigmax)
	y<- rnorm(n,muy,sigmay)
	breaks<- range(c(x,y))
	breaks[1]<- breaks[1]*ifelse(breaks[1]<0, 3, 0.8)
	breaks[2]<- breaks[2]*ifelse(breaks[2]>0, 3, 0.8)
	breaks<- seq(breaks[1],to=breaks[2], by= diff(breaks)/150)

	if(0)
	{
		pdf(paste(dir.name,"/talk_nABC1a.pdf",sep=''),version="1.4",width=4,height=1.5)
		par(mar=c(2,0,0,0))
		plot(1,1,xlab='', xlim= range(c(x,y)),type='n',ylim=c(0,0.15),ylab='',yaxt='n',bty='n',main="")
		lines(c(mean(x),mean(x)),c(0,0.1),col="black",lwd=2)
		lines(c(mean(y),mean(y)),c(0,0.1),col="blue",lwd=2)
		points(mean(x),0.1,pch=19)
		points(mean(y),0.1,pch=19,col=my.fade.col("#0080FFFF",1))
		p.arrows((mean(x)+mean(y))/2,0.05,mean(y),0.05,fill="red",col="red")
		p.arrows((mean(x)+mean(y))/2,0.05,mean(x),0.05,fill="red",col="red")
		legend("topleft",expression(paste("sim | ",theta)),fill=my.fade.col("#0080FFFF",0.6),bty='n')
		legend("topright","obs",fill="gray60",bty='n')
		dev.off()
	}
	if(0)	#compare summaries with nABC
	{
		pdf(paste(dir.name,"/talk_nABC1.pdf",sep=''),version="1.4",width=4,height=2.5)
		hx<- hist(x,	breaks= breaks,	plot= F)
		hy<- hist(y,	breaks= breaks,	plot= F)
		ylim<- c(0,max(c(hx$intensities,hy$intensities),na.rm=TRUE))
		xlim<- range(c(x,y))
		par(mar=c(0,0,0,0))
		plot(1,1,xlab='', xlim= xlim,type='n',ylim=ylim,ylab='',main="", xaxt='n', yaxt='n', bty='n')
		plot(hx,freq=F,col="gray60",border=NA,add=TRUE)
		plot(hy,freq=F,col=my.fade.col("#0080FFFF",0.6),border=NA,add=TRUE)
		lines(c(mean(x),mean(x)),c(0,0.4),col="black",lwd=2)
		lines(c(mean(y),mean(y)),c(0,0.4),col="blue",lwd=2)
		p.arrows((mean(x)+mean(y))/2,0.1,mean(y),0.1,fill="red",col="red")
		p.arrows((mean(x)+mean(y))/2,0.1,mean(x),0.1,fill="red",col="red")
		legend("topleft",expression(paste("sim | ",theta)),fill=my.fade.col("#0080FFFF",0.6),bty='n')
		legend("topright","obs",fill="gray60",bty='n')
		dev.off()
	}
	if(0)	#simple auxiliary model of summaries
	{
		pdf(paste(dir.name,"/talk_nABC2.pdf",sep=''),version="1.4",width=4,height=4)
		hx<- hist(x,	breaks= breaks,	plot= F)
		hy<- hist(y,	breaks= breaks,	plot= F)
		ylim<- c(0,max(c(hx$intensities,hy$intensities),na.rm=TRUE))
		xlim<- range(c(x,y))
		par(mar=c(1,0,0,0))
		plot(1,1,xlab='', xlim= xlim,type='n',ylim=ylim,ylab='',main="", xaxt='n', yaxt='n', bty='n')
		plot(hx,freq=F,col="gray60",border=NA,add=TRUE)
		plot(hy,freq=F,col=my.fade.col("#0080FFFF",0.6),border=NA,add=TRUE)
		lines(c(mean(x),mean(x)),c(0,0.4),col="black",lwd=1)
		lines(c(mean(y),mean(y)),c(0,0.4),col="blue",lwd=1)
		legend("topleft",expression(paste("sim | ",theta)),fill=my.fade.col("#0080FFFF",0.6),bty='n')
		legend("topright","obs",fill="gray60",bty='n')
		mtext(expression(nu[k](theta)),side=1,line=-0.5,at=mean(y),col="blue")
		mtext(expression(nu[k](x)),side=1,line=-0.5,at=mean(x))
		xlim<- range(x)
		lines(seq( xlim[1], xlim[2], by= diff(xlim)/100 ), dnorm( seq( xlim[1], xlim[2], by= diff(xlim)/100 ), mean(x), sd(x) ),lwd=2)
		xlim<- range(y)
		lines(seq( xlim[1], xlim[2], by= diff(xlim)/100 ), dnorm( seq( xlim[1], xlim[2], by= diff(xlim)/100 ), mean(y), sd(y) ), col="blue",lwd=2)
		dev.off()
	}
	if(0)
	{
		pdf(paste(dir.name,"/talk_nABC3.pdf",sep=''),version="1.4",width=4,height=4)
		par(mar=c(4.5,4,1,1))
		t<- seq(0,5,0.01)
		plot(t, -(t-1.5)^2 , type='l', xlab= expression(paste("original parameter ",theta)), ylab= expression(paste("auxiliary parameter  ",nu[k])),las=0,bty='n', col="blue")
		lines(c(0,5),c(-6,-6))
		legend(x=3.8,y=-6,legend=expression(nu[k](x)),bty='n')
		legend(x=3.2,y=-3,legend=expression(nu[k](theta)),bty='n',text.col="blue")
		dev.off()
	}
	if(0)
	{
		n<- 10
		t<- seq(-5,5,0.01)
		alpha<- 0.05
		pdf(paste(dir.name,"/talk_nABC4.pdf",sep=''),version="1.4",width=4,height=4)
		par(mar=c(4.5,0,0,0))
		plot(t,dt(t, df=2*n-2),type='l', xlab="Student T-test",yaxt='n',bty='n')
		polygon(  c(seq(t[1],qt(alpha/2, df=2*n-2),0.01),qt(alpha/2, df=2*n-2),t[1]), c(dt(seq(t[1],qt(alpha/2, df=2*n-2),0.01), df=2*n-2),0,0), border=NA, col=my.fade.col("red",0.6)  )
		polygon(  c(seq(qt(1-alpha/2, df=2*n-2),t[length(t)],0.01),t[length(t)],qt(1-alpha/2, df=2*n-2)), c(dt(seq(qt(1-alpha/2, df=2*n-2),t[length(t)],0.01), df=2*n-2),0,0), border=NA, col=my.fade.col("red",0.6)  )
		legend("topleft",legend=expression( paste("P( R | ",H[0]," )") ),fill="red",bty='n')
		legend("topright",legend=expression( paste(H[0]," : ",nu[k](theta)-nu[k](x)," = 0") ),bty='n')
		dev.off()
	}
	if(0)
	{
		n<- 10
		t<- seq(0,5,0.01)
		alpha<- 0.05
		tau<-  3.5
		pdf(paste(dir.name,"/talk_nABC5.pdf",sep=''),version="1.4",width=4,height=4)
		par(mar=c(4.5,0,0,0))
		plot(-t, dt(-t+tau, df=2*n-2),type='l', xlab="Schuirmann T-test",xlim=c(-t[length(t)],t[length(t)]),ylim=c(0,1.1*max(dt(-t, df=2*n-2))),yaxt='n',ylab='',bty='n')
		lines(c(-tau,-tau),c(0,dt(0, df=2*n-2)), lty=2)
		polygon(  c(seq(qt(1-alpha, df=2*n-2)-tau,0,0.01),0,qt(1-alpha, df=2*n-2)-tau),
						c(dt(seq(qt(1-alpha, df=2*n-2),tau,0.01), df=2*n-2),0,0), border=NA, col=my.fade.col("red",0.6)  )
		polygon(  c(seq(0,qt(alpha, df=2*n-2)+tau,0.01),qt(alpha, df=2*n-2)+tau,0),
						c(dt(seq(-tau,qt(alpha, df=2*n-2),0.01), df=2*n-2),0,0), border=NA, col=my.fade.col("red",0.6)  )
		lines(t, dt(t-tau, df=2*n-2) )
		lines(c(tau,tau),c(0,dt(0, df=2*n-2)), lty=2)
		legend("topleft",legend=expression( paste("P( R | ",H[0]," )") ),fill="red",bty='n')
		legend("topright",legend=expression( paste(H[0]," : |",nu[k](theta)-nu[k](x),"| > ",tau) ),bty='n')
		legend(y=0.03,x=-tau-0.7,legend=expression(-tau),bty='n')
		legend(y=0.03,x=tau-0.7,legend=expression(tau),bty='n')
		dev.off()
	}
	if(0)
	{
		pdf(paste(dir.name,"/talk_nABC6.pdf",sep=''),version="1.4",width=6,height=4)
		par(mar=c(2,2,2,2))
		plot(1,1,type='n',xlim=c(0,4),ylim=c(0,4),xlab='',ylab='',xaxt='n',yaxt='n',bty='n')
		c11<- 1.25; c12<- 2.75
		polygon(c(-1,c11,c11,-1),c(-1,-1,5,5),bty='n',col=my.fade.col("red",0.4),border=NA)
		polygon(c(5,c12,c12,5),c(-1,-1,5,5),bty='n',col=my.fade.col("red",0.4),border=NA)
		c21<- 1.5; c22<- 2.5
		polygon(c(-1,-1,5,5),c(c21,-1,-1,c21),bty='n',col=my.fade.col("#0080FFFF",0.4),border=NA)
		polygon(c(-1,-1,5,5),c(c22,5,5,c22),bty='n',col=my.fade.col("#0080FFFF",0.4),border=NA)
		#abline(v=2,lty=3)
		#abline(h=2,lty=3)
		points(2,2,pch=19,cex=0.5)
		legend(x=2-0.2,y=2+0.2,legend="c(0,0)",bty='n')
		legend(x=c11-0.2,y=c22,legend=expression(paste("reject ",T[1]," and ",T[2])),bty='n')
		legend(x=-0.2,y=c22,legend=expression(paste("accept ",T[1])),text.col="red",bty='n')
		legend(x=c11-0.2,y=3.75,legend=expression(paste("accept ",T[2])),text.col="blue",bty='n')
		mtext(expression(c[2]^L),at=c21,side=2,line=0.5,col="blue")
		mtext(expression(c[2]^U),at=c22,side=2,line=0.5,col="blue")
		mtext(expression(c[1]^L),at=c11,side=1,line=0.5,col="red")
		mtext(expression(c[1]^U),at=c12,side=1,line=0.5,col="red")
		dev.off()
	}
	if(0)
	{
		pdf(paste(dir.name,"/talk_nABC7.pdf",sep=''),version="1.4",width=6,height=4)
		par(mar=c(4,4,1,0))
		acc<- seq(0.01,0.4,0.01)
		alpha<- c(0.05,0.01)
		tp<- sapply(alpha, function(x){ tmp<- 1-x/acc; tmp[tmp<0]<-0; tmp })
		plot(1,1,type='n',bty='n',xlim=range(acc),ylim=range(tp),xlab="ABC acceptance %",ylab=expression(paste("P(",H[1],"| R ) is at least ")) )
		lines(acc,tp[,1],lty= 1)
		lines(acc,tp[,2],lty= 2,col="blue")
		legend(x=0.15,y=0.75,expression(paste(alpha,"= 0.05")),bty='n')
		legend(x=0.03,y=0.85,expression(paste(alpha,"= 0.01")),bty='n',text.col="blue")
		dev.off()
		stop()
	}
	if(1)
	{
		pdf(paste(dir.name,"/talk_nABC8.pdf",sep=''),version="1.4",width=4,height=4)
		n<- 10
		t<- seq(-5,5,0.005)
		s<- dt(t, df=2*n-2)
		tau<- 1.75
		alpha<- dt(-tau, df=2*n-2)
		par(mar=c(1,0.5,0.5,0.5))
		plot(1,1,type='n',xlim=range(t),ylim=range(s),xaxt='n',yaxt='n',bty='n')
		polygon(  c(t[1],-tau,-tau,t[1]), c(0,0,alpha,alpha), border=NA, col=my.fade.col("red",0.4)  )
		polygon(  c(tau,t[length(t)],t[length(t)],tau), c(0,0,alpha,alpha), border=NA, col=my.fade.col("red",0.4)  )
		polygon(  c(-tau,tau,tau,-tau), c(1,1,alpha,alpha), border=NA, col=my.fade.col("blue",0.4)  )
		lines(c(0,0),c(0,dt(0, df=2*n-2)),lty=4,col="gray40")
		lines(c(-5,-tau),c(0,0),col="red")
		lines(c(-tau,tau),c(0,0),col="blue")
		lines(c(tau,5),c(0,0),col="red")
		legend(x=-5.5,y=0.07,legend=expression(paste("p( R | ",H[0]," )")),text.col="red",bty='n')
		legend(x=-4,y=0.3,legend=expression(paste("p( R | ",H[1]," )")),text.col="blue",bty='n')
		mtext(expression(alpha),side=2,line=-0.5,at=alpha,col="red")
		mtext(c(expression(-tau),0,expression(tau)),side=1,line=-0.5,at=c(-tau,0,tau))
		mtext(c(expression(H[0]),expression(H[0])),side=1,line=-0.5,at=c(-4,3),col="red")
		mtext(expression(H[1]),side=1,line=-0.5,at=-1,col="blue")
		lines(t,s,lty=2)
		dev.off()
	}
	if(0)
	{
		pdf(paste(dir.name,"/talk_nABC3.pdf",sep=''),version="1.4",width=6,height=6)
		y2<- rnorm(n,mux+8,sigmay*1.2)
		hy2<- hist(y2,	breaks= breaks,	plot= F)
		y3<- rnorm(n,mux+2,sigmax*0.8)
		hy3<- hist(y3,	breaks= breaks,	plot= F)
		ylim<- c(0,max(c(hx$intensities,hy$intensities,hy2$intensities,hy3$intensities),na.rm=TRUE))
		xlim<- range(c(x,y))
		par(mar=c(0,0,0,0))
		plot(1,1,xlab='', xlim= xlim,type='n',ylim=ylim,ylab='',main="", xaxt='n', yaxt='n', bty='n')

		plot(hy,freq=F,col=my.fade.col("#0080FFFF",0.4),border=NA,add=TRUE)
		plot(hy2,freq=F,col=my.fade.col("red",0.4),border=NA,add=TRUE)
		plot(hy3,freq=F,col=my.fade.col("green",0.4),border=NA,add=TRUE)
		plot(hx,freq=F,col="gray60",border=NA,add=TRUE)
		lines(c(mean(x),mean(x)),c(0,0.4),col="black",lwd=2)
		lines(c(mean(y),mean(y)),c(0,0.4),col="blue",lwd=2)
		lines(c(mean(y2),mean(y2)),c(0,0.4),col="red",lwd=2)
		lines(c(mean(y3),mean(y3)),c(0,0.4),col="green",lwd=2)
		legend("topleft",c(expression(paste("sim | ",theta,"1")),expression(paste("sim | ",theta,"2")),expression(paste("sim | ",theta,"3"))),fill=c(my.fade.col("#0080FFFF",0.6),my.fade.col("red",0.6),my.fade.col("green",0.6)),bty='n')
		legend("topright","obs",fill="gray60",bty='n')
		dev.off()
	}
	stop()
}
#------------------------------------------------------------------------------------------------------------------------
project.nABC.binom<- function()
{
	#simulated data: compute p-values for one sided test
	ny<- 30
	alpha<- 0.05
	alternative<- "greater"
	dir.name<- "/Users/olli0601/duke/2012_frequencyABC/sim.data"

	project.nABC.binom.fbinom<- function(px, ny, py, alternative)
	{
		y<- rbinom(ny,1,py)
#print(y)
		t.h0<- binom.test(		length(which(as.logical(y))),	 length(y), p=px, alternative= alternative)
		t.h1<- binom.test(		length(which(as.logical(y))),	 length(y), p=px, alternative= ifelse(alternative=="greater","less","greater"))	#flipped H0
		ans<- c(	length(which(as.logical(y))),  py, ny,  t.h0$p.value, t.h1$p.value, px	)
		names(ans)<- c("1s","py","ny","H0.pval","H1.pval","px")
		ans
	}
	project.nABC.binom.tostbinom<- function(px, tau, ny, py)
	{
		y<- rbinom(ny,1,py)
		tl<- binom.test(		length(which(as.logical(y))),	 length(y), p=px-tau, alternative= "greater")		#reject if py>=px-tau
		tr<- binom.test(		length(which(as.logical(y))),	 length(y), p=px+tau, alternative= "less")			#reject if py<=px+tau
		te<- binom.test(		length(which(as.logical(y))),	 length(y), p=px, alternative= "two.sided")			#reject if py neq px
		ans<- c(	length(which(as.logical(y))),  py, ny,  max(tl$p.value, tr$p.value), te$p.value, tl$p.value, tr$p.value, px	)
		names(ans)<- c("1s","py","ny","tost.pval","te.pval","tl.pval","tr.pval","px")
		ans
	}

	if(0)
	{
		px<- 0.4
		py<- 0.4
		tau<- 0.2
		#ny<- 1e3
		#px and py equal
		pdf(paste(dir.name,"/binom_tost_px_py_equal.pdf",sep=''),version="1.4",width=6,height=6)
		tmp<- replicate(  1000, project.nABC.binom.tostbinom(px, tau, ny, py) )
		def.par <- par(no.readonly = TRUE)
		layout.m<- matrix(data= c(1,1,2,3,3,4),ncol=3,nrow=2,byrow=1)
		layout(layout.m)
		plot(	1:ncol(tmp), tmp["tost.pval",], type='n',  xlab= "replicates", ylab="TOST H0.pval"	)
		polygon(c(1,ncol(tmp),ncol(tmp),1),c(alpha,alpha,0,0), col= "grey80", border=NA)
		points(	1:ncol(tmp), tmp["tost.pval",], type='s', col="blue")
		qq<- qqplot((1:ncol(tmp)-0.5) / ncol(tmp), tmp["tost.pval",], plot.it= 0)
		plot(qq,type='p',col="blue", pch=19, xlab="expected quantiles under H0",ylab="quantiles of TOST pval", xlim=c(0,1),ylim=c(0,1))
		abline(a=0,b=1, lty= 2)

		plot(	1:ncol(tmp), tmp["te.pval",], type='n', xlab= "replicates", ylab="equal H0.pval"	)
		#polygon(c(1,ncol(tmp),ncol(tmp),1),c(alpha,alpha,0,0), col= "grey80", border=NA)
		points(	1:ncol(tmp), tmp["te.pval",], type='s', col="red")
		qq<- qqplot((1:ncol(tmp)-0.5) / ncol(tmp), tmp["te.pval",], plot.it= 0)
		plot(qq,type='p',col="red",pch=19, xlab="expected quantiles under H1",ylab="quantiles of equal H0 pval", xlim=c(0,1),ylim=c(0,1))
		abline(a=0,b=1, lty= 2)
		par(def.par)
		dev.off()
	}
	if(0)
	{
		px<- 0.4
		py<- 0.2
		tau<- 0.2
		#py on edge
		pdf(paste(dir.name,"/binom_tost_py_onedge.pdf",sep=''),version="1.4",width=6,height=6)
		tmp<- replicate(  1000, project.nABC.binom.tostbinom(px, tau, ny, py) )
		def.par <- par(no.readonly = TRUE)
		layout.m<- matrix(data= c(1,1,2,3,3,4),ncol=3,nrow=2,byrow=1)
		layout(layout.m)
		plot(	1:ncol(tmp), tmp["tost.pval",], type='n',  xlab= "replicates", ylab="TOST H0.pval"	)
		polygon(c(1,ncol(tmp),ncol(tmp),1),c(alpha,alpha,0,0), col= "grey80", border=NA)
		points(	1:ncol(tmp), tmp["tost.pval",], type='s', col="blue")
		qq<- qqplot((1:ncol(tmp)-0.5) / ncol(tmp), tmp["tost.pval",], plot.it= 0)
		plot(qq,type='p',col="blue", pch=19, xlab="expected quantiles under H0",ylab="quantiles of TOST pval", xlim=c(0,1),ylim=c(0,1))
		abline(a=0,b=1, lty= 2)

		plot(	1:ncol(tmp), tmp["te.pval",], type='n', xlab= "replicates", ylab="equal H0.pval"	)
		#polygon(c(1,ncol(tmp),ncol(tmp),1),c(alpha,alpha,0,0), col= "grey80", border=NA)
		points(	1:ncol(tmp), tmp["te.pval",], type='s', col="red")
		qq<- qqplot((1:ncol(tmp)-0.5) / ncol(tmp), tmp["te.pval",], plot.it= 0)
		plot(qq,type='p',col="red",pch=19, xlab="expected quantiles under H1",ylab="quantiles of equal H0 pval", xlim=c(0,1),ylim=c(0,1))
		abline(a=0,b=1, lty= 2)
		par(def.par)
		dev.off()
	}
	if(1)
	{
		px<- 0.4
		py<- 0.25
		tau<- 0.2
		ny<- seq(30,2e3,1)
		#py on edge, and n goes large
		pdf(paste(dir.name,"/binom_tost_py_onedge_nlarge.pdf",sep=''),version="1.4",width=6,height=6)
		tmp<- sapply( seq_along(ny),function(i){			project.nABC.binom.tostbinom(px, tau, ny[i], py)				})
		plot(	ny, tmp["tost.pval",], type='n',  xlab= "ny", ylab="TOST H0.pval"	)
		polygon(c(ny[1],ny[length(ny)],ny[length(ny)],ny[1]),c(alpha,alpha,0,0), col= "grey80", border=NA)
		points(	ny, tmp["tost.pval",], type='s', col="blue")
		dev.off()
	}
	if(0)
	{
		px<- 0.4
		py<- 0.4
		#px and py equal
		pdf(paste(dir.name,"/binom_px_py_equal.pdf",sep=''),version="1.4",width=6,height=6)
		tmp<- replicate(  1000, project.nABC.binom.fbinom(px,ny,py,alternative) )
		def.par <- par(no.readonly = TRUE)
		layout.m<- matrix(data= c(1,1,2,3,3,4),ncol=3,nrow=2,byrow=1)
		layout(layout.m)
		plot(	1:ncol(tmp), tmp["H0.pval",], type='n',  xlab= "replicates", ylab="H0.pval"	)
		polygon(c(1,ncol(tmp),ncol(tmp),1),c(1-alpha,1-alpha,1,1), col= "grey80", border=NA)
		points(	1:ncol(tmp), tmp["H0.pval",], type='s', col="blue")
		qq<- qqplot((1:ncol(tmp)-0.5) / ncol(tmp), tmp["H0.pval",], plot.it= 0)
		plot(qq,type='p',col="blue", pch=19, xlab="expected quantiles under H0",ylab="quantiles of H0 pval", xlim=c(0,1),ylim=c(0,1))
		abline(a=0,b=1, lty= 2)
		plot(	1:ncol(tmp), tmp["H1.pval",], type='n', xlab= "replicates", ylab="flipped H0.pval"	)
		polygon(c(1,ncol(tmp),ncol(tmp),1),c(alpha,alpha,0,0), col= "grey80", border=NA)
		points(	1:ncol(tmp), tmp["H1.pval",], type='s', col="red")
		qq<- qqplot((1:ncol(tmp)-0.5) / ncol(tmp), tmp["H1.pval",], plot.it= 0)
		plot(qq,type='p',col="red",pch=19, xlab="expected quantiles under H1",ylab="quantiles of flipped H0 pval", xlim=c(0,1),ylim=c(0,1))
		abline(a=0,b=1, lty= 2)
		par(def.par)
		dev.off()
	}
	if(0)
	{
		px<- 0.4
		py<- 0.1
		#py smaller than px, this is still H0.		the pval is not uniform when py is fixed to a specific value < px
		pdf(paste(dir.name,"/binom_py_smaller_px.pdf",sep=''),version="1.4",width=6,height=6)
		tmp<- replicate(  1000, project.nABC.binom.fbinom(px,ny,py,alternative) )
		def.par <- par(no.readonly = TRUE)
		layout.m<- matrix(data= c(1,1,2,3,3,4),ncol=3,nrow=2,byrow=1)
		layout(layout.m)
		plot(	1:ncol(tmp), tmp["H0.pval",], type='n',  xlab= "replicates", ylab="H0.pval"	)
		polygon(c(1,ncol(tmp),ncol(tmp),1),c(1-alpha,1-alpha,1,1), col= "grey80", border=NA)
		points(	1:ncol(tmp), tmp["H0.pval",], type='s', col="blue")
		qq<- qqplot((1:ncol(tmp)-0.5) / ncol(tmp), tmp["H0.pval",], plot.it= 0)
		plot(qq,type='p',col="blue", pch=19, xlab="expected quantiles under H0",ylab="quantiles of H0 pval", xlim=c(0,1),ylim=c(0,1))
		abline(a=0,b=1, lty= 2)
		plot(	1:ncol(tmp), tmp["H1.pval",], type='n', xlab= "replicates", ylab="flipped H0.pval"	)
		polygon(c(1,ncol(tmp),ncol(tmp),1),c(alpha,alpha,0,0), col= "grey80", border=NA)
		points(	1:ncol(tmp), tmp["H1.pval",], type='s', col="red")
		qq<- qqplot((1:ncol(tmp)-0.5) / ncol(tmp), tmp["H1.pval",], plot.it= 0)
		plot(qq,type='p',col="red",pch=19, xlab="expected quantiles under H1",ylab="quantiles of flipped H0 pval", xlim=c(0,1),ylim=c(0,1))
		abline(a=0,b=1, lty= 2)
		par(def.par)
		dev.off()
	}
	if(0)
	{
		px<- 0.4
		py<- runif(1e4, 0, px)
		#py smaller than px, this is still H0.		the pval is not uniform when py is fixed to a specific value < px
		pdf(paste(dir.name,"/binom_py_in_U(0,px).pdf",sep=''),version="1.4",width=6,height=6)
		tmp<- sapply(seq_along(py), function(i){  	project.nABC.binom.fbinom(px,ny,py[i],alternative) 		})
		def.par <- par(no.readonly = TRUE)
		layout.m<- matrix(data= c(1,1,2,3,3,4),ncol=3,nrow=2,byrow=1)
		layout(layout.m)
		plot(	1:ncol(tmp), tmp["H0.pval",], type='n',  xlab= "replicates", ylab="H0.pval"	)
		polygon(c(1,ncol(tmp),ncol(tmp),1),c(1-alpha,1-alpha,1,1), col= "grey80", border=NA)
		points(	1:ncol(tmp), tmp["H0.pval",], type='s', col="blue")
		qq<- qqplot((1:ncol(tmp)-0.5) / ncol(tmp), tmp["H0.pval",], plot.it= 0)
		plot(qq,type='p',col="blue", pch=19, xlab="expected quantiles under H0",ylab="quantiles of H0 pval", xlim=c(0,1),ylim=c(0,1))
		abline(a=0,b=1, lty= 2)
		plot(	1:ncol(tmp), tmp["H1.pval",], type='n', xlab= "replicates", ylab="flipped H0.pval"	)
		polygon(c(1,ncol(tmp),ncol(tmp),1),c(alpha,alpha,0,0), col= "grey80", border=NA)
		points(	1:ncol(tmp), tmp["H1.pval",], type='s', col="red")
		qq<- qqplot((1:ncol(tmp)-0.5) / ncol(tmp), tmp["H1.pval",], plot.it= 0)
		plot(qq,type='p',col="red",pch=19, xlab="expected quantiles under H1",ylab="quantiles of flipped H0 pval", xlim=c(0,1),ylim=c(0,1))
		abline(a=0,b=1, lty= 2)
		par(def.par)
		dev.off()
	}
	#print(tmp)
	stop()
}
#------------------------------------------------------------------------------------------------------------------------
project.nABC<- function()
{
	#more power packages		http://wiki.math.yorku.ca/index.php/R:_Power
	#binomial http://www.stat.ucl.ac.be/ISdidactique/Rhelp/library/Hmisc/html/bpower.html
	#power.t.test has df=n-1 only when variance is equal in both samples, so we cannot use this here
	UPPER<- 4
	#n<- round( 10^seq(0.8, UPPER, by=0.003) )
	n<- round( seq(10, 500, by=1) )
	RESUME<- 1
	dir.name<- "/Users/olli0601/duke/2012_frequencyABC/sim.data"

	project.nABC.interval.test.obsphylo<- function(n, what, mu, delta)
	{
		fname<- paste(CODE.HOME,"examples/NBH3N2_NL_EU101027_manyobservedsus.R",sep='')
		load(fname)
		os<- list()
		os[["DIST2ROOT"]]<- many.obs.epdf[["DIST2ROOT.slope"]]
		os[["MED.BR.LEN"]]<- many.obs.epdf[["MED.BR.LEN.swL1"]]
		os[["LINEAGE"]]<- many.obs.epdf[["ME.LINEAGE.1991"]]
		os[["MX.PH.TMRCA"]]<- many.obs.epdf[["MX.PH.TMRCA.1991"]]

		sapply(seq_along(n),function(i)
				{
					x<- sample(os[[what]],n[i])
					#mean is around 0.0057
					y<- rnorm( n[i], mu, sd= sd(x) )

					moments<- matrix(c( length(x), mean(x), var(x),  length(y), mean(y), var(y) ), 2, 3, byrow=1)
					tmp<- moments[,3]/moments[,1]										#temporarily store 		var(sim)/sim.n, var(obs)/obs.n
					tmp<- c(	diff( moments[, 2] ) / sqrt( sum(tmp) ),				#[1]  test statistic
									sqrt( sum(tmp) ),												#[2]	estimate of common standard deviation
									sum(tmp)^2 / (		tmp[2]^2/(moments[2,1]-1)	+ tmp[1]^2/(moments[1,1]-1)		)		)		#[3]	Welch Satterthwaite approximation to degrees of freedom
					ans<- c(	n[i], tmp[1], tmp[3], tmp[2],
									1-pt( abs(tmp[1])-delta/tmp[2], df= tmp[3] )+pt( -abs(tmp[1])-delta/tmp[2], df= tmp[3] )	)
					names(ans)<- c("n","statistic","df","sd","p.val")
					ans
				})
	}

	project.nABC.t.test<- function(n,mu, sd)
	{
		cat(paste("project.nABC.t.test: mu",mu,"sd",sd))
		sapply(seq_along(n),function(i)
				{
					x<- rnorm(n[i],0,sd)
					y<- rnorm(n[i],mu,sd)
					ans<- t.test(x,y)
					c(n[i], ans$statistic, ans$parameter, sqrt( var(x)/length(x) + var(y)/length(y) ), ans$p.value)
				})
	}
	project.nABC.type2error<- function(alpha,t.test.df,t.test.sd,t.test.H1,t.test.type2error)
	{
		pt( qt( 1-alpha/2, df=t.test.df ) - t.test.H1 / t.test.sd,  df= t.test.df)-pt( qt( alpha/2, df=t.test.df ) - t.test.H1 / t.test.sd,  df= t.test.df) - t.test.type2error
	}
	project.nABC.plot.type2error<- function(t.test.n,t.test.sd, t.test.df, t.test.H1)
	{

		alpha<- seq(0,1,1e-3)
		plot(1,1,type='n',xlim=c(0,1),ylim=c(0,1),xlab="alpha",ylab="beta",main=bquote(n == .(t.test.n)))
		sapply(seq_along(t.test.H1),function(i)
				{
					lines(alpha, project.nABC.type2error(alpha,t.test.df,t.test.sd,t.test.H1[i],0), lty= i )
				})
		legend("topright", legend= t.test.H1, lty= seq_along(t.test.H1), bty= 'n')
	}
	project.nABC.get.criticalvalue<- function(n, t.test.sd, t.test.df, t.test.type2error, t.test.H1)
	{
		if(length(n)!=length(t.test.sd))	stop("project.nABC.getalpha: 	error at 1a")
		if(length(n)!=length(t.test.df))	stop("project.nABC.getalpha: 	error at 1b")
		if(t.test.H1>0)
			alpha<- sapply(seq_along(n), function(i)
					{
						tmp<- uniroot(	project.nABC.type2error,	c(0,1), tol=.Machine$double.eps^0.5, t.test.df=t.test.df[i],	t.test.sd=t.test.sd[i],	t.test.H1=t.test.H1,	t.test.type2error=t.test.type2error	)
						c(tmp$root,tmp$f.root)
					})
		else
			alpha<- rbind(rep(1-t.test.type2error,length(n)),rep(0,length(n)))

		alpha[1, abs(alpha[2,])>.Machine$double.eps^0.5]<- NA		#cannot compute alpha reliably when tau is too large so that the beta~alpha curve approaches a discontinuity in (alpha=0,beta=1)
		alpha<- rbind(n,alpha, qt(alpha[1,]/2, df= t.test.df) )
		alpha<- rbind(alpha, -alpha[4,] )
		rownames(alpha)<- c("n","alpha","uniroot.tol","cil","cir")
		alpha
	}

	what<- "DIST2ROOT"
	beta<- c(0.8, 0.5, 0.05)
	mu<- c(0.0035, 0.0037, 0.0039)
	m.DISTROOT<- lapply(mu, function(x){		project.nABC.interval.test.obsphylo(n, what, x, 0.002)		})
	names(m.DISTROOT)<- mu

	sapply(seq_along(m.DISTROOT),function(i)
			{
				m<- m.DISTROOT[[i]]
				pdf(paste(dir.name,paste("intervaltest_",what,"_mu",names(m.DISTROOT)[i],".pdf",sep=''),sep='/'),version="1.4",width=6,height=6)
				plot(m["n",],m["p.val",],type='l', xlab="sample size",ylab="p-value",main=bquote(mu==.(as.numeric(names(m.DISTROOT))[i] )))
				dev.off()
			})

stop()
	m<- m.DISTROOT[["0.004"]]

	t.test.H1<- c(0,1e-5, 3e-5, 6e-5, 1e-4)
	project.nABC.plot.type2error(m[1,491], m[4,491], m[3,491], t.test.H1)
	project.nABC.plot.type2error(m[1,91], m[4,91], m[3,91], t.test.H1)





	tmp<- project.nABC.get.criticalvalue(m[1,], m[4,], m[3,], 0.8, 1e-4)
	print(tmp)

	plot(m[1,],m[2,],type='s', xlab="sample size", ylab="t-test")
	sapply(seq_along(beta), function(i){			lines(m[1,], project.nABC.get.criticalvalue(m[1,], m[4,], m[3,], beta[i], 1e-5)["cir",], col='blue', lty=i,lwd=1.25)				})
	#sapply(seq_along(beta), function(i){			lines(m[1,], project.nABC.get.criticalvalue(m[1,], m[4,], m[3,], beta[i], 0.01)["cil",], col='red', lty=i,lwd=1.25)				})
	#sapply(seq_along(beta), function(i){			lines(m[1,], project.nABC.get.criticalvalue(m[1,], m[4,], m[3,], beta[i], 0.02)["cil",], col='green', lty=i,lwd=1.25)				})
	legend("topleft", legend= beta, lty= seq_along(beta), bty= 'n')
	legend("bottomright",legend=c(0,0.01,0.02),fill=c("blue","red","green"), bty='n')



	stop()

	if(RESUME)		#//load if there is R file
	{
		options(show.error.messages = FALSE)
		readAttempt<-try(suppressWarnings(load(		paste(dir.name,paste("ttestasymptotics_",UPPER,".R",sep=''),sep='/')					)))
		options(show.error.messages = TRUE)
	}
	if(!RESUME || inherits(readAttempt, "try-error"))
	{
		m.0<- project.nABC.t.test(n, 0, 1)
		m.0.01<- project.nABC.t.test(n, 1e-2, 1)
		m.0.02<- project.nABC.t.test(n, 2e-2, 1)
		m.0.1<- project.nABC.t.test(n, 0.1, 1)
		save(m.0, m.0.01, m.0.02, m.0.1, file=paste("/Users/olli0601/duke/2012_frequencyABC/sim.data/ttestasymptotics_",UPPER,".R",sep=''))
	}
	#cat(print.m(ans, print.char=0,as.R=1))

	#dependence on power and ABC threshold, case mu=0.01
	m<- m.0.01
	#power for n= 100
	pdf(paste(dir.name,"ttest_mu0.01_power_n100.pdf",sep='/'),version="1.4",width=6,height=6)
	t.test.H1<- c(0,0.04,0.08, 0.16,0.32)
	project.nABC.plot.type2error(m[1,121], m[4,121], m[3,121], t.test.H1)
	dev.off()
	#power for n= 1000
	pdf(paste(dir.name,"ttest_mu0.01_power_n1000.pdf",sep='/'),version="1.4",width=6,height=6)
	t.test.H1<- c(0,0.04,0.08,0.16,0.32)
	project.nABC.plot.type2error(m[1,221], m[4,221], m[3,221], t.test.H1)
	dev.off()

	beta<- c(0.8, 0.5, 0.05)

	pdf(paste(dir.name,paste("ttest_",UPPER,"_mu0.01.pdf",sep=''),sep='/'),version="1.4",width=6,height=6)
	plot(m[1,],m[2,],type='s', xlab="sample size", ylab="t-test")
	sapply(seq_along(beta), function(i){			lines(m[1,], project.nABC.get.criticalvalue(m[1,], m[4,], m[3,], beta[i], 0)["cil",], col='blue', lty=i,lwd=1.25)				})
	sapply(seq_along(beta), function(i){			lines(m[1,], project.nABC.get.criticalvalue(m[1,], m[4,], m[3,], beta[i], 0.01)["cil",], col='red', lty=i,lwd=1.25)				})
	sapply(seq_along(beta), function(i){			lines(m[1,], project.nABC.get.criticalvalue(m[1,], m[4,], m[3,], beta[i], 0.02)["cil",], col='green', lty=i,lwd=1.25)				})
	legend("bottomleft", legend= beta, lty= seq_along(beta), bty= 'n')
	legend("topright",legend=c(0,0.01,0.02),fill=c("blue","red","green"), bty='n')
	dev.off()

	pdf(paste(dir.name,paste("ttest_",UPPER,"_mu0.02.pdf",sep=''),sep='/'),version="1.4",width=6,height=6)
	m<- m.0.02
	plot(m[1,],m[2,],type='s', xlab="sample size", ylab="t-test")
	sapply(seq_along(beta), function(i){			lines(m[1,], project.nABC.get.criticalvalue(m[1,], m[4,], m[3,], beta[i], 0)["cil",], col='blue', lty=i,lwd=1.25)				})
	sapply(seq_along(beta), function(i){			lines(m[1,], project.nABC.get.criticalvalue(m[1,], m[4,], m[3,], beta[i], 0.02)["cil",], col='red', lty=i,lwd=1.25)				})
	sapply(seq_along(beta), function(i){			lines(m[1,], project.nABC.get.criticalvalue(m[1,], m[4,], m[3,], beta[i], 0.04)["cil",], col='green', lty=i,lwd=1.25)				})
	legend("bottomleft", legend= beta, lty= seq_along(beta), bty= 'n')
	legend("topright",legend=c(0,0.02,0.04),fill=c("blue","red","green"), bty='n')
	dev.off()
	stop()
	pdf(paste(dir.name,"ttest_mu0.1.pdf",sep='/'),version="1.4",width=6,height=6)
	m<- m.0.1
	plot(m[1,],m[2,],type='s', xlab="sample size", ylab="t-test")
	sapply(seq_along(beta), function(i){			lines(m[1,], project.nABC.getalpha(m[1,], m[4,], m[3,], beta[i], 0.05), col='blue', lty=i)				})
	sapply(seq_along(beta), function(i){			lines(m[1,], project.nABC.getalpha(m[1,], m[4,], m[3,], beta[i], 0.1), col='red', lty=i)				})
	sapply(seq_along(beta), function(i){			lines(m[1,], project.nABC.getalpha(m[1,], m[4,], m[3,], beta[i], 0.2), col='green', lty=i)				})
	legend("bottomleft", legend= beta, lty= seq_along(beta), bty= 'n')
	legend("topright",legend=c(0.05,0.1,0.2),fill=c("blue","red","green"), bty='n')
	dev.off()
	stop()
}


