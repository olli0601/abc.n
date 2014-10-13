#
#	ABC* calibrations on the Moving average (order 1) example
#
xa				<- 0.1
xsigma2			<- 1	
xn				<- 150
leave.out.a		<- 2
leave.out.sig2	<- 1
#	load precomputed MCMC output targeting the exact posterior for a=0.1, sig2=1, n=150. The MCMC chain is thinned.
#	in object 'ma.exact'
#	see also 'data(package='abc.star')'
data('ma_mcmc_a=0.1')
x					<- ma.exact$data$x
#	exact MAP 
x.map				<- ma.get.2D.mode(ma.exact$posterior[,a],ma.exact$posterior[,sig2], xlim= c(-0.4,0.4),ylim=c(0.6,1/0.6),plot=0, nbin=10,  method="ash")
#	exact MAP on auxiliary space
x.map.on.rho		<- ma.rho2a( ma.exact$data$unthinned$s_stat$autocorr )
x.map.on.rho		<- c( x.map.on.rho, ma.rho2sig2( ma.exact$data$unthinned$s_stat$variance, x.map.on.rho ) )					
#	load precomputed ABC* output targeting the ABC* posterior for a=0.1, sig2=1, n=150. 
#	This is only the first 1e6 samples. For publication, we used 5e6 samples so results may differ.
#	in object 'ma.abc.star'
data('ma_abc.star_a=0.1')
#	convert data back to double
ma.abc.star$data	<- matrix(ma.abc.star$data, nrow=9, dimnames=list(c("th.a","rho.a", "T.a", "T.a2", "T.a3", "th.s2", "rho.s2",  "T.s2",  "T.s22"),c()) ) / 1e4
#	calibrate ABC*
#	5 sets of summary values: 2 for the variance test and 3 for the autocorrelation test
#	abc.param.sig2.2 is the same as abc.param.sig2.1 because the length of the two sets of summary values is equal; only included for illustration
#	abc.param.a.2	is the same as abc.param.a.1 for the same reason
suppressWarnings({
	zx					<- ma.cor(x, leave.out=leave.out.a)
	abc.param.a.1		<- corrz.calibrate(zx["n"], mx.pw=0.9, alpha=0.01, max.it=100, pow_scale=2, debug=FALSE, plot=FALSE)					
	zx					<- ma.cor(x[-1], leave.out=leave.out.a)
	abc.param.a.2		<- corrz.calibrate(zx["n"], mx.pw=0.9, alpha=0.01, max.it=100, pow_scale=2, debug=FALSE, plot=FALSE)						
	zx					<- ma.cor(x[-c(1,2)], leave.out=leave.out.a)
	abc.param.a.3		<- corrz.calibrate(zx["n"], mx.pw=0.9, alpha=0.01, max.it=100, pow_scale=2, debug=FALSE, plot=FALSE)								
	vx					<- x[seq.int(1,length(x),by=1+leave.out.sig2)]
	abc.param.sig2.1	<- chisqstretch.calibrate(length(vx), sd(vx), mx.pw=0.9, alpha=0.01, max.it=100, debug=FALSE, plot=FALSE)
	vx					<- x[seq.int(2,length(x),by=1+leave.out.sig2)]
	abc.param.sig2.2	<- chisqstretch.calibrate(length(vx), sd(vx), mx.pw=0.9, alpha=0.01, max.it=100, debug=FALSE, plot=FALSE)
		})	
#
#	run using all 5 sets of summary values
#
acc					<- which( 	ma.abc.star[["data"]]["T.s2",]>=abc.param.sig2.1["cl"]  &  ma.abc.star[["data"]]["T.s2",]<=abc.param.sig2.1["cu"]	&
								ma.abc.star[["data"]]["T.s22",]>=abc.param.sig2.2["cl"]  &  ma.abc.star[["data"]]["T.s22",]<=abc.param.sig2.2["cu"]	&
								ma.abc.star[["data"]]["T.a",]*sqrt(abc.param.a.1["n.of.y"]-3)>=abc.param.a.1["cl"]  &  ma.abc.star[["data"]]["T.a",]*sqrt(abc.param.a.1["n.of.y"]-3)<=abc.param.a.1["cu"]	&
								ma.abc.star[["data"]]["T.a2",]*sqrt(abc.param.a.2["n.of.y"]-3)>=abc.param.a.2["cl"]  &  ma.abc.star[["data"]]["T.a2",]*sqrt(abc.param.a.2["n.of.y"]-3)<=abc.param.a.2["cu"] &
								ma.abc.star[["data"]]["T.a3",]*sqrt(abc.param.a.3["n.of.y"]-3)>=abc.param.a.3["cl"]  &  ma.abc.star[["data"]]["T.a3",]*sqrt(abc.param.a.3["n.of.y"]-3)<=abc.param.a.3["cu"]
							)
acc.prob			<- length(acc)/ncol(ma.abc.star[["data"]])					
tmp					<- ma.get.2D.mode(ma.abc.star[["data"]]["th.a",acc],ma.abc.star[["data"]]["th.s2",acc], xlim= c(-0.5,0.5),ylim=c(0.5,2),plot=1, nbin=10, levels=c(1,3,5,10), method="ash", xlab="a", ylab=expression(sigma^2), cols=head( gray(seq(.3,.7,len=50)), 50))
abline(h=xsigma2, lty=2)
abline(v=xa, lty=2)
dist.MAP			<- sqrt(sum(c(tmp-x.map)^2))
dist.MAP.on.rho		<- sqrt(sum(c(tmp-x.map.on.rho)^2))
ma.add.contour(ma.exact$posterior[,a], ma.exact$posterior[,sig2], levels=c(1,3,5,10), contour.col="white")
points(x.map, pch=18, col="white")										
df1					<- data.table(	th1=ma.abc.star[["data"]]["th.a",acc],	th2=ma.abc.star[["data"]]["th.s2",acc]	)			
df2					<- data.table(	th1=ma.exact$posterior[,a], th2=ma.exact$posterior[,sig2]			)
kl					<- kl.2D(df1, df2, nbin=100)$two
print(c(acc.prob, kl, dist.MAP, dist.MAP.on.rho))
#
#	run using only variance test 
#	(link function not bijective)
#
acc					<- which( 	ma.abc.star[["data"]]["T.s2",]>=abc.param.sig2.1["cl"]  &  ma.abc.star[["data"]]["T.s2",]<=abc.param.sig2.1["cu"]	&
								ma.abc.star[["data"]]["T.s22",]>=abc.param.sig2.2["cl"]  &  ma.abc.star[["data"]]["T.s22",]<=abc.param.sig2.2["cu"]					
							)
acc.prob			<- length(acc)/ncol(ma.abc.star[["data"]])					
tmp					<- ma.get.2D.mode(ma.abc.star[["data"]]["th.a",acc],ma.abc.star[["data"]]["th.s2",acc], xlim= c(-0.5,0.5),ylim=c(0.5,2),plot=1, nbin=10, levels=c(1,3,5,10), method="ash", xlab="a", ylab=expression(sigma^2), cols=head( gray(seq(.3,.7,len=50)), 50))
abline(h=xsigma2, lty=2)
abline(v=xa, lty=2)
dist.MAP			<- sqrt(sum(c(tmp-x.map)^2))
dist.MAP.on.rho		<- sqrt(sum(c(tmp-x.map.on.rho)^2))
ma.add.contour(ma.exact$posterior[,a], ma.exact$posterior[,sig2], levels=c(1,3,5,10), contour.col="white")
points(x.map, pch=18, col="white")										
df1					<- data.table(	th1=ma.abc.star[["data"]]["th.a",acc],	th2=ma.abc.star[["data"]]["th.s2",acc]	)			
df2					<- data.table(	th1=ma.exact$posterior[,a], th2=ma.exact$posterior[,sig2]			)
kl					<- kl.2D(df1, df2, nbin=100)$two
print(c(acc.prob, kl, dist.MAP, dist.MAP.on.rho))
#
#	run using only 2 sets of summary values, effectively only using a subset of the data 
#	(summary values not sufficient for theta)
#
acc					<- which( 	ma.abc.star[["data"]]["T.s2",]>=abc.param.sig2.1["cl"]  &  ma.abc.star[["data"]]["T.s2",]<=abc.param.sig2.1["cu"]	&								
								ma.abc.star[["data"]]["T.a",]*sqrt(abc.param.a.1["n.of.y"]-3)>=abc.param.a.1["cl"]  &  ma.abc.star[["data"]]["T.a",]*sqrt(abc.param.a.1["n.of.y"]-3)<=abc.param.a.1["cu"]								
							)
acc.prob			<- length(acc)/ncol(ma.abc.star[["data"]])					
tmp					<- ma.get.2D.mode(ma.abc.star[["data"]]["th.a",acc],ma.abc.star[["data"]]["th.s2",acc], xlim= c(-0.5,0.5),ylim=c(0.5,2),plot=1, nbin=10, levels=c(1,3,5,10), method="ash", xlab="a", ylab=expression(sigma^2), cols=head( gray(seq(.3,.7,len=50)), 50))
abline(h=xsigma2, lty=2)
abline(v=xa, lty=2)
dist.MAP			<- sqrt(sum(c(tmp-x.map)^2))
dist.MAP.on.rho		<- sqrt(sum(c(tmp-x.map.on.rho)^2))
ma.add.contour(ma.exact$posterior[,a], ma.exact$posterior[,sig2], levels=c(1,3,5,10), contour.col="white")
points(x.map, pch=18, col="white")										
df1					<- data.table(	th1=ma.abc.star[["data"]]["th.a",acc],	th2=ma.abc.star[["data"]]["th.s2",acc]	)			
df2					<- data.table(	th1=ma.exact$posterior[,a], th2=ma.exact$posterior[,sig2]			)
kl					<- kl.2D(df1, df2, nbin=100)$two
print(c(acc.prob, kl, dist.MAP, dist.MAP.on.rho))
