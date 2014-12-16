alpha	<- 0.01
n.of.x	<- 60
n.of.y	<- 60

chisqstretch.calibrate(n.of.x=n.of.x, n.of.y=n.of.y, tau.u=2, what='MXPW_AT_EQU', alpha=0.01, plot=FALSE, verbose=FALSE)

scale		<- n.of.x<- 60
df			<- 299
tmp			<- chisqstretch.calibrate.taulow(2, scale, df, alpha, rho.star=1, tol= 1e-5, max.it=100, verbose=1)
rho			<- seq(0.1, 5, len=1024)
pw			<- chisqstretch.pow(rho, n.of.x, df, tmp["cl"], tmp["cu"])
plot(rho,pw,col="black", type='l')
tmp			<- chisqstretch.calibrate.taulow(1.8, scale, df, alpha, rho.star=1, tol= 1e-5, max.it=100, verbose=1)
pw			<- chisqstretch.pow(rho, n.of.x, df, tmp["cl"], tmp["cu"])
lines(rho,pw,col="blue")
tmp			<- chisqstretch.calibrate.taulow(1.5, scale, df, alpha, rho.star=1, tol= 1e-5, max.it=100, verbose=1)
pw			<- chisqstretch.pow(rho, n.of.x, df, tmp["cl"], tmp["cu"])
lines(rho,pw,col="red")
tmp			<- chisqstretch.calibrate.taulow(1.3, scale, df, alpha, rho.star=1, tol= 1e-5, max.it=100, verbose=1)
pw			<- chisqstretch.pow(rho, n.of.x, df, tmp["cl"], tmp["cu"])
lines(rho,pw,col="green")	

#
#	illustrate calibration of tau.l, tau.u and m 
#
n.of.x		<- 60
x			<- rnorm(n.of.x,0,1)
s.of.x		<- sd(x)
abc.param	<- chisqstretch.calibrate(n.of.x, s.of.x, scale=n.of.x, n.of.y=n.of.x, mx.pw=0.9, alpha=0.01, max.it=100, debug=FALSE, plot=TRUE)
abc.param
