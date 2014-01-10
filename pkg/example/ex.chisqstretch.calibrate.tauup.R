#
#	illustrate calibration of tau.l and tau.u for given target maximum power at rho.star
#
alpha		<- 0.01
scale		<- n.of.x<- 60
n.of.y		<- 299
df			<- n.of.y-1
x			<- rnorm(n.of.x,0,1)
s.of.x		<- sd(x)
tmp			<- chisqstretch.calibrate.tauup(0.9, 3*s.of.x, scale, df, alpha, rho.star=1, tol= 1e-5, max.it=100, verbose=1)
rho			<- seq(0.4, 2, len=1024)
pw			<- chisqstretch.pow(rho, n.of.x, df, tmp["cl"], tmp["cu"])		
plot(rho,pw,col="black", type='l')
