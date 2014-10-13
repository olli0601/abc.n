#
#	illustrate calibration of tau.l for given tau.u such that the mode of the
#	power function is at rho.star
#
alpha		<- 0.01
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
