#
#	illustrate power function of the variance test for normal summary values
#
alpha		<- 0.01
n.of.x		<- 60 
n.of.y		<- 60
df			<- n.of.y-1
rho			<- seq(0.1, 5, len=1024)
#	compute tolerances
tmp			<- .Call("abcScaledChiSq",	c(n.of.x, df, 1/2.2, 2.2, alpha, 1e-10, 100, 0.05)	)
#	compute power function for given tolerances
pw			<- chisqstretch.pow(rho, n.of.x, df, tmp[1], tmp[2])
plot(rho,pw,col="red", type='l')
lines(rho,pw,col="red")
#	repeat for more simulated summary values
n.of.y		<- 110
df			<- n.of.y-1
tmp			<- .Call("abcScaledChiSq",	c(n.of.x, df, 1/2.2, 2.2, alpha, 1e-10, 100, 0.05)	)
pw			<- chisqstretch.pow(rho, n.of.x, df, tmp[1], tmp[2])
plot(rho,pw,col="red", type='l')
lines(rho,pw,col="red")
