#
#	illustrate power function of the autocorrelation test for normal summary values
#
alpha	<- 0.01
tau.u	<- 0.09
tau.l	<- -tau.u
sim.n	<-	5e3
rho		<- seq(tau.l,tau.u,0.001)
pw		<- corrz.pow(rho, tau.u, alpha, 1/sqrt(floor(sim.n/3)-3))
plot(rho,pw,col="red", ylim=c(0,1), type='l')
#	flat power function when the number of simulated summary values is much lower	
sim.n	<-	5e2
pw		<- corrz.pow(rho, tau.u, alpha, 1/sqrt(floor(sim.n/3)-3))
lines(rho,pw,col="blue")
#	flat power function when the number of simulated summary values is much higher	
sim.n	<-	2e4
pw		<- corrz.pow(rho, tau.u, alpha, 1/sqrt(floor(sim.n/3)-3))
lines(rho,pw,col="green")
