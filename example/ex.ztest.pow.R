# power function of the Z-test,
# to test the equivalence of autocorrelations of normal summary values in the MA(1) model

n2s		<- function(n){ 1/sqrt(floor(n)-3) }	#need formula to convert n.of.y into s.of.T
s2n		<- function(s){ (1/s)^2+3 }				#need formula to convert s.of.T into n.of.y
n.of.y	<- 1500		#number of independent pairs of the form (y_t, y_{t+1}) for simulated summary values y_t, t=1, ...

#	reasonable power properties of the ABC accept/reject step

rho		<- seq(-0.3,0.3,0.001)
tau.u	<- 0.09
tmp		<- data.table(rho=rho, power=ztest.pow(rho, tau.u, alpha=0.01, n2s(n.of.y)))
p		<- ggplot(tmp,aes(x=rho,y=power)) + geom_line() + labs(y='Power\n(ABC acceptance probability)')
print(p)

# power increases with tau.u and becomes flat

tmp <- lapply(c(0.06,0.07,0.08,0.1,0.15,0.2),function(tau.u)
		{
			data.table(tau.u=as.factor(tau.u),rho=rho,power=ztest.pow(rho, tau.u, alpha=0.01, n2s(n.of.y)), tau.u=tau.u)
		})
tmp	<- do.call('rbind', tmp)
p <- ggplot(tmp,aes(x=rho,y=power,colour=tau.u,group=tau.u)) + geom_line() + labs(y='Power\n(ABC acceptance probability)')
print(p)

# power increases with n.of.y and becomes flat

tau.u <- 0.1
tmp <- lapply(c(750,1000,1500,3000),function(n.of.y)
		{
			data.table(n.of.y=as.factor(n.of.y),rho=rho,power=ztest.pow(rho, tau.u, alpha=0.01, n2s(n.of.y)),n.of.y=n.of.y)
		})
tmp	<- do.call('rbind', tmp)
p <- ggplot(tmp,aes(x=rho,y=power,colour=n.of.y,group=n.of.y)) + geom_line() + labs(y='Power\n(ABC acceptance probability)')
print(p)
