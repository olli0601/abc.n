n.of.y 	<- 30
s.of.T	<- 0.2		#this is s.of.y / squrt(n.of.y)

# useful equivalence region for given stochasticity in T: 
# power is not zero, does not plateau around 1, and still high around rho=0

rho		<- seq(-2,2,0.01)
tmp	<- data.table(rho=rho, power=mutost.pow(rho, df=n.of.y-1, s.of.T, tau.u=0.8, alpha=0.01))
p	<- ggplot(tmp,aes(x=rho,y=power)) + geom_line() + labs(y='Power\n(ABC acceptance probability)')
print(p)

# power increases with tau.u and becomes flat

tmp <- lapply(seq(0.5,1.6,0.2),function(tau.u)
		{
			data.table(tau.u=as.factor(tau.u),rho=rho,power=mutost.pow(rho, df=n.of.y-1, s.of.T, tau.u=tau.u, alpha=0.01))
		})
tmp	<- do.call('rbind', tmp)
p <- ggplot(tmp,aes(x=rho,y=power,colour=tau.u,group=tau.u)) + geom_line() + labs(y='Power\n(ABC acceptance probability)')
print(p)

# power increases with n.of.y and becomes flat

tau.u 	<- 0.8
s.of.y	<- 0.8
tmp <- lapply(c(10,15,20,25,50,100),function(n.of.y)
		{
			data.table(n.of.y=as.factor(n.of.y),rho=rho,power=mutost.pow(rho, df=n.of.y-1, s.of.y/sqrt(n.of.y), tau.u=tau.u, alpha=0.01))
		})
tmp	<- do.call('rbind', tmp)
p <- ggplot(tmp,aes(x=rho,y=power,colour=n.of.y,group=n.of.y)) + geom_line() + labs(y='Power\n(ABC acceptance probability)')
print(p)
