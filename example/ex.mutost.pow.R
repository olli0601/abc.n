require(plyr)
require(ggplot2)

rho 	<- seq(-2,2,0.01)
n.of.y 	<- 30
s.of.T 	<- 0.2
alpha 	<- 0.01

# useful equivalence region for given stochasticity in T: 
# power is not zero, does not plateau around 1, and still high around rho=0
tmp		<- data.frame(rho=rho, power=mutost.pow(rho, df=n.of.y-1, s.of.T, tau.u=0.8, alpha))
p <- ggplot(tmp,aes(x=rho,y=power)) + geom_line() + labs(y='Power\n(ABC acceptance probability)')
print(p)

# power increases with tau.u and becomes flat
tmp <- ldply(seq(0.2,1.6,0.2),function(tau.u)
		{
			return(data.frame(tau.u=as.factor(tau.u),rho=rho,power=mutost.pow(rho, df=n.of.y-1, s.of.T, tau.u=tau.u, alpha)))
		})

p <- ggplot(tmp,aes(x=rho,y=power,colour=tau.u)) + geom_line() + labs(y='Power\n(ABC acceptance probability)')
print(p)


# power increases with n.of.y and becomes flat
tau.u <- 1.5
tmp <- ldply(c(seq(1,5,1),seq(20,100,20)),function(n.of.y)
		{
			return(data.frame(n.of.y=as.factor(n.of.y),rho=rho,power=mutost.pow(rho, df=n.of.y-1, s.of.T, tau.u, alpha)))
		})

p <- ggplot(tmp,aes(x=rho,y=power,colour=n.of.y)) + geom_line() + labs(y='Power\n(ABC acceptance probability)')
print(p)
