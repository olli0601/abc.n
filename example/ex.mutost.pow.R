require(plyr)
require(ggplot2)

rho <- seq(-2,2,0.01)
n.of.y <- 30
s.of.T <- 0.2
alpha <- 0.01

# power increases with tau.u and becomes flat
tmp <- ldply(seq(0.2,1.6,0.2),function(tau.u){
	return(data.frame(tau.u=as.factor(tau.u),rho=rho,power=mutost.pow(rho, df=n.of.y-1, s.of.T, tau.u=tau.u, alpha)))
})

p <- ggplot(tmp,aes(x=rho,y=power,colour=tau.u))
p <- p+geom_line()
print(p)


# power increases with n.of.y and becomes flat
tau.u <- 1.5
tmp <- ldply(c(seq(1,5,1),seq(20,100,20)),function(n.of.y){
	return(data.frame(n.of.y=as.factor(n.of.y),rho=rho,power=mutost.pow(rho, df=n.of.y-1, s.of.T, tau.u, alpha)))
})

p <- ggplot(tmp,aes(x=rho,y=power,colour=n.of.y))
p <- p+geom_line()
print(p)