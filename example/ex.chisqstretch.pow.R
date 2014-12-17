n.of.x		<- 60 
n.of.y		<- 60

#	compute ABC tolerances

cali	<- vartest.calibrate(n.of.x=n.of.x, n.of.y=n.of.y, tau.l=1/2, tau.u=2, what='CR', alpha=0.01)

# problematic ABC tolerances for ABC inference: 
# although power is not zero, does not plateau around 1, and still high around rho=1 (desirable),
# the power is not maximised at the point of equality (rho=1).

rho		<- seq(0.1, 3, len=1024)
tmp		<- data.frame(rho=rho, power=vartest.pow(rho, n.of.x, n.of.y-1, cali['c.l'], cali['c.u']))
p 		<- ggplot(tmp,aes(x=rho,y=power)) + geom_line() + labs(y='Power\n(ABC acceptance probability)')
print(p)

# power increases with tau.u and becomes flat

tmp	<- lapply(seq(1.2,2.8,0.4),function(tau.u)
		{
			cali	<- vartest.calibrate(n.of.x=n.of.x, n.of.y=n.of.y, tau.l=1/tau.u, tau.u=tau.u, what='CR', alpha=0.01)
			data.table(rho=rho, power=vartest.pow(rho, n.of.x, n.of.y-1, cali['c.l'], cali['c.u']), tau.u=tau.u)			
		})
tmp	<- do.call('rbind', tmp)
p	<- ggplot(tmp,aes(x=rho,y=power,colour=tau.u, group=tau.u)) + geom_line() + labs(x=expression(rho), y='Power\n(ABC acceptance probability)')
print(p)
