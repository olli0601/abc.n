n.of.y <- 40

# compute ABC tolerances
cali	<- ratetest.calibrate(tau.l=1/2, tau.u=2, n.of.y=n.of.y, what='CR', alpha=0.01)
# compute the power for the range (0.1, 3)
rho  	<- seq(0.1, 3, len=1024)
tmp		<- data.frame(rho=rho, power=ratetest.pow(rho, cali['c.l'], cali['c.u'], m=n.of.y))

library(ggplot2)
p 		<- ggplot(tmp,aes(x=rho,y=power)) + geom_line() + labs(y='Power\n(ABC acceptance probability)')
print(p)



