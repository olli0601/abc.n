# power function of the F-test, to test equality of means for multivariate
# normal samples with unknown covariance matrix

#set number of variables (i.e. summary statistics)
p <- 3
#set number of simulations
n <- 100

#calculate power for fixed equivalence value
tau <- 1.2
rho <- seq(0, .5, length = 1024)
ftest.pow(rho, tau, n = n, p = p)

# power increases as size of equivalence region increases but power function
# flattens out as equivalence region gets large
tmp <- lapply(c(0.05, 0.1, 0.2, 0.3), function(tau)
		{
			data.table(tau = as.factor(tau), rho = rho, power = ftest.pow(rho, tau, n, p, alpha = 0.01))
		})
tmp	<- do.call('rbind', tmp)
pp <- ggplot(tmp, aes(x = rho, y = power, colour = tau, group = tau)) + geom_line() + labs(y = 'Power\n(ABC acceptance probability)')
print(pp)

# power increases as number of simulations increase
tau 	<- 0.2
rho		<- seq(0, .3, length = 1024)
tmp  	<- lapply(c(25, 50, 100, 200, 400), function(n)
		{
			data.table(n = as.factor(n), rho = rho, y = ftest.pow(rho, tau, n, p, alpha = 0.01), d='power')
		})
tmp		<- do.call('rbind', tmp)
pp 		<- ggplot(tmp, aes(x = rho, y = y, colour = n, group = n)) + geom_line() + labs(y = 'Power\n(ABC acceptance probability)')
print(pp)

# add likelihood density to last power plot
t2.x	<- 0.25
tmp		<- rbind(tmp, data.table(n=n, rho=rho, y=ftest.sulkl(rho, t2.x, n, p, norm = 1, support= c(0,Inf), log=FALSE), d='prtl.lkl'))
pp 		<- ggplot(tmp, aes(x = rho, y = y, colour = n, linetype=d, group = interaction(n,d))) + geom_line() + labs(y = 'Power\n(ABC acceptance probability)')
print(pp)

tmp		<- rbind(tmp, data.table(n=n, rho=rho, y=ftestz.sulkl(rho, n, p, norm = 1, support= c(0,Inf), log=FALSE), d='prtl.lkl'))
pp 		<- ggplot(tmp, aes(x = rho, y = y, colour = n, linetype=d, group = interaction(n,d))) + geom_line() + labs(y = 'Power\n(ABC acceptance probability)')
print(pp)
