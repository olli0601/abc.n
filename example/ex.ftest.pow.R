# power function of the F-test, to test equality of means for multivariate
# normal samples with unknown covariance matrix

#set number of variables (i.e. summary statistics)
p <- 3
#set number of simulations
n <- 100

#calculate power for fixed equivalence value
tau2 <- 1.2
rho2 <- seq(0, .5, length = 1024)
ftest.pow(rho2, tau2, n = n, p = p)

# power increases as size of equivalence region increases but power function
# flattens out as equivalence region gets large
tmp <- lapply(c(0.05, 0.1, 0.2, 0.3), function(tau2)
		{
			data.table(tau2 = as.factor(tau2), rho2 = rho2, power = ftest.pow(rho2, tau2, alpha = 0.01, n, p))
		})
tmp	<- do.call('rbind', tmp)
pp <- ggplot(tmp, aes(x = rho2, y = power, colour = tau2, group = tau2)) + geom_line() + labs(y = 'Power\n(ABC acceptance probability)')
print(pp)

# power increases as number of simulations increase
tau2 <- 1.2
rho2 <- seq(0, 1.5, length = 1024)
tmp  <- lapply(c(5, 10, 20, 30, 40, 50), function(n)
		{
			data.table(n = as.factor(n), rho2 = rho2, power = ftest.pow(rho2, tau2, alpha = 0.01, n, p))
		})
tmp	<- do.call('rbind', tmp)
pp <- ggplot(tmp, aes(x = rho2, y = power, colour = n, group = n)) + geom_line() + labs(y = 'Power\n(ABC acceptance probability)')
print(pp)
