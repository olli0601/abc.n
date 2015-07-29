# power function of the two-sided F-test, to test equality of means for multivariate
# normal samples with unknown covariance matrix

n		<- 10
p		<- 3
cu		<- 8
rho 	<- seq(0, 2, length = 1024)
tmp <- lapply(c(0.01, 0.5, 1, 2, 3, 4, 5, 6, 7), function(cl)
		{
			data.table(cl = as.factor(cl), cu=cu, rho = rho, power = tsftest.pow(rho, cl, cu, n, p))
		})
tmp	<- do.call('rbind', tmp)
pp <- ggplot(tmp, aes(x = rho, y = power, colour = cl, group = cl)) + geom_line() + labs(y = 'Power\n(ABC acceptance probability)')
print(pp)
