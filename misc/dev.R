start_me <- function() {

	dir_pkg <<- path.expand("~/work/projects/abc_star/git/abc.n")

}

dev <- function() {

	document(dir_pkg)
	# check(dir_pkg)

}

test <- function() {

	?ftest.calibrate
# set number of variables (i.e. summary statistics)
	p <- 3
# set number of observations
	n.of.x <- 100
# set number of simulations
	n.of.y <- 100
# set T2 calculated on the observed data
	t2.x <- 0.25
	ftest.calibrate(n.of.x=n.of.x, t2.x=t2.x, p=p, what='KL', mx.pw=0.9, alpha=0.01, use.R= TRUE, debug=TRUE, plot=TRUE, verbose=FALSE)

	ftest.calibrate(n.of.x=n.of.x, t2.x=t2.x, p=p, what='KL', mx.pw=0.9, alpha=0.01, use.R= FALSE, plot=FALSE, verbose=FALSE)

}

fix_ratetest <- function() {

	ans <- ratetest.calibrate(n.of.x=30, mean.x=2, n.of.y=30, what='KL', plot=TRUE, debug = TRUE, verbose = TRUE)
	print(ans)
	# ratetest.getkl(mean.x = 2, n.of.y=40, tau.u = 2*2, mx.pw = 0.9, alpha = 0.01, pow_scale=1.5, calibrate.tau.u=TRUE, plot=TRUE)

}

main <- function() {

	start_me()
	dev()
	fix_ratetest()

}

main()