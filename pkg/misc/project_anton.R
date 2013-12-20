
nabc_MA1_compute_rho_bounds <- function(a_bounds, sig2_bounds, variance, autocorr) {

	#compute bounds on rho1 and rho2 based on those of a and sig2
	rho_1_bounds <- nabc.acf.sig22rho(sort(sig2_bounds), a = c(ifelse(prod(a_bounds) <= 0, 0, min(abs(a_bounds))), 
		max(abs(a_bounds))), vx = variance)
	rho_2_bounds <- nabc.acf.a2rho(x = a_bounds, vx = autocorr)

	return(list(rho_1 = rho_1_bounds, rho_2 = rho_2_bounds))
}



nabc_MA1_rprior <- function(sample_size=1, a_bounds = c(-0.3, 0.3), sig2_bounds = c(0.5, 2), prior_dist=c("uniform","uniform_on_rho"), variance=NULL,autocorr=NULL){
	
	prior_dist <- match.arg(prior_dist)
	
	if(prior_dist=="uniform"){
		a <- runif(sample_size,min=min(a_bounds),max=max(a_bounds))
		sig2 <- runif(sample_size,min=min(sig2_bounds),max=max(sig2_bounds))

	}
	
	if(prior_dist=="uniform_on_rho"){
		
		#estimate of nu must be provided to compute rho
		stopifnot(!is.null(variance),!is.null(autocorr))
		
		#compute bounds on rho1 and rho2 based on those of a and sig2
		rho_bounds <- nabc_MA1_compute_rho_bounds(a_bounds, sig2_bounds, variance, autocorr)
		
		#sample a from prior induced by uniform prior on rho
		rho_2 <- runif(sample_size, rho_bounds$rho_2[1], rho_bounds$rho_2[2]) #uniform on rho
		a <- nabc.acf.rho2a(rho_2,autocorr)
		
		#sample sig2 from uniform
		rho_1 <- runif(sample_size, rho_bounds$rho_1[1], rho_bounds$rho_1[2]) #uniform on rho
		sig2 <- nabc.acf.rho2sig2(rho_1, a = a, vx=  variance)
		
	}

	ans <- data.frame(a=a,sig2=sig2)

	return(ans)	
	
}

is_within_bounds <- function(value,bounds){
	stopifnot(length(bounds)==2)
	return(value>=min(bounds) & value<=max(bounds))
}

nabc_MA1_is_within_prior_support <- function(a,sig2,a_bounds = c(-0.3, 0.3), sig2_bounds = c(0.5, 2), prior_dist=c("uniform","uniform_on_rho"),variance=NULL,autocorr=NULL){
	
	prior_dist <- match.arg(prior_dist)
	
	if(prior_dist=="uniform"){
		if(!is_within_bounds(a, a_bounds) || !is_within_bounds(sig2, sig2_bounds)){
			return(FALSE)
		}
		
	}
	
	if(prior_dist=="uniform_on_rho"){
		
		#estimate of nu must be provided to compute rho
		stopifnot(!is.null(variance),!is.null(autocorr))
		
		#compute bounds on rho1 and rho2 based on those of a and sig2
		rho_bounds <- nabc_MA1_compute_rho_bounds(a_bounds, sig2_bounds, variance, autocorr)

		#test if a is within the prior (equivalent to test rho_2)
		rho_2 <- nabc.acf.a2rho(a,autocorr)
		if(!is_within_bounds(rho_2,rho_bounds$rho_2)){
			return(FALSE)
		}
		
		#test if sig2 is within the prior (equivalent to test rho_1)
		rho_1 <- nabc.acf.sig22rho(sig2,a,variance)
		if(!is_within_bounds(rho_1,rho_bounds$rho_1)){
			return(FALSE)
		}
		
	}

	return(TRUE)		
	
}

nabc_MA1_dprior_a <- function(a,a_bounds= c(-0.3, 0.3),prior_dist=c("uniform","uniform_on_rho"),autocorr){
	
	prior_dist <- match.arg(prior_dist)
	
	if(prior_dist=="uniform"){
		a <- rep(1/abs(diff(a_bounds)),length(a))
	}
	
	if(prior_dist=="uniform_on_rho"){
		
		#compute bounds on rho2 based on those of a
		rho_2_bounds <- nabc.acf.a2rho(x = a_bounds, vx = autocorr)
		
		#test which a are within the prior (equivalent to test rho_2)
		rho_2 <- nabc.acf.a2rho(a,autocorr)
		ind <- which(is_within_bounds(rho_2, rho_2_bounds))
		
		x <- a[ind]
		dx <- (1-x^2)/(1+x^2+x^4)/abs(diff(rho_2_bounds))
		
		a[ind] <- dx
		a[-ind] <- 0
		
	}
	
	return(a)
	
}

nabc_MA1_dprior <- function(a,sig2,a_bounds = c(-0.3, 0.3), sig2_bounds = c(0.5, 2), prior_dist=c("uniform","uniform_on_rho"),variance=NULL,autocorr=NULL,give_log=FALSE,test_support=TRUE){
	
	prior_dist <- match.arg(prior_dist)
	stopifnot(length(a)==1,length(sig2)==1)
	
	if(test_support && !nabc_MA1_is_within_prior_support(a,sig2,a_bounds, sig2_bounds, prior_dist,variance,autocorr)){
		return(ifelse(give_log,log(0),0))
	}
	
	if(prior_dist=="uniform"){
		dprior <- 1/abs(diff(a_bounds)*diff(sig2_bounds))	
	}
	
	if(prior_dist=="uniform_on_rho"){
		
		#estimate of nu must be provided to compute rho
		stopifnot(!is.null(variance),!is.null(autocorr))

		#compute bounds on rho
		rho_bounds <- nabc_MA1_compute_rho_bounds(a_bounds, sig2_bounds, variance, autocorr)

		#compute dprior
		dprior <- (1-a^4)/(1+a^2+a^4)/abs(prod(sapply(rho_bounds,diff)))/variance
		
	}
	
	return(ifelse(give_log,log(dprior),dprior))

	
}



nabc_plot_density2d <- function(data=NULL,var_names=NULL,x_lab=NULL,y_lab=NULL, contour=TRUE, mode=FALSE, smoothing = c("none","ash", "kde"), ash_smooth = c(5, 5), ash_kopt = c(2, 2), kde_width_infl = 0.25, grid_size = NULL, plot=TRUE) {

	
	smoothing <- match.arg(smoothing)
	
	stopifnot(!is.null(data),length(var_names)==2,all(var_names%in%names(data)),smoothing!="none" | "density"%in%names(data))
	
	
	require(ggplot2)
	require(reshape2)
	require(RColorBrewer)
	require(plyr)

	
	#2D
	x <- data[[var_names[1]]]
	y <- data[[var_names[2]]]

	xlim <- range(x) 
	ylim <- range(y)

	if(is.null(grid_size)){
		grid_size=c(nclass.Sturges(x), nclass.Sturges(y))
	}
	
	if (smoothing == "ash") {
		
		require(ash)
		bins <- bin2(cbind(x, y), ab = rbind(xlim, ylim), nbin = grid_size)
		f <- ash2(bins, ash_smooth,ash_kopt)
		f <- rename(f, c(z = "density"))

	}

	if (smoothing == "kde") {
		require(KernSmooth)

		#fit kde
		x_bw <- kde_width_infl * diff(summary(x)[c(2, 5)])
		y_bw <- kde_width_infl * diff(summary(y)[c(2, 5)])
		f <- bkde2D(cbind(x, y), range.x = list(xlim, ylim), bandwidth = c(x_bw, y_bw), gridsize = grid_size)
		f <- rename(f, c(fhat = "density", x1 = "x", x2 = "y"))

	}

	if(smoothing == "none"){
		
		df_density2d <- data[c(var_names,"density")]
		names(df_density2d)[1:2] <- c("x","y")
		
	}else{

		df_density2d <- data.frame(x = rep(f$x, length(f$y)), y = rep(f$y, each = length(f$x)), density = as.vector(f$density))
		
	}
	
	
	if(mode){
		if(smoothing=="none"){
			df_summary <- subset(df_density2d,density==max(density))
			df_summary$summary <- "mode"
			
		}else{
			#get mode
			mxidx <- c((which.max(f$density) - 1)%%nrow(f$density) + 1, (which.max(f$density) - 1)%/%ncol(f$density) + 1) #row, col
			mx <- c(mean(f$x[c(mxidx[1], ifelse(mxidx[1] < length(f$x), mxidx[1] + 1, mxidx[1]))]), mean(f$y[c(mxidx[2], ifelse(mxidx[2] < 
				length(f$y), mxidx[2] + 1, mxidx[2]))]))	
			df_summary <- data.frame(x = mx[1], y = mx[2], summary = "mode")
			
		}


		
	}	
	
	#my plot	
	p <- ggplot(df_density2d, aes(x = x, y = y))
	p <- p + geom_tile(aes(fill = density), alpha = 0.85)
	if(contour){
		p <- p + geom_contour(aes(z = density, linetype = factor(..level..)), colour = "black")
		p <- p + scale_linetype("contour")	
	}
	if(mode){
		p <- p + geom_point(data = df_summary, aes(x = x, y = y, shape = summary))		
	}
	p <- p + scale_fill_gradientn("density", colours = rev(brewer.pal(11, "Spectral")))
	p <- p + guides(fill = guide_colourbar(order = 1), linetype = guide_legend(reverse = T, 
		keywidth = 2, order = 2))
	p <- p + xlab(ifelse(is.null(x_lab), var_names[1], x_lab)) + ylab(ifelse(is.null(y_lab), var_names[2], y_lab))
	
	if(plot){
		print(p)
	}else{
		return(p)
	}
}


nabc_MA1_plot_prior <- function(a_bounds = c(-0.3, 0.3), sig2_bounds = c(0.5, 2), prior_dist = c("uniform", "uniform_on_rho"), 
	variance = NULL, autocorr = NULL, method = c("analytic", "monte-carlo"), sample_size = 1e+05, smoothing = c("ash", 
		"kde"), ash_smooth = c(5, 5), ash_kopt = c(2, 2), kde_width_infl = 0.25, grid_size = c(100, 100)) {

	prior_dist <- match.arg(prior_dist)
	method <- match.arg(method)

	require(ggplot2)
	require(reshape2)
	require(RColorBrewer)
	require(plyr)

	if (prior_dist == "uniform_on_rho") {
		rho_1_bounds <- nabc.acf.sig22rho(sort(sig2_bounds), a = c(ifelse(prod(a_bounds) <= 0, 0, min(abs(a_bounds))), 
			max(abs(a_bounds))), vx = variance)
		rho_2_bounds <- nabc.acf.a2rho(x = a_bounds, vx = autocorr)
		new_sig2_bounds <- sort(rho_1_bounds) * variance/(1 + c(max(abs(a_bounds)), ifelse(prod(a_bounds) <= 0, 0, 
			min(abs(a_bounds))))^2)

	}

	if (prior_dist == "uniform") {
		new_sig2_bounds <- sig2_bounds
	}

	
	if (method == "monte-carlo") {

		df_prior <- nabc_MA1_rprior(sample_size, a_bounds, sig2_bounds, prior_dist, variance, autocorr)
		p <- nabc_plot_density2d(data = df_prior, var_names = c("a", "sig2"), y_lab = expression(sigma^2), mode = T, 
			smoothing = smoothing, ash_smooth = ash_smooth, ash_kopt = ash_kopt, kde_width_infl = kde_width_infl, grid_size = grid_size, 
			plot = F)


	}

	if (method == "analytic") {


		a <- seq(a_bounds[1], a_bounds[2], length = grid_size[1])
		sig2 <- seq(new_sig2_bounds[1], new_sig2_bounds[2], length = grid_size[2])

		if (prior_dist == "uniform") {
			
			d_a <- 1/abs(diff(a_bounds))
			d_sig2 <- 1/abs(diff(sig2_bounds))
			f <- list()
			f$a <- a
			f$sig2 <- sig2
			f$density <- rep(d_sig2 * d_a, length(a)*length(sig2))

		}

		if (prior_dist == "uniform_on_rho") {

			d_a <- (1 - a^2)/(1 + a^2 + a^4)/abs(diff(rho_2_bounds))
			d_sig2_a <- (1 + a^2)/(abs(diff(rho_1_bounds)) * variance)

			f <- list()
			f$a <- a
			f$sig2 <- sig2
			f$density <- rep(d_sig2_a * d_a, length(sig2))
		}

		#prior
		df_prior <- data.frame(a = rep(f$a, length(f$sig2)), sig2 = rep(f$sig2, each = length(f$a)), density = as.vector(f$density))

		if (prior_dist == "uniform_on_rho") {
			#limits on sig2 for given a
			df_prior <- ddply(df_prior, "a", function(df) {

				a <- df$a[1]
				ind <- which(df$sig2 < min(rho_1_bounds) * variance/(1 + a^2) | df$sig2 > max(rho_1_bounds) * variance/(1 + 
					a^2))
				df$density[ind] <- NA
				return(df)
			})
			df_prior <- na.omit(df_prior)

			#limits on sig2 for given a
			df_limits <- data.frame(a = f$a, ymin = min(rho_1_bounds) * variance/(1 + (f$a)^2), ymax = max(rho_1_bounds) * 
				variance/(1 + (f$a)^2))
			df_limits <- melt(df_limits, id.vars = "a", value.name = "sig2")

		}
		p <- nabc_plot_density2d(data = df_prior, var_names = c("a", "sig2"), y_lab = expression(sigma^2), mode =F, 
			smoothing = "none",contour=T, plot = F)
		if (prior_dist == "uniform_on_rho") {
			p <- p + geom_line(data = df_limits, aes(x = a, y = sig2, group = variable), colour = "white")
		}
	}


	print(p)

}


is.diag <- function(m) {
	if (!is.matrix(m)) {
		stop("is.diag: exception 1")
	}
	n <- ncol(m)
	if (n != nrow(m)) {
		return(FALSE)
	} else {
		return(all(m == (diag(rep(1, n)) * diag(m))))
	}
}


nabc_dratio_propwgausskernel <- function(support, ctheta, ptheta, covmat, give_log = F) {
	require(mvtnorm)
	
	#normalizing factor of N_{|[0,1]}(\mu,\sigma) is \Phi((1-\mu)/\sigma)-\Phi(-\mu/\sigma)
	#note that due to small variance some of the factors for each parameter will be 1
	#thus, since the support of many model parameters is restricted, the gaussian proposal kernel is not symmetric and we need to divide by the appropriate normalizing constants
	
	if (is.diag(covmat)) {
		#circumvent Cholesky decomposition if diagonal covariance matrix
		sigma <- sqrt(diag(covmat))
		names.sd.pos <- names(sigma)[sigma > 0]
		if (give_log) {

			r <- sum(log(pnorm(support["upper", names.sd.pos], mean = ctheta[names.sd.pos], sd = sigma[names.sd.pos]) - pnorm(support["lower", names.sd.pos], mean = ctheta[names.sd.pos], sd = sigma[names.sd.pos])) - log(pnorm(support["upper", names.sd.pos], mean = ptheta[names.sd.pos], sd = sigma[names.sd.pos]) - pnorm(support["lower", names.sd.pos], mean = ptheta[names.sd.pos], sd = sigma[names.sd.pos])))
		} else {
			#take ratios for each model parameter to help with numerical accuracy
			r <- prod((pnorm(support["upper", names.sd.pos], mean = ctheta[names.sd.pos], sd = sigma[names.sd.pos]) - pnorm(support["lower", names.sd.pos], mean = ctheta[names.sd.pos], sd = sigma[names.sd.pos]))/(pnorm(support["upper", names.sd.pos], mean = ptheta[names.sd.pos], sd = sigma[names.sd.pos]) - pnorm(support["lower", names.sd.pos], mean = ptheta[names.sd.pos], sd = sigma[names.sd.pos])))
		}
	} else {
		names.posdef <- rownames(covmat)[diag(covmat) > 0]

		if (give_log) {
			r <- log(pmvnorm(lower = support["lower", names.posdef], upper = support["upper", names.posdef], mean = ctheta[names.posdef],sigma = covmat[names.posdef, names.posdef])) - log(pmvnorm(lower = support["lower", names.posdef], upper = support["upper", names.posdef], mean = ptheta[names.posdef], sigma = covmat[names.posdef, names.posdef]))
		} else {
			r <- pmvnorm(lower = support["lower", names.posdef], upper = support["upper", names.posdef], mean = ctheta[names.posdef], sigma = covmat[names.posdef, names.posdef])/pmvnorm(lower = support["lower", names.posdef], upper = support["upper", names.posdef], mean = ptheta[names.posdef], sigma = covmat[names.posdef, names.posdef])
		}
	}
	return(r)
}

nabc_rpropwgausskernel<- function(theta, covmat, support)
{
	require(mvtnorm)
	ntrials<- 50
	itrials<- 1
	maxtrials<- 100

	#if(length(theta[["v"]])!=ncol(covmat))	stop("ABC.rpropwgausskernel: error at 1a")
	#if(!all(names(theta[["v"]])==colnames(covmat)))	stop("ABC.rpropwgausskernel: error at 1b")

	posdef_names<- names(theta)[diag(covmat)>0]
	prop<-	matrix(data= theta, ncol= length(theta), nrow= ntrials, dimnames= list(c(), names(theta)), byrow= TRUE)

	valid<-FALSE			#// the support of model parameters is often restricted, so draws from the normal distribution may not fall into that support. to speed up computation, make 'ntrials' attempts and choose the first successfull one, otherwise repeat

	while(all(valid==FALSE) && itrials<maxtrials)
	{
		if(length(posdef_names)==1){
			prop[,posdef_names]<- rnorm(ntrials, theta[posdef_names], sqrt(covmat[posdef_names,posdef_names]))
		}else{
			prop[,posdef_names]<- rmvnorm(ntrials, theta[posdef_names], covmat[posdef_names,posdef_names])			
		}		

		valid<- apply(prop,1, function(row){
			all(row[posdef_names]<=support["upper",posdef_names]	& row[posdef_names]>=support["lower",posdef_names])
		})
		itrials<- itrials+1
	}


	if(itrials>=maxtrials){	print(theta);	print(support); stop("ABC.rpropwgausskernel_AC: error at 1c");		}

	theta<- prop[which(valid)[1],]

	return(theta)
}

nabc_MA1_MLE <- function(x, variance_thin=0, autocorr_thin=0){

	n <- length(x)	
	x_thin <- x[seq.int(1,n,by=1+variance_thin)]
	n_thin <- length(x_thin)
	x_cor <- nabc.acf.equivalence.cor(x,leave=autocorr_thin)[["cor"]]
	x_var <- var(x_thin)*(n_thin-1)/n_thin

	a_MLE <- nabc.acf.nu2a(x_cor)
	sig2_MLE <- x_var/(1+ a_MLE^2)

	return(list(MLE=data.frame(a=a_MLE,sig2=sig2_MLE),s_stat=data.frame(variance= x_var,autocorr= x_cor)))
}

nabc_MA1_simulate <- function(n=1000,a=0.1,sig2=1,match_MLE=F,tol=c(a=1e-3,sig2=1e-3),variance_thin=0,autocorr_thin=0,plot=F){

	true_value <- c(a=a,sig2=sig2)

	if(!match_MLE){
		tol <- c(a=Inf,sig2=Inf)
	}

	i <- 0
	while(1){

		i <- i+1
		if(i%%1000==0){print(i)}

		eps <- rnorm(n+1,sd=sqrt(sig2))	
		eps_0 <- eps[1]
		x <- eps[-1] + a*eps[-(n+1)]

		#compute MLE	 unthinned
		unthinned <- nabc_MA1_MLE(x)

		#compute MLE thinned 
		thinned <- nabc_MA1_MLE(x, variance_thin, autocorr_thin)

		#compute errors
		error <- abs(unlist(c(unthinned$MLE-true_value, thinned$MLE-true_value, unthinned$MLE-thinned$MLE)))

		#break condition
		if(all(error<=tol)){break}
	}

	if(plot){
		df <- data.frame(t=1:n,x=x)
		p <- ggplot(data=df,aes(x=t,y=x))+geom_line()+scale_y_continuous(expression(x[t]))
		print(p)
	}


	return(list(param=data.frame(a=a,sig2=sig2),eps_0=eps_0,n=n,thin=data.frame(variance=variance_thin,autocorr= autocorr_thin), unthinned=unthinned, thinned=thinned, precision=list(tol=tol,error=matrix(error,ncol=2,byrow=T)),x=x))
}

check_MA1_simulator <- function(n_replicate = 10000, n_x = 5000, a = 0.1, sig2 = 1, dir_save = ".", RDS_file=NULL, dir_pdf= dir_save, grid_size=NULL, contour=FALSE,mode=FALSE) {

	if (!file.exists(dir_save)) {
		dir.create(dir_save)
	}

	if(is.null(RDS_file)){
		RDS_file <- file.path(dir_save, paste0("check_MA1_simulator_nx=", n_x, "_a=", a, "_sig2=", sig2, "_nrep=", n_replicate, 
			".rds"))
	}

	if (!file.exists(RDS_file)) {

		df_replicate <- ldply(1:n_replicate, function(i) {

			simu <- nabc_MA1_simulate(n = n_x, a = a, sig2 = sig2)
			return(simu$unthinned$MLE)

		}, .progress = "text")

		saveRDS(df_replicate, file = RDS_file)
	} else {
		df_replicate <- readRDS(file = RDS_file)
	}

	p <- nabc_plot_density2d(data = df_replicate, var_names = c("a", "sig2"), y_lab = expression(sigma^2), contour = contour, 
		mode = mode, smoothing = "ash", ash_smooth = c(5, 5), plot = F,grid_size=grid_size)
	p <- p + geom_hline(yintercept = sig2) + geom_vline(xintercept = a)

	pdf_file <- file.path(dir_pdf, paste0("hist2D_MLE_nx=", n_x, "_a=", a, "_sig2=", sig2, "_nrep=", n_replicate, ".pdf"))

	cairo_pdf(pdf_file, width = 7, height = 6)	
	print(p)
	dev.off()

	return(df_replicate)
}


nabc_MA1_conditional_loglikelihood <- function(a, sig2, x, eps_0) {

	#compute eps_hat_0:n
	eps_hat <- c(eps_0, x)
	for (i in 2:length(eps_hat)) {
		eps_hat[i] <- eps_hat[i] - a * eps_hat[i - 1]
	}
	eps_hat <- eps_hat[-1]
	loglike <- -(sum(eps_hat^2)/sig2 + length(x) * log(sig2))/2
	names(loglike) <- "loglike"

	return(loglike)
}


nabc_MA1_MCMC_MH <- function(data=NULL,theta_init=NULL,covmat_mvn_proposal=NULL,mcmc=NULL,a_true=0.1,sig2_true=1,eps_0_true=NULL,n_x=1000,a_bounds = c(-0.3, 0.3), sig2_bounds = c(0.5, 2), prior_dist=c("uniform","uniform_on_rho"),n_iter=1000,iter_adapt=n_iter,beta_adapt=0.5,plot=FALSE,dir_pdf=NULL,sample_from_prior= FALSE){

	prior_dist <- match.arg(prior_dist)

	if(!is.null(mcmc)){
		#go on mcmc
		data <- mcmc$data
		posterior <- mcmc$posterior
		theta_init <- posterior[[length(posterior)]]
		theta_init <- theta_init[setdiff(names(theta_init),"weight")]
		covmat_mvn_proposal <- mcmc$covmat_mvn_proposal
		a_bounds <- mcmc$bounds$a
		sig2_bounds <- mcmc$bounds$sig2
	}

	if(is.null(data)){
		#simulate data and compute variance and autocorr
		data <- nabc_MA1_simulate(n=n_x,a=a_true,sig2=sig2_true)	
	}

	if(sample_from_prior){
		#fix variance and autocorr to their true values
		a_true <- data$param$a 
		sig2_true <- data$param$sig2 
		data$unthinned$s_stat$variance <- (1 + a_true^2) * sig2_true
		data$unthinned$s_stat$autocorr <- a_true/(1 + a_true^2)

	}

	stopifnot(data$param$a>min(a_bounds),data$param$a<max(a_bounds), data$param$sig2>min(sig2_bounds), data$param$sig2<max(sig2_bounds))				

	if(plot){
		stopifnot(!is.null(dir_pdf))
		#plot prior
		pdf(file = file.path(dir_pdf, paste0("full_prior_analytic.pdf")), 6, 6)
		nabc_MA1_plot_prior(a_bounds, sig2_bounds, prior_dist, variance = data$unthinned$s_stat$variance, 
			autocorr = data$unthinned$s_stat$autocorr, method = "analytic",grid_size = c(100, 100))
		dev.off()
	}

	##choose init state
	if(is.null(theta_init)){
		#true value
		theta_init <- c(a=data$param$a,sig2=data$param$sig2,eps_0=data$eps_0)
		#sample from prior
		#theta_init <- c(unlist(nabc_MA1_rprior(1, a_bounds, sig2_bounds, prior_dist,data$s_stat$variance,data$s_stat$autocorr)),eps_0=data$eps_0)	
	}


	##specify support for parameters
	support <- matrix(c(-0.5, 0, -Inf, 0.5, Inf,Inf), ncol = length(theta_init), byrow = T, dimnames = list(c("lower","upper"), names(theta_init)))

	if(is.null(covmat_mvn_proposal)){
		##specify default covmat for gaussian proposal (eps_0 fixed to true value)
		covmat_mvn_proposal <- matrix(c(0.1, 0, 0, 0, 0.1, 0, 0, 0, 0), nrow = length(theta_init), byrow = T, dimnames = list(names(theta_init), names(theta_init)))	
	}


	#initial theta + loglike + dprior
	theta_curr <- theta_init
	#theta_curr <- nabc_rpropwgausskernel(theta_init, covmat_mvn_proposal, support)
	#theta_curr <- c(unlist(nabc_MA1_rprior(1, a_bounds, sig2_bounds, prior_dist,data$s_stat$variance,data$s_stat$autocorr)),eps_0=data$eps_0)	

	ll_curr <-  ifelse(sample_from_prior,0,nabc_MA1_conditional_loglikelihood(theta_curr["a"], theta_curr["sig2"],data$x,theta_curr["eps_0"]))
	log_dprior_curr <- nabc_MA1_dprior(a=theta_curr["a"],sig2=theta_curr["sig2"],a_bounds, sig2_bounds, prior_dist,variance=data$unthinned$s_stat$variance,autocorr=data$unthinned$s_stat$autocorr,give_log=T,test_support=FALSE)

	#run n_iterations
	progress_bar <- txtProgressBar(min = 1, max = n_iter, style = 3)

	if(is.null(mcmc)){
		posterior <- list(c(theta_curr,weight=0))		
	}

	time_start <- proc.time()[3]

	for (i_iter in seq_len(n_iter)) {

		if(i_iter > iter_adapt){
			df_posterior <- ldply(posterior)
			df_posterior <-as.data.frame(apply(df_posterior,2,rep,times=df_posterior$weight))
			df_posterior["weight"] <- NULL

			empirical_covmat <- cov(df_posterior)
			#param_fixed <- names(which(diag(covmat_mvn_proposal)==0))
			#empirical_covmat[param_fixed,] <- empirical_covmat[,param_fixed] <- 0
			covmat_mvn_proposal_adapted	<- (1-beta_adapt)*(2.38^2/length(sum(diag(covmat_mvn_proposal)>0)))*empirical_covmat + beta_adapt*covmat_mvn_proposal		
		}else{
			covmat_mvn_proposal_adapted <- covmat_mvn_proposal
		}

		#make a proposal + test if on prior's support
		theta_prop <- nabc_rpropwgausskernel(theta_curr, covmat_mvn_proposal_adapted, support)

		if(!nabc_MA1_is_within_prior_support(a=theta_prop["a"],sig2=theta_prop["sig2"],a_bounds, sig2_bounds, prior_dist,variance=data$unthinned$s_stat$variance,autocorr=data$unthinned$s_stat$autocorr)){
			#reject without computing likelihood
			posterior[[length(posterior)]]["weight"] <- posterior[[length(posterior)]]["weight"]+1 

		}else{
			#compute loglikelihood
			ll_prop <- ifelse(sample_from_prior,0,nabc_MA1_conditional_loglikelihood(theta_prop["a"], theta_prop["sig2"],data$x, theta_prop["eps_0"]))
			#compute dprior
			log_dprior_prop <- nabc_MA1_dprior(a= theta_prop["a"],sig2= theta_prop["sig2"],a_bounds, sig2_bounds, prior_dist,variance=data$unthinned$s_stat$variance,autocorr=data$unthinned$s_stat$autocorr,give_log=T,test_support=FALSE)

			#compute acceptance ratio:
			log_accept_ratio <- ll_prop + log_dprior_prop - ll_curr - log_dprior_curr + nabc_dratio_propwgausskernel(support, theta_curr, theta_prop, covmat_mvn_proposal_adapted, give_log = T)
			#log_accept_ratio <- ll_prop - ll_curr #- log_dprior_curr + nabc_dratio_propwgausskernel(support, theta_curr, theta_prop, covmat_mvn_proposal_adapted, give_log = T)

			if(log(runif(1))<log_accept_ratio){
				#accept
				theta_curr <- theta_prop
				ll_curr <- ll_prop
				log_dprior_curr <- log_dprior_prop
				posterior[[length(posterior)+1]] <- c(theta_curr,weight=1)

			}else{
				#reject
				posterior[[length(posterior)]]["weight"] <- posterior[[length(posterior)]]["weight"]+1 

			}

		}

		setTxtProgressBar(progress_bar, i_iter)

	}
	close(progress_bar)

	time_end <- proc.time()[3]
	cat("MCMC with",n_iter,"iterations lasted", round((time_end-time_start)/60,digits=2),"min\n")

	return(list(data=data,theta_init= theta_init,bounds=list(a=a_bounds,sig2=sig2_bounds),prior_dist= prior_dist,posterior= posterior,covmat_mvn_proposal= covmat_mvn_proposal,covmat_mvn_proposal_adapted= covmat_mvn_proposal_adapted))	

}

analyse_MCMC_MA1 <- function(mcmc,dir_pdf, smoothing=c("ash","kde"), ash_smooth=c(5,5),thin_every=0,burn=0,grid_size=NULL){

	require(coda)
	smoothing <- match.arg(smoothing)

	if(!file.exists(dir_pdf)){
		dir.create(dir_pdf,rec=T)
	}

	#true value and MLE	
	df_estimate <- rbind(mcmc$data$param,mcmc$data$unthinned$MLE,mcmc$data$thinned$MLE)
	df_estimate$type <- c("true value","MLE","MLE_thinned")

	cat("mvn_proposal_init:\n")
	print(mcmc$covmat_mvn_proposal)
	cat("\nmvn_proposal_adapted:\n")
	print(mcmc$covmat_mvn_proposal_adapted)

	df_posterior <- ldply(mcmc$posterior)		
	cat("\nmcmc acceptance rate:",100*round(nrow(df_posterior)/sum(df_posterior$weight),3),"%\n")
	df_posterior <- as.data.frame(apply(df_posterior, 2, rep, times = df_posterior$weight))
	df_posterior["weight"] <- NULL

	#remove fixed param
	ind <- names(which(sapply(df_posterior,function(x){length(unique(x))!=1})))
	df_posterior <- df_posterior[,ind,drop=F]

	#acf
	mc<-as.mcmc(df_posterior)
	cairo_pdf(file.path(dir_pdf,"autocorr.pdf"), width = 8, height = 6)	
	autocorr.plot(mc,lag.max=20)
	dev.off()

	#thin and burn
	x <- df_posterior
	x<-x[seq(1,nrow(x),thin_every),,drop=F]
	x<-x[seq(round(burn*nrow(x)),nrow(x),1),,drop=F]
	df_posterior <- x

	if(length(ind)>1){
		cat("\ncovariance matrix of thin and burned posterior:\n")
		print(cov(df_posterior[,ind,drop=F]))
	}

	#melt
	df_posterior$sample <- 1:nrow(df_posterior)
	gdf_posterior <- melt(df_posterior,id.vars="sample")	

	if(0){
		#trace
		p <- ggplot(gdf_posterior,aes(x=sample,y=value))+facet_wrap(~variable,scales="free_y")
		p <- p+geom_line()
		pdf(file = file.path(dir_pdf, "trace.pdf"), 6, 3)	
		print(p)
		dev.off()
	}

	#marginal posterior for a + prior
	gdf_a <- subset(gdf_posterior,variable=="a")
	p <- ggplot(gdf_a,aes(x=value))
	p <- p+geom_histogram(aes(y=..density..),binwidth=0.005,alpha=0.5)
	p <- p+geom_density()
	p <- p+stat_function(fun=nabc_MA1_dprior_a,args=list(a_bounds= mcmc$bounds$a,prior_dist=mcmc$prior_dist,autocorr= mcmc$data$unthinned$s_stat$autocorr),col="red")
	p <- p+geom_vline(data= df_estimate,aes(xintercept=a,linetype=type))
	pdf(file = file.path(dir_pdf, "marginal_posterior_a.pdf"), 4, 4)	
	print(p)
	dev.off()

	#plot posterior 2d
	pdf(file = file.path(dir_pdf, "full_posterior.pdf"), 6, 6)
	p <- nabc_plot_density2d(data = df_posterior, var_names = c("a", "sig2"), y_lab = expression(sigma^2), smoothing = smoothing, 
		ash_smooth = ash_smooth,plot=F,grid_size= grid_size)
	#add MLE and true value
	p <- p+geom_point(data=df_estimate,aes(x=a,y=sig2,shape=type))
	p <- p+scale_shape("estimates")
	print(p)	
	dev.off()


}	

check_MCMC_sampler <- function(a_true=0.1, sig2_true=1, a_bounds=c(-0.3, 0.3),sig2_bounds=c(0.9, 1.1),n_iter=100000,dir_pdf) {

	if(0){
		a_true <- 0.1
		sig2_true <- 1
		variance <- (1 + a_true^2) * sig2_true
		autocorr <- a_true/(1 + a_true^2)

		a_bounds <- c(-0.3, 0.3)
		sig2_bounds <- c(0.9, 1.1)
		n_iter <- 100000
	}

	smoothing <- "ash"
	ash_smooth <- c(5,1)

	mcmc <- nabc_MA1_MCMC_MH(a_true = a_true, sig2_true = sig2_true, n_x = 1, a_bounds = a_bounds, sig2_bounds = sig2_bounds, prior_dist = "uniform_on_rho", n_iter = n_iter, iter_adapt = n_iter/1, plot = T, dir_pdf=dir_pdf, sample_from_prior = T)
	saveRDS(mcmc,file=file.path(dir_pdf,"mcmc.rds"))

	analyse_MCMC_MA1(mcmc, dir_pdf, smoothing, ash_smooth)

	#compare with monte carlo estimate of the prior with same sample size
	pdf(file = file.path(dir_pdf, "full_prior_MC.pdf"), 6, 6)
	nabc_MA1_plot_prior(a_bounds = a_bounds, sig2_bounds = sig2_bounds, prior_dist = "uniform_on_rho", variance = variance, 
		autocorr = autocorr, method = "monte-carlo", sample_size = n_iter, smoothing = smoothing, ash_smooth = ash_smooth)
	dev.off()

}


run_MCMC_MA1 <- function(data=NULL,n_iter=1000,a_true=0.1,sig2_true=1,n_x=2000,a_bounds=c(-0.3, 0.3),sig2_bounds=c(0.5, 2),variance_thin=0,autocorr_thin=0){

	if(0){
		a_true <- 0.1
		sig2_true <- 1
		n_x <- 2000

		a_bounds <- c(-0.3, 0.3)
		sig2_bounds <- c(0.5, 2)

	}

	##
	dir_mcmc <- paste0("~/Documents/GitProjects/nABC/pdf/mcmc_MA1_a=",a_true,"_sig2=",sig2_true,"_nx=",n_x,"_nIter=",n_iter,"_thinVar=", variance_thin,"_thinCor=", autocorr_thin)
	dir.create(dir_mcmc)

	if(is.null(data)){
		#create data
		data <- nabc_MA1_simulate(n=n_x,a=a_true,sig2=sig2_true,match_MLE=T,tol=c(a=1e-3,sig2=1e-3),variance_thin= variance_thin,autocorr_thin= autocorr_thin)			
	}
	#theta_init
	theta_init <- c(a=data$param$a,sig2=data$param$sig2,eps_0=data$eps_0)		

	#default proposal
	covmat_mvn_proposal <- matrix(c(1e-3, 1e-5, 0, 1e-5, 1e-3, 0, 0, 0, 0), nrow = length(theta_init), byrow = T, dimnames = list(names(theta_init), names(theta_init)))	
	#covmat of a and sig2 based on the likelihood surface
	tmp <- check_MA1_simulator(n_replicate=10000,n_x=n_x,a=a_true,sig2=sig2_true,dir_save="~/Documents/GitProjects/nABC/pdf/covmat_MLE",dir_pdf= dir_mcmc,RDS_file=NULL,contour=T,grid_size=c(50,50))
	cov_a_sig2_MLE <- cov(tmp)
	covmat_mvn_proposal[rownames(cov_a_sig2_MLE),colnames(cov_a_sig2_MLE)] <- cov_a_sig2_MLE
	iter_adapt <- n_iter/10

	mcmc <- nabc_MA1_MCMC_MH(data= data, theta_init= theta_init, covmat_mvn_proposal= covmat_mvn_proposal, a_bounds = a_bounds, sig2_bounds = sig2_bounds, prior_dist = "uniform_on_rho", n_iter = n_iter, iter_adapt = iter_adapt, plot = T, dir_pdf= dir_mcmc)
	saveRDS(mcmc,file=file.path(dir_mcmc,"mcmc.rds"))

	#mcmc <- readRDS(file=file.path(dir_pdf,"mcmc.rds"))

	analyse_MCMC_MA1(mcmc, dir_mcmc,smoothing="ash",ash_smooth=c(5,5),thin_every=10,burn=0)


} 

continue_MCMC_MA1 <- function(mcmc,n_iter_more=1000){

	a_true <- mcmc$data$param$a
	sig2_true <- mcmc$data$param$sig2
	n_x <- mcmc$data$n
	df_posterior <- ldply(mcmc$posterior)
	n_iter <- sum(df_posterior$weight)+n_iter_more
	variance_thin <- mcmc$thin$variance
	autocorr_thin <- mcmc$thin$autocorr

	dir_pdf <- paste0("~/Documents/GitProjects/nABC/pdf/mcmc_MA1_a=",a_true,"_sig2=",sig2_true,"_nx=",n_x,"_nIter=",n_iter,"_thinVar=", variance_thin,"_thinCor=", autocorr_thin)
	dir.create(dir_pdf)

	mcmc <- nabc_MA1_MCMC_MH(mcmc= mcmc, prior_dist = "uniform_on_rho", n_iter = n_iter_more, iter_adapt = 0, beta_adapt=0.8, plot = T, dir_pdf=dir_pdf)
	saveRDS(mcmc,file=file.path(dir_pdf,"mcmc.rds"))

	#mcmc <- readRDS(file=file.path(dir_pdf,"mcmc.rds"))

	analyse_MCMC_MA1(mcmc, dir_pdf,smoothing="ash",ash_smooth=c(5,5),thin_every=10,burn=0)


}


run_parallel_MCMC_MA1 <- function(i_process, n_CPU, stream_names, data=NULL, n_iter=1000, iter_adapt= n_iter,a_true=0.1,sig2_true=1,n_x=2000,a_bounds=c(-0.3, 0.3),sig2_bounds=c(0.5, 2), prior_dist=c("uniform","uniform_on_rho"),variance_thin=0,autocorr_thin=0, dir_pdf){

	library(multicore)
	prior_dist <- match.arg(prior_dist)

	if(0){
		a_true <- 0.1
		sig2_true <- 1
		n_x <- 2000

		a_bounds <- c(-0.3, 0.3)
		sig2_bounds <- c(0.5, 2)

	}

	if(is.null(data)){
		#create data
		data <- nabc_MA1_simulate(n=n_x,a=a_true,sig2=sig2_true,match_MLE=T,tol=c(a=1e-2,sig2=1e-2),variance_thin= variance_thin,autocorr_thin= autocorr_thin)			
	}

	##
	dir_mcmc <- file.path(dir_pdf,paste0(i_process,"_mcmc_MA1_a=",data$param$a,"_sig2=",data$param$sig2,"_prior=",prior_dist,"_nx=",data$n,"_nIter=",n_iter,"_thinVar=", data$thin$variance,"_thinCor=", data$thin$autocorr,"_nChains=",n_CPU))
	dir.create(dir_mcmc,rec=T)


	#theta_init (the theta init of each chain will be a gaussian perturbation of the true value)
	theta_init <- c(a=data$param$a,sig2=data$param$sig2,eps_0=data$eps_0)		
	#theta_init <- c(a=0,sig2=0.9,eps_0=data$eps_0)		

	#default proposal
	covmat_mvn_proposal <- 50*matrix(c(1e-3, 1e-5, 0, 1e-5, 1e-3, 0, 0, 0, 0), nrow = length(theta_init), byrow = T, dimnames = list(names(theta_init), names(theta_init)))	
	#covmat of a and sig2 based on the likelihood surface
	#tmp <- check_MA1_simulator(n_replicate=20000,n_x=data$n,a=data$param$a,sig2=data$param$sig2,dir_save=file.path(dir_pdf,"covmat_MLE"),dir_pdf= dir_mcmc,RDS_file=NULL,contour=T,grid_size=c(50,50))
	#cov_a_sig2_MLE <- cov(tmp)
	#covmat_mvn_proposal[rownames(cov_a_sig2_MLE),colnames(cov_a_sig2_MLE)] <- cov_a_sig2_MLE


	if(prior_dist!="uniform"){
		#plot prior
		pdf(file = file.path(dir_mcmc, paste0("full_prior_analytic.pdf")), 6, 6)
		nabc_MA1_plot_prior(a_bounds, sig2_bounds, prior_dist= prior_dist, variance = data$unthinned$s_stat$variance, 
			autocorr = data$unthinned$s_stat$autocorr, method = "analytic",grid_size = c(100, 100))
		dev.off()		
	}


	all_chains<-mclapply(stream_names,FUN=run_foo_on_RNGstream, foo_name="nabc_MA1_MCMC_MH", data= data, theta_init= theta_init, covmat_mvn_proposal= covmat_mvn_proposal, a_bounds = a_bounds, sig2_bounds = sig2_bounds, prior_dist = prior_dist, n_iter = n_iter, iter_adapt = iter_adapt, plot = F)

	#combine

	saveRDS(all_chains,file=file.path(dir_mcmc,"all_chains.rds"))

	#combine posterior
	mcmc <- all_chains[[1]]	
	if(length(all_chains)>1){
		for(i in 2:length(all_chains)){	
			mcmc$posterior <- c(mcmc$posterior,all_chains[[i]]$posterior)
		}
	}
	saveRDS(mcmc,file=file.path(dir_mcmc,"mcmc_combined.rds"))

	#mcmc <- readRDS(file=file.path(dir_pdf,"mcmc_combined.rds"))

	analyse_MCMC_MA1(mcmc, dir_mcmc,smoothing="kde",ash_smooth=c(5,5),thin_every=20,burn=0,grid_size=c(100,100))


} 




run_foo_on_nCPU <- function(foo_name, n_CPU, use_cluster, ...) {

	require(rlecuyer)

	#read argument (Process between 0 and n_machines-1), add 1 to avoid 0, so that different machines have different seed
	if (use_cluster) {
		i_process <- as.numeric(Sys.getenv("ARG1")) + 1
	} else {
		i_process <- 1
	}

	cat("i_process=", i_process, "\n")
	#set seed multiplicator (time difference in second since 01-12-2012) so that simulations at different time can be combined (different parameter)
	seed_mult <- as.numeric(Sys.time() - ISOdate(2012, 12, 1)) * 24 * 3600
	cat("seed_mult=", seed_mult, "\n")
	cat("i_process*seed_mult=", i_process * seed_mult, "\n")

	#set seed to have different parameter set accross machines
	set.seed(i_process * seed_mult)

	#set seed to have different RNG stream accross CPU
	.lec.SetPackageSeed(round(runif(6, 1, 10000)))
	stream_names <- paste("rng_stream",1:n_CPU,sep="_")
	.lec.CreateStream(stream_names)

	time1 <- proc.time()[3]

	#start foo
	cat("Call", foo_name, "on", n_CPU, "CPUs with the following arguments:\n")
	arg_list <- c(i_process = i_process, n_CPU = n_CPU, stream_names = list(stream_names), list(...))
	print(arg_list)
	Sys.sleep(0.1)

	foo_res <- do.call(foo_name, arg_list)

	time2 <- proc.time()[3]
	cat("total simulation time: ", time2 - time1)

	return(foo_res)
}

run_foo_on_RNGstream <- function(stream_name, foo_name, ...){

	require(rlecuyer)

	.lec.CurrentStream(stream_name)

	foo_res <- do.call(foo_name, list(...))

	.lec.CurrentStreamEnd()

	return(foo_res)
}



main <- function() {

	require(ggplot2)
	require(reshape2)
	require(RColorBrewer)
	require(plyr)
	require(devtools)
	require(data.table)

	USE_CLUSTER <- T

	dev_mode()
	NABC_PKG <- ifelse(USE_CLUSTER,"/users/ecologie/camacho/GitProjects/abc.n/pkg","~/Documents/GitProjects/nABC/git_abc.n/pkg")
	load_all(NABC_PKG)

	#source Olli's prjct:
	source(file.path(NABC_PKG,"misc","nabc.prjcts.R"))

	dir_pdf <- ifelse(USE_CLUSTER,"/users/ecologie/camacho/nABC/MA1_a_0_to_0.3_tol=5e-3","~/Documents/GitProjects/nABC/pdf")
	dir.create(dir_pdf,rec=T)

	if(0){
		#plot prior
		a_0 <- 0
		sig2_0 <- 1
		variance <- (1 + a_0^2) * sig2_0
		autocorr <- a_0/(1 + a_0^2)

		pdf(file = file.path(dir_pdf, paste0("prior_induced_on_a=",a_0,"_sig2=",sig2_0,"_analytic.pdf")), 5, 6)
		nabc_MA1_plot_prior(a_bounds = c(-0.4, 0.4), sig2_bounds = c(0.8, 1.5), prior_dist = "uniform_on_rho", variance = variance, 
			autocorr = autocorr, method = "analytic", sample_size = 1e+05, smoothing = "ash", ash_smooth = c(5, 1), ash_kopt = c(2, 
				2), kde_width_infl = 0.25, grid_size = c(100, 100))
		dev.off()

		pdf(file = file.path(dir_pdf, paste0("prior_induced_on_a=",a_0,"_sig2=",sig2_0,"_numerical.pdf")), 6, 6)
		nabc_MA1_plot_prior(a_bounds = c(-0.3, 0.3), sig2_bounds = c(0.5, 2), prior_dist = "uniform_on_rho", variance = variance, 
			autocorr = autocorr, method = "monte-carlo", sample_size = 1000, smoothing = "ash", ash_smooth = c(5, 1), 
			ash_kopt = c(2, 2), kde_width_infl = 0.25, grid_size = c(100, 100))
		dev.off()

		#check sampler
		dir_pdf_check <- file.path(dir_pdf,"check_sampler")
		dir.create(dir_pdf_check)
		check_MCMC_sampler(a_true=0.1, sig2_true=1, a_bounds=c(-0.3, 0.3),sig2_bounds=c(0.9, 1.1),n_iter=10000,dir_pdf= dir_pdf_check) 
	}

	#check simulator of MA1
	#check_MA1_simulator(n_replicate=100000,n_x=300,a=0.1,sig2=1,dir_save=file.path(dir_pdf,"check_simulator_MA1"),contour=T,grid_size=c(50,50))

	#run sampler
	#run_MCMC_MA1(n_iter=1000,a_true=0.1,sig2_true=1,n_x=300,a_bounds=c(-0.3, 0.3),sig2_bounds=c(0.5, 2),leave.out.a=2,leave.out.s2=1)

	#continue
	#mcmc <- readRDS(file=file.path(dir_pdf,"mcmc_MA1_a=0.1_sig2=1_nx=300_nIter=21000_louta=0_louts2=0","mcmc.rds"))
	#continue_MCMC_MA1(mcmc,50000)

	#mcmc <- readRDS(file=file.path(dir_pdf,"mcmc_MA1_a=0.1_sig2=1_nx=300_nIter=10000_louta=2_louts2=1","mcmc.rds"))
	#continue_MCMC_MA1(mcmc,50000)

	n_x <- 150
	index_a	<- as.numeric(Sys.getenv("ARG1")) + 1
	a_true <- 0.1 # seq(0, 0.3, 0.025)[index_a]
	sig2_true <- 1
	tol <- 5e-3
	variance_thin <- 1
	autocorr_thin <- 2
	n_iter <- 2000000
	iter_adapt <- n_iter
	a_bounds <- c(-0.45, 0.45)
	sig2_bounds <- c(0.3, 1.7)
	prior_dist <- "uniform_on_rho"


	if(1){
		#foo n CPU
		file_data <- file.path(dir_pdf,paste0("data_with_a=",a_true,"_nx=",n_x,"_tol=", tol,"_varThin=", variance_thin,"_corThin=", autocorr_thin,".rds"))
		if(!file.exists(file_data)){
			data <- project.nABC.movingavg.get.fixed.ts(n=n_x, mu=0, a=a_true, sd= sqrt(sig2_true), leave.out.a= autocorr_thin, leave.out.s2= variance_thin, verbose=1, tol=tol, return_eps_0=T)
			#data <- nabc_MA1_simulate(n=n_x,a=a_true,sig2=sig2_true,match_MLE=T,tol=c(a= a_tol,sig2= sig2_tol),variance_thin=1,autocorr_thin= 2)				
			saveRDS(data,file= file_data)
		}else{
			data <- readRDS(file= file_data)
		}
	}


	#parallel
	run_foo_on_nCPU(foo_name="run_parallel_MCMC_MA1", n_CPU=ifelse(USE_CLUSTER,1,2), use_cluster= USE_CLUSTER, data=data, n_iter= n_iter, iter_adapt= iter_adapt,a_bounds= a_bounds,sig2_bounds= sig2_bounds,prior_dist= prior_dist, dir_pdf=dir_pdf) 

	if(0){

		dir_mcmc <- file.path(dir_pdf,"1_mcmc_MA1_a=0.1_sig2=1_prior=uniform_nx=150_nIter=5e+05_thinVar=1_thinCor=2_nChains=12")
		mcmc <- readRDS(file=file.path(dir_mcmc,"mcmc_combined.rds"))
		analyse_MCMC_MA1(mcmc, dir_pdf=dir_mcmc,smoothing="kde",ash_smooth=c(5,5),thin_every=20,burn=0,grid_size=c(50,50))
		#			analyse_MCMC_MA1(mcmc, dir_pdf=dir_mcmc,smoothing="ash",ash_smooth=c(5,5),thin_every=20,burn=0,grid_size=c(50,50))

	}

	if(0){
		dir_mcmc1 <- file.path(dir_pdf,"1_mcmc_MA1_a=0.1_sig2=1_nx=300_nIter=50000_thinVar=1_thinCor=2_nChains=10")
		mcmc1 <- readRDS(file=file.path(dir_mcmc1,"mcmc.rds"))
		dir_mcmc2 <- file.path(dir_pdf,"1_mcmc_MA1_a=0.1_sig2=1_nx=300_nIter=1e+05_thinVar=1_thinCor=2_nChains=12")
		mcmc <- readRDS(file=file.path(dir_mcmc2,"mcmc_combined.rds"))
		mcmc$posterior <- c(mcmc$posterior,mcmc1$posterior)
		saveRDS(mcmc,file=file.path(dir_mcmc2,"mcmc_combined_all.rds"))
		analyse_MCMC_MA1(mcmc, dir_pdf=file.path(dir_mcmc2,"all_combined"),smoothing="ash",ash_smooth=c(5,5),thin_every=10,burn=0)
		analyse_MCMC_MA1(mcmc, dir_pdf=file.path(dir_mcmc2,"big_grid"),smoothing="ash",ash_smooth=c(5,5),thin_every=10,burn=0,grid_size=c(100,100))
		analyse_MCMC_MA1(mcmc, dir_pdf=file.path(dir_mcmc2,"no_trim"),smoothing="ash",ash_smooth=c(5,5),thin_every=1,burn=0,grid_size=c(100,100))
	}


	if(0){
		dir_mcmc <- file.path(dir_pdf,"1_mcmc_MA1_a=0.1_sig2=1_prior=uniform_nx=150_nIter=1e+05_thinVar=1_thinCor=2_nChains=12 7")
		#combine mcmc
		all_chains <- readRDS(file.path(dir_mcmc,"all_chains.rds"))

		#combine posterior
		mcmc <- all_chains[[1]]
		n_iter <- length(mcmc$posterior)
		n_burn <- 0.1*n_iter
		mcmc$posterior <- mcmc$posterior[-(1:n_burn)]
		for(i in 2:length(all_chains)){	
			mcmc$posterior <- c(mcmc$posterior,all_chains[[i]]$posterior[-(1:n_burn)])
		}
		saveRDS(mcmc,file=file.path(dir_mcmc,"mcmc_burned.rds"))
		#
		analyse_MCMC_MA1(mcmc, dir_pdf=dir_mcmc,smoothing="ash",ash_smooth=c(5,5),thin_every=10,burn=0)
		analyse_MCMC_MA1(mcmc, dir_pdf=file.path(dir_mcmc,"no_trim_big_grid"),smoothing="ash",ash_smooth=c(5,5),thin_every=1,burn=0,grid_size=c(100,100))
	}
}

main()





