
nabc_MA_sample_from_prior <- function(sample_size, a_bounds = c(-0.3, 0.3), sig2_bounds = c(0.5, 2), nu_1=NULL,nu_2=NULL, constraint=c("uniform_prior","uniform_prior_on_rho")){
	
	constraint <- match.arg(constraint)
	stopifnot(a_bounds[2]>= a_bounds[1],sig2_bounds[2]>= sig2_bounds[1])
	
	
	if(constraint=="uniform_prior_on_rho"){
		
		#estimate of nu must be provided to compute rho
		stopifnot(!is.null(nu_1),!is.null(nu_2))
		
		#compute bounds on rho1 and rho2 based on those of a and sig2
		rho_2_bounds <- nabc.acf.a2rho(min(a_prior))
		max_rho2 <- nabc.acf.a2rho(max(a_prior))	

		
		#first sample a from	
		a <- 
		
	}
	
	
}

nabc_MA_induced_prior_numerical <- function(a_prior = c(-0.3, 0.3), sig2_prior = c(0.5, 2), sample_size = 1e+05, method = c("ash", 
	"kde"), ash_smooth = c(5, 5), ash_kopt = c(2, 2), width_infl = 0.25, grid_size = c(100, 100)) {

	method <- match.arg(method)
	#		method <- method[1]
	require(ggplot2)
	require(reshape2)
	require(RColorBrewer)
	require(plyr)

	sample_size <- 1e6	
	min_rho2 <- nabc.acf.a2rho(min(a_prior))
	max_rho2 <- nabc.acf.a2rho(max(a_prior))	
	xrho2 <- runif(sample_size, min_rho2, max_rho2 ) #uniform on rho
	xa <- nabc.acf.rho2a(xrho2)
	
	#x <- seq(a_prior[1], a_prior[2],length=1000)
	#y <- (1-x^2)/(1+x^2+x^4)/(max_rho2-min_rho2)
	
	hist(xa,freq=FALSE)
	lines(x,y,col="red")
	abline(v=c(min_rho2,max_rho2))


	min_rho1 <- nabc.acf.sig22rho(min(sig2_prior),a=0)
	max_rho1 <- nabc.acf.sig22rho(max(sig2_prior),a=max(a_prior))
	xrho1 <- runif(sample_size, min_rho1, max_rho1) #uniform on rho
	xsig2 <- nabc.acf.rho2sig2(xrho1, a = xa)
	
	#2D
	x <- xa
	y <- xsig2
	xlim <- range(x) #* 1.05
	ylim <- range(y) #* 1.05

	if (method == "ash") {
		require(ash)

		bins <- bin2(cbind(x, y), ab = rbind(xlim, ylim), nbin = c(nclass.Sturges(x), nclass.Sturges(y)))
		f <- ash2(bins, ash_smooth,ash_kopt)

	}

	if (method == "kde") {
		require(KernSmooth)

		#fit kde
		x_bw <- width_infl * diff(summary(x)[c(2, 5)])
		y_bw <- width_infl * diff(summary(y)[c(2, 5)])
		f <- bkde2D(cbind(x, y), range.x = list(xlim, ylim), bandwidth = c(x_bw, y_bw), gridsize = grid_size)
		f <- rename(f, c(fhat = "z", x1 = "x", x2 = "y"))

	}

	
	#prior
	df_prior <- data.frame(x = rep(f$x, length(f$y)), y = rep(f$y, each = length(f$x)), density = as.vector(f$z))
	#limits on sig2
	df_prior <- ddply(df_prior,"x",function(df){
		
		x <- df$x[1]
		ind <- which(df$y<=min_rho1/(1+x^2) | df$y>=max_rho1/(1+x^2))
		df$density[ind] <- NA
		return(df)
	})
	df_prior <- na.omit(df_prior)
	#limits on sig2 for given a
	df_limits <-  data.frame(x=f$x,ymin=min_rho1/(1+(f$x)^2),ymax=max_rho1/(1+(f$x)^2))
	df_limits <- melt(df_limits,id.vars="x",value.name="y")
	#get mode
	mxidx <- c((which.max(f$z) - 1)%%nrow(f$z) + 1, (which.max(f$z) - 1)%/%ncol(f$z) + 1) #row, col
	mx <- c(mean(f$x[c(mxidx[1], ifelse(mxidx[1] < length(f$x), mxidx[1] + 1, mxidx[1]))]), mean(f$y[c(mxidx[2], ifelse(mxidx[2] < 
		length(f$y), mxidx[2] + 1, mxidx[2]))]))


	df_summary <- data.frame(x = mx[1], y = mx[2], summary = "mode")

	#my plot	
	p <- ggplot(df_prior, aes(x = x, y = y))
	p <- p + geom_tile(aes(fill = density), alpha = 0.85)
	p <- p + geom_contour(aes(z = density, linetype = factor(..level..)), colour = "black")
	p <- p + geom_point(data = df_summary, aes(x = x, y = y, shape = summary))
	p <- p + geom_line(data = df_limits, aes(group=variable), colour="white")
	p <- p + scale_fill_gradientn("density", colours = rev(brewer.pal(11, "Spectral")))
	p <- p + scale_linetype("contour") + guides(fill = guide_colourbar(order = 1), linetype = guide_legend(reverse = T, 
		keywidth = 2, order = 2))
	p <- p + xlab("a") + ylab(expression(sigma^2))
	print(p)

}

main <- function(){
	dev_mode()
	NABC_PKG <- "~/Documents/GitProjects/nABC/git_abc.n/pkg/"
	load_all(NABC_PKG)
	
	dir_pdf <- "~/Documents/GitProjects/nABC/pdf/"
	dir.create(dir_pdf)
	
	pdf(file=file.path(dir_pdf,"prior_induced_on_a_sig2_large.pdf"),5,6)
	induced_prior_MA(a_prior=c(-1,1),sig2_prior=c(0.1,3),sample_size=10000000,method="ash",ash_smooth=c(5,1),ash_kopt=c(2,2),grid_size=c(100,100))	
	dev.off()

	pdf(file=file.path(dir_pdf,"prior_induced_on_a_sig2_zoom.pdf"),5,6)
	induced_prior_MA(a_prior=c(-0.3,0.3),sig2_prior=c(0.5,2),sample_size=1000000,method="ash",ash_smooth=c(5,1),ash_kopt=c(2,2),grid_size=c(100,100))	
	dev.off()


}

#main()





