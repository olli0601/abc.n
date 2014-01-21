\name{corrz.calibrate}
\alias{corrz.calibrate}
\title{Calibrate the power function of \code{corrz}}
\usage{
corrz.calibrate(n.of.x, n.of.y = n.of.x, n2s = function(n) {     1/sqrt(n -
  3) }, s2n = function(s) {     (1/s)^2 + 3 }, mx.pw = 0.9, alpha = 0.01,
  max.it = 100, pow_scale = 2, debug = F, plot = F)
}
\arguments{
  \item{n.of.x}{Number of observed summary values. Each
  summary data point is a pair of two values.}

  \item{n.of.y}{Number of simulated summary values. Each
  summary data point is a pair of two values.}

  \item{n2s}{Function to transform the number of data
  points into the standard deviation of the test
  statistic.}

  \item{s2n}{Function to transform the standard deviation
  of the test statistic into the number of data points.}

  \item{plot}{Logical. If \code{plot==TRUE}, the calibrated
  power function is plotted along with the summary
  likelihood.}

  \item{debug}{Logical. If \code{debug==TRUE}, detailed
  optimization output for the number of simulated summary
  values is printed to the console.}

  \item{max.it}{this algorithm stops prematurely when the
  number of iterations to calibrate the number of simulated
  data points exceeds 'max.it'}

  \item{alpha}{Level of the equivalence test.}

  \item{mx.pw}{maximum power at the point of reference
  (rho.star).}

  \item{pow_scale}{Used to set the support of the pdf
  associated to the power function. The power is truncated
  between \code{[tau.l/pow_scale,tau.u*pow_scale]} and then
  standardized.}
}
\value{
vector of length 7 \item{n.of.y}{number of simulated
summary values} \item{tau.l}{lower tolerance of the
equivalence region} \item{tau.u}{upper tolerance of the
equivalence region} \item{pw.cmax}{actual power at
rho.star} \item{KL_div}{actual KL divergence between the
ABC approximation and the summary likelihood}
\item{cl}{lower point of the critical region, i.e. lower
standard ABC tolerance} \item{cu}{upper point of the
critical region, i.e. upper standard ABC tolerance}
}
\description{
This function calibrates the power function of \code{corrz}
so that its mode coincides with the mode of the summary
likelihood and so that its KL divergence to the summary
likelihood is minimized. The function minimizes the KL
divergences and includes recursive calls to re-calibrate
the upper and lower tolerance regions for every new
proposed number of simulated summary values.
}
\examples{
#
#	illustrate calibration of tau.l, tau.u and m 
#
n.of.x		<- 60
x			<- rnorm(n.of.x,0,1)
s.of.x		<- sd(x)
abc.param	<- chisqstretch.calibrate(n.of.x, s.of.x, scale=n.of.x, n.of.y=n.of.x, mx.pw=0.9, alpha=0.01, max.it=100, debug=FALSE, plot=TRUE)
abc.param
#
#	ABC* rejection sampler for the normal variance example 
#	
xn			<- yn	<- 60
df			<- yn-1
ymu			<- 0 
xmu			<- 0
xsigma2		<- 1
prior.u		<- 4
prior.l		<- 0.2
#	NOTE set this to 1e6; only low to pass R CMD check
N			<- 1e3
#	function to precompute ABC* simulations
chi2stretch.simu<- function(N, prior.l, prior.u, x, yn, ymu)		
{		
	if(prior.u<1)	stop("project.nABC.StretchedChi2.fix.x.uprior.ysig2: error at 1c")
	if(prior.l>1)	stop("project.nABC.StretchedChi2.fix.x.uprior.ysig2: error at 1d")
	ans						<- vector("list",3)
	names(ans)				<- c("x","xsigma2","data")
	ans[["x"]]				<- x
	ans[["xsigma2"]]		<- (length(x)-1)*var(x)/length(x)			#MLE of sig2
	ans[["data"]]			<- sapply(1:N,function(i)
			{					
				ysigma2		<- runif(1,prior.l,prior.u)
				y			<- rnorm(yn,ymu,sd=sqrt(ysigma2))
				tmp			<- c( ysigma2, (var(y)*(yn-1))/(var(x)*(length(x)-1)), log( var(y)/var(x) ), var(y)-var(x) )									
				tmp					
			})								
	rownames(ans[["data"]])	<- c("ysigma2","T","log.sy2/sx2","sy2-sx2")
	ans
}
#	function to produce histogram
chi2stretch.hist<- function(x, theta, nbreaks= 20, breaks= NULL, width= 0.5, plot=0, rtn.dens=0,...)
{
	EPS	<- 1e-12
	#compute break points sth theta is in the middle
	if(is.null(breaks))
	{
		breaks<- max(abs( theta - x ))*1.1 / nbreaks								
		breaks<- c( rev(seq(from= theta-breaks/2, by=-breaks, length.out= nbreaks )), seq(from= theta+breaks/2, by=breaks, length.out= nbreaks ) )		
	}
	ans.h<- hist(x, breaks=breaks, plot= 0)
	ans.h[["mean"]]<- mean(x)		
	ans.h[["hmode"]]<- mean(ans.h[["breaks"]][seq(which.max(ans.h[["counts"]]),length.out=2)])	
	tmp<- density(x, kernel="biweight",from=breaks[1],to=breaks[length(breaks)],width = max(EPS,width*diff(summary(x)[c(2,5)])))
	ans.h[["dmode"]]<- tmp[["x"]][which.max( tmp[["y"]])]
	if(rtn.dens)
		ans.h[["dens"]]<- tmp
	if(plot)
	{
		plot(ans.h, freq=0,...)
		lines(tmp)
	}
	ans.h
}

#	pseudo data
x				<- rnorm(xn,xmu,sd=sqrt(xsigma2))
x 				<- (x - mean(x))/sd(x) * sqrt(xsigma2) + xmu
#	calibrated ABC* parameters
abc.param		<- chisqstretch.calibrate(length(x), sd(x), mx.pw=0.9, alpha=0.01, max.it=100, debug=FALSE, plot=FALSE)
#	get ABC* simulations
simu			<- chi2stretch.simu(N, prior.l, prior.u, x, abc.param['n.of.y'], ymu)
#	evaluate ABC* simulations
tstat			<- simu[["data"]]["T",] 
acc				<- which( tstat>=abc.param["cl"]  &  tstat<=abc.param["cu"] )
chi2stretch.hist(simu[["data"]]["ysigma2",acc], simu[["xsigma2"]], nbreaks= 50, width= 0.5, plot=0, ylim=c(0,2.25))
}
\seealso{
\code{\link{corrz.pow}}, \code{\link{corrz.sulkl}}
}

