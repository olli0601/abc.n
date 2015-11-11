##--------------------------------------------------------------------------------------------------------
##
##
##	stuff to run on HPC
##
##
##--------------------------------------------------------------------------------------------------------
#PR.STARTME					<- '/Users/Oliver/git/HPTN071sim/source/rPANGEAHIVsim/misc/rPANGEAHIV.startme.R'
PR.STARTME					<- '/work/or105/libs/abc.star/misc/nabc.startme.R'
PR.VARIOUS					<- paste(PR.STARTME," -exe=VARIOUS",sep='')
HPC.CX1.IMPERIAL			<- "cx1.hpc.ic.ac.uk"		#this is set to system('domainname',intern=T) for the hpc cluster of choice
HPC.CX1.IMPERIAL.LOAD		<- "module load intel-suite R/3.2.0"
##--------------------------------------------------------------------------------------------------------
cmd.hpcsys<- function()
{
	tmp<- system('domainname',intern=T)
	if(!nchar(tmp))	tmp<- "debug"
	tmp
}
##--------------------------------------------------------------------------------------------------------
cmd.hpcwrapper<- function(cmd, hpcsys= cmd.hpcsys(), hpc.walltime=24, hpc.mem="1750mb", hpc.nproc=1, hpc.q=NA)
{
	wrap<- "#!/bin/sh"
	#hpcsys<- HPC.CX1.IMPERIAL
	if(hpcsys%in%c(HPC.CX1.IMPERIAL,'(none)'))
	{				
		tmp	<- paste("#PBS -l walltime=",hpc.walltime,":59:59,pcput=",hpc.walltime,":45:00",sep='')
		wrap<- paste(wrap, tmp, sep='\n')		
		tmp	<- paste("#PBS -l select=1:ncpus=",hpc.nproc,":mem=",hpc.mem,sep='')
		wrap<- paste(wrap, tmp, sep='\n')
		wrap<- paste(wrap, "#PBS -j oe", sep='\n')
		if(!is.na(hpc.q))
			wrap<- paste(wrap, paste("#PBS -q",hpc.q), sep='\n\n')
		wrap<- paste(wrap, HPC.CX1.IMPERIAL.LOAD, sep='\n')
	}
	else if(hpcsys=='debug')
		cat(paste("\ndetected no HPC system and no hpcwrapper generated, domain name is",hpcsys))
	else
		stop(paste("unknown hpc system with domain name",hpcsys))
	
	cmd<- lapply(seq_along(cmd),function(i){	paste(wrap,cmd[[i]],sep='\n')	})
	if(length(cmd)==1)
		cmd<- unlist(cmd)
	cmd	
}
##--------------------------------------------------------------------------------------------------------
cmd.hpccaller<- function(outdir, outfile, cmd)
{
	if( nchar( Sys.which("qsub") ) )
	{
		file	<- paste(outdir,'/',gsub(':','',outfile),'.qsub',sep='')
		cat(paste("\nwrite HPC script to",file,"\n"))
		cat(cmd,file=file)
		cmd		<- paste("qsub",file)
		cat( cmd )
		cat( system(cmd, intern=TRUE) )
		Sys.sleep(1)
	}
	else
	{
		file	<- paste(outdir,'/',gsub(':','',outfile),'.sh',sep='')
		cat(paste("\nwrite Shell script to\n",file,"\nStart this shell file manually\n"))
		cat(cmd,file=file)
		Sys.chmod(file, mode = "777")	
		Sys.sleep(1)
	}
	file
}
##--------------------------------------------------------------------------------------------------------
cmd.various<- function(prog= PR.VARIOUS)
{
	cmd		<- "#######################################################
# start: run VARIOUS
#######################################################"
	cmd		<- paste(cmd, '\n', prog, '\n', sep='')
	cmd		<- paste(cmd,"#######################################################
# end: run VARIOUS
#######################################################\n",sep='')
	cmd
}
##--------------------------------------------------------------------------------------------------------
gof.pipeline<- function()
{
	outdir		<- getwd()
	cmd			<- cmd.various()
	#cmd			<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeelab', hpc.walltime=71, hpc.mem="5000mb")
	cmd			<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q=NA, hpc.walltime=3, hpc.mem="1890mb")
	cat(cmd)			
	outfile		<- paste("gof",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
	cmd.hpccaller(outdir, outfile, cmd)
	quit("no")	
}
##--------------------------------------------------------------------------------------------------------
##
##
##	ABC and ABCSTAR pre-simulations
##
##
##--------------------------------------------------------------------------------------------------------
abc.presim.uprior.mu<- function(abc.nit, xn, xmean, xsigma, prior.l, prior.u, ysigma, yn=NA )		
{		
	ans					<- vector("list",5)
	names(ans)			<- c("x","xn","xmean","xsigma","sim")
	obs 				<- rnorm(xn, xmean, xsigma)
	obs 				<- (obs - mean(obs))/sd(obs) * xsigma + xmean
	ans[["x"]]			<- obs
	ans[["xn"]]			<- xn
	ans[["xmean"]]		<- xmean
	ans[["xsigma"]]		<- xsigma
	
	ans[["sim"]]		<- sapply(1:abc.nit, function(i)
			{					
				ymu		<- runif(1, prior.l, prior.u)
				y		<- rnorm(yn, ymu, sd=ysigma)
				tmp		<- c(yn, ymu, ysigma, mean(y), sd(y), quantile(y, prob=0.25), quantile(y, prob=0.75), max(y) )									
				tmp					
			})					
	rownames(ans[["sim"]])	<- toupper(c('ym','ymu','ysigma','ysmean','yssd','ysq25','ysq75','ysmx'))
	ans[["sim"]]		<- as.data.table(t(ans$sim))
	ans[["sim"]][, IT:=seq_len(nrow(ans[["sim"]]))]
	ans
}
##--------------------------------------------------------------------------------------------------------
abcstar.presim.uprior.mu<- function(abc.nit, xn, xmean, xsigma, prior.l, prior.u )		
{		
	SIG		<<- xsigma								#sig assumed known
	n2s		<- function(n){ SIG/sqrt(floor(n)) }	#need formula to convert n.of.y into s.of.T, depends on application
	s2n		<- function(s){ (SIG/s)^2 }				#need formula to convert s.of.T into n.of.y, depends on application
	ztest.calibrate(n.of.x=xn, n2s=n2s, s2n=s2n, mx.pw=0.9, alpha=0.01, what='KL', plot=TRUE)
	
	ans					<- vector("list",5)
	names(ans)			<- c("x","xn","xmean","xsigma","sim")
	obs 				<- rnorm(xn, xmean, xsigma)
	obs 				<- (obs - mean(obs))/sd(obs) * xsigma + xmean
	ans[["x"]]			<- obs
	ans[["xn"]]			<- xn
	ans[["xmean"]]		<- xmean
	ans[["xsigma"]]		<- xsigma
	
	ans[["sim"]]		<- sapply(1:abc.nit, function(i)
			{					
				ymu		<- runif(1, prior.l, prior.u)
				
				y		<- rnorm(yn, ymu, sd=ysigma)
				tmp		<- c(yn, ymu, ysigma, mean(y), sd(y), quantile(y, prob=0.25), quantile(y, prob=0.75), max(y) )									
				tmp					
			})					
	rownames(ans[["sim"]])	<- toupper(c('ym','ymu','ysigma','ysmean','yssd','ysq25','ysq75','ysmx'))
	ans[["sim"]]		<- as.data.table(t(ans$sim))
	ans[["sim"]][, IT:=seq_len(nrow(ans[["sim"]]))]
	ans
}


##--------------------------------------------------------------------------------------------------------
gof.mutostabc.presim.mu<- function(outdir, outfile, rep=10)
{
	for(i in seq_len(rep))
	{
		dt	<- abc.presim.uprior.mu( 	abc.nit=1e7, xn=60, xmean=1.34, xsigma=1.4, 
										prior.l=1.34-5, prior.u=1.34+5, ysigma=1.4, yn=60 )
		file<- paste(outdir, '/', gsub('\\.rda',paste('_R',i,'.rda',sep=''), outfile), sep='')
		cat('save to', file)
		save(dt, file=file)		
		dt	<- NULL
		gc()
	}	
}
##--------------------------------------------------------------------------------------------------------
gof.mutostabc.presim.mu2<- function(outdir, outfile, rep=10)
{
	for(i in seq_len(rep))
	{
		dt	<- abcstar.presim.uprior.mu(abc.nit=1e7, xn=60, xmean=1.34, xsigma=1.4, prior.l=1.34-5, prior.u=1.34+5 )		
		file<- paste(outdir, '/', gsub('\\.rda',paste('_R',i,'.rda',sep=''), outfile), sep='')
		cat('save to', file)
		save(dt, file=file)		
		dt	<- NULL
		gc()
	}	
}
##--------------------------------------------------------------------------------------------------------
gof.mutostabc.main<- function()
{
	require(data.table)
	require(abc.star)
	#outdir	<- '~/Dropbox (Infectious Disease)/gof-abc/calc/example-paper'
	outdir	<- paste(HOME, '/gof', sep='')
	if(1)
	{
		outfile	<- 'Normal-ME-OR151111.rda'
		gof.mutostabc.presim.mu(outdir, outfile, rep=10)
	}
	if(0)
	{
		library(devtools)
		code.dir	<- "/Users/Oliver/git/abc.star"
		devtools::install(code.dir)
		
		require(abc.star)
		
		
		tmp			<- abc.df14[,  as.list( mutost.calibrate(	n.of.y=m, s.of.y=yssd, tau.u.ub=1, what='MXPW', mx.pw=0.9, alpha=0.01)[1:4] ), by='it']
		abc.df14		<- merge(abc.df14, tmp, by='it')								
		
		
		
		
		abc.presim14<- abc.presim.uprior.mu( 	abc.nit=1e7, xn=60, xmean=1.34, xsigma=1.4, 
				prior.l=1.34-5, prior.u=1.34+5, ysigma=1.4, yn=60 )
		abc.df14	<- as.data.table(t(abc.presim14$sim))								
		abc.df14[, it:=seq_len(nrow(abc.df14))]								
		tmp			<- abc.df14[,  as.list( mutost.calibrate(	n.of.y=m, s.of.y=yssd, tau.u.ub=1, what='MXPW', mx.pw=0.9, alpha=0.01)[1:4] ), by='it']
		abc.df14		<- merge(abc.df14, tmp, by='it')								
		save(file="/Users/Oliver/duke/2014_ABCBookChapter/data/abc.presim.mutost.sigma14.Rdata", abc.df14)
	}
	
	
}