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
	if(1)
	{
		outdir		<- getwd()
		cmd			<- cmd.various()
		#cmd			<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeelab', hpc.walltime=71, hpc.mem="5000mb")
		cmd			<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeelab', hpc.walltime=200, hpc.mem="5000mb")
		#cmd			<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeelab', hpc.walltime=200, hpc.mem="17000mb")
		cat(cmd)			
		outfile		<- paste("gof",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
		cmd.hpccaller(outdir, outfile, cmd)
		quit("no")	
	}		
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
abc.presim.uprior.musig<- function(abc.nit, xn, xmean, xsigma, mu.prior.l, mu.prior.u, sig.prior.l, sig.prior.u )		
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
				ymu		<- runif(1, mu.prior.l, mu.prior.u)
				ysigma	<- exp(runif(1, log(sig.prior.l), log(sig.prior.u)))
				y		<- rnorm(xn, ymu, sd=ysigma)
				tmp		<- c(xn, ymu, ysigma, mean(y), sd(y), quantile(y, prob=0.25), quantile(y, prob=0.75), max(y) )									
				tmp					
			})					
	rownames(ans[["sim"]])	<- toupper(c('ym','ymu','ysigma','ysmean','yssd','ysq25','ysq75','ysmx'))
	ans[["sim"]]		<- as.data.table(t(ans$sim))
	ans[["sim"]][, IT:=seq_len(nrow(ans[["sim"]]))]
	ans
}
##--------------------------------------------------------------------------------------------------------
abcstar.presim.uprior.musig<- function(abc.nit, xn, xmean, xsigma, mu.prior.l, mu.prior.u, sig.prior.l, sig.prior.u )		
{		
	ans					<- vector("list",5)
	names(ans)			<- c("x","xn","xmean","xsigma","sim")
	obs 				<- rnorm(xn, xmean, xsigma)
	obs 				<- (obs - mean(obs))/sd(obs) * xsigma + xmean
	ans[["x"]]			<- obs
	ans[["xn"]]			<- xn
	ans[["xmean"]]		<- xmean
	ans[["xsigma"]]		<- xsigma
	#	upper limit for simulations
	yn					<- ceiling( 1.5*mutost.calibrate(n.of.x=length(obs), s.of.x= sd(obs), s.of.y=sig.prior.u, what='KL', mx.pw=0.9, alpha=0.01, plot=FALSE)['n.of.y'] )				
	#	always the same
	sig.cali			<- vartest.calibrate(n.of.x=xn, s.of.x=sd(obs), what='KL', mx.pw=0.9, alpha=0.01, plot=FALSE, verbose=FALSE)		
	
	ans[["sim"]]		<- sapply(1:abc.nit, function(i)
			{					
				ymu		<- runif(1, mu.prior.l, mu.prior.u)
				ysigma	<- exp(runif(1, log(sig.prior.l), log(sig.prior.u)))
				sim		<- rnorm(yn, ymu, sd=ysigma)			
				mu.cali <- mutost.calibrate(n.of.x=length(obs), s.of.x= sd(obs), s.of.y=sd(sim), what='KL', mx.pw=0.9, alpha=0.01, plot=FALSE)
				stopifnot(mu.cali['n.of.y']<=yn)
				tmp		<- unname(c(mu.cali['n.of.y'], sig.cali['n.of.y'], ymu, ysigma, mean(sim[seq_len(mu.cali['n.of.y'])]), sd(sim[seq_len(mu.cali['n.of.y'])]), sd(sim[seq_len(sig.cali['n.of.y'])]), quantile(sim, prob=0.25), quantile(sim, prob=0.75), max(sim[seq_len(xn)]) ))									
				tmp					
			})					
	rownames(ans[["sim"]])	<- toupper(c('ym.mu','ym.sigma','ymu','ysigma','ysmean','yssd.mu','yssd.sig','ysq25','ysq75','ysmx'))
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
	tmp		<- ztest.calibrate(n.of.x=xn, n2s=n2s, s2n=s2n, mx.pw=0.9, alpha=0.01, what='KL', plot=FALSE)
	yn		<- tmp['n.of.y']
	ysigma	<- xsigma
	
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
				tmp		<- c(yn, ymu, ysigma, mean(y), sd(y), quantile(y, prob=0.25), quantile(y, prob=0.75), max(y[seq_len(xn)]) )									
				tmp					
			})					
	rownames(ans[["sim"]])	<- toupper(c('ym','ymu','ysigma','ysmean','yssd','ysq25','ysq75','ysmx'))
	ans[["sim"]]		<- as.data.table(t(ans$sim))
	ans[["sim"]][, IT:=seq_len(nrow(ans[["sim"]]))]
	ans
}
##--------------------------------------------------------------------------------------------------------
gof.mutostabc.presim.mu<- function(outdir, outfile, n.rep=10)
{
	for(i in seq_len(n.rep))
	{
		dt	<- abc.presim.uprior.mu( 	abc.nit=1e6, xn=60, xmean=1.34, xsigma=1.4, 
										prior.l=1.34-5, prior.u=1.34+5, ysigma=1.4, yn=60 )
		file<- paste(outdir, '/', gsub('\\.rda',paste('_R',i,'.rda',sep=''), outfile), sep='')
		cat('save to', file)
		save(dt, file=file)		
		dt	<- NULL
		gc()
	}	
}
##--------------------------------------------------------------------------------------------------------
gof.mutostabc.presim.musig<- function(outdir, outfile, n.rep=10)
{
	for(i in seq_len(n.rep))
	{		
		dt	<- abc.presim.uprior.musig(abc.nit=1e6, xn=60, xmean=1.34, xsigma=1.4, mu.prior.l=1.34-2, mu.prior.u=1.34+2, sig.prior.l=1/3, sig.prior.u=3)
		file<- paste(outdir, '/', gsub('\\.rda',paste('_R',i,'.rda',sep=''), outfile), sep='')
		cat('save to', file)
		save(dt, file=file)		
		dt	<- NULL
		gc()
	}	
}
##--------------------------------------------------------------------------------------------------------
gof.mutostabc.presim.musig.ABCstar<- function(outdir, outfile, n.rep=10)
{
	for(i in seq_len(n.rep))
	{		
		dt	<- abcstar.presim.uprior.musig(abc.nit=1e6, xn=60, xmean=1.34, xsigma=1.4, mu.prior.l=1.34-2, mu.prior.u=1.34+2, sig.prior.l=1/3, sig.prior.u=3)
		file<- paste(outdir, '/', gsub('\\.rda',paste('_R',i,'.rda',sep=''), outfile), sep='')
		cat('save to', file)
		save(dt, file=file)		
		dt	<- NULL
		gc()
	}	
}
##--------------------------------------------------------------------------------------------------------
gof.mutostabc.MX.mu.evalcpp<- function(indir='~/Dropbox (Infectious Disease)/gof-abc/calc/example-paper')
{
	#	get p-val
	infiles	<- data.table(FILE=list.files(indir, pattern='CPP\\.rda$'))
	set(infiles, NULL, 'TYPE', infiles[, gsub('-OR.*','',FILE)])
	load( paste(indir, subset(infiles, TYPE=='Normal-ME')[, FILE], sep='/') )
	cpps.std	<- copy(cpps)
	cpps.std[, TYPE:='standard ABC']
	load( paste(indir, subset(infiles, TYPE=='Normal-ME-MforZTEST')[, FILE], sep='/') )
	cpps[, TYPE:='calibrated ABC']
	cpps		<- rbind(cpps, cpps.std)
	set(cpps, NULL, 'TOLM', cpps[, factor(as.character(TOLM), levels=sort(as.numeric(as.character(unique(TOLM)))), labels=sort(as.numeric(as.character(unique(TOLM)))))])
	set(cpps, NULL, 'TYPE', cpps[, factor(TYPE, levels=c('standard ABC','calibrated ABC'), labels=c('standard ABC','calibrated ABC'))])
	
	ggplot(cpps, aes(x=CPP)) + 
			geom_histogram(aes(fill=TOLM), colour='black', position="identity", binwidth=0.05) +
			coord_cartesian(xlim=c(0,1)) +
			labs(title='ABC model correct\n', x='\nconditional predictive p-value', fill='tolerance multiplier\nrelative to\ncalibrated value') +
			facet_grid(TOLM~TYPE) + theme_bw() + theme(panel.margin=unit(2, "lines"))
	ggsave(file=paste(indir,'/','Normal-ME_cpp_pdf.pdf',sep=''), w=7, h=6)
	
	ggplot(cpps, aes(x=CPP)) + stat_ecdf(aes(colour=TOLM)) +
			geom_abline(intercept=0, slope=1, colour='black') +
			coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
			facet_grid(~TYPE) +
			labs(y='empirical c. d. f.\n', title='ABC model correct\n', x='\nconditional predictive p-value', colour='tolerance multiplier\nrelative to\ncalibrated value') +
			theme_bw() + theme(panel.margin=unit(2, "lines"))
	ggsave(file=paste(indir,'/','Normal-ME_cpp_cdf.pdf',sep=''), w=7, h=4)
	
}
##--------------------------------------------------------------------------------------------------------
gof.mutostabc.MX.musig.evalcpp<- function(indir='~/Dropbox (Infectious Disease)/gof-abc/calc/example-paper')
{
	#	get p-val
	infiles	<- data.table(FILE=list.files(indir, pattern='CPP\\.rda$'))
	set(infiles, NULL, 'TYPE', infiles[, gsub('-OR.*','',FILE)])
	set(infiles, NULL, 'INSUFF', infiles[, grepl('insuff',FILE)])
	
	load( paste(indir, subset(infiles, TYPE=='Normal-MESIG' & !INSUFF)[, FILE], sep='/') )
	cpps.std	<- copy(cpps)
	cpps.std[, TYPE:='standard ABC']
	cpps.std[, SUS:='summaries sufficient']
	load( paste(indir, subset(infiles, TYPE=='Normal-MESIG' & INSUFF)[, FILE], sep='/') )
	cpps[, TYPE:='standard ABC']
	cpps[, SUS:='summaries not sufficient']
	cpps		<- rbind(cpps, cpps.std)
	set(cpps, NULL, 'TOLMU', cpps[, factor(substring(regmatches(TOL, regexpr('^[^;]*',TOL)), 4))])
	cpps[, TOLSIG:='None']
	tmp			<- cpps[, which(SUS=='summaries sufficient')]
	set(cpps, tmp, 'TOLSIG', cpps[tmp, substring(regmatches(TOL, regexpr(';.*',TOL)), 7)])	
	set(cpps, NULL, 'TOLSIG', cpps[, factor(TOLSIG, levels=cpps[, unique(TOLSIG)], labels=cpps[, unique(TOLSIG)])] )
	set(cpps, NULL, 'TOLLEG', cpps[, paste('tolerances\n',gsub(':','= ',gsub('; ','\n',TOL)),sep='')])
	#set(cpps, NULL, 'TOL', cpps[, factor(as.character(TOLM), levels=sort(as.numeric(as.character(unique(TOLM)))), labels=sort(as.numeric(as.character(unique(TOLM)))))])
	#set(cpps, NULL, 'TYPE', cpps[, factor(TYPE, levels=c('standard ABC','calibrated ABC'), labels=c('standard ABC','calibrated ABC'))])
	
	ggplot(cpps, aes(x=CPP)) + 
			geom_histogram(aes(fill=TOLMU), colour='black', position="identity", binwidth=0.05) +
			coord_cartesian(xlim=c(0,1)) +
			scale_fill_brewer(palette='Set1', guide=FALSE) +
			labs(title='ABC model correct\n', x='\nconditional predictive p-value', fill='tolerances') +
			facet_grid(TOLLEG~SUS) + theme_bw() + theme(panel.margin=unit(2, "lines"))
	ggsave(file=paste(indir,'/','Normal-MESIG_cpp_pdf_insuff_vs_suff.pdf',sep=''), w=8, h=13)
	
	ggplot(cpps, aes(x=CPP)) + stat_ecdf(aes(colour=TOLMU, linetype=TOLSIG)) +
			geom_abline(intercept=0, slope=1, colour='black') +
			coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
			scale_colour_brewer(palette='Set1') +
			facet_grid(~SUS) +
			labs(y='empirical c. d. f.\n', title='ABC model correct\n', x='\nconditional predictive p-value', colour='tolerances\nfor comparing means', linetype='tolerances\nfor comparing std devs') +
			theme_bw() + theme(panel.margin=unit(2, "lines"))
	ggsave(file=paste(indir,'/','Normal-MESIG_cpp_cdf_insuff_vs_suff.pdf',sep=''), w=7, h=4)
	
	load( paste(indir, subset(infiles, TYPE=='Normal-MESIG')[, FILE], sep='/') )
	cpps[, TYPE:='standard ABC']
	#cpps		<- rbind(cpps, cpps.std)
	#set(cpps, NULL, 'TOL', cpps[, factor(as.character(TOL), levels=sort(as.numeric(as.character(unique(TOL)))), labels=sort(as.numeric(as.character(unique(TOL)))))])
	#set(cpps, NULL, 'TYPE', cpps[, factor(TYPE, levels=c('standard ABC','calibrated ABC'), labels=c('standard ABC','calibrated ABC'))])
	ggplot(cpps, aes(x=CPP)) + stat_ecdf(aes(colour=TOL)) +
			geom_abline(intercept=0, slope=1, colour='black') +
			coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
			#facet_grid(~TYPE) +
			labs(y='empirical c. d. f.\n', title='ABC model correct\n', x='\nconditional predictive p-value', colour='tolerances') +
			theme_bw() + theme(panel.margin=unit(2, "lines"))
	ggsave(file=paste(indir,'/','Normal-MESIG_cpp_cdf.pdf',sep=''), w=5, h=4)
}
##--------------------------------------------------------------------------------------------------------
gof.mutostabc.MX.mu<- function(indir='~/Dropbox (Infectious Disease)/gof-abc/calc/example-paper')
{
	infiles	<- data.table(FILE=list.files(indir, pattern='rda$'))
	set(infiles, NULL, 'TYPE', infiles[, gsub('-OR.*','',FILE)])
	set(infiles, NULL, 'REP', infiles[, as.integer(substring(regmatches(FILE, regexpr('_R[0-9]+',FILE)),3))])
	infiles	<- subset(infiles, REP>0L & TYPE=='Normal-ME')
	cpps	<- infiles[, {
				file	<- paste(indir, FILE,sep= '/')
				cat('\n',file)
				load(file)
				dtc		<- copy(dt$sim)
				#	exact posterior density
				de		<- data.table(YMU= seq(min(dtc[, YMU]), max(dtc[, YMU]), len=2048) )
				de[, DENS:= de[, dnorm(YMU, dt$xmean, dt$xsigma/sqrt(dt$xn))]]
				#	get Hyp test stat H so thresholds will be comparable to ABCSTAR version
				dtc[, HSTAT:= dtc[, (dt$xmean-YSMEAN)/YSIGMA*sqrt(YM)]]
				#	upper lower tolerance will be +-1.645, set ABC tolerances around that
				tols	<- 1.645*c(0.1, 0.5, 1, 2, 4, 6)
				#	get accepted iterations for each tolerance
				abca	<- do.call('rbind',lapply(tols, function(tol)
								{
									tmp	<- subset(dtc, abs(HSTAT)<tol)
									tmp[, TOL:=tol]
									tmp
								}))
				abca[, TOLM:= factor(TOL/1.645)]
				#	plot abc posterior density
				ggplot(abca, aes(x=YMU)) + 
						geom_density(aes(colour=TOLM, group=TOLM)) + 
						geom_line(data=de, aes(y=DENS), colour='black') +
						theme_bw() + theme(legend.position='bottom') + labs(title='ABC posterior density\n', x='mu', y='', colour='tolerance multiplier\nrelative to\nABC* calibrated tolerance')
				ggsave(file=gsub('\\.rda','_ABCposterior.pdf',file), w=6, h=5)
				#	acceptance prob
				acc.prob	<- abca[, list(PERC_ACC= length(YMU)/nrow(dt$sim)), by='TOLM']
				#	CPP
				cpp			<- abca[, list( CPP= mean(YSMX>=max(dt$x)) ),by='TOLM']
				ans			<- merge(cpp, acc.prob, by='TOLM')
				ans
			},by='FILE']
	cpps	<- merge(cpps, infiles, by='FILE')
	file	<- paste(indir, infiles[1, gsub('_R[0-9]+\\.rda','_CPP.rda',FILE)], sep='/')
	save(cpps, file=file)	
}
##--------------------------------------------------------------------------------------------------------
gof.mutostabc.MX.musig.insuff<- function(indir='~/Dropbox (Infectious Disease)/gof-abc/calc/example-paper')
{
	infiles	<- data.table(FILE=list.files(indir, pattern='rda$'))
	infiles	<- subset(infiles, !grepl('CPP', FILE))
	set(infiles, NULL, 'TYPE', infiles[, gsub('-OR.*','',FILE)])
	set(infiles, NULL, 'REP', infiles[, as.integer(substring(regmatches(FILE, regexpr('_R[0-9]+',FILE)),3))])
	infiles	<- subset(infiles, REP>0L & TYPE=='Normal-MESIG')
	cpps	<- infiles[, {
				file	<- paste(indir, FILE,sep= '/')
				cat('\n',file)
				load(file)
				dt$sim	<- copy(dt$sim)
				#	exact posterior density
				scale	<- dt$xsigma*dt$xsigma*(dt$xn-1)
				de		<- data.table(	YMU= rt(1e6, dt$xn-1)*(dt$xsigma/sqrt(dt$xn))+dt$xmean, 
						YSIG2=rigamma(1e6, (dt$xn-2)/2, dt$xsigma*dt$xsigma*(dt$xn-1)/2)	
				)					
				#	get Hyp test stat H so thresholds will be comparable to ABCSTAR version
				set(dt$sim, NULL, 'HSTAT_MU', dt$sim[, (dt$xmean-YSMEAN)/YSSD*sqrt(YM)])
				set(dt$sim, NULL, 'HSTAT_SIG', dt$sim[, (YSSD*YSSD)*(YM-1)/(dt$xsigma*dt$xsigma*(dt$xn-1))])
				#	upper lower tolerance will be +-1.645, set ABC tolerances around that
				#	vartest.calibrate(n.of.x=dt$xn, s.of.x=sd(dt$x), what='KL', mx.pw=0.9, alpha=0.01, plot=TRUE, verbose=FALSE)
				mu.tols		<- 1.645*c(0.25, 0.5, 1, 2, 4)
				#s.tols		<- c(1.3, 1.5, 1.7, 2.3, 2.9)				
				tols		<- data.table(T_MU=mu.tols)
				abca		<- tols[,	{
											subset(dt$sim, abs(HSTAT_MU)<T_MU)				
										}, by=c('T_MU')]
				abca[, TOL:= abca[, factor(paste('mu:',round(T_MU,d=2),sep=''))]]
				#	plot abc posterior density
				ggplot(abca, aes(x=YMU)) + 			 
						geom_vline(xintercept= dt$xmean, colour='grey80', size=1) +
						geom_hline(yintercept= dt$xsigma*dt$xsigma, colour='grey80', size=1) +
						#geom_density2d(data=de, aes(y=YSIG2), colour='black') +
						geom_density2d(aes(y=YSIGMA*YSIGMA, colour=TOL, group=TOL)) +
						facet_wrap(~TOL, ncol=3) +
						theme_bw() + theme(legend.position='bottom') + labs(title='ABC posterior density\n', x='mu', y='sigma2', colour='tolerances')			
				ggsave(file=gsub('\\.rda','_insuff_ABCposterior.pdf',file), w=10, h=8)
				#	acceptance prob
				acc.prob	<- abca[, list(PERC_ACC= length(YMU)/nrow(dt$sim)), by='TOL']
				#	CPP
				cpp			<- abca[, list( CPP= mean(YSMX>=max(dt$x)) ),by='TOL']
				ans			<- merge(cpp, acc.prob, by='TOL')
				ans				
			},by='FILE']
	cpps	<- merge(cpps, infiles, by='FILE')
	file	<- paste(indir, infiles[1, gsub('_R[0-9]+\\.rda','_insuff_CPP.rda',FILE)], sep='/')
	save(cpps, file=file)		
}
##--------------------------------------------------------------------------------------------------------
gof.mutostabc.MX.musig<- function(indir='~/Dropbox (Infectious Disease)/gof-abc/calc/example-paper')
{
	infiles	<- data.table(FILE=list.files(indir, pattern='rda$'))
	infiles	<- subset(infiles, !grepl('CPP', FILE))
	set(infiles, NULL, 'TYPE', infiles[, gsub('-OR.*','',FILE)])
	set(infiles, NULL, 'REP', infiles[, as.integer(substring(regmatches(FILE, regexpr('_R[0-9]+',FILE)),3))])
	infiles	<- subset(infiles, REP>0L & TYPE=='Normal-MESIG')
	cpps	<- infiles[, {
				file	<- paste(indir, FILE,sep= '/')
				cat('\n',file)
				load(file)
				dt$sim	<- copy(dt$sim)
				#	exact posterior density
				scale	<- dt$xsigma*dt$xsigma*(dt$xn-1)
				de		<- data.table(	YMU= rt(1e6, dt$xn-1)*(dt$xsigma/sqrt(dt$xn))+dt$xmean, 
										YSIG2=rigamma(1e6, (dt$xn-2)/2, dt$xsigma*dt$xsigma*(dt$xn-1)/2)	
										)					
				#	get Hyp test stat H so thresholds will be comparable to ABCSTAR version
				set(dt$sim, NULL, 'HSTAT_MU', dt$sim[, (dt$xmean-YSMEAN)/YSSD*sqrt(YM)])
				set(dt$sim, NULL, 'HSTAT_SIG', dt$sim[, (YSSD*YSSD)*(YM-1)/(dt$xsigma*dt$xsigma*(dt$xn-1))])
				#	upper lower tolerance will be +-1.645, set ABC tolerances around that
				#	vartest.calibrate(n.of.x=dt$xn, s.of.x=sd(dt$x), what='KL', mx.pw=0.9, alpha=0.01, plot=TRUE, verbose=FALSE)
				mu.tols		<- 1.645*c(0.25, 0.5, 1, 2, 4)
				s.tols		<- c(1.3, 1.5, 1.7, 2.3, 2.9)
				#tols		<- as.data.table(expand.grid(T_MU=mu.tols, T_SIG=s.tols))
				tols		<- data.table(T_MU=mu.tols, T_SIG=s.tols)
				abca		<- tols[,	{
							subset(dt$sim, abs(HSTAT_MU)<T_MU & HSTAT_SIG<T_SIG & HSTAT_SIG>1/T_SIG)				
						}, by=c('T_MU','T_SIG')]
				abca[, TOL:= abca[, factor(paste('mu:',round(T_MU,d=2),'; sig:',round(T_SIG,d=2),sep=''))]]
				#	plot abc posterior density
				ggplot(abca, aes(x=YMU)) + 			 
						geom_vline(xintercept= dt$xmean, colour='grey80', size=1) +
						geom_hline(yintercept= dt$xsigma*dt$xsigma, colour='grey80', size=1) +
						geom_density2d(data=de, aes(y=YSIG2), colour='black') +
						geom_density2d(aes(y=YSIGMA*YSIGMA, colour=TOL, group=TOL)) +
						facet_wrap(~TOL, ncol=3) +
						theme_bw() + theme(legend.position='bottom') + labs(title='ABC posterior density\n', x='mu', y='sigma2', colour='tolerances')			
				ggsave(file=gsub('\\.rda','_ABCposterior.pdf',file), w=10, h=8)
				#	acceptance prob
				acc.prob	<- abca[, list(PERC_ACC= length(YMU)/nrow(dt$sim)), by='TOL']
				#	CPP
				cpp			<- abca[, list( CPP= mean(YSMX>=max(dt$x)) ),by='TOL']
				ans			<- merge(cpp, acc.prob, by='TOL')
				abca	<- de	<- df<- NULL
				gc()
				ans				
			},by='FILE']
	cpps	<- merge(cpps, infiles, by='FILE')
	file	<- paste(indir, infiles[1, gsub('_R[0-9]+\\.rda','_CPP.rda',FILE)], sep='/')
	save(cpps, file=file)		
}
##--------------------------------------------------------------------------------------------------------
gof.mutostabc.MX.mu.ABCstar<- function(indir='~/Dropbox (Infectious Disease)/gof-abc/calc/example-paper')
{
	infiles	<- data.table(FILE=list.files(indir, pattern='rda$'))
	set(infiles, NULL, 'TYPE', infiles[, gsub('-OR.*','',FILE)])
	infiles	<- subset(infiles, !grepl('CPP', FILE))
	set(infiles, NULL, 'REP', infiles[, as.integer(substring(regmatches(FILE, regexpr('_R[0-9]+',FILE)),3))])
	infiles	<- subset(infiles, REP>0L & TYPE=='Normal-ME-MforZTEST')
	cpps	<- infiles[, {
				file	<- paste(indir, FILE,sep= '/')
				cat('\n',file)
				load(file)
				dtc		<- copy(dt$sim)
				#	exact posterior density
				de		<- data.table(YMU= seq(min(dtc[, YMU]), max(dtc[, YMU]), len=2048) )
				de[, DENS:= de[, dnorm(YMU, dt$xmean, dt$xsigma/sqrt(dt$xn))]]				
				#	get Hyp test stat H so thresholds will be comparable to ABCSTAR version
				set(dtc, NULL, 'HSTAT', dtc[, (dt$xmean-YSMEAN)/YSIGMA*sqrt(YM)])
				#	calibrated tolerance
				SIG		<<- dt$xsigma							#sig assumed known
				n2s		<- function(n){ SIG/sqrt(floor(n)) }	#need formula to convert n.of.y into s.of.T, depends on application
				s2n		<- function(s){ (SIG/s)^2 }				#need formula to convert s.of.T into n.of.y, depends on application
				tmp		<- ztest.calibrate(dt$xn, n2s=n2s, s2n=s2n, mx.pw=0.9, alpha=0.01, what='KL', plot=FALSE)
				#	get accepted iterations for each tolerance
				abca	<- subset(dtc, abs(HSTAT)<tmp['c.u'])
				abca[, TOL:=tmp['c.u']]
				abca[, TOLM:= 1L]
				#	plot abc posterior density
				ggplot(abca, aes(x=YMU)) + 
						geom_density(aes(fill=TOLM, colour=TOLM, group=TOLM), alpha=0.6) + 
						geom_line(data=de, aes(y=DENS), colour='black') +
						scale_fill_continuous(guide = FALSE) + scale_colour_continuous(guide = FALSE) +
						theme_bw() + theme(legend.position='bottom') + labs(title='calibrated ABC posterior density\n', x='mu', y='')
				ggsave(file=gsub('\\.rda','_ABCStarposterior.pdf',file), w=6, h=5)
				#	acceptance prob
				acc.prob	<- abca[, list(PERC_ACC= length(YMU)/nrow(dt$sim)), by='TOLM']
				#	CPP
				cpp			<- abca[, list( CPP= mean(YSMX>=max(dt$x)) ),by='TOLM']
				ans			<- merge(cpp, acc.prob, by='TOLM')
				ans
			},by='FILE']
	cpps	<- merge(cpps, infiles, by='FILE')
	file	<- paste(indir, infiles[1, gsub('_R[0-9]+\\.rda','_CPP.rda',FILE)], sep='/')
	save(cpps, file=file)	
}
##--------------------------------------------------------------------------------------------------------
gof.mutostabc.presim.mu.ABCstar<- function(outdir, outfile, n.rep=10)
{
	for(i in seq_len(n.rep))
	{
		dt	<- abcstar.presim.uprior.mu(abc.nit=1e6, xn=60, xmean=1.34, xsigma=1.4, prior.l=1.34-5, prior.u=1.34+5 )		
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
	require(pscl)
	require(grid)
	require(ggplot2)	
	require(abc.star)
	if(1)
	{
		#outdir	<- '~/Dropbox (Infectious Disease)/gof-abc/calc/example-paper'
		outdir	<- paste(HOME, '/data/gof', sep='')		
		#outfile	<- 'Normal-ME-OR151111.rda'
		#gof.mutostabc.presim.mu(outdir, outfile, n.rep=200)
		outfile	<- 'Normal-ME-MforZTEST-OR151111.rda'
		gof.mutostabc.presim.mu.ABCstar(outdir, outfile, n.rep=200)
		#outfile	<- 'Normal-MESIG-OR151111.rda'
		#gof.mutostabc.presim.musig(outdir, outfile, n.rep=200)
		#outfile	<- 'Normal-MESIG-MforMUTOST-OR151111.rda'
		#gof.mutostabc.presim.musig.ABCstar(outdir, outfile, n.rep=200)		
	}
	if(0)
	{
		indir	<- '~/Dropbox (Infectious Disease)/gof-abc/calc/example-paper'
		indir	<- paste(HOME, '/data/gof', sep='')
		#gof.mutostabc.MX.mu(indir=indir)
		#gof.mutostabc.MX.mu.ABCstar(indir=indir)
		#gof.mutostabc.MX.musig(indir=indir)
		# insufficient stat, then compare test --> aim to show that test not distr uniformly
		# idea: power can still drive vary small MAX
		gof.mutostabc.MX.musig.insuff(indir=indir)
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