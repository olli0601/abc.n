#! /Library/Frameworks/R.framework/Versions/2.15/Resources/bin/Rscript
##
##	first line in shell script starts with #! and points to absolute path to Rscript
##	CHANGE  as needed
##
##! /apps/R/2.15/lib64/R/bin/Rscript
###############################################################################
#
#	project scripts that can be run from command line, without re-building the package all the time,
# 	because the R files are re-loaded below
#
# usage from R:
#> setwd("/Users/Oliver/git/abc.n/pkg")
#> source("misc/nabc.startme.R")
# usage from bash:
#> misc/abc-n.startme.R --help
#
###############################################################################
args <- commandArgs()
if(!any(args=='--args'))
	args<- vector("numeric",0)
if(any(args=='--args'))
	args<- args[-(1:match("--args", args)) ]
###############################################################################
CODE.HOME	<<- "/Users/Oliver/git/abc.n/pkg"
#CODE.HOME	<<- "/home/koelle/or7/utils/abc.n/pkg"
#CODE.HOME	<<- "/work/or105/libs/abc.n/pkg"
HOME		<<- "/Users/Oliver/workspace_sandbox/phylody"
#HOME		<<- "/home/koelle/or7/phylody"
#HOME		<<- "/work/or105/phylody"
DATA		<<- paste(HOME,"data",sep='/')
NABC.DEBUG	<<- 0
LIB.LOC		<<- NULL
#LIB.LOC		<<- paste(CODE.HOME,"../",sep='')
EPS			<<- 1e-12
###############################################################################
#	the default script to be called if -exe is not specified on the command line
default.fun	<- "my.make.documentation"
default.fun	<- "project.nABC.TOST"
default.fun <- "project.nABC.compareSEIRS"
default.fun	<- "project.nABC.StretchedChi2"
#default.fun	<- "nabc.test.mutost.calibrate"
#default.fun	<- "nabc.test.chi2stretch.calibrate"
#default.fun	<- "nabc.test.chi2stretch.montecarlo.calibrated.tau.and.m"
#default.fun	<- "nabc.test.chi2stretch.montecarlo.calibrated.tau.and.increasing.m"
#default.fun	<- "nabc.test.acf.montecarlo.calibrated.tau.and.m"
default.fun	<- "nabc.test.acf.montecarlo.vary.a"
#default.fun<- "project.nABC.movingavg"
###############################################################################
#	select script specified with -exe on the command line. If missing, start default script 'default.fun'.
argv<- list()
if(length(args))
{
	tmp<- na.omit(sapply(args,function(arg)
					{
						switch(substr(arg,2,4),
								exe= return(substr(arg,6,nchar(arg))),
								NA)
					}))
	if(length(tmp)!=0)
	{
		if(length(tmp)>1) stop("abc-n.startme.R: duplicate -exe")
		else default.fun<- switch(tmp[1],
					ROXYGENIZE				 = "package.roxygenize",
					MUTOST					 = "project.nABC.TOST",
					CHISQU					 = "project.nABC.StretchedChi2",
					ACFTOST					 = "project.nABC.movingavg"
			)
	}
	tmp<- na.omit(sapply(args,function(arg)
					{
						switch(substr(arg,2,10),
								code.home= return(substr(arg,12,nchar(arg))),
								NA)
					}))	
	if(length(tmp)!=0)	CODE.HOME<<- tmp[1]	
	tmp<- na.omit(sapply(args,function(arg)
					{
						switch(substr(arg,2,6),
								debug= 1,
								NA)
					}))		
	if(length(tmp)!=0)	NABC.DEBUG<<- tmp[1]	
	argv<<- args
}
###############################################################################
#	re-load all R files
require(data.table)
function.list<-c(list.files(path= paste(CODE.HOME,"R",sep='/'), pattern = ".R$", all.files = FALSE,
		full.names = TRUE, recursive = FALSE),paste(CODE.HOME,"misc","nabc.prjcts.R",sep='/'),paste(CODE.HOME,"misc","profiling.R",sep='/'))
sapply(function.list,function(x) source(x,echo=FALSE,print.eval=FALSE, verbose=FALSE))
###############################################################################
#	run script
if(NABC.DEBUG)	options(error= package.dumpframes)	
cat(paste("\nabc-n.startme.R: ",ifelse(NABC.DEBUG,"debug",""),"call",default.fun))
do.call(default.fun,list()) 	
cat("\nabc-n: ",ifelse(NABC.DEBUG,"debug","")," end\n")
