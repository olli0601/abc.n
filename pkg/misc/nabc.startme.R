#! /Library/Frameworks/R.framework/Versions/2.15/Resources/bin/Rscript
#--DSCR------- #! /opt/apps/R-2.15.1/lib64/R/bin/Rscript
#--CX1-------- #! /apps/R/2.13.0/lib64/R/bin/Rscript
###############################################################################
#
#
# Author: oliver ratmann
# file: phdes.startme.R
#
# usage from R:
#> setwd("/Users/Oliver/git/abc.n/pkg")
#> source("misc/nabc.startme.R")
# usage from bash:
#> misc/abc-n.startme.R --help
#
#
# Installation instructions:
#	1 make sure the first line in this file points to your Rscript
#	2.1	create a directory CODE.HOME and set the CODE.HOME variable below to this path
#	2.2	create CODE.HOME/src_tipclust and copy all the R files into this directory
#	2.3 create a directory HOME and set the HOME variable below to this path
#
#
###############################################################################
args <- commandArgs()
if(!any(args=='--args'))
	args<- vector("numeric",0)
if(any(args=='--args'))
	args<- args[-(1:match("--args", args)) ]

CODE.HOME	<<- "/Users/Oliver/git/abc.n/pkg"
#CODE.HOME	<<- "/work/or105/libs/abc.n/pkg"
HOME		<<- "/Users/Oliver/workspace_sandbox/phylody"
#HOME		<<- "/work/or105/phylody"
DATA		<<- paste(HOME,"data",sep='/')
NABC.DEBUG	<<- 0
LIB.LOC		<<- NULL
#LIB.LOC		<<- paste(CODE.HOME,"../",sep='')
EPS			<<- 1e-12
default.fun	<- "my.make.documentation"
default.fun	<- "project.nABC.TOST"
default.fun	<- "project.nABC.StretchedChi2"
#default.fun<- "project.nABC.movingavg"
###############################################################################
#if(length(args) && !is.loaded("tipc_tabulate_after_sample"))
#{
#	file<- paste(CODE.HOME,"src_tipclust",paste("libtipcr",.Platform$dynlib.ext,sep=''),sep='/')
#	cat(paste("\nloading",file,'\n',sep=' '))
#	dyn.load(file)
#}
#cat(paste("is.loaded('tipcr')->",is.loaded("tipc_tabulate_after_sample"),'\n'))
###############################################################################
function.list<-c(list.files(path= paste(CODE.HOME,"R",sep='/'), pattern = ".R$", all.files = FALSE,
		full.names = TRUE, recursive = FALSE),paste(CODE.HOME,"misc","nabc.prjcts.R",sep='/'))
sapply(function.list,function(x) source(x,echo=FALSE,print.eval=FALSE, verbose=FALSE))
###############################################################################
my.mkdir<-function(root,data.name)
{
	if(length(dir(root,pattern=paste('^',data.name,'$',sep='')))==0)
		system(paste("mkdir ",paste(root,data.name,sep='/'),sep=''))
}

my.make.documentation<- function()
{
	require(roxygen2)		
	roxygenize(CODE.HOME)
}

my.fade.col<-function(col,alpha=0.5)
{
	return(rgb(col2rgb(col)[1]/255,col2rgb(col)[2]/255,col2rgb(col)[3]/255,alpha))
}
###############################################################################
my.mkdir(HOME,"data")
my.mkdir(HOME,"pdf")
my.mkdir(HOME,"script")
argv<- list()
if(length(args))
{
	tmp<- na.omit(sapply(args,function(arg)
					{
						switch(substr(arg,2,4),
								exe= return(substr(arg,5,nchar(arg))),
								NA)
					}))
	if(length(tmp)!=0)
	{
		if(length(tmp)>1) stop("abc-n.startme.R: duplicate -exe")
		else default.fun<- switch(tmp[1],
					MAKE.DOCUMENTATION		 = "my.make.documentation",
					MUTOST					 = "project.nABC.TOST",
					CHISQU					 = "project.nABC.StretchedChi2",
					ACFTOST					 = "project.nABC.movingavg"
					)
	}
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
if(NABC.DEBUG)	options(error= my.dumpframes)	
require(abc.n)
cat(paste("\nabc-n.startme.R: ",ifelse(NABC.DEBUG,"debug",""),"call",default.fun))
do.call(default.fun,list()) 	
cat("\nabc-n: ",ifelse(NABC.DEBUG,"debug","")," end\n")
