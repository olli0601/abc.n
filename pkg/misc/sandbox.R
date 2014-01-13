NABC.DEFAULT.ANS<- {tmp<- c(0, 50, 1, NA, NA, NA, 0, 0, 0, 0, 0, 1, 1, 1, NA, NA, NA); names(tmp)<- c("lkl", "error", "pval","link.mc.obs","link.mc.sim", "rho.mc", "cil", "cir","tl","tr","al","ar","pfam.pval","ts.pval","nsim","mx.pow","rho.pow"); tmp}

#------------------------------------------------------------------------------------------------------------------------
#' Estimate summary parameter errors rho from unbiased Monte Carlo estimates rho.mc for all proposed theta including rejections 
#' @export
#' @param df			data frame with all proposed theta and corresponding rho.mc for each summary of interest 
#' @param theta.names	vector of theta names (columns in df)
#' @param rho.names		vector of rho names (columns in df)
#' @param thin			thinning factor in case there are many rows in df
#' @return	matrix containing the estimated rho (per column). The ith row corresponds to the ith theta in df.
nabc.exprho.at.theta<- function(df, theta.names, rho.names, thin=1)
{
	require(locfit)
	if(any(is.na(df)))	stop("unexpected NA in df")
	links.exp	<- sapply(rho.names,function(rho)
					{
						tmp		<- paste("locfit(",rho,'~',paste(theta.names,collapse=':',sep=''),", data=df, maxk=400)",sep='')
						lnk.fit	<- eval(parse(text=tmp))
						tmp		<- locfit:::preplot.locfit(lnk.fit, newdata= NULL, where="data", band = "none", tr = NULL, what = "coef", get.data = 0, f3d = 0)
						tmp$fit
					})	
	if(!is.matrix(links.exp))	
		links.exp<- as.matrix(links.exp)			
	colnames(links.exp)<- rho.names
	links.exp
}
#------------------------------------------------------------------------------------------------------------------------
nabc.get.locfit.links<- function(th.d, m, th.thin= 1, th.sep=100)
{
	require(locfit)
	th.names	<- colnames(m)[1:th.d]
	su.idx		<- (th.d+1):ncol(m)
	su.names	<- colnames(m)[su.idx]
	colnames(m)[1:th.d]<- paste("th",1:th.d,sep='')
	colnames(m)[su.idx]<- paste("su",seq_along(su.idx),sep='')
	links	<- lapply(seq_along(su.idx),function(k)
			{
				colnames(m)[su.idx]<- paste("su",seq_along(su.idx),sep='')
				colnames(m)[th.d+k]<- "rho"				
				lnk<-try(locfit(rho~th1:th2, data=m[seq.int(1,nrow(m),by=th.thin),],mint=20,deg=2,alpha=2,kern="epan"))
				if(inherits(lnk, "try-error"))
					lnk<- locfit(rho~th1:th2, data=m[seq.int(1,nrow(m),by=th.thin),],maxk=200,mint=20,deg=2,alpha=2,kern="epan")	
				#do locpoly regression on loc poly regression to save memory space
				theta<- expand.grid(lfmarg(lnk$box, rep(th.sep, lnk$mi["d"])))
				colnames(theta)<- paste("th",1:th.d,sep='')
				rhok<- predict(lnk, theta)
				lnk<-try(locfit(rho~th1:th2, data=cbind(theta,rho=rhok),mint=20,deg=2,alpha=2,kern="epan"))
				if(inherits(lnk, "try-error"))
					lnk<- locfit(rho~th1:th2, data=cbind(theta,rho=rhok),maxk=200,mint=20,deg=2,alpha=2,kern="epan")
				lnk	
			})
	links
}
#------------------------------------------------------------------------------------------------------------------------
nabc.get.locfit.jac<- function(th.d, m, th.thin= 1, th.sep=100)
{
	require(locfit)
	th.idx		<- 1:th.d
	th.names	<- colnames(m)[th.idx]
	su.idx		<- (th.d+1):ncol(m)
	su.names	<- colnames(m)[su.idx]
	colnames(m)[1:th.d]<- paste("th",th.idx,sep='')
	colnames(m)[su.idx]<- paste("su",seq_along(su.idx),sep='')
	jac	<- lapply(seq_along(su.idx),function(k)
			{
				jacrow<- lapply(th.idx,function(d)
						{
							colnames(m)[su.idx]<- paste("su",seq_along(su.idx),sep='')
							colnames(m)[th.d+k]<- "rho"				
							partialder<-try(locfit(rho~th1:th2, data=m[seq.int(1,nrow(m),by=th.thin),],mint=20,deg=2,alpha=2,kern="epan",deriv=paste("th",d,sep='')))
							if(inherits(partialder, "try-error"))
								partialder<- locfit(rho~th1:th2, data=m[seq.int(1,nrow(m),by=th.thin),],maxk=200,mint=20,deg=2,alpha=2,kern="epan",deriv=paste("th",d,sep=''))	
							#do locpoly regression on loc poly regression to save memory space
							theta<- expand.grid(lfmarg(partialder$box, rep(th.sep, partialder$mi["d"])))
							colnames(theta)<- paste("th",1:th.d,sep='')
							rhok<- predict(partialder, theta)
							partialder<-try(locfit(rho~th1:th2, data=cbind(theta,rho=rhok),mint=20,deg=2,alpha=2,kern="epan",deriv=paste("th",d,sep='')))
							if(inherits(partialder, "try-error"))
								partialder<- locfit(rho~th1:th2, data=cbind(theta,rho=rhok),maxk=200,mint=20,deg=2,alpha=2,kern="epan",deriv=paste("th",d,sep=''))
							partialder	
						})		
				names(jacrow)<- paste("dth",th.idx,sep='')
				jacrow
			})
	names(jac)<- paste("dL",su.idx,sep='')
	jac
}
#------------------------------------------------------------------------------------------------------------------------
checkjac<- function(a,sig2,ax,vx)
{
	if(length(sig2)>1 || length(ax)>1 || length(vx)>1)	stop("checkjac not vectorized")
	c( 2*a*sig2/vx, (1+a*a)/vx, (1-a*a)/(1+a*a+a^4), rep(0,length(a)))	
}
#------------------------------------------------------------------------------------------------------------------------
checklink<- function(a,sig2,ax,vx)
{
	c((1+a*a)*sig2/vx, project.nABC.movingavg.a2rho(a)-project.nABC.movingavg.a2rho(ax) )
}
#------------------------------------------------------------------------------------------------------------------------
nabc.estimate.jac<- function( links, th.eval, th.dh, ax, vx )
{
	th.d	<- ncol(th.eval)
	#print(th.d)
	if(length(th.dh)==1)	th.dh<- rep(th.dh,ncol(th.eval))
	if(length(th.dh)!=th.d)	stop("nabc.estimate.Jac: incorrect length of dh")	
	if(th.d!=2)	stop("nabc.estimate.Jac: incorrect th.d")		
	th.dh	<- matrix(c(th.dh[1],0, -th.dh[1],0, 0,th.dh[2], 0,-th.dh[2]), nrow=4, byrow=1)	#TODO th.d for higher th.d
#	print(th.dh); #plot(links[[2]]); print(th.eval)
	#th.eval<- th.eval[150:160,]
	jac		<- apply(th.eval,1,function(theta)
			{
#print(theta)
				if(0)			#exact to double check
				{
					diffq<- t(sapply(seq_along(links),function(k)
						{				
							tmp<- matrix(rep(theta,length(nrow(th.dh))),nrow=nrow(th.dh),ncol=th.d,byrow=1)+th.dh
#print(tmp)							
							rhok<- sapply(seq_len(nrow(tmp)),function(i)	checklink(tmp[i,1],tmp[i,2],ax, vx)		)[k,]
#print(rhok)
#stop()
							tmp2<- matrix(links[[k]]$box,ncol=th.d,byrow=1)
							tmp2<- apply(tmp,1,function(row){		all(row>=links[[k]]$box[1:th.d]	& row<=links[[k]]$box[(th.d+1):(2*th.d)]) 	})
							rhok[!tmp2]<- NA
							#print(tmp); print( tmp2 ); print(rhok);  							
							if(length(rhok)!=nrow(th.dh))	stop("nabc.estimate.Jac: unexpected behavior of predict.locfit - inappropriate links?")							
							rhok<- -apply(matrix(rhok,nrow=2), 2, diff)					#L_k(theta+h_d) - L_k(theta-h_d) for all d
							tmp<- rhok/apply(matrix(abs(th.dh[,1]+th.dh[,2]),nrow=2),2,sum)	#difference quotient for L_k and all d
							#print(tmp)
							tmp
						}))
#print(diffq); print(det(diffq))		
				}
#	
#stop()
				if(1)
				{
				diffq<- t(sapply(seq_along(links),function(k)
						{				
							tmp<- matrix(rep(theta,length(nrow(th.dh))),nrow=nrow(th.dh),ncol=th.d,byrow=1)+th.dh
							rhok<- predict(links[[k]], tmp )
#print(rhok)										
							tmp2<- matrix(links[[k]]$box,ncol=th.d,byrow=1)
							tmp2<- apply(tmp,1,function(row){		all(row>=links[[k]]$box[1:th.d]	& row<=links[[k]]$box[(th.d+1):(2*th.d)]) 	})
							rhok[!tmp2]<- NA
							#print(tmp); print( tmp2 ); print(rhok);  							
							if(length(rhok)!=nrow(th.dh))	stop("nabc.estimate.Jac: unexpected behavior of predict.locfit - inappropriate links?")							
							rhok<- -apply(matrix(rhok,nrow=2), 2, diff)					#L_k(theta+h_d) - L_k(theta-h_d) for all d
							tmp<- rhok/apply(matrix(abs(th.dh[,1]+th.dh[,2]),nrow=2),2,sum)	#difference quotient for L_k and all d
							#print(tmp)
							tmp
						}))
#print(diffq); print(det(diffq))		
				}

#stop()
				rownames(diffq)<- paste("L",seq_along(links),sep='')
				colnames(diffq)<- paste("th",1:th.d,sep='')
				#diffq[2,2]<- 0
				abs(det(diffq))								
			})
	jac
}
#------------------------------------------------------------------------------------------------------------------------
nabc.get.jacobian.2d<- function(m, th.mode, th.sep= rep(10,2), th.thin= 2000)
{
	require(locfit)
	th.d		<- 2	
	th.names	<- colnames(m)[1:th.d]
	su.idx		<- (th.d+1):ncol(m)
	su.names	<- colnames(m)[su.idx]
	colnames(m)[1:th.d]<- paste("th",1:th.d,sep='')
	colnames(m)[su.idx]<- paste("su",seq_along(su.idx),sep='')
	
	#print(colnames(m))
	th.range<- apply(m[,1:th.d],2,range)
	#print(th.range)
	th.dh	<- apply(th.range,2,diff) / th.sep
	#print(th.mode)
	th.eval	<- (th.range-rep(th.mode,each=2))/2		#+rep(th.mode,each=2)	
	th.eval	<- cbind(	th.mode+c(th.eval[1,1],0),th.mode+c(th.eval[2,1],0),
						th.mode+c(0,th.eval[1,2]),th.mode+c(0,th.eval[2,2]), th.mode		)
	th.dh	<- matrix(c(th.dh[1],0, -th.dh[1],0, 0,th.dh[2], 0,-th.dh[2]), ncol=4)				
	#print(th.dh)
	#print(th.eval)
	#remains to compute difference quotient to all rows in 'th.eval', using the dh in 'th.dh'
	#now first construct approximations to L_k
	links	<- lapply(seq_along(su.idx),function(k)
			{
				colnames(m)[su.idx]<- paste("su",seq_along(su.idx),sep='')
				colnames(m)[th.d+k]<- "rho"
				locfit(rho~th1:th2, data=m[seq.int(1,nrow(m),by=th.thin),])				
			})
	jac		<- apply(th.eval,2,function(theta)
			{
				#print(theta)
				diffq<- sapply(seq_along(links),function(k)
					{
						rhok<- apply(th.dh,2,function(h)		#the rho_k's for each of theta+h_1, theta-h_1, theta+h_2, theta-h_2, ...
							{																
								predict(links[[k]], matrix(theta+h, nrow=1, dimnames=list(c(),paste("th",1:th.d,sep=''))))								
							})					
						rhok<- -apply(matrix(rhok,nrow=2), 2, diff)					#L_k(theta+h_d) - L_k(theta-h_d) for all d
						rhok/apply(matrix(abs(th.dh[1,]+th.dh[2,]),nrow=2),2,sum)	#difference quotient for L_k and all d												
					})
				colnames(diffq)<- paste("L",seq_along(links),sep='')
				rownames(diffq)<- paste("th",1:th.d,sep='')
				abs(det(t(diffq)))				
			})
	ans		<- rbind(th.eval, jac)
	rownames(ans)<- c(th.names,"jac")
	ans
}
#------------------------------------------------------------------------------------------------------------------------
get.dist.fstretch<- function(sim, obs, args=NA, verbose= FALSE, alpha=0, tau.l=1, tau.u=1, plot=FALSE, xlab= NA, nbreaks=40, normal.test= "sf.test")
{
	#verbose<- 1
	#sim<- rnorm(100, 8.1, 1.1)
	if(any(is.na(sim)))	stop("get.dist.fstretch: error at 1a")
	if(any(is.na(obs)))	stop("get.dist.fstretch: error at 1b")	
	if(length(obs)!=length(sim))		stop("get.dist.fstretch: error at 1c")
	if(!is.na(args))
	{
		args<- strsplit(args,'/')[[1]]
		if(length(args)==3)	
		{
			tau.u<- as.numeric( args[2] )
			tau.l<- 1/tau.u
			alpha<- as.numeric( args[3] )
		}
		else if(length(args)==4)
		{
			tau.l<- as.numeric( args[2] )
			tau.u<- as.numeric( args[3] )
			alpha<- as.numeric( args[4] )
		}
		else 
			stop("get.dist.fstretch: error at 1A")
		args<- args[1]
	}
	if(alpha<0 || alpha>1)		stop("get.dist.fstretch: error at 1e")
	if(tau.u<1 )		stop("get.dist.fstretch: error at 1f")
	if(tau.l>1 )		stop("get.dist.fstretch: error at 1g")
	if(!normal.test%in%c("shapiro.test","lillie.test","cvm.test","ad.test","pearson.test","sf.test"))		stop("get.dist.fstretch: error at 1h")
	if(normal.test%in%c("lillie.test","cvm.test","ad.test","pearson.test","sf.test"))		require(nortest, lib.loc=LIB.LOC)
	df.sim<- length(sim) - 1
	df.obs<- length(obs) - 1 
	ans<- NABC.DEFAULT.ANS
	
	ans["pfam.pval"]<- ifelse(any(diff(sim)>0), ifelse(normal.test%in%c("shapiro.test","sf.test") && length(sim)>5000, eval(call(normal.test,sim[1:5000]))$p.value, eval(call(normal.test,sim))$p.value), 0)

	#get confidence intervals by numerical approximation; this is the fstretch R command by Wellek
	if(tau.l==1/tau.u && df.sim==df.obs)
	{
		fstretch.solvefor.cir<- function(cir, tau.u, n) pf( cir / tau.u,  n, n) - pf( 1 / (cir*tau.u),  n, n) - alpha
		catch<- try({
			ans["cir"]<- uniroot(	fstretch.solvefor.cir,	c(1,100), tol=.Machine$double.eps^0.5, tau.u=tau.u, n=df.sim)$root
			ans["cil"]<- 1/ans["cir"]
		})
		if(inherits(catch, "try-error"))
		{
			geterrmessage()						#the upper and lower bound of 'fstretch.solvefor.cir' may not cross 0 if tau.u is too large
			ans[c("cil","cir")]<- c(0,1e5)
		}
		ans["mx.pow"]<- pf(ans["cir"],df.obs, df.sim) - pf(ans["cil"],df.obs,df.sim)
	}
	else
	{
		tmp<- .Call("abcScaledF",	c(df.obs,df.sim,tau.l,tau.u,alpha,1e-10,100,0.05)	)
		if(tmp[4]>1e-10)	stop("get.dist.fstretch: error at 3a")
		ans[c("cil","cir","mx.pow")]<- tmp[1:3]
		
	}
	ans["error"]<- 			var(obs) / var(sim)
	ans["lkl"]<- 			df(ans["error"],df.obs,df.sim)	
	ans["pval"]<- 			pf(ans["error"],df.obs,df.sim)
	ans[c("al","ar")]<- 	c(0, 1 - diff( pf(ans[c("cil","cir")],df.obs, df.sim) ) )
	ans["pval"]<-			( ans["pval"] - ans["ar"]/2 ) / ( 1 - ans["ar"] )
	ans["link.mc.sim"]<- 	var(sim)
	ans["link.mc.obs"]<- 	var(obs)
	ans["rho.mc"]<- log(var(obs) / var(sim))
	if(verbose)	cat(paste(paste("\n{",args,"<-list(sim.var=",var(sim) ," , obs.var=",var(obs)," , alpha=",alpha," , tau.l=",tau.l," , tau.u=",tau.u,", log.ciu=",log(ans["cir"]),", ",sep=''),paste(names(ans), ans, collapse=', ', sep='='),")}",sep=''))
	ans
}
#------------------------------------------------------------------------------------------------------------------------
get.dist.mwu.equivalence<- function(sim, obs, args= NA, verbose= FALSE, alpha= 0.05, tau= 0, tau.translate= 0, plot= FALSE, xlab= NA, nbreaks= 40)
{
	verbose<- 1
	get.dist.mwu.equivalence.psd<- function(sim,obs)
	{
		if(length(sim)==1) sd(obs)
		else if(length(obs)==1)	sd(sim)
		else	sqrt(1/length(sim)+1/length(obs)) * sqrt(		(  (length(sim)-1)*var(sim)+(length(obs)-1)*var(obs)  )		/		(length(sim)+length(obs)-2)				)
	}
	get.dist.mwu.equivalence.tau<- function(x,c=8)
	{
		ans<- numeric(length(x));
		tmp<- x<c;
		ans[tmp]<- pnorm(x[tmp])-.5;
		ans[!tmp]<- (x[!tmp]-(c-4))/8;
		ans
	}
	#compute two sample t-test on either z-scores or untransformed data points
	if(any(is.na(sim)))	stop("get.dist.mwu.equivalence: error at 1a")
	if(any(is.na(obs)))	stop("get.dist.mwu.equivalence: error at 1b")
	if(length(sim)<1)		stop("get.dist.mwu.equivalence: error at 1c")
	if(length(obs)<1)		stop("get.dist.mwu.equivalence: error at 1d")
	if(length(obs)<2 && length(sim)<2)		stop("get.dist.mwu.equivalence: error at 1e")
	#if not missing, arguments 'args' always overwrite 'alpha' and 'nc'
	#expect args string of the form studenttXXX/beta/non-centrality
	if(!is.na(args))
	{
		args<- strsplit(args,'/')[[1]]
		if(length(args)!=5)	stop("get.dist.mwu.equivalence: error at 2a")
		standardize<- as.numeric( args[2] )
		tau.translate<- as.numeric( args[3] )
		tau<- as.numeric( args[4] )
		alpha<- as.numeric( args[5] )
		args<- args[1]
	}
	if(!standardize%in%c(0,1))	stop("get.dist.mwu.equivalence: error at 2b")
	if(alpha<0 || alpha>1)		stop("get.dist.mwu.equivalence: error at 2c")
	if(tau<=0 )		stop("get.dist.mwu.equivalence: error at 2d")
	if(tau.translate && standardize)	stop("get.dist.mwu.equivalence: error at 2e")
	ans<- NABC.DEFAULT.ANS
	#transform tau from  cdf tolerance of pi+ into location tolerance based on 6.10 and 6.2 in Wellek 2010
	#we want tau for |nu_k(x)-nu_k(theta)|<tau_k. so need an estimate of sigma. we pool over sim and obs, assuming equal sample sizes.
	if(!tau.translate && !standardize)
	{
		if(tau>.5)	tau<- 0.5
		if(tau<=0)					stop("get.dist.mwu.equivalence: error at 3a")
		tau.loc<- qnorm(.5+tau)*get.dist.mwu.equivalence.psd(sim,obs)
	}
	else if(tau.translate)
	{
		tau.loc<- tau
		tau<- tau.loc / get.dist.mwu.equivalence.psd(sim,obs)
		if(tau>8)		#allow for tau>0.5 -- this does not make much sense for MWU but can be useful for annealing purposes. in any case raise a warning
		{
			tau<- get.dist.mwu.equivalence.tau(tau)
			options(warn=1)
			warning(paste("get.dist.mwu.equivalence: tau is larger than 0.5:",tau))
			options(warn=2)
		}
		else
			tau<- pnorm( tau ) - 0.5
	}
	else
		tau.loc<- tau
	if(length(sim)==1)	sim<- obs-median(obs)+sim		#apply shift hypothesis
	if(length(obs)==1)	obs<- sim-median(sim)+obs		#apply shift hypothesis
#print(obs); print(sim)
	if(0)
	{
		#compute mwu test statistic
		w<- sum( sapply(obs, function(x){			length(which(x>=sim))		}) )
		#compute PIxxy
		tmp<-   obs[-1]
		tmp2<- matrix(ncol= length(obs)-1, nrow= length(obs)-1)
		tmp2<- col(tmp2) >= row(tmp2)
		tmp2<-	c(  sapply(seq_along(tmp),function(i){     tmp[tmp2[i,]]      }), recursive=1 )	#elements of upper triangular matrix with entries x[2] .. x[m],       x[3] .. x[m], 	..., x[m,m]		in this order
		tmp<- rep(obs[-length(obs)], times= seq(length(obs)-1,1,-1))										#m-1 times x[1], m-2 times x[2], ..., 1 times x[m-1]	in this order
		#compute PIxxy and append to W
		w<- c(w, sum(sapply(seq_along(tmp), function(i){		length(which(tmp[i]>=sim & tmp2[i]>=sim))		})))
		#compute PIxyy
		tmp<-   sim[-1]
		tmp2<- matrix(ncol= length(sim)-1, nrow= length(sim)-1)
		tmp2<- col(tmp2) >= row(tmp2)
		tmp2<-	c(  sapply(seq_along(tmp),function(i){     tmp[tmp2[i,]]      }), recursive=1 )	#elements of upper triangular matrix with entries x[2] .. x[m],       x[3] .. x[m], 	..., x[m,m]		in this order
		tmp<- rep(sim[-length(sim)], times= seq(length(sim)-1,1,-1))										#m-1 times x[1], m-2 times x[2], ..., 1 times x[m-1]	in this order
		#compute PIxyy and append to W
		w<- c(w, sum(sapply(seq_along(tmp), function(i){		length(which(obs>=tmp[i] & obs>=tmp2[i]))		})))
		#check with simple implementation taken from mawi.R in Wellek 2010
		#m<- length(obs); n<- length(obs); x<- obs; y<- sim; wxy<- pihxxy<- pihxyy<- 0
		#for (i in 1:m) for (j in 1:n) wxy <- wxy + trunc(0.5*(sign(x[i] - y[j]) + 1))
		#for (i in 1:m) for (j1 in 1:(n-1)) for (j2 in (j1+1):n) pihxyy <- pihxyy + trunc(0.5*(sign(x[i] - max(y[j1],y[j2])) + 1))
		#for (i1 in 1:(m-1)) for (i2 in (i1+1):m) for (j in 1:n) pihxxy <- pihxxy + trunc(0.5*(sign(min(x[i1],x[i2]) - y[j]) + 1))
		#normalize
		tmp<- length(obs)
		tmp2<- length(sim)
		w<- w * c(  1/(tmp*tmp2),    2/(tmp*(tmp-1)*tmp2),    	2/(tmp2*(tmp2-1)*tmp)   )
		#estimate sd of mwu statistic
		w<- c(w,     	sqrt( (w[1] - (tmp+tmp2-1)*w[1]*w[1] + (tmp-1)*w[2]  + (tmp2-1)*w[3])/(tmp*tmp2) )    )
	}
	else		#faster than vectorized version. w[1]= W+  w[2]= (m-1)PIxxy   w[3]= (n-1)PIxyy  w[4]= sd(W+)   w[5]=W-  w[6]=sd(W-)
		w<- .Call("abcMWUE",obs,sim)
	#can compute lkl and pval at tau=0 explicitly for small sample sizes
	if(0)
	{	#this is slow and currently not needed
		options(warn=1)
		require(coin, warn.conflicts = FALSE)
		options(warn=2)

		tmp<- data.frame( val= c(obs,sim), group= c(rep("1",length(obs)), rep("2",length(sim))) )
		tmp<- wilcox_test(val ~ group, data = tmp,    distribution = "exact" )		#this assumes unpaired data
		idx<-  which(support(tmp)<=as.numeric(statistic(tmp,"standardized")))
		ans[c("lkl","pval")]<-	c( dperm( tmp, support(tmp)[	ifelse(!length(idx), 1, idx[length(idx)])	] ), pvalue(tmp) )
	}
#print(w)
	#compute test statistic and critical bound, assuming symmetric tolerance.
	#the CI regions in (6.19) are numerically not stable (large ncp). We use (4.5) and construct a TOST.
	#this leads to a numerically stable procedure
	if(is.na(w[4]) && sd(sim)==sd(obs))							#if variance in obs & sim equal, then sigma(w-)=sigma(w+), and sigma(w-) might not be NA
		w[4]<- w[6]
	if(is.na(w[4]) || !w[4] )
	{
		tmp<- rep(NA,2)
		ans["error"]<- ifelse(tau>=0.5, 1, 2)*alpha			#always accept when tau >= 0.5
	}
	else if(standardize)
	{
		tmp<- c( (w[1] - 0.5) / w[4] + tau, (w[1] - 0.5) / w[4] - tau)		#here, tau is the same as tau/sigma(W_+) below, and hence independent of y 
		which.not.reject<- which( c(pnorm( tmp[1] )<1-alpha, pnorm( tmp[2] )>alpha) )		
		if(length(which.not.reject)%in%c(0,2))
			which.not.reject<- which.min(abs(c(1-alpha-pnorm( tmp[1] ), pnorm( tmp[2] )-alpha )))
		ans["error"]<-	ifelse(which.not.reject==1, 		1-pnorm( tmp[which.not.reject] ),		pnorm( tmp[which.not.reject] ) )
	}
	else
	{
		tmp<- c( (w[1] - 0.5 + tau) / w[4], (w[1] - 0.5 - tau) / w[4] )		#compute ZL= ( W_+ - 0.5 + tau ) / sigma(W_+) and ZU for TOST; 		
		which.not.reject<- which( c(pnorm( tmp[1] )<1-alpha, pnorm( tmp[2] )>alpha) )
		#only if length==0, TOST will be rejected ie ABC accept
		if(length(which.not.reject)%in%c(0,2))
		{
			#decide which of the two pvalues are to be reported --> figure out which one is closer to boundary
			#both upper and lower test statistics indicate mean difference < -tau and mean difference > tau  	OR 		mean difference >= -tau and mean difference <= tau
			which.not.reject<- which.min(abs(c(1-alpha-pnorm( tmp[1] ), pnorm( tmp[2] )-alpha )))
		}
		#the pvalue of ZU is the lower tail, but for ZL it is the upper tail, so..
		ans["error"]<-	ifelse(which.not.reject==1, 		1-pnorm( tmp[which.not.reject] ),		pnorm( tmp[which.not.reject] ) )
	}
	ans[c("al","ar")]<- ans[c("cil","cir")]<-	c(0,alpha)
	ans["mx.pow"]<- tau.loc			#just for now
	ans["link.mc.sim"]<- 	-w[1]
	ans["link.mc.obs"]<- 	w[1]
	ans["rho.mc"]<- 		w[1] - 0.5
	
	if(verbose)	cat(paste(paste("\n{",args,"<-list(sim.median=",median(sim) ,", obs.median=",median(obs) ,", ZL=",tmp[1] ,", ZU=",tmp[2] ,", CU=",qnorm( alpha ) ,", alpha=",alpha,", tau.cdf=",tau,", tau.loc=",tau.loc,", ",sep=''),paste(names(ans), ans, collapse=', ', sep='='),")}",sep=''))
	if(plot)
	{
		breaks<-  range(c(sim,obs))
		breaks[1]<- breaks[1]*ifelse(breaks[1]<0, 1.2, 0.8)
		breaks[2]<- breaks[2]*ifelse(breaks[2]>0, 1.2, 0.8)
		breaks<- seq(from= breaks[1], to= breaks[2], by= (breaks[2]-breaks[1])/nbreaks)
		xlim<- range(c(sim,obs))
		sim.h<- hist(sim,	breaks= breaks,	plot= F)
		obs.h<- hist(obs,	breaks= breaks,	plot= F)
		ylim<- c(0,max(c(obs.h$intensities, sim.h$intensities),na.rm=TRUE)*1.1)
		plot(1,1,xlab=xlab, xlim= xlim,type='n',ylim=ylim,ylab="probability",main="")
		plot(sim.h,freq=F,add=TRUE,col=myFadeCol("#0080FFFF",0.6),border=myFadeCol("#0080FFFF",0.6))
		plot(obs.h,freq=F,add=TRUE)
		points( sim, rep(ylim[2],length(sim)),col=myFadeCol("#0080FFFF",0.6),pch=20,cex=2.5 )
		points( obs, rep(ylim[2],length(obs)),pch=2,cex=1.5 )
	}
	ans
}

#---------------------------------------------------------------------------------------------------------------------------
nabc.getlevelset.2d<- function(df, lnk.name, theta.names, rho.eq=0, rho.eq.sep=100, rho.eq.q= 0.075, theta.sep= 100, method="quantile",plot=0, verbose=0)
{
	if(!method%in%c("quantile","fixed"))	stop("nabc.getlevelset.2d: error at 1a")
	if(method=="quantile")
	{
		rho.eq.sep.u<- length(df[,lnk.name])
		rho.eq.sep.l<- 1		
	}
	require(locfit)
	tmp<- paste("locfit(",lnk.name,'~',paste(theta.names,collapse=':',sep=''),", data=df,maxk=200)",sep='')
	lnk.locfit<- eval(parse(text=tmp))
	if(plot)
	{
		plot(lnk.locfit)
		stop()
	}
	lnk.xrange<- lfmarg(lnk.locfit$box, rep(ifelse(lnk.locfit$mi["d"]==1, 100, theta.sep), lnk.locfit$mi["d"]))
	names(lnk.xrange)<- theta.names
	lnk.pred<- locfit:::preplot.locfit(lnk.locfit, lnk.xrange, band = "none", tr = NULL, what = "coef", get.data = 0, f3d = 0)
	lnk.len<- sapply(lnk.xrange,length)
	lnk.pred<- array(lnk.pred$fit, dim=lnk.len)
	lnk.xsep<- sapply(lnk.xrange, function(x)	diff(x[1:2])	)
	if(method=="quantile")
	{	
		while(1)
		{		
			rho.eq.sep<- mean(c(rho.eq.sep.l,rho.eq.sep.u))
			rho.eq.eps<- range(df[,lnk.name]) /  rho.eq.sep
			lnk.eqidx<- which( lnk.pred<(rho.eq+rho.eq.eps) & lnk.pred>(rho.eq-rho.eq.eps) )
			if(verbose) cat(paste("\nrho.eq.seps are",rho.eq.sep.l,rho.eq.sep.u,"\tlevel set length is",length(lnk.eqidx)))
			if(	abs(length(lnk.eqidx)/length(df[,lnk.name])-rho.eq.q)<rho.eq.q/20 ||
					rho.eq.sep.l==rho.eq.sep.u	)
				break
			if(length(lnk.eqidx)/length(df[,lnk.name])<rho.eq.q)
				rho.eq.sep.u<- rho.eq.sep
			else
				rho.eq.sep.l<- rho.eq.sep
		}
	}
	else
	{
		rho.eq.eps<- range(df[,lnk.name]) / rho.eq.sep
		lnk.eqidx<- which( lnk.pred<(rho.eq+rho.eq.eps) & lnk.pred>(rho.eq-rho.eq.eps) )
	}	
	#print(lnk.locfit)
	#print(lnk.len)
	#print(lnk.xsep)
	#print(dim(lnk.pred))
	#print(lnk.eqidx)	
	lnk.eqx<- matrix(NA,nrow=length(lnk.len),ncol=length(lnk.eqidx),dimnames=list(names(lnk.xrange),c()))
	for(i in rev(seq_along(lnk.len)))
	{
		lnk.eqx[i,]<- (lnk.eqidx-1) %/% prod(lnk.len[-length(lnk.len)]) + 1		#idx wrt ith theta
		lnk.eqidx<- (lnk.eqidx-1) %% prod(lnk.len[-length(lnk.len)]) + 1
		lnk.len<- lnk.len[-length(lnk.len)]
	}
	lnk.eqx<- sapply(seq_len(ncol(lnk.eqx)),function(j)
			{
				sapply(seq_len(nrow(lnk.eqx)),function(i)		lnk.xrange[[i]][ lnk.eqx[i,j] ]	)
			}) 	
	rownames(lnk.eqx)<- names(lnk.xrange)
	#print(lnk.eqx)
	lnk.eqx	
}
#---------------------------------------------------------------------------------------------------------------------------
nabc.getlevelsetintersection.2d<- function(lsets,theta.names,rho.eq.sep= 25,plot=0,verbose=0)
{
	theta.range<- sapply(seq_len(nrow(lsets[[1]])),function(d)
			{
				range(c(sapply(seq_along(lsets),function(k){		range(lsets[[k]][d,])	}),recursive=1))		
			})
	colnames(theta.range)<- theta.names
	theta.sep<- apply(theta.range,2,diff)/rho.eq.sep 
	
	theta.intersection<- lsets[[1]]
	lsets.intersection<- vector("list",length(lsets)-1)
	for(k in seq_along(lsets)[-1])	
	{
		if(verbose)	cat(paste("\nnumber of theta in intersection",ncol(theta.intersection),"\nnow intersect w summary",names(lsets)[k],"\n"))		
		tmp<- .Call("abcIntersectLevelSets", lsets[[k]], theta.intersection, theta.sep)				
		tmp<- which(tmp<sum(theta.sep*theta.sep))				
		theta.close<- matrix(	c( (tmp-1) %/% ncol(lsets[[k]]) + 1, (tmp-1) %% ncol(lsets[[k]]) + 1 ), nrow=2,ncol=length(tmp),byrow=1,dimnames=list(c("1","k"),c()) )		
		tmp<- unique(theta.close["1",])
		lsets.intersection[[k-1]]<- theta.intersection[,unique(theta.close["1",]),drop=0]				
		if(!length(tmp))
			break
		theta.intersection<- theta.intersection[,tmp,drop=0]
	}	
	if(plot)
	{
		plot(	theta.intersection[1,],theta.intersection[2,],
				xlim=range(theta.intersection[1,]),ylim=range(theta.intersection[2,]),xlab=rownames(theta.intersection)[1],ylab=rownames(theta.intersection)[2],
				type='p',pch=19,col=myFadeCol("black",0.3))
	}
	if(verbose) cat(paste("\nfinal number of theta in intersection",ncol(theta.intersection),"\n"))
	theta.intersection
}
#---------------------------------------------------------------------------------------------------------------------------
nabc.getlevelset.3d<- function(df, lnk.name, theta.names, rho.eq=0, rho.eq.sep=100, rho.eq.q= 0.075, theta.sep= 100, method="quantile",plot=0)
{	
	if(!method%in%c("quantile","fixed"))	stop("nabc.getlevelset.3d: error at 1a")
	if(method=="quantile")
	{
		rho.eq.sep.u<- length(df[,lnk.name])
		rho.eq.sep.l<- 1		
	}			
	tmp<- paste("locfit(",lnk.name,'~',paste(theta.names,collapse=':',sep=''),", data=df)",sep='')								
	lnk.locfit<- eval(parse(text=tmp))					
	lnk.xrange<- lfmarg(lnk.locfit$box, rep(ifelse(lnk.locfit$mi["d"]==1, 100, theta.sep), lnk.locfit$mi["d"]))
	names(lnk.xrange)<- theta.names
	lnk.pred<- locfit:::preplot.locfit(lnk.locfit, lnk.xrange, band = "none", tr = NULL, what = "coef", get.data = 0, f3d = 0)
	lnk.len<- sapply(lnk.xrange,length)
	lnk.pred<- array(lnk.pred$fit, dim=lnk.len)
	lnk.xsep<- sapply(lnk.xrange, function(x)	diff(x[1:2])	)
	
	if(method=="quantile")
	{	
		while(1)
		{		
			rho.eq.sep<- mean(c(rho.eq.sep.l,rho.eq.sep.u))
			rho.eq.eps<- range(df[,lnk.name]) /  rho.eq.sep
			lnk.eqidx<- which( lnk.pred<(rho.eq+rho.eq.eps) & lnk.pred>(rho.eq-rho.eq.eps) )
			cat(paste("\nrho.eq.seps are",rho.eq.sep.l,rho.eq.sep.u,"\tlevel set length is",length(lnk.eqidx)))
			if(	abs(length(lnk.eqidx)/length(df[,lnk.name])-rho.eq.q)<rho.eq.q/20 ||
					rho.eq.sep.l==rho.eq.sep.u	)
				break
			if(length(lnk.eqidx)/length(df[,lnk.name])<rho.eq.q)
				rho.eq.sep.u<- rho.eq.sep
			else
				rho.eq.sep.l<- rho.eq.sep
		}
	}
	else
	{
		rho.eq.eps<- range(df[,lnk.name]) / rho.eq.sep
		lnk.eqidx<- which( lnk.pred<(rho.eq+rho.eq.eps) & lnk.pred>(rho.eq-rho.eq.eps) )
	}
	lnk.eqx<- matrix(NA,nrow=length(lnk.len),ncol=length(lnk.eqidx),dimnames=list(names(lnk.xrange),c()))
	for(i in rev(seq_along(lnk.len)))
	{
		lnk.eqx[i,]<- (lnk.eqidx-1) %/% prod(lnk.len[-length(lnk.len)]) + 1		#idx wrt ith theta
		lnk.eqidx<- (lnk.eqidx-1) %% prod(lnk.len[-length(lnk.len)]) + 1
		lnk.len<- lnk.len[-length(lnk.len)]
	}
	lnk.eqx<- sapply(seq_len(ncol(lnk.eqx)),function(j)
			{
				sapply(seq_len(nrow(lnk.eqx)),function(i)		lnk.xrange[[i]][ lnk.eqx[i,j] ]	)
			}) 	
	rownames(lnk.eqx)<- names(lnk.xrange)
	if(plot)
	{
		ABC.CI.MMCMC.plot.trellis.levelset(lnk.locfit, zlab=theta.names[length(theta.names)], ysep=rho.eq.eps)		
	}
	lnk.eqx
}
###############################################################################
my.qqline<- function (y, datax = FALSE, u=1.95, ...) 
{
	y 	<- quantile(y[!is.na(y)], c(0.25, 0.75))
	x 	<- qnorm(c(0.25, 0.75))
	if (datax) 
	{
		slope 	<- diff(x)/diff(y)
		int 	<- x[1L] - slope * y[1L]
	}
	else 
	{
		slope 	<- diff(y)/diff(x)
		int 	<- y[1L] - slope * x[1L]
	}
	x	<- qnorm(c(0.01, 0.99))
	x	<- seq(x[1],x[2],len=1e3)
	z	<- int + slope*x
	lines(x,z, ...)
	z	<- u*sqrt( pnorm(x)*(1-pnorm(x)) ) / ( dnorm(x)*sqrt(length(y)) )
	y95u<- int + slope*(x+z)
	y95l<- int + slope*(x-z)
	lines(x,y95u, ...)
	lines(x,y95l, ...)
}
###############################################################################
plot.2D.dens<- function(x,y,xlab,ylab,xlim=NA,ylim=NA,nbin=NA,width.infl=2,n.hists=5,method="gauss", palette= "topo", persp.theta= -30, persp.phi= 30, zero.abline=TRUE, ...)
{
	if(!method%in%c("gauss","ash","persp"))	stop("plot.2D.dens: exception 1a")
	if(!palette%in%c("topo","heat","gray"))	stop("plot.2D.dens: exception 1b")
	switch(method,
			gauss=
					{
						require(KernSmooth)
						require(fields)
						x.bw<- width.infl*diff(summary(x)[c(2,5)])
						y.bw<- width.infl*diff(summary(y)[c(2,5)])
						if(!x.bw) x.bw<- EPS
						if(!y.bw) y.bw<- EPS
						if(any(is.na(xlim)))	xlim<- range(x)+c(-1.5,1.5)*x.bw
						if(any(is.na(ylim)))	ylim<- range(y)+c(-1.5,1.5)*y.bw
						dens <- bkde2D(cbind(x, y), range.x=list(xlim,ylim),bandwidth=c(x.bw,y.bw))
						contour(dens$x1, dens$x2, dens$fhat,xlab=xlab,ylab=ylab)
						if(zero.abline) abline(v=0,col="black",lty=3,lwd=1.5)
						if(zero.abline) abline(h=0,col="black",lty=3,lwd=1.5)
					},
			ash={
				require(ash)
				if(any(is.na(xlim))) xlim<- range(x)*1.05
				if(any(is.na(ylim))) ylim<- range(y)*1.05
				if(any(is.na(nbin))) nbin<- 2*c(nclass.Sturges(x),nclass.Sturges(y))
				bins<- bin2(cbind(x, y), ab=rbind(xlim,ylim),nbin=nbin)
				f <- ash2(bins,rep(n.hists,2))
				#image(f$x,f$y,f$z, col=rainbow(50,start= 4/6,end=0),xlab=xlab,ylab=ylab)
				if(palette=="topo")					image(f$x,f$y,f$z, col=tail( topo.colors(trunc(50*1.4)), 50 ),xlab=xlab,ylab=ylab,...)
				else if(palette=="gray")			image(f$x,f$y,f$z, col=head( rev(gray(seq(0,.95,len=trunc(50*1.4)))), 50),xlab=xlab,ylab=ylab,...)					
				else								image(f$x,f$y,f$z, col=heat.colors( 50 ),xlab=xlab,ylab=ylab,...)
				contour(f$x,f$y,f$z,add=TRUE, nlevels= 5)
				if(zero.abline) abline(v=0,col="black",lty=2,lwd=2)
				if(zero.abline) abline(h=0,col="black",lty=2,lwd=2)
			},
			persp={
				require(ash)
				require(MASS)
				if(any(is.na(xlim))) xlim<- range(x)*1.05
				if(any(is.na(ylim))) ylim<- range(y)*1.05
				bins<- bin2(cbind(x, y), ab=rbind(xlim,ylim),nbin=2*c(nclass.Sturges(x),nclass.Sturges(y)))
				f <- ash2(bins,rep(n.hists,2))
				
				nrz <- nrow(f$z)
				ncz <- ncol(f$z)
				col<- tail(topo.colors(trunc(1.4 * 50)),50)
				fcol      <- col[trunc(f$z / max(f$z)*(50-1))+1]
				dim(fcol) <- c(nrz,ncz)
				fcol      <- fcol[-nrz,-ncz]
				par(mar=c(1/2,1/2,1/2,1/2))
				persp(x=f$x,y=f$y,z=f$z,col= fcol,theta=persp.theta,phi=persp.phi,xlab=xlab,ylab=ylab,zlab='', ticktype= "detailed" )
			})
}
###############################################################################
plot.persplocfit<- function(x, pv, theta= 30, phi= 20, palette= "gray",tcl=-0.05,...)
{	
	d <- x$mi["d"]
	ev <- x$mi["ev"]
	where <- "grid"
	
	pv <- match(pv, x$vnames)
	tv <- (1:d)[-pv]
	vrs <- c(pv, tv)
	if (any(duplicated(vrs))) 
		warning("Duplicated variables in pv, tv")
	if (any((vrs <= 0) | (vrs > d))) 
		stop("Invalid variable numbers in pv, tv")
	m <- ifelse(d == 1, 100, 40) 
	m <- rep(m, d)
	m[tv] <- mtv<- 6		
	xl <- x$box
	marg <- lfmarg(xl, m)
	pred <- locfit:::preplot.locfit(x, marg, band = "none", tr = NULL, what = "coef", get.data = 0, f3d = 0)
	z<- matrix(pred$fit, nrow = length(marg[[1]]))
	nbcol <- 100
	if(palette=="gray")	
		color<- head( rev(gray(seq(0,0.95,len=trunc(nbcol*1.4)))), nbcol)
	else
	{
		jet.colors <- colorRampPalette( c("blue", "green") )
		color <- jet.colors(nbcol)
	}
	# Compute the z-value at the facet centres
	nrz <- nrow(z)
	ncz <- ncol(z)
	zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
	# Recode facet z-values into color indices
	facetcol <- cut(zfacet, nbcol)		
	par(mar=c(0,2.5,0,0), tcl=tcl)
	pmat<- persp(marg[[1]], marg[[2]], z, zlim= range(z)*1.1, col = color[facetcol], shade= 0.1, border=NA, ticktype = "detailed", ltheta = 120,theta = theta, phi = phi, expand = 0.75, box=1, ... )
	
	list(pmat=pmat, x= marg[[1]], y= marg[[2]], z= z)
}