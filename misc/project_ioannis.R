yia.check.normality<- function()
{
	indir	<- '~/Dropbox (Personal)/abc/Ioannis_LondonSchool'
	infile	<- '141003_ibmoutput.RData'
	infile	<- 'Data.RData'
	file	<- paste(indir, infile, sep='/')
	z		<- load(file)
	
	#	Outputs is an 10000 x 14 x 18 array, containing the 10000 reps, 14 different design point and 18 outputs. 
	#	X is a 22 x 14 matrix containing the 14 22-dimensional design points
	
	#	convert into data.table
	df	<- lapply( 1:dim(Outputs)[3], function(s)
				{
					tmp				<- Outputs[,,s]
					colnames(tmp)	<- paste('design',seq_len(ncol(tmp)),sep='')
					tmp				<- as.data.table(tmp)
					tmp[, rep:=seq.int(1,nrow(tmp))]
					tmp[, summary:=paste('summary',s,sep='')]
					tmp
				})
	df	<- do.call('rbind', df)	
	df	<- melt(df, id.vars=c('summary','rep'), variable.name='design_pt', value.name='simu_data')
	set(df, NULL, 'summary', df[, factor(summary)])
	
	#	
	df.normal	<- df[, list(mu= mean(simu_data), sd=sd(simu_data), min=min(simu_data), max=max(simu_data)), by=c('summary', 'design_pt')]
	df.normal	<- df.normal[, list(x= seq(min, max, by=(max-min)/1e3), y=dnorm(seq(min, max, by=(max-min)/1e3), mean=mu, sd=sd) ), by=c('summary', 'design_pt')]	
	#dfe			<- subset(df, design_pt=='design1' & summary%in%c('summary1','summary2'))
	#df.normal2	<- subset(df.normal, design_pt=='design1' & summary%in%c('summary1','summary2'))
	ggplot(df, aes(x=simu_data, group=summary)) + 
			geom_line(aes(y = ..density.., colour = 'Empirical'), stat = 'density') +  
			geom_line(data=df.normal, aes(x=x, y=y, colour='Normal MLE')) + geom_histogram(aes(y=..density..), alpha=0.4) +
			scale_colour_brewer(palette='Set1', name='Density fits') + 
			facet_grid(design_pt~summary, scales='free')		
	outfile	<- '141007_ibmoutput_normal1e5.pdf'
	file	<- paste(indir, outfile, sep='/')
	ggsave(file=file, w=18*4, h=14*4, limitsize=FALSE)	
}

yia.check.correlation<- function()
{
	indir	<- '~/Dropbox (Personal)/abc/Ioannis_LondonSchool'
	infile	<- '141003_ibmoutput.RData'
	file	<- paste(indir, infile, sep='/')
	z		<- load(file)
	
	#	Outputs is an 10000 x 14 x 18 array, containing the 10000 reps, 14 different design point and 18 outputs. 
	#	X is a 22 x 14 matrix containing the 14 22-dimensional design points
	
	#	convert into data.table
	df	<- lapply( 1:dim(Outputs)[3], function(s)
			{
				tmp				<- Outputs[,,s]
				colnames(tmp)	<- paste('design',seq_len(ncol(tmp)),sep='')
				tmp				<- as.data.table(tmp)
				tmp[, rep:=seq.int(1,nrow(tmp))]
				tmp[, summary:=paste('summary',s,sep='')]
				tmp
			})
	df	<- do.call('rbind', df)	
	df	<- melt(df, id.vars=c('summary','rep'), variable.name='design_pt', value.name='simu_data')
	set(df, NULL, 'summary', df[, factor(summary, levels=paste('summary',1:18,sep=''), labels=paste('summary',1:18,sep=''))])
	
	require(GGally)
	for(d in df[, unique(as.character(design_pt))])
	{
		#print(d)
		df2		<- subset(df, design_pt==d)
		df2		<- dcast.data.table(df2, design_pt+rep~summary, value.var='simu_data')
		#print(df2)
		outfile	<- paste('141007_ibmoutput_corr1e5_',d,'.pdf', sep='')
		file	<- paste(indir, outfile, sep='/')
		print(file)
		pdf(file=file, w=18*2, h=18*2)
		print( ggpairs(data=df2, columns=3:ncol(df2), title=d, lower = list(continuous = "density", combo = "box"), upper=list(params=list(size=7)), params=list(labelSize=5, gridLabelSize=5)) )
		dev.off()		
	}
}

yia.check.information_vs_noise<- function()
{
	indir	<- '~/Dropbox (Personal)/abc/Ioannis_LondonSchool'
	infile	<- '141003_ibmoutput.RData'
	file	<- paste(indir, infile, sep='/')
	z		<- load(file)
	
	#	Outputs is an 10000 x 14 x 18 array, containing the 10000 reps, 14 different design point and 18 outputs. 
	#	X is a 22 x 14 matrix containing the 14 22-dimensional design points
	
	#	convert into data.table
	df	<- lapply( 1:dim(Outputs)[3], function(s)
			{
				tmp				<- Outputs[,,s]
				colnames(tmp)	<- paste('design',seq_len(ncol(tmp)),sep='')
				tmp				<- as.data.table(tmp)
				tmp[, rep:=seq.int(1,nrow(tmp))]
				tmp[, summary:=paste('summary',s,sep='')]
				tmp
			})
	df	<- do.call('rbind', df)	
	df	<- melt(df, id.vars=c('summary','rep'), variable.name='design_pt', value.name='simu_data')
	set(df, NULL, 'summary', df[, factor(summary, levels=paste('summary',1:18,sep=''), labels=paste('summary',1:18,sep=''))])
	
	ggplot(df, aes(x=design_pt, y=simu_data, colour=summary)) + geom_boxplot() + facet_grid(summary~., scales='free') 
	outfile	<- paste('141007_ibmoutput_change1e5.pdf', sep='')
	file	<- paste(indir, outfile, sep='/')
	ggsave(file=file, w=10, h=14*3)
	
}