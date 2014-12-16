########################################################################################

	Calibration procedures for accurate ABC

 	The corresponding paper is now on arxiv http://arxiv.org/abs/1305.4283	

 	authors 	Oliver Ratmann oliver.ratmann@ic.ac.uk
      			Anton Camacho ntncmch@gmail.com	
			Sen Hu ethansen.hu@gmail.com
			Caroline Colijn c.colijn@imperial.ac.uk

########################################################################################

	Contributors:

	If you d like to add calibration routines for different equivalence tests to this library,
	please let us know and we are happy to give you write access to the github repository

#########################################################################################

	Installation instructions for UNIX / MACOS:

	To clone this R package from the git repository:
	git clone https://github.com/olli0601/abc.star.git
	If there are issues like "Permission denied (publickey)" , try 
	git clone http://github.com/olli0601/abc.star.git

	To build this R package:
	Go to the folder containing abc.star and type into iterm or a similar bash:
	R CMD build abc.star

	To install this R package, type then
	R CMD INSTALL abc.star_1.0-0.tar.gz

	To work with this R package:
	fire up R, and type library(help=abc.star)

#########################################################################################

	Content:

	Calibration routines for dispersion equivalence of normally distributed summaries
	start with
	chisqstretch.*

	Calibration routines for location equivalence of normally distributed summaries
	start with
	nabc.mutost.onesample.*

	Calibration routines for asymptotic equivalence of autocorrelation
	start with
	ma.*

	Various helper functions to reconstruct and display link functions and level sets
	
 
#########################################################################################
