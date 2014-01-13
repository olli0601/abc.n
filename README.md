########################################################################################

	Calibration procedures for accurate ABC

 	The corresponding paper is now on arxiv http://arxiv.org/abs/1305.4283
	and is under review and not yet accepted

 	authors 	Oliver Ratmann oliver.ratmann@ic.ac.uk
      			Anton Camacho ntncmch@gmail.com	

########################################################################################

	Contributors:

	If you d like to add calibration routines for different equivalence tests to this library,
	please let us know and we are happy to give you write access to the github repository

#########################################################################################

	Installation instructions:


	To clone this R package from the git repository:
	git clone https://github.com/olli0601/abc.n.git
	If there are issues like "Permission denied (publickey)" , try 
	git clone http://github.com/olli0601/abc.n.git

  To build the R documentation: 
  run roxygen in the roxygen R package on the abc.n directory

	To build this R package:
	R CMD build pkg

	To install this R package:
	R CMD INSTALL abc.star_1.0-0.tar.gz


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
