abc.n
=====
test procedures for ABC based on n summary values

to build this package:

* R CMD build pkg: builds the package (generates an archive called foo.tar.gz).
* R CMD build pkg --compact-vignettes --resave-data: same, making the package as small as possible.
* R CMD check abc.n_1.0-0.gz --as-cran: runs package quality checks similar to the checks made by CRAN; must be passed without error/warning before.

to install this R package:
* R CMD INSTALL abc.n_1.0-0.tar.gz

to re-build the R package:
* R CMD build pkg

to re-install the R package:
* R CMD INSTALL  abc.n_1.0-0.tar.gz
