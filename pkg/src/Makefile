CWD= $$(pwd)

all:
	@echo "*********************************************** "
	@echo "********* MAKE DOCUMENTATION ****************** "
	@echo "*********************************************** "	
	make roxygen
	@echo "*********************************************** "
	@echo "********* BUILD SOURCE FILES ****************** "
	@echo "*********************************************** "		
	MAKEFLAGS="SHLIB_CXXLD=g++ CXX=g++" R CMD SHLIB -o nabc.so *.cpp
	@echo "*********************************************** "
	@echo "********* CLEANUP ***************************** "
	@echo "*********************************************** "			
	make cleanup
	
make roxygen:
	echo $(CWD)
	../misc/nabc.startme.R -exe=ROXYGENIZE -code.home=$(CWD)/..
	
make cleanup:
	rm -f *.o		