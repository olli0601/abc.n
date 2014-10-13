.First.lib <- function(lib, pkg)  
{
	library.dynam("abc.star", pkg, lib)
	packageStartupMessage("This is abc.star ", utils::packageDescription("abc.star", field="Version"), appendLF = TRUE)
}