.onLoad <- function(libname, pkgname) { }

.onAttach <- function(...) {
  theLib <- dirname(system.file(package = "rethinking"))
  pkgdesc <- packageDescription("rethinking", lib.loc = theLib)
  builddate <- gsub(';.*$', '', pkgdesc$Packaged)
  packageStartupMessage(paste("rethinking (Version ", pkgdesc$Version, ", packaged: ", builddate, ")", sep = ""))
} 