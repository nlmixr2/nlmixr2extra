
.iniS3 <- function() {
  rxode2::.s3register("rxode2::rxUiGet", "linearizeError")
}
.onLoad <- function(libname, pkgname) {
  .iniS3()
}

.onAttach <- function(libname, pkgname) {
  .iniS3()
}
