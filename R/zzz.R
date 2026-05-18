.onLoad <- function(libname, pkgname) {
  Rcpp::loadModule("stan_fit4barcode_mod", what = TRUE)
  Rcpp::loadModule("stan_fit4confluency_mod", what = TRUE)
  Rcpp::loadModule("stan_fit4tumor_volume_mod", what = TRUE)
}