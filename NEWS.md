# nlmixr2extra 2.0.8

* Use `assignInMyNamespace()` instead of using the global assignment
  operator for the horseshoe prior
  
* Be specific in version requirements (as requested by CRAN checks)

# nlmixr2extra 2.0.7

* Fix `cli` issues with the new `cli` 3.4+ release that will allow
  bootstrapping to run again (before `cli` would error, this fixes the
  `donttest` issues on CRAN).
  
* Fixed step-wise covariate selection to work a bit better with the
  updated UI, thanks to Vishal Sarsani
  
* Added lasso covariate selection (thanks to Vishal Sarsani)

* Added horseshoe prior covarite selecion (thanks to Vishal Sarsani)

* Added a `NEWS.md` file to track changes to the package.
