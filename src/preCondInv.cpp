#define STRICT_R_HEADER
#include <RcppArmadillo.h>

using namespace Rcpp;

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("nlmixr2extra", String)
/* replace pkg as appropriate */
#else
#define _(String) (String)
#endif

//' Calculate the inverse preconditioning matrix
//'
//' @param Rin The R matrix input
//'
//' @return The inverse preconditioning matrix
//'
//[[Rcpp::export]]
SEXP preCondInv(SEXP Rin) {
  // Assumes Rin is symmetric
  arma::vec eigval;
  arma::mat eigvec;
  arma::mat R = as<arma::mat>(Rin);
  bool success = eig_sym(eigval, eigvec, R);
  if (success){
    // Now calculate the norm
    arma::mat eignorm = normalise(eigvec);
    arma::mat v12 = diagmat(1/(abs(eigval)));
    // This is equation 15 in the Aoki 2016 paper
    R = eignorm*v12;
    SEXP out = wrap(R);
    Rf_setAttrib(out, R_DimNamesSymbol, Rf_getAttrib(Rin, R_DimNamesSymbol));
    return out;
  }
  Rcpp::stop("cannot calculate the eigenvectors/eigenvalues required for preconditioning");
  return R_NilValue;
}
