#define STRICT_R_HEADER
#include <RcppArmadillo.h>
// NLS not available on all platforms; no translations provided
#define _(String) (String)

using namespace Rcpp;
// #ifdef ENABLE_NLS
// #include <libintl.h>
// #define _(String) dgettext ("nlmixr2extra", String)
// /* replace pkg as appropriate */
// #else
#define _(String) (String)
// #endif

//' Calculate the inverse preconditioning matrix
//'
//' @param Rin The R matrix input
//'
//' @return The inverse preconditioning matrix
//' @noRd
//[[Rcpp::export]]
SEXP preCondInv(SEXP Rin) {
  // Assumes Rin is symmetric
  arma::vec eigval;
  arma::mat eigvec;
  arma::mat R = as<arma::mat>(Rin);
  bool success = eig_sym(eigval, eigvec, R);
  if (success){
    // Guard against singular/near-singular matrices: a zero (or near-zero)
    // eigenvalue causes 1/|lambda| = Inf which silently corrupts the
    // preconditioner (Aoki 2016 eq. 15).
    double eigTol = 1e-10;
    arma::uvec zeroEigs = arma::find(arma::abs(eigval) < eigTol);
    if (zeroEigs.n_elem > 0) {
      Rcpp::stop(_("matrix is singular or near-singular (has zero or near-zero eigenvalues); cannot calculate the preconditioning matrix"));
    }
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
