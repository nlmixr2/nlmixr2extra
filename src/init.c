#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

SEXP _nlmixr2extra_preCondInv(SEXP Rin);

static const R_CallMethodDef CallEntries[] = {
  {"_nlmixr2extra_preCondInv", (DL_FUNC) &_nlmixr2extra_preCondInv, 1},
  {NULL, NULL, 0}
};

void R_init_nlmixr2extra(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, TRUE);
  R_forceSymbols(dll,FALSE);
}

void R_unload_nlmixr2extra(DllInfo *info){
}

