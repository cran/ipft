#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
  Check these declarations against the C/Fortran source code.
*/

  /* .Call calls */
extern SEXP ipft_ipfEuclidean(SEXP, SEXP);
extern SEXP ipft_ipfLGD(SEXP, SEXP, SEXP, SEXP);
extern SEXP ipft_ipfManhattan(SEXP, SEXP);
extern SEXP ipft_ipfNormDistance(SEXP, SEXP, SEXP);
extern SEXP ipft_ipfPLGD(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"ipft_ipfEuclidean",    (DL_FUNC) &ipft_ipfEuclidean,    2},
  {"ipft_ipfLGD",          (DL_FUNC) &ipft_ipfLGD,          4},
  {"ipft_ipfManhattan",    (DL_FUNC) &ipft_ipfManhattan,    2},
  {"ipft_ipfNormDistance", (DL_FUNC) &ipft_ipfNormDistance, 3},
  {"ipft_ipfPLGD",         (DL_FUNC) &ipft_ipfPLGD,         6},
  {NULL, NULL, 0}
};

void R_init_ipft(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
