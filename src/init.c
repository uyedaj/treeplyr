#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP geiger_descendants(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"geiger_descendants", (DL_FUNC) &geiger_descendants, 1},
    {NULL, NULL, 0}
};

void R_init_treeplyr(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
