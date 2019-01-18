#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP branchingTimesCall(SEXP edge, SEXP edgelength)
{
    double *p, *el;
    int *e;
    SEXP xx;
    PROTECT(edge = coerceVector(edge, INTSXP));
    PROTECT(edgelength = coerceVector(edgelength, REALSXP));
    int N = LENGTH(edge)/2; /* number of edges */
    int n = N/2 + 1; /* number of tips */
    PROTECT(xx = allocVector(REALSXP, n - 1));
    e = INTEGER(edge);
    el = REAL(edgelength);
    p = REAL(xx);
    memset(p, 0, (n - 1) * sizeof(double));

    for (int i = 0; i < N; i++) {
	if (e[i + N] <= n) continue;
	p[e[i + N] - n - 1] = p[e[i] - n - 1] + el[i];
    }
    double depth = p[e[N - 1] - n - 1] + el[N - 1];
    for (int i = 0; i < n - 1; i++) {
	p[i] *= -1;
	p[i] += depth;
    }
    UNPROTECT(3);
    return xx;
}

static R_CallMethodDef Call_entries[] = {
    {"branchingTimesCall", (DL_FUNC) &branchingTimesCall, 2},
    {NULL, NULL, 0}
};

void R_init_coalescentMCMC(DllInfo *info)
{
    R_registerRoutines(info, NULL, Call_entries, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}
