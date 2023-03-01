#include "header.h"
#include <R_ext/Rdynload.h>

static const R_CMethodDef cMethods[] = {
    {"gevlik",  (DL_FUNC) &gevlik,  4},
    {"gevlikt", (DL_FUNC) &gevlikt, 5},
    {"gpdlik",  (DL_FUNC) &gpdlik,  4},
    {"gpdlikt", (DL_FUNC) &gpdlikt, 5},
    {"pplik",  (DL_FUNC) &pplik,  7},
    {"pplikt", (DL_FUNC) &pplikt, 9},
    {"oslik",  (DL_FUNC) &oslik,  6},
    {"oslikt", (DL_FUNC) &oslikt, 8},

    {"dprior_norm",  (DL_FUNC) &dprior_norm,  5},
    {"dprior_loglognorm", (DL_FUNC) &dprior_loglognorm, 5},
    {"dprior_prob",  (DL_FUNC) &dprior_prob,  5},
    {"dprior_quant", (DL_FUNC) &dprior_quant, 6},
	
    {NULL, NULL, 0}
};	   

static const R_CallMethodDef callMethods[]  = {
  {"gibbs", (DL_FUNC) &gibbs, 7},
  {"gibbsmix", (DL_FUNC) &gibbs, 11},

  {NULL, NULL, 0}
};

void R_init_evdbayes(DllInfo *dll)
{
    R_registerRoutines(dll, cMethods, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
	R_forceSymbols(dll, TRUE);
}


