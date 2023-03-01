#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

#define RANDIN GetRNGstate()
#define RANDOUT PutRNGstate()

void gevlik(double *data, int *n, double *par, double *dns);
void gevlikt(double *data, int *n, double *par, double *trend, 
             double *dns);
void gpdlik(double *data, int *n, double *par, double *dns);
void gpdlikt(double *data, int *n, double *par, double *trend, 
             double *dns);
void pplik(double *data, int *nh, double *par, double *thresh, 
           int *n, double *noy, double *dns);
void pplikt(double *data, int *nh, double *par, double *thresh, 
            int *n, double *noy, double *trend, double *htrend, 
            double *dns);
void oslik(double *data, double *thresh, int *n, int *m, double *par, 
           double *dns);
void oslikt(double *data, double *thresh, int *n, int *m, int *r, 
            double *par, double *trend, double *dns);
void dprior_norm(double *par, double *mean, double *icov, 
                 double *trendsd, double *dns);
void dprior_loglognorm(double *par, double *mean, double *icov, 
		       double *trendsd, double *dns);
void dprior_prob(double *par, double *quant, double *alpha, 
                 double *trendsd, double *dns);
void dprior_quant(double *par, double *prob, double *shape, 
                  double *scale, double *trendsd, double *dns);
SEXP gibbs(SEXP n, SEXP np, SEXP thin,  
           SEXP prow, SEXP propsd, SEXP f, SEXP rho);

