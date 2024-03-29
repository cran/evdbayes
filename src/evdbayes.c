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

void gevlik(double *data, int *n, double *par, double *dns)
{
  int i;
  double *dvec, eps;

  dvec = (double *)R_alloc(*n, sizeof(double));
  eps = R_pow(DBL_EPSILON, 0.3);


  for(i=0;i<*n;i++)  {
    data[i] = (data[i] - par[0]) / par[1];
    if(fabs(par[2]) <= eps)
      dvec[i] = log(1 / par[1]) - data[i] - exp(-data[i]);
    else {
      data[i] = 1 + par[2] * data[i];
      if(data[i] <= 0) {
        *dns = R_NegInf;
        return;
      }
      dvec[i] = log(1 / par[1]) - R_pow(data[i], -1 / par[2]) -
                (1 / par[2] + 1) * log(data[i]);
    }
  }

  for(i=0;i<*n;i++)
    *dns = *dns + dvec[i];
}

void gevlikt(double *data, int *n, double *par, double *trend, 
             double *dns)
{
  int i;
  double *loc, *dvec, eps;

  loc = (double *)R_alloc(*n, sizeof(double));
  dvec = (double *)R_alloc(*n, sizeof(double));
  eps = R_pow(DBL_EPSILON, 0.3);

  for(i=0;i<*n;i++) loc[i] = par[0] + trend[i] * par[3];

  for(i=0;i<*n;i++)  {
    data[i] = (data[i] - loc[i]) / par[1];
    if(fabs(par[2]) <= eps)
      dvec[i] = log(1 / par[1]) - data[i] - exp(-data[i]);
    else {
      data[i] = 1 + par[2] * data[i];
      if(data[i] <= 0) {
        *dns = R_NegInf;
        return;
      }
      dvec[i] = log(1 / par[1]) - R_pow(data[i], -1 / par[2]) -
                (1 / par[2] + 1) * log(data[i]);
    }
  }

  for(i=0;i<*n;i++)
    *dns = *dns + dvec[i];
}

void gpdlik(double *data, int *n, double *par, double *dns)
{
  int i;
  double *dvec, eps;
  
  dvec = (double *)R_alloc(*n, sizeof(double));
  eps = R_pow(DBL_EPSILON, 0.3);

  if(par[1] <= 0) {
     *dns = R_NegInf;
     return;
  }

  for(i=0;i<*n;i++)  {
    data[i] = (data[i] - par[0]) / par[1];
    if (data[i] <= 0) {
      *dns = R_NegInf;
      return;
    }
    if(fabs(par[2]) <= eps) 
      dvec[i] = log(1 / par[1]) - data[i];
    else {
      data[i] = 1 + par[2] * data[i];
      if(data[i] <= 0) {
	*dns = R_NegInf;
	return;
      }
      dvec[i] = log(1 / par[1]) - (1 / par[2] + 1) * log(data[i]);
    }
  }
  
  for(i=0;i<*n;i++) 
    *dns = *dns + dvec[i];
}

void gpdlikt(double *data, int *n, double *par, double *trend,
	     double *dns)
{
  int i;
  double *loc, *dvec, eps;
  
  loc = (double *)R_alloc(*n, sizeof(double));
  dvec = (double *)R_alloc(*n, sizeof(double));
  eps = R_pow(DBL_EPSILON, 0.3);

  for(i=0;i<*n;i++) loc[i] = par[0] + trend[i] * par[3];

  if(par[1] <= 0) {
     *dns = R_NegInf;
     return;
  }

  for(i=0;i<*n;i++)  {
    data[i] = (data[i] - loc[i]) / par[1];
    if (data[i] <= 0) {
      *dns = R_NegInf;
      return;
    }
    if(fabs(par[2]) <= eps) 
      dvec[i] = log(1 / par[1]) - data[i];
    else {
      data[i] = 1 + par[2] * data[i];
      if(data[i] <= 0) {
	*dns = R_NegInf;
	return;
      }
      dvec[i] = log(1 / par[1]) - (1 / par[2] + 1) * log(data[i]);
    }
  }
  
  for(i=0;i<*n;i++) 
    *dns = *dns + dvec[i];
}

void oslik(double *data, double *thresh, int *n, int *m, 
           double *par, double *dns)
{
  int i;
  double *dvec1, *dvec2, eps;

  dvec1 = (double *)R_alloc(*n, sizeof(double));
  dvec2 = (double *)R_alloc(*m, sizeof(double));
  eps = R_pow(DBL_EPSILON, 0.3);

  for(i=0;i<*n;i++)  {
    data[i] = (data[i] - par[0]) / par[1];
    if(fabs(par[2]) <= eps)
      dvec1[i] = - data[i];
    else {
      data[i] = 1 + par[2] * data[i];
      if(data[i] <= 0) {
        *dns = R_NegInf;
        return;
      }
      dvec1[i] = - (1 / par[2] + 1) * log(data[i]);
    }
  }

  for(i=0;i<*n;i++)
    *dns = *dns + dvec1[i];

  for(i=0;i<*m;i++)  {
    thresh[i] = (thresh[i] - par[0]) / par[1];
    if(fabs(par[2]) <= eps)
      dvec2[i] = - exp(-thresh[i]);
    else {
      thresh[i] = 1 + par[2] * thresh[i];
      dvec2[i] = - R_pow(thresh[i], -1 / par[2]);
    }
  }

  for(i=0;i<*m;i++)
    *dns = *dns + dvec2[i];

  *dns = *dns - *n * log(par[1]);
}

void oslikt(double *data, double *thresh, int *n, int *m, int *r,
            double *par, double *trend, double *dns)
{
  int i, k=0;
  double *loc, *dvec1, *dvec2, eps;

  dvec1 = (double *)R_alloc(*n, sizeof(double));
  dvec2 = (double *)R_alloc(*m, sizeof(double));
  loc = (double *)R_alloc(*m, sizeof(double));
  eps = R_pow(DBL_EPSILON, 0.3);

  for(i=0;i<*m;i++) loc[i] = par[0] + trend[i] * par[3];  

  for(i=0;i<*n;i++)  {
    data[i] = (data[i] - loc[k]) / par[1];
    if(fabs(par[2]) <= eps)
      dvec1[i] = - data[i];
    else {
      data[i] = 1 + par[2] * data[i];
      if(data[i] <= 0) {
        *dns = R_NegInf;
        return;
      }
      dvec1[i] = - (1 / par[2] + 1) * log(data[i]);
    }
    if(i == r[k]-1) k = k + 1;
  }

  for(i=0;i<*n;i++)
    *dns = *dns + dvec1[i];

  for(i=0;i<*m;i++)  {
    thresh[i] = (thresh[i] - loc[i]) / par[1];
    if(fabs(par[2]) <= eps)
      dvec2[i] = - exp(-thresh[i]);
    else {
      thresh[i] = 1 + par[2] * thresh[i];
      dvec2[i] = - R_pow(thresh[i], -1 / par[2]);
    }
  }

  for(i=0;i<*m;i++)
    *dns = *dns + dvec2[i];

  *dns = *dns - *n * log(par[1]);
}

void pplik(double *data, int *nh, double *par, double *thresh, 
           int *n, double *noy, double *dns)
{
  int i;
  double *dvec1, *dvec2, eps;

  dvec1 = (double *)R_alloc(*nh, sizeof(double));
  dvec2 = (double *)R_alloc(*n, sizeof(double));
  eps = R_pow(DBL_EPSILON, 0.3);
    
  for(i=0;i<*nh;i++) {
    data[i] = (data[i] - par[0]) / par[1];
    if(fabs(par[2]) <= eps) 
      dvec1[i] = log(1 / par[1]) - data[i];
    else {
      data[i] = 1 + par[2] * data[i];
      if(data[i] <= 0) {
        *dns = R_NegInf;
        return;
      }
      dvec1[i] = log(1 / par[1])  - (1 / par[2] + 1) * log(data[i]);
    }
  }

  for(i=0;i<*n;i++) {
    thresh[i] = (thresh[i] - par[0]) / par[1];
    if(fabs(par[2]) <= eps)
      dvec2[i] = - *noy / *n * exp(-thresh[i]);
    else {
      thresh[i] = 1 + par[2] * thresh[i];
      if(thresh[i] <= 0 && par[2] > 0) {
        *dns = R_NegInf;
        return;
      }
      if(thresh[i] <= 0 && par[2] < 0) dvec2[i] = 0;
      else dvec2[i] = - *noy / *n * R_pow(thresh[i], -1 / par[2]);
    }
  }

  for(i=0;i<*nh;i++)
    *dns = *dns + dvec1[i];
  for(i=0;i<*n;i++)
    *dns = *dns + dvec2[i];
}

void pplikt(double *data, int *nh, double *par, double *thresh, 
            int *n, double *noy,  double *trend, double *htrend,
            double *dns)
{
  int i;
  double *dvec1, *dvec2, eps, *loc, *hloc;

  dvec1 = (double *)R_alloc(*nh, sizeof(double));
  dvec2 = (double *)R_alloc(*n, sizeof(double));
  hloc = (double *)R_alloc(*nh, sizeof(double));
  loc = (double *)R_alloc(*n, sizeof(double));
  eps = R_pow(DBL_EPSILON, 0.3);
  
  for(i=0;i<*nh;i++) {
    hloc[i] = par[0] + htrend[i] * par[3];
    data[i] = (data[i] - hloc[i]) / par[1];
    if(fabs(par[2]) <= eps) 
      dvec1[i] = log(1 / par[1]) - data[i];
    else {
      data[i] = 1 + par[2] * data[i];
      if(data[i] <= 0) {
        *dns = R_NegInf;
        return;
      }
      dvec1[i] = log(1 / par[1])  - (1 / par[2] + 1) * log(data[i]);
    }
  }

  for(i=0;i<*n;i++) {
    loc[i] = par[0] + trend[i] * par[3];
    thresh[i] = (thresh[i] - loc[i]) / par[1];
    if(fabs(par[2]) <= eps)
      dvec2[i] = - *noy / *n * exp(-thresh[i]);
    else {
      thresh[i] = 1 + par[2] * thresh[i];
      if(thresh[i] <= 0) {
        *dns = R_NegInf;
        return;
      }
      dvec2[i] = - *noy / *n * R_pow(thresh[i], -1 / par[2]);
    }
  }

  for(i=0;i<*nh;i++)
    *dns = *dns + dvec1[i];
  for(i=0;i<*n;i++)
    *dns = *dns + dvec2[i];
}

void dprior_norm(double *par, double *mean, double *icov, 
                 double *trendsd, double *dns)
{
  double cpar[3], ld;
  int i;

  par[1] = log(par[1]);
  for(i=0;i<3;i++) cpar[i] = par[i] - mean[i];
  ld = icov[0] * R_pow_di(cpar[0], 2) + icov[3] * R_pow_di(cpar[1], 2) +
    icov[5] * R_pow_di(cpar[2], 2) + 
    2 * icov[1] * cpar[0] * cpar[1] + 2 * icov[2] * cpar[0] * cpar[2] +
    2 * icov[4] * cpar[1] * cpar[2];
  ld = -ld / 2 - par[1];
  if(*trendsd != 0) ld = ld - R_pow_di(par[3] / *trendsd, 2) / 2;
  *dns = ld;
}

void dprior_loglognorm(double *par, double *mean, double *icov, 
		       double *trendsd, double *dns)
{
  double cpar[3], ld;
  int i;
  
  par[0] = log(par[0]);
  par[1] = log(par[1]);
  for(i=0;i<3;i++) cpar[i] = par[i] - mean[i];
  ld = icov[0] * R_pow_di(cpar[0], 2) + icov[3] * R_pow_di(cpar[1], 2) +
    icov[5] * R_pow_di(cpar[2], 2) + 
    2 * icov[1] * cpar[0] * cpar[1] + 2 * icov[2] * cpar[0] * cpar[2] +
    2 * icov[4] * cpar[1] * cpar[2];
  ld = -ld / 2 - par[1] - par[0];
  if(*trendsd != 0) ld = ld - R_pow_di(par[3] / *trendsd, 2) / 2;
  *dns = ld;
}

void dprior_prob(double *par, double *quant, double *alpha, 
                 double *trendsd, double *dns)
{
  int i;
  double z[3], pr[5], pd[4], eps, lj, ld;

  eps = R_pow(DBL_EPSILON, 0.3);
  pr[0] = 1 ; pr[4] = 0; 
  
  for(i=0;i<3;i++) {
    if(fabs(par[2]) <= eps) 
      z[i] = exp(-(quant[i] - par[0])/par[1]);
    else {
      z[i] = 1 + par[2] * (quant[i] - par[0]) / par[1];
      if(z[i] <= 0) {
        *dns = R_NegInf;
        return;
      }
      z[i] = R_pow(z[i], -1 / par[2]);
    }
    pr[i+1] = z[i];
    if(pr[i+1] > 1e-7) pr[i+1] = 1 - exp(-pr[i+1]);
  }

  for(i=0;i<4;i++) {
    pd[i] = pr[i] - pr[i+1];
    if(pd[i] <= 0) {
      *dns = R_NegInf;
      return;
    }
  }
  
  if(fabs(par[2]) <= eps) {
    lj = quant[0] * quant[1] * (quant[0] - quant[1]) -
      quant[0] * quant[2] * (quant[0] - quant[2]) +
      quant[1] * quant[2] * (quant[1] - quant[2]);
    lj = log(fabs(lj)) - 5 * log(par[1]) - log(2);
  }
  else {
    lj = R_pow(z[0]*z[1], -par[2]) * log(z[1]/z[0]) -
      R_pow(z[0]*z[2], -par[2]) * log(z[2]/z[0]) +
      R_pow(z[1]*z[2], -par[2]) * log(z[2]/z[1]);
    lj = log(fabs(lj)) - 2 * log(par[1]) - log(R_pow_di(par[2], 2));
  }
  lj = lj + (1 + par[2]) * (log(z[0]) + log(z[1]) + log(z[2])) -
    (z[0] + z[1] + z[2]);

  ld = (alpha[0]-1) * log(pd[0]) + (alpha[1]-1) * log(pd[1]) +
       (alpha[2]-1) * log(pd[2]) + (alpha[3]-1) * log(pd[3]);
  ld = ld + lj;
  if(*trendsd != 0) ld = ld - R_pow_di(par[3] / *trendsd, 2) / 2;
  *dns = ld;
}

void dprior_quant(double *par, double *prob, double *shape, 
                  double *scale, double *trendsd, double *dns)
{
  int i;
  double z[3], quant[4], qd[4], eps, lj, ld;

  eps = R_pow(DBL_EPSILON, 0.3);
  quant[0] = 0; 
  
  for(i=0;i<3;i++) {
    if(fabs(par[2]) <= eps) {
      z[i] = log(-log(1 - prob[i])); 
      quant[i+1] = par[0] - par[1] * z[i];
    }
    else {
      z[i] = -log(1 - prob[i]);
      quant[i+1] = par[0] + par[1] * (R_pow(z[i], -par[2]) - 1) / par[2];
    }
  }

  for(i=0;i<3;i++) {
    qd[i] = quant[i+1] - quant[i];
    if(qd[i] <= 0) {
      *dns = R_NegInf;
      return;
    }
  }
  
  if(fabs(par[2]) <= eps) {
    lj = z[0] * z[1] * (z[1] - z[0]) -
      z[0] * z[2] * (z[2] - z[0]) +
      z[1] * z[2] * (z[2] - z[1]);
    lj = log(fabs(lj)) + log(par[1]) - log(2);
  }
  else {
    lj = R_pow(z[0]*z[1], -par[2]) * log(z[1]/z[0]) -
      R_pow(z[0]*z[2], -par[2]) * log(z[2]/z[0]) +
      R_pow(z[1]*z[2], -par[2]) * log(z[2]/z[1]);
    lj = log(fabs(lj)) + log(par[1]) - log(R_pow_di(par[2], 2));
  }

  ld = (shape[0]-1) * log(qd[0]) - qd[0]/scale[0] +
    (shape[1]-1) * log(qd[1]) - qd[1]/scale[1] +
    (shape[2]-1) * log(qd[2]) - qd[2]/scale[2];

  ld = ld + lj;
  if(*trendsd != 0) ld = ld - R_pow_di(par[3] / *trendsd, 2) / 2;
  *dns = ld;
}

SEXP gibbs(SEXP n, SEXP np, SEXP thin,  
           SEXP init, SEXP propsd, SEXP f, SEXP rho)
{
  int i,j,k,nr;
  int nn = INTEGER(n)[0], nnp = INTEGER(np)[0], thinn = INTEGER(thin)[0];
  double prop, prop_ratio, acc_prob, post_ratio;
  double *crow, *prow;
  SEXP ans, nacc, nex, mc, current, dpst_lower, dpst_upper;

  nr = 1 + ftrunc(nn/thinn);
  crow = (double *)R_alloc(nnp, sizeof(double));
  prow = (double *)R_alloc(nnp, sizeof(double));
  PROTECT(current = allocVector(REALSXP, nnp));
  PROTECT(nacc = allocVector(REALSXP, nnp));
  PROTECT(nex = allocVector(REALSXP, nnp));
  PROTECT(mc = allocVector(REALSXP, nr * nnp));
  PROTECT(ans = allocVector(VECSXP, 3));
  PROTECT(dpst_lower = allocVector(REALSXP, 1));
  PROTECT(dpst_upper = allocVector(REALSXP, 1));

  for(i=0;i<nnp;i++) {
    prow[i] = REAL(init)[i];
    REAL(mc)[i] = REAL(init)[i];
    REAL(nex)[i] = REAL(nacc)[i] = 0.0;
  }

  RANDIN;
  for(i=0;i<nn;i++) {
    for(j=0;j<nnp;j++) {
      if(j==1) {
        prop = rlnorm(log(prow[1]), REAL(propsd)[1]);
        prop_ratio = prop / prow[1];
      }
      else {
        prop = rnorm(prow[j], REAL(propsd)[j]);
        prop_ratio = 1;
      }
      for(k=0;k<nnp;k++) {
	if(k < j) REAL(current)[k] = crow[k];
        else REAL(current)[k] = prow[k];
      }
      defineVar(install("x"), current, rho);
      dpst_lower = PROTECT(eval(f, rho));

      if (TYPEOF(dpst_lower) != REALSXP)
	error("non-numeric result");

      REAL(current)[j] = prop;

      defineVar(install("x"), current, rho);
      dpst_upper = PROTECT(eval(f, rho));

      if (TYPEOF(dpst_upper) != REALSXP)
	error("non-numeric result");

      post_ratio = exp(REAL(dpst_upper)[0] - REAL(dpst_lower)[0]);

      if(!R_FINITE(REAL(dpst_upper)[0]))
        REAL(nex)[j] = REAL(nex)[j] + 1;

      UNPROTECT(2);

      acc_prob = fmin2(1, prop_ratio * post_ratio);
      if(R_IsNaN(acc_prob)) {
        acc_prob = 0;
        warning("NaN returned for posterior density");
      }
      if(runif(0, 1) < acc_prob) {
        crow[j] = prop;
        REAL(nacc)[j] = REAL(nacc)[j] + 1;
      }
      else crow[j] = prow[j];
    }
    if(((i+1) % thinn) == 0)
      for(j=0;j<nnp;j++) REAL(mc)[(i+1)/thinn * nnp + j] = crow[j];
    for(j=0;j<nnp;j++) prow[j] = crow[j];
  }

  RANDOUT;
  SET_VECTOR_ELT(ans, 0, mc);
  SET_VECTOR_ELT(ans, 1, nacc);
  SET_VECTOR_ELT(ans, 2, nex);
  UNPROTECT(7);
  return(ans);
}



