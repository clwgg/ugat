#include <math.h>
#include <stdio.h>

#include "gsl/gsl_vector.h"
#include "gsl/gsl_multifit_nlin.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_cdf.h"

typedef struct my_data {

  double *x;
  double *y;
  int n;

} my_data;


int my_exp(const gsl_vector *x, void *data, gsl_vector *f)
{

  int i;
  double *xval, *yval;

  /* need to type-cast from void to avoid warning
   * of incompatible types since its defined as void
   * in gsl. */
  xval = ((my_data*)data)->x;
  yval = ((my_data*)data)->y;
  int n = ((my_data*)data)->n;

  double lambda = gsl_vector_get(x, 0);
  double N = gsl_vector_get(x, 1);

  for (i = 0; i < n; i++) {

    double Yi = N * exp(-lambda*xval[i]);
    gsl_vector_set(f, i, Yi-yval[i]);

  }

  return 0;
}


int my_dfexp(const gsl_vector *x, void *data, gsl_matrix *J)
{

  int i;
  double *xval;

  /* need to type-cast from void to avoid warning
   * of incompatible types since its defined as void
   * in gsl. */
  xval = ((my_data*)data)->x;
  int n = ((my_data*)data)->n;

  double lambda = gsl_vector_get(x, 0);
  double N = gsl_vector_get(x, 1);

  for (i = 0; i < n; i++) {

    double e = exp(-lambda*xval[i]);
    gsl_matrix_set(J, i, 0, -xval[i] * N * e);
    gsl_matrix_set(J, i, 1, e);

  }

  return 0;
}


double run_fit(int freq[20][6])
{

  my_data data;
  double tmpx[20];
  double tmpy[20];
  data.x = tmpx;
  data.y = tmpy;
  data.n = 20;

  int i;
  for (i = 0; i < 20; i++) {
    data.x[i] = i + 1;
    data.y[i] = (double)freq[i][5] / freq[i][1];
  }

  double param[2] = {1.0, 0.1};
  int n_param = 2;
  gsl_vector_view x = gsl_vector_view_array(param, n_param);

  const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
  gsl_multifit_fdfsolver *s;
  s = gsl_multifit_fdfsolver_alloc(T, data.n, n_param);

  gsl_multifit_function_fdf f;

  f.f = &my_exp;
  f.df = my_dfexp;
  f.n = data.n;
  f.p = n_param;
  f.params = &data;

  gsl_multifit_fdfsolver_wset(s, &f, &x.vector, NULL);

  gsl_vector *res_f = gsl_multifit_fdfsolver_residual(s);

  int info;

  gsl_multifit_fdfsolver_driver(s, 100, 1e-8, 1e-8, 0.0, &info);

  double chi = gsl_blas_dnrm2(res_f);

  gsl_matrix *J = gsl_matrix_alloc(data.n, n_param);
  gsl_matrix *covar = gsl_matrix_alloc(n_param, n_param);

  gsl_multifit_fdfsolver_jac(s, J);
  gsl_multifit_covar(J, 0.0, covar);

  // calculate sigma for scaling of error
  gsl_vector *res_sq = gsl_vector_alloc(data.n);
  gsl_vector_memcpy(res_sq, res_f);
  gsl_vector_mul(res_sq, res_f);
  double ssr = gsl_blas_dasum(res_sq);
  gsl_vector_free(res_sq);

  double dof = data.n - n_param;
  double sigma = sqrt(ssr/dof);

  #define FIT(i) gsl_vector_get(s->x, i)
  #define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

  double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 

  // get t-value and one-sided p-value
  double lam_t = FIT(0) / (c*ERR(0) * sigma);
  double lam_p = gsl_cdf_tdist_Q(lam_t, dof);

  gsl_multifit_fdfsolver_free (s);
  gsl_matrix_free (covar);
  gsl_matrix_free (J);

  #undef FIT
  #undef ERR

  return lam_p;
}

