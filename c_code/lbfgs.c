/* #########################################
 * Luis Manuel Román García
 * luis.roangarci@gmail.com
 * #########################################
 *
 * -----------------------------------------
 * Functions for applying LBFGS
 * algorithm. Ideas from
 * José Luis Morales's, ITAM 2015 LBFGS
 * implementation have been taken
 * as a guide for this script.
 * -----------------------------------------

 */

#include "trunc_newton.c"


/* -------------------------------------
 * FindH
 * IN
 * func:  Function to be minimized.
 * nRow:  Number of rows (length) of x.
 * N_max: Maximum number of CG iterates.
 * TOL:   Minimum size for gradient of f.
 * OUT
 * x: Local minimum of func.
 * -------------------------------------
 */
double * findH(double* grad, double** s, double** y,
               int nRow, int m, int k){
  double *q, *r, *alpha, *rho;
  double  constant, beta;
  int i, state;
  q = grad;
  // Initialize variables
  alpha = (double*) malloc(nRow * sizeof(double));
  rho   = (double*) malloc(nRow * sizeof(double));
  state = min(k, m);
  // Fill in rho
  for(i = 0; i < state; i++){
    rho[i] = 1 / (dotProd(y[i], s[i], nRow));
  }
  // First Loop
  for(i = (state - 1); i > 0; i--){
    alpha[i] = rho[i] * dotProd(s[i], q, nRow);
    q        = vSum(q, vProd(y[i], -alpha[i], nRow), nRow);
  }
  // r = H0
  if(k == 0){
    r = q;
  }else{
    constant = dotProd(s[state - 1],
                       y[state - 1], nRow) / dotProd(y[state - 1],
                                                     y[state - 1], nRow);
    r        = vProd(q, constant, nRow);
  }
  // Second Loop
  for(i = 0; i < state; i ++){
    beta  = rho[i] * dotProd(y[i], r, nRow);
    r     = vSum(r, vProd(s[i], (alpha[i] - beta), nRow), nRow);
  }
  // Memory release.
  free(alpha);
  free(rho);
  // Return result.
  return r;
}


/* -------------------------------------
 * UpdateSY
 * IN
 * func:  Function to be minimized.
 * nRow:  Number of rows (length) of x.
 * N_max: Maximum number of CG iterates.
 * TOL:   Minimum size for gradient of f.
 * OUT
 * x: Local minimum of func.
 * -------------------------------------
 */
void updateSY(double** s, double** y, double * s_new,
              double* y_new, int m, int k){
  int i;
  if(k < m){
    s[k] = s_new;
    y[k] = y_new;
  }else{
    for(i = 0; i < (m - 1); i++){
      s[i] = s[i + 1];
      y[i] = y[i + 1];
    }
    s[m - 1] = s_new;
    y[m - 1] = y_new;
  }
}


/* -------------------------------------
 * LBFGS
 * IN
 * func:  Function to be minimized.
 * nRow:  Number of rows (length) of x.
 * N_max: Maximum number of CG iterates.
 * TOL:   Minimum size for gradient of f.
 * OUT
 * x: Local minimum of func.
 * -------------------------------------
 */
double * LBFGS(double (* func)(double*, int),
               int nRow, int m, double TOL, int verbose){
  // Variable declaration.
  double **s, **y;
  double *x, *grad, *p, *x_new, *grad_new;
  double alpha, norm_grad0;
  int i, k, MAX_ITER;
  // Space allocation.
  x = (double *)malloc(nRow * sizeof(double));
  s = (double **)malloc((nRow * nRow) * sizeof(double));
  y = (double **)malloc((nRow * nRow) * sizeof(double));
  // Initialize x.
  for(i = 0; i < nRow; i++){
    x[i] = ((double) rand() / INT_MAX) ;
  }
  // Until Convergence or MAX_ITER.
  MAX_ITER = 1.5e4;
  grad     = gradCentralDiff(func, x, nRow);
  // Update s, y.
  k = 0;
  // Initial norm of gradient.
  norm_grad0 = norm(grad, nRow);
  while(norm(grad, nRow) > TOL*(1 + norm_grad0) && k < MAX_ITER){
    // p = -Hgrad(f)
    p        = vProd(findH(grad, s, y, nRow, m, k), -1, nRow);
    // Alpha that statifies Wolfe conditions.
    alpha    = backTrack(func, x, p, nRow, verbose);
    x_new    = vSum(x, vProd(p, alpha, nRow), nRow);
    grad_new = gradCentralDiff(func, x_new, nRow);
    // Update s, y.
    updateSY(s, y, vProd(p, alpha, nRow),
             vSum(grad_new, vProd(grad, -1, nRow), nRow), m, k);

    // ---------------- PRINT ------------------- //
    if(verbose){
      printf("\n ITER = %d; f(x) = %.10e ; "
             "||x|| = %.10e ; ||grad|| =  %.10e ; "
             "||p|| =  %.10e ; sTy =  %.10e ; "
             "alpha = %.10e",
             k,
             func(x, nRow),
             norm(x, nRow),
             norm(grad, nRow),
             norm(p, nRow),
             dotProd(s[(int)min(k , (m - 1))],
                     y[(int)min(k , (m - 1))], nRow),
             alpha);
    }
    // ---------------- PRINT ------------------- //y
    // Update k, x, grad.
    x    = x_new;
    grad = grad_new;
    k    = k + 1;
  }
  free(grad);
  free(s);
  free(y);
  return x;
}
