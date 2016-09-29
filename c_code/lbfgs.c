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
 *
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
  double *q, *r;
  double rho, alpha, constant, beta, denom;
  int i, state;
  // Copy grad
  q = grad;
  // State
  state = min(k, m);
  /* -------------------------
   * Firts Loop
   * -------------------------
  */
  for(i = (state - 1); i >= 0; i--){
    rho   = 1/(1e-6 + dotProd(y[i], s[i], nRow)); // avoid dividing by numerical zero.
    alpha = rho * dotProd(s[i], q, nRow);
    q     = vSum(q, vProd(y[i], -alpha, nRow), nRow);
    //printf("\n||grad|| = %f     sTy = %f \n", norm(grad, nRow), dotProd(s[i], y[i], nRow));
  }
  // r = H'q
  if(k == 0){
    r = mProd(identity(nRow), q, nRow, nRow);
  }else{
    constant = dotProd(s[state - 1], y[state - 1], nRow) /
      dotProd(y[state - 1], y[state - 1], nRow);
    r        = mProd(vProd(identity(nRow),
                           constant,
                           nRow * nRow), q, nRow, nRow);
  }
  /* -------------------------
   * Second Loop
   * -------------------------
  */
  for(i = 0; i < state; i ++){
    rho  = 1 /(1e-6 + dotProd(y[i], s[i], nRow)); // avoid dividing by numerical zero.
   beta = rho * dotProd(y[i], r, nRow);
   r    = vSum(r, vProd(s[i], (alpha - beta), nRow), nRow);
  }
  // Memory release.
  free(q);
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
               int nRow, int m, double TOL){
  // Variable declaration.
  double **s, **y;
  double *x, *grad, *p, *x_new, *grad_new;
  double alpha;
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
  MAX_ITER = 1.5e3;
  grad     = gradCentralDiff(func, x, nRow);
  // Update s, y.
  k = 0;
  updateSY(s, y, x, grad, m, k); // With k = 0; s = x, y = grad(f)
  while(norm(grad, nRow) > TOL && k < MAX_ITER){
    // p = -Hgrad(f)
    p        = vProd(findH(grad, s, y, nRow, m, k), -1, nRow);
    // Alpha that statifies Wolfe conditions.
    alpha    = backTrack(func, x, p, nRow);
    // ---------------- PRINT ------------------- //
    printf("\n ITER = %d; f(x) = %f ; ||grad|| =  %.10e ; ||p|| =  %.10e ; sTy =  %.10e ; alpha = %f",
           k,
           func(x, nRow),
           norm(grad, nRow),
           norm(p, nRow),
           dotProd(s[k % m], y[k % m], nRow),
           alpha);
    // ---------------- PRINT ------------------- //y
    x_new    = vSum(x, vProd(p, alpha, nRow), nRow);
    grad_new = gradCentralDiff(func, x_new, nRow);
    // Update k.
    k = k + 1;
    // Update s, y.
    updateSY(s, y, vSum(x_new, vProd(x, -1, nRow), nRow),
             vSum(grad_new, vProd(grad, -1, nRow), nRow), m, k);
    // Update x, grad.
    x    = x_new;
    grad = grad_new;
  }
  return x;
}
