/* #########################################
 * Luis Manuel Román García
 * luis.roangarci@gmail.com
 * #########################################
 *
 * -----------------------------------------
 * Functions for implementing SLM-LBFGS
 * -----------------------------------------

 */

#include "lbfgs.c"

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
double * findHSLM(double (*func)(double*, int), double *x,
                  double** s, double** y,
                  int nRow, int m, int k, int N_max){
  double *q, *r, *alpha, *rho, *Bd, *r_cg, *d, *z, *r_new;
  double beta, beta_cg, alpha_cg, epsilon;
  int i, state, j;
  // Calculate gradient.
  r = q = gradCentralDiff(func, x, nRow);
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
  /*
   * -----------------------------------
   * ########### CG Iteration ##########
   * -----------------------------------
   * Outputs: r
   */
  // Initialize: epsilon, d, r_cg, z
  epsilon = min(.5, sqrt(norm(q, nRow))) * norm(q, nRow);
  d    = vProd(q,  1, nRow);
  r_cg = vProd(q,  1, nRow);
  z    = vProd(q,  0, nRow);
  for(j = 0; j < N_max; j++){
    Bd = hessCentralDiff(func, x, d, nRow);
    // Check if d'Bd <= 0 i.e. d is a descent direction.
    if(dotProd(d, Bd, nRow) <= 0){
      if(j == 0){
        r = d;
        break;
      }else{
        r = z;
        break;
      }
    }
    // alpha_j = rj'rj/d_j'Bd_j
    alpha_cg = dotProd(r_cg, r_cg, nRow) / dotProd(d, Bd, nRow);
    // z_{j+1} = z_j + alpha_j*d_j
    z        = vSum(z, vProd(d, alpha_cg, nRow), nRow);
    // r_{j+1} = r_j + alpha_j*B_kd_j
    r_new = vSum(r_cg, vProd(Bd, alpha_cg, nRow), nRow);
    if(norm(r_new, nRow) < epsilon){
      r = z;
      break;
    }
    // Update beta, d, r_cg.
    beta_cg = dotProd(r_new, r_new, nRow) / dotProd(r_cg, r_cg, nRow);
    d       = vSum(vProd(r_new, -1, nRow),
                   vProd(d, beta_cg, nRow), nRow);
    r_cg    = r_new;
  }
  /* -----------------------------------
   * ######### CG Iteration End ########
   * -----------------------------------
   */

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
double * SLM_LBFGS(double (* func)(double*, int),
                   int nRow, int m, double TOL, int N_max, int verbose){
  // Variable declaration.
  double **s, **y;
  double *x, *grad, *p, *x_new, *grad_new;
  double alpha, norm_grad0;
  int i, k, MAX_ITER;
  // Space allocation.
  x = (double *)malloc(nRow * sizeof(double));
  s = (double **)malloc((nRow*m) * sizeof(double));
  y = (double **)malloc((nRow*m) * sizeof(double));
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
    if(k > 0){
      p = vProd(findHSLM(func, x, s,
                         y, nRow, m,
                         k, N_max),
                -1, nRow);
    }else{
      p = vProd(grad, -1, nRow);
    }
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
