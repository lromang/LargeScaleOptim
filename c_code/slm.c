#include "lbfgs.c"


double * twoLoop_r(double (*func)(double*, int), int nRow, int N_max, double TOL,
                   double * x, double** s, double** y, int m, int k){
  // Variable declaration.
  double *grad, *grad_cg, *z, *d, *Bd, *p, *grad_new;
  double rho, alpha, alpha_loop, epsilon, beta, beta_loop;
  int i, j, state;

  // Variable initialization.
  grad  = gradCentralDiff(func, x, nRow);
  state = (int) min(k, m);
  z     = (double*) malloc(nRow * sizeof(double));

  // ####  First Loop  ####
  for(i = (state - 1); i > 0; i--){
    rho        = 1 / (dotProd(y[i], s[i], nRow));
    alpha_loop = rho * dotProd(s[i], grad, nRow);
    grad  =     vSum(grad, vProd(y[i], -alpha_loop, nRow), nRow);
  }
  // ####  CG VARIABLES ####
  // Initialize d, z, eta, r_cg, epsilon.
  d       = vProd(grad, -1, nRow);
  grad_cg = vProd(grad, 1, nRow);
  epsilon = min(.5, sqrt(norm(grad, nRow))) * norm(grad, nRow);
  for(i = 0; i < nRow; i++){
    z[i] = 0;
  }

  /* -----------------------------------
   * ########### CG Iteration ##########
   * -----------------------------------
   */
  for(j = 0; j < N_max; j++){
    Bd = hessCentralDiff(func, x, d, nRow);
    // Check if d'Bd <= 0 i.e. d is a descent direction.
    if(dotProd(d, Bd, nRow) <= 0){
      if(j == 0){
        p = d;
        break;
      }else{
        p = z;
        break;
      }
    }
    alpha    = dotProd(grad_cg, grad_cg, nRow) / dotProd(d, Bd, nRow);  // alpha_j = rj'rj/d_j'Bd_j
    z        = vSum(z, vProd(d, alpha, nRow), nRow);        // z_{j+1} = z_j + alpha_j*d_j
    grad_new = vSum(grad_cg, vProd(Bd, alpha, nRow), nRow);       // r_{j+1} = r_j + alpha_j*B_kd_j
    if(norm(grad_new, nRow) < epsilon){
      p = z;
      break;
    }
    // Update beta, d, r_cg.
    beta    = dotProd(grad_new, grad_new, nRow) / dotProd(grad_cg, grad_cg, nRow);
    d       = vSum(vProd(grad_new, -1, nRow), vProd(d, beta, nRow), nRow);
    grad_cg = grad_new;
  }
  /* -----------------------------------
   * ########### CG Iteration ##########
   * -----------------------------------
   */

  // ####  First Loop  ####
  for(i = 0; i < state; i ++){
    rho       = 1 / (dotProd(y[i], s[i], nRow));
    beta_loop = rho * dotProd(y[i], p, nRow);
    p         = vSum(p, vProd(s[i], (alpha_loop - beta_loop), nRow), nRow);
  }

}
