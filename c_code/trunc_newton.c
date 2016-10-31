/* #########################################
 *
 * Luis Manuel Román García
 * luis.roangarci@gmail.com
 *
 * #########################################
 *
 * -----------------------------------------
 * Functions for the Truncated Newton
 * algorithm.
 * -----------------------------------------
 *
 */

#include "gen_optim.c"

/* -------------------------------------
 * Conjugate gradient.
 * Applies conjugate gradient method to
 * solve problems of the form: Ax = b
 * where A is symmetric positive definite
 * matrix.
 * IN
 * A:    A spd matrix
 * b:    The vector to be generated by the
 *       columns of A.
 * nRow: Number of rows of A.
 * OUT
 * x: Solution of the system.
 * -------------------------------------
 */
double* GC(double* A, double* b, int nRow){
  double *r, *p, *x, *x_new, *r_new, *p_new;
  double alpha, beta, Tol;
  int k, i, N_max;
  // Maximum number of iterations.
  N_max = 1e3;
  Tol   = 1e-3;
  // Space allocation.
  x     = (double*)malloc(sizeof(double)*nRow);
  x_new = (double*)malloc(sizeof(double)*nRow);
  r_new = (double*)malloc(sizeof(double)*nRow);
  p_new = (double*)malloc(sizeof(double)*nRow);
  // Set initial point.
  for(i = 0; i < nRow; i++){
    x[i]    = rand() % 10;
  }
  // Initialize r, p y k.
  // r = Ax0 - b
  r = vSum(mProd(A, x, nRow, nRow), vProd(b, -1, nRow), nRow);
  // p = -r
  p = vProd(r, -1, nRow);
  k = 0;
  // Loop
  // While ||r|| > Tol && k < N_max
  while((sqrt(dotProd(r, r, nRow)) > Tol) && (k < N_max)){
    // alpha = r'r / p'Ap
    alpha = dotProd(r, r, nRow) / dotProd(p, mProd(A, p, nRow, nRow), nRow);
    // x_new = x + alpha * p
    x_new = vSum(x, vProd(p, alpha, nRow), nRow);
    x     = x_new;
    // r_new = r + alpha * A * p
    r_new = vSum(r, vProd(mProd(A, p, nRow, nRow), alpha, nRow), nRow);
    // beta  = r_new'r_new / r'r
    beta  = dotProd(r_new, r_new, nRow)/dotProd(r, r, nRow);
    r = r_new;
    // p_new = -r + beta * p
    p_new = vSum(vProd(r, -1, nRow), vProd(p, beta, nRow), nRow);
    p     = p_new;
    k++;
  }
  // Memory release.
  free(r);
  free(p);
  // Return result.
  return x;
}


/* -------------------------------------
 * Truncated Newton with Conjugate Gradient
 * IN
 * func:  Function to be minimized.
 * x:     Initial point.
 * nRow:  Number of rows (length) of x.
 * N_max: Maximum number of CG iterates.
 * TOL:   Minimum size for gradient of f.
 * OUT
 * x: Local minimum of func.
 * -------------------------------------
 */
double* NGC(double (*func)(double*, int), int nRow,
            int N_max, double TOL, int verbose,
            int gradN, int gradSamp){
  // Variable declaration.
  double *r, *d, *z, *r_new,
    *Bd, *p, *x_new, *x, *r_cg;
  int k, i, j, stop, wolf_cond, grad_iter;
  double epsilon, alpha, beta, step, eta, sg_alpha;

  // Space allocation.
  z = (double*) malloc(nRow * sizeof(double));
  x = (double*) malloc(nRow * sizeof(double));

  /*
   * Initial point x
   * Alternative: Some iterations of Stochastic
   * Gradient Descent for a better initial point.
   */
  for(i = 0; i < nRow; i++){
    // Improve initial point with Stochastic Gradient Descent!
    x[i] = ((double) rand() / INT_MAX);
  }
  /*
   * STOCHASTIC GRADIENT DESCENT BEGIN
   * Improve Loop Conditions...
   * Parameters Tolerance, size of sample.
   */
  if(run_logistic){
    for(grad_iter = 0; grad_iter < gradN; grad_iter++){
    // Choose random observation
      SAMPLE = gradSamp;
      create_sample(0);
      r        = gradCentralDiff(func, x, nRow);
      sg_alpha = .0001; //backTrack(func, x, r, nRow, verbose); // beware
      x        = vSum(x, vProd(r, -sg_alpha, nRow), nRow);
      if(verbose && (grad_iter % 10 == 0)){
        printf("\n ITER = %d; f(x) = %.10e;  ||grad|| =  %.10e ; "
               " alpha =  %.10e;",
               grad_iter,
               func(x, nRow),
               norm(r, nRow),
               sg_alpha);
      }
    }
  }
  /*
   * STOCHASTIC GRADIENT DESCENT END
   */
  // Sample mode
  if(stocMode){
      printf("\nRUNNING STOCASTIC MODE\n");
      SAMPLE      = rand() % (int)(MAX_FILE_ROWS * sampProp);
      create_sample(verbose);
  }
  // Calculate the gradient.
  r    = gradCentralDiff(func, x, nRow);
  stop = 8;
  // Outer loop, this modifies x!
  for(k = 0; (norm(r, nRow) >= TOL) && (k < stop); k++){
    // Stochastic mode.
    if(stocMode && k){
      printf("\nRUNNING STOCASTIC MODE\n");
      SAMPLE = rand() % (int)(MAX_FILE_ROWS * sampProp);
      create_sample(verbose);
    }
    // ############ CG Variables ###########
    // Set tolerance.
    epsilon = min(.5, sqrt(norm(r, nRow))) * norm(r, nRow);
    // Initialize d, z, eta, r_cg.
    d    = vProd(r, -1, nRow);
    r_cg = r; // testing
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
      alpha = dotProd(r_cg, r_cg, nRow) / dotProd(d, Bd, nRow);  // alpha_j = rj'rj/d_j'Bd_j
      z     = vSum(z, vProd(d, alpha, nRow), nRow);        // z_{j+1} = z_j + alpha_j*d_j
      r_new = vSum(r_cg, vProd(Bd, alpha, nRow), nRow);       // r_{j+1} = r_j + alpha_j*B_kd_j
      if(norm(r_new, nRow) < epsilon){
        p = z;
        break;
      }
      // Update beta, d, r_cg.
      beta = dotProd(r_new, r_new, nRow) / dotProd(r_cg, r_cg, nRow);
      d    = vSum(vProd(r_new, -1, nRow), vProd(d, beta, nRow), nRow);
      r_cg = r_new;
    }
    /* -----------------------------------
     * ######### CG Iteration End ########
     * -----------------------------------
     */

    // Choose step via backtracking
    step = 1;
    x_new = vSum(x, vProd(p, step, nRow), nRow);
    // Backtrack loop
    eta = 1e-4;
    wolf_cond = 0;
    while(func(x_new, nRow) > func(x, nRow) + (eta * step *  dotProd(r, p, nRow))){
      // Update x
      x_new = vSum(x, vProd(p, step, nRow), nRow);
      step  = step / 2;
      wolf_cond = wolf_cond + 1;
    }
    x = x_new;
    // Update r
    r = gradCentralDiff(func, x, nRow);
    // ---------------- PRINT ------------------- //
    if(verbose){
      printf("\n ITER = %d; f(x) = %.10e;  ||grad|| =  %.10e ; "
             "||p|| =  %.10e ; alpha =  %.10e; backtrack iters = %d",
             k,
             func(x, nRow),
             norm(r, nRow),
             norm(p, nRow),
             step,
             wolf_cond);
    }
    // ---------------- PRINT ------------------- //y
  } // Outer loop, modifies x!
  // Memory release.
  free(z);
  // Return result.
  return x;
}
