/*
 * Luis Manuel Román García
 *
 * -------------------------------------
 * Funciones para calcular Newton truncado
 * -------------------------------------
 */
#include "gen_optim.c"

/* -------------------------------------
 * Gradiente Conjugado
 * Resuelve equaciones de la forma
 * Ax = b
 * IN
 * A: Matríz simétrica positiva definida
 * b: Vector
 * length: longitud de b.
 * OUT
 * x
 * -------------------------------------
 */
double* GC(double* A, double* b, int nrow){
  double *r, *p, *x, *x_new, *r_new, *p_new;
  double alpha, beta, Tol;
  int k, i, N_max;

  // Número máximo de ciclos.
  N_max = 1e3;
  Tol   = 1e-3;

  // Alocar espacio para vectores.
  x     = (double*)malloc(sizeof(double)*nrow);
  x_new = (double*)malloc(sizeof(double)*nrow);
  r_new = (double*)malloc(sizeof(double)*nrow);
  p_new = (double*)malloc(sizeof(double)*nrow);

  // Punto inicial.
  for(i = 0; i < nrow; i++){
    x[i]    = rand() % 10;
  }

  // Inicializar r, p y k.
  // r = Ax0 - b
  r = vSum(mProd(A, x, nrow, nrow), vProd(b, -1, nrow), nrow);
  // p = -r
  p = vProd(r, -1, nrow);
  k = 0;

  // Comenzar loop.
  // Mientras ||r|| > Tol && k < N_max
  while((sqrt(dotProd(r, r, nrow)) > Tol) && (k < N_max)){
    // alpha = r'r/p'Ap
    alpha = dotProd(r, r, nrow) / dotProd(p, mProd(A, p, nrow, nrow), nrow);
    // x_new = x + alpha*p
    x_new = vSum(x, vProd(p, alpha, nrow), nrow);
    x     = x_new;
    // r_new = r + alpha*A*p
    r_new = vSum(r, vProd(mProd(A, p, nrow, nrow), alpha, nrow), nrow);
    // beta  = r_new'r_new/r'r
    beta  = dotProd(r_new, r_new, nrow)/dotProd(r, r, nrow);
    r = r_new;
    // p_new = -r + beta*p
    p_new = vSum(vProd(r, -1, nrow), vProd(p, beta, nrow), nrow);
    p     = p_new;
    k++;
  }

  // Liberar Espacio.
  free(r);

  return x;
}


/* -------------------------------------
 * Newton Gradiente Conjugado
 * IN
 * A: Hessian
 * b: -Gradient
 * OUT
 * p: Dirección de descenso que hace
 * que ||grad(f(x))|| < epsilon||f(x)||
 * -------------------------------------
 */
double* NGC(double (*func)(double*, int), double* b, double* x, int nrow, int N_max, double TOL){
  // Variable declaration.
  double *r, *d, *z, *r_new, *Bd, *p;
  int k, i, j;
  double epsilon, alpha, beta, step;
  // Space allocation.
  z = (double*) malloc(nrow * sizeof(double));
  // Calculate gradient.
  r = gradCentralDiff(func, x, nrow);
  // Outer loop, this modifies x! (stop criteria is missing)
  for(k = 0; norm(r, nrow) >= TOL; k++){
    // Set tolerance.
    epsilon = min(.5, sqrt(norm(r, nrow)) * norm(r, nrow));
    // Initialize z, r, d
    d = vProd(r, -1, nrow);
    for(i = 0; i < nrow; i++){
      z[i] = 0;
    }
    /* -----------------------------------
     * ########### CG Iteration ##########
     * -----------------------------------
     */
    for(j = 0; j < N_max; j++){
      Bd = hessCentralDiff(func, x, d, nrow);
      // Check if d'Bd <= 0
      if(dotProd(d, Bd, nrow) <= 0){
        if(j == 0){
          p = d;
          break;
        }else{
          p = z;
          break;
        }
      }
      alpha = dotProd(r, r, nrow) / dotProd(d, Bd, nrow);  // alpha_j = rj'rj/d_j'Bd_j
      z     = vSum(z, vProd(d, alpha, nrow), nrow);        // z_{j+1} = z_j + alpha_j*d_j
      r_new = vSum(r, vProd(Bd, alpha, nrow), nrow);       // r_{j+1} = r_j + alpha_j*B_kd_j
      if(norm(r_new, nrow) < epsilon){
        p = z;
        break;
      }
      beta = dotProd(r_new, r_new, nrow) / dotProd(r, r, nrow);
      d    = vSum(vProd(r_new, -1, nrow), vProd(d, beta, nrow), nrow);
      /* -----------------------------------
       * ###################################
       * -----------------------------------
       */
    }
    // Choose step via backtracking
    alpha = backTrack(func, x, p, nrow);
    // Update x
    x = vSum(x, vProd(p, alpha, nrow), nrow);
    // Update r
    r = gradCentralDiff(func, x, nrow);
  } // Outer loop, modifies x!
  return x;
}
