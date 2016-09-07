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
double* NGC(double* A, double* b, double* x, int nrow, int N_max){
  double *r, *p, *z, *z_new, *r_new, *d, *d_new;
  double alpha, beta, Tol, epsilon, norm;
  int k, i;

  // Número máximo de ciclos.
  Tol   = 1e-3;

  // Alocar espacio para vectores.
  z_new = (double*)malloc(sizeof(double)*nrow);
  r_new = (double*)malloc(sizeof(double)*nrow);
  d     = (double*)malloc(sizeof(double)*nrow);
  d_new = (double*)malloc(sizeof(double)*nrow);
  z     = (double*)malloc(sizeof(double)*nrow);

  // Inicializar x y z.
  for(i = 0; i < nrow; i++){
    z[i]    = 0;
  }

  // Inicializar r y epsilon.
  // r = grad(fx) = Ax-b
  r       = vSum(mProd(A, x, nrow, nrow), vProd(b, -1, nrow), nrow);
  norm    = sqrt(dotProd(r, r, nrow));
  epsilon = min(.5, sqrt(norm)) * norm;

  // d = -r
  d = vProd(r, -1, nrow);
  k = 0;

  // Comenzar loop.
  // Mientras ||r|| > Tol && k < N_max
  while((sqrt(dotProd(r, r, nrow)) > Tol) && (k < N_max)){
    // Verificar si estamos en una dirección de descenso.
    if(dotProd(d, mProd(A, d, nrow, nrow), nrow) <= 0){
      if(k == 0){
        p = d;
        break; // Break sale de while.
      }else{
        p = z;
        break; // Break sale de while.
      }
    }
    // alpha = r'r/d'Ad
    alpha = dotProd(r, r, nrow) / dotProd(d, mProd(A, d, nrow, nrow), nrow);
    // x_new = x + alpha*d
    z_new = vSum(z, vProd(d, alpha, nrow), nrow);
    z     = z_new;
    // r_new = r + alpha*A*d
    r_new = vSum(r, vProd(mProd(A, d, nrow, nrow), alpha, nrow), nrow);
    // Truncamiento de GC.
    if(sqrt(dotProd(r_new, r_new, norm)) < epsilon){
      p = z;
      break; // Break sale de while.
    }
    // beta  = r_new'r_new/r'r
    beta  = dotProd(r_new, r_new, nrow) / dotProd(r, r, nrow);
    r = r_new;
    // p_new = -r + beta*p
    d_new = vSum(vProd(r, -1, nrow), vProd(d, beta, nrow), nrow);
    d     = d_new;
    k++;
  }

  // Liberar Espacio.
  free(r);

  return p;
}
