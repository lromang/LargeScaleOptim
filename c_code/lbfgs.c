/*
 * Luis Manuel Román García
 *
 * -------------------------------------
 * Funciones para calcular Newton truncado
 * -------------------------------------
 */
#include "gen_optim.c"

/* -------------------------------------
 * Método para obtener una aproximación
 * de H en LBFGS.
 * IN
 * fun: Apuntador a una función que recibe
 * un apuntador a un double y un entero.
 * x: Apuntador a un vector sobre el cual
 * se quiere aproximar el gradiente
 * length: Longitud del vector.
 * m: memoria
 * k: iteración actual.
 * s: apuntador a iteraciones pasadas de s
 * y: apuntador a iteraciones pasadas de y
 * OUT
 * res: Gradiente de fun evaluado en x.
 * -------------------------------------
 */
double* findH(double (*func)(double*, int), double* x, double ** s, double ** y, int length, int m, int k){
  // Variable declaration
  int i, constant;
  double alpha, beta, rho;
  double *q, *r;
  // Variable initialization
  q = gradCentralDiff(func, x, length);
  for(i = 0; i < m; i++){
    rho   = 1 / dotProd(y[k - i], s[k - i], length);
    alpha = rho * dotProd(s[k - i], q, length);
    q     = vSum(q , vProd(y[k - i], -alpha, length), length);
  }
  constant = dotProd(s[k - 1], y[k - 1], length) / dotProd(y[k - 1], y[k - 1], length);
  r        = mProd(vProd(identity(length), constant, length * length), q, length, length);
  for(i = 0; i < m; i ++){
    rho  = 1 / dotProd(y[k - m + i], s[k - m + i], length);
    beta = rho * dotProd(y[k - m + i], r, length);
    r    = vSum(r, vProd(s[k - m + i], alpha - beta, length), length);
  }
  return r;
}


/* -------------------------------------
 * UpdateSY
 * IN
 * s: apuntador a iteraciones pasadas de s
 * y: apuntador a iteraciones pasadas de y
 * s_new: new element of s
 * y_new: new element of y
 * k: iteration
 * m: memory.
 * OUT
 * s
 * y
 * -------------------------------------
 */
void updateSY(double ** s, double ** y, double * s_new, double * y_new, int k,  int m){
  int i;
  // k == m
  if(k >= m){
    for(i = 1; i < m; i++){
      s[m - i - 1] = s[m - i];
      y[m - i - 1] = y[m - i];
    }
    s[m - 1] = s_new;
    y[m - 1] = y_new;
  }else{
    s[k] = s_new;
    y[k] = y_new;
  }
}

/* -------------------------------------
 * LBFGS.
 * IN
 * fun: Apuntador a una función que recibe
 * un apuntador a un double y un entero.
 * x: Apuntador a un vector sobre el cual
 * se quiere aproximar el gradiente
 * length: Longitud del vector.
 * m: memoria
 * k: iteración actual.
 * s: apuntador a iteraciones pasadas de s
 * y: apuntador a iteraciones pasadas de y
 * OUT
 * res: Gradiente de fun evaluado en x.
 * -------------------------------------
 */
double* LBFGS(double (*func)(double*, int), double* x, double ** y, double **s, int m, int nrow, int N_max, double TOL){
  // Declare variables
  double *grad, *H0, *p, *x_new, *grad_new;
  double constant, alpha;
  int k;
  // Initialize variables
  grad = gradCentralDiff(func, x, nrow);
  k    = 0;
  while(norm(grad, nrow) > TOL){
    constant = dotProd(s[k - 1], y[k - 1], nrow) / dotProd(y[k - 1], y[k - 1], nrow);
    H0       = vProd(identity(nrow), constant, nrow * nrow);
    p        = findH(func, x, s, y, nrow, m, k);
    alpha    = backTrack(func, x, p, nrow);
    x_new    = vSum(x, vProd(p, alpha, nrow), nrow);
    grad_new = gradCentralDiff(func, x_new, nrow);
    // Update s y.
    updateSY(s, y, vSum(x_new, vProd(x, -1, nrow), nrow), vSum(grad_new, vProd(grad, -1, nrow), nrow), k, m);
    // Update grad, x and k.
    grad = grad_new;
    x    = x_new;
    k ++;
  }
  return x;
}
