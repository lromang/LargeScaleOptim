/*
 * Luis Manuel Román García
 *
 * -------------------------------------
 * Funciones para calcular Newton truncado
 * -------------------------------------
 */
#include "trunc_newton.c"

/* -------------------------------------
 * Approximate gf'H in LBFGS method
 * IN
 * fun:    Function of which the Hessian is
 *         going to be approximated.
 * x:      Point where the Hessian is going
 *         to be approximated.
 * length: x's length.
 * m:      Size of 'memory'
 * k:      Number of iteration.
 * s:      Pointer to las 'm' iterations of s.
 * y:      Pointer to las 'm' iterations of y.
 * OUT
 * res: gf'H evaluated on x.
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
  // Memory release.
  free(q);
  // Return result.
  return r;
}


/* -------------------------------------
 * Update s and y
 * Keeps track of last 'm' instances of
 * 's' and 'y'
 * IN
 * s:     Pointer to las 'm' iterations of s.
 * y:     Pointer to las 'm' iterations of y.
 * s_new: new element of s
 * y_new: new element of y
 * k:     Iteration number.
 * m:     Size of 'memory'.
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
double* LBFGS(double (*func)(double*, int), int m, int nRow, double TOL){
  // Declare variables
  double **s, **y;
  double *grad, *p, *x_new, *grad_new, *x;
  double  alpha;
  int k, i;
  // Space allocation.
  // x
  x = (double*) malloc(nRow * sizeof(double));
  // s
  s = (double**) malloc((m * nRow) * sizeof(double));
  // y
  y = (double**) malloc((m * nRow) * sizeof(double));
  // Initial x.
  for(i = 0; i < nRow; i++){
    x[i] = ((double) rand()/INT_MAX) + 1;
  }
  // Initialize variables
  grad = gradCentralDiff(func, x, nRow);
  k    = 0;
  // Initial values s y.
  updateSY(s, y, x, grad, k, m);
  while(norm(grad, nRow) > TOL){
    p        = findH(func, x, s, y, nRow, m, k); // caso k = 0
    alpha    = backTrack(func, x, p, nRow);
    x_new    = vSum(x, vProd(p, alpha, nRow), nRow);
    grad_new = gradCentralDiff(func, x_new, nRow);
    // Update s y.
    //updateSY(s, y, vSum(x_new, vProd(x, -1, nRow), nRow), vSum(grad_new, vProd(grad, -1, nRow), nRow), k, m);
    // Update grad, x and k.
    grad = grad_new;
    x    = x_new;
    k ++;
  }
 /*// Memory release.
  free(grad);
  free(p);
  free(x_new);
  free(grad_new);*/
  // Return results.
  return x;
}
