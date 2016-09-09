/*
 * Luis Manuel Román García
 *
 * -------------------------------------
 * Rutinas para llevar a cabo
 * operaciones generales para algoritmos
 * de optimización.
 * -------------------------------------
 */

#include "line_alg.c"
#include "utileries.c"

/* -------------------------------------
 * Método de diferencias centrales para
 * aproximar el gradiente de una función
 * IN
 * fun: Apuntador a una función que recibe
 * un apuntador a un double y un entero.
 * x: Apuntador a un vector sobre el cual
 * se quiere aproximar el gradiente
 * length: Longitud del vector.
 * OUT
 * res: Gradiente de fun evaluado en x.
 * -------------------------------------
 */
double* gradCentralDiff(double (*func)(double*, int), double* x, int length){
  // Declaración de variables.
  int i;
  double epsilon, u;
  double *res, *laux, *raux;

  // Inicialización de variables.
  u       = 1.1e-16;
  epsilon = sqrt(u);
  res     = (double*)malloc(length * sizeof(double));
  laux     = (double*)malloc(length * sizeof(double));
  raux     = (double*)malloc(length * sizeof(double));

  // Llenado de vector auxiliar.
  for(i = 0; i < length; i++){
    laux[i] = x[i];
    raux[i] = x[i];
  }

  // Evaluación de derivada.
  for(i = 0; i < length; i++){
    laux[i] = laux[i] + epsilon;
    raux[i] = raux[i] - epsilon;
    res[i]  = (func(laux, length) - func(raux, length)) / (2*epsilon);
    raux[i] = x[i];
    laux[i] = x[i];
  }

  // Liberar memoria
  // free(aux);
  return res;
}

/* -------------------------------------
 * Método de diferencias centrales para
 * aproximar el producto de la hessiana
 * de una función por un vector.
 * IN
 * fun: Apuntador a una función que recibe
 * un apuntador a un double y un entero.
 * x: Apuntador a un vector sobre el cual
 * se quiere aproximar el gradiente
 * length: Longitud del vector.
 * OUT
 * res: Gradiente de fun evaluado en x.
 * -------------------------------------
 */
double* hessCentralDiff(double (*func)(double*, int), double* x, double* p, int length){
  // Declaración de variables.
  double *grad, *aux_grad;
  double epsilon;

  // Initializar epsilon
  epsilon = sqrt(1.1e-6);

  // Cálculo del gradiente y el gradiente auxiliar.
  grad     = gradCentralDiff(func, x, length);
  aux_grad = gradCentralDiff(func, vSum(x, vProd(p, epsilon, length), length), length);

  // Approximación final.
  return vProd(vSum(aux_grad, vProd(grad, -1, length), length), 1/epsilon, length);
}
