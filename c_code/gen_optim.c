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
double* centralDiff(double (*func)(double*, int), double* x, int length){
  // Declaración de variables.
  int i;
  double epsilon, u;
  double *res, *aux;

  // Inicialización de variables.
  u       = 1.1e-16;
  epsilon = sqrt(u);
  res     = (double*)malloc(length * sizeof(double));
  aux     = (double*)malloc(length * sizeof(double));

  // Llenado de vector auxiliar.
  for(i = 0; i < length; i++){
    aux[i] = x[i];
  }

  // Evaluación de derivada.
  for(i = 0; i < length; i++){
    aux[i] = aux[i] + epsilon;
    res[i] = (func(aux, length) - func(x, length)) / epsilon;
    aux[i] = x[i];
  }

  // Liberar memoria
  free(aux);
  return res;
}

/* -------------------------------------
 * Método para buscar tamaño de paso
 * en la dirección de descenso.
 * IN
 * func: Apuntador a una función que recibe
 * cómo parámetros un apuntador a un double
 * y un entero.
 * x: Vector donde se lleva a cabo la
 * evaluación.
 * p: Dirección de descenso.
 * length: Longitud de los vectores.
 * OUT
 * alpha: Tamaño del paso en la dirección
 * de descenso.
 * -------------------------------------
 */
double backTrack(double (*func)(double*, int),
                 double* x,
                 double* p,
                 int length){

  // Decalaración de variables
  double alpha, rho, c, fx, fx_p;
  double *gradx, *aux_vprod, *aux_vsum;

  // Inicialización de variables
  alpha = rand() % 10;
  rho   = rand() % 1;
  c     = rand() % 1;
  aux_vprod = (double*)malloc(length * sizeof(double));
  aux_vsum  = (double*)malloc(length * sizeof(double));
  gradx     = (double*)malloc(length * sizeof(double));

  // Condiciones del Loop
  gradx = centralDiff(func, x, length);
  do{
    aux_vprod = vProd(p, alpha, length);
    aux_vsum  = vSum(x, aux_vprod, length);
    fx_p      = func(aux_vsum, length);
    fx        = func(x, length) + c * alpha * dotProd(gradx, p, length);
    alpha     = alpha * rho;
  } while(fx_p <= fx);

  // Liberar espacio
  free(aux_vprod);
  free(aux_vsum);
  free(gradx);

  return alpha;
}
