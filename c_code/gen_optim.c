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
 * p: El vector por el cual se multiplica
 * la hessiana.
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


/* -------------------------------------
 * Auxiliar function for line search
 * algorithm.
 * -------------------------------------
 */
double phi(double (*func)(double*, int), double* x, double* p, double alpha, int length){
  // phi(a) = f(x + ap)
  return func(vSum(x, vProd(p, alpha, length), length), length);
}

/* -------------------------------------
 * Método auxiliar para búsqueda de línea.
 * IN
 * fun: Apuntador a una función que recibe
 * un apuntador a un double y un entero.
 * x: Apuntador a un vector sobre el cual
 * se quiere aproximar el gradiente
 * p: Una dirección de descenso.
 * length: Longitud del vector.
 * OUT
 * res: Gradiente de fun evaluado en x.
 * -------------------------------------
 */
double zoom(double (*func)(double*, int), double alpha_lo, double alpha_hi, double* x, double*p,  int length){
  double phi_new, alpha, c1, c2, grad_phi;
  c1 = 1; // modify
  c2 = 1; // modify
  while(1){
    // Obtain alpha
    alpha   = (alpha_lo + alpha_hi)/2; // modify
    phi_new = phi(func, x, p, alpha, length);
    if((phi_new > func(x, length) + c1*alpha*dotProd(gradCentralDiff(func, x, length), p, length)) || (phi_new >= phi(func, x, p, alpha_lo, length))){
      alpha_hi = alpha;
    }else{
      grad_phi = dotProd(gradCentralDiff(func, vSum(x, vProd(p, alpha, length), length), length), p, length);
      if(abs(grad_phi) <= -c2 * dotProd(gradCentralDiff(func, x, length), p, length))
        return alpha;
      if(grad_phi * (alpha_hi - alpha_lo))
        alpha_hi = alpha_lo;
      alpha_lo = alpha;
    }
  }
}


/* -------------------------------------
 * Método de búsqueda de línea para tamaño
 * de paso.
 * IN
 * fun: Apuntador a una función que recibe
 * un apuntador a un double y un entero.
 * x: Apuntador a un vector sobre el cual
 * se quiere aproximar el gradiente
 * p: Una dirección de descenso.
 * length: Longitud del vector.
 * OUT
 * res: Gradiente de fun evaluado en x.
 * -------------------------------------
 */
double lineSearch(double (*func)(double*, int), double* x, double* p, int length){
  // Declare variables.
  double alpha_old, alphaM, alpha_new, phi_new, phi_old, c1, c2, grad_phi, alpha_aux;
  int i;
  // Initialize variables.
  alpha_old = 0;
  alphaM    = rand() % 1000;
  alpha_new = rand() % ((int) alphaM);
  i  = 0;
  c1 = 1;  // Modify
  c2 = 1;  // Modify
  while(1){
    phi_new = phi(func, x, p, alpha_new, length);
    if((phi_new > func(x, length) + dotProd(gradCentralDiff(func, x, length), p, length) * (c1 * alpha_new)) ||
       ((phi_new >= phi_old) && i > 1))
      return zoom(func, alpha_old, alpha_new, x, p, length);
    grad_phi = dotProd(gradCentralDiff(func, vSum(x, vProd(p, alpha_new, length), length), length), p, length);
    if(abs(grad_phi) <= -c2*dotProd(gradCentralDiff(func, x, length), p, length))
      return alpha_new;
    if(grad_phi >= 0)
      return zoom(func, alpha_new, alpha_old, x, p, length);
    alpha_aux = alpha_new;
    alpha_new = alpha_new + (rand() % (int)(alpha_new - alpha_old));
    alpha_old = alpha_aux;
  }
}


/* -------------------------------------
 * Backtrack
 * IN
 * fun: Apuntador a una función que recibe
 * un apuntador a un double y un entero.
 * x: Apuntador a un vector sobre el cual
 * se quiere aproximar el gradiente
 * p: Una dirección de descenso.
 * length: Longitud del vector.
 * OUT
 * res: Gradiente de fun evaluado en x.
 * -------------------------------------
 */
double backTrack(double (*func)(double*, int), double* x, double* p, int length){
  // Declare variables.
  double *x_new;
  double alpha, rho, c;
  // Initialize variables.
  alpha = 1; // Return alpha=1 whenever possible.
  rho   = rand() % 1;
  c     = rand() % 1;
  // Iterate
  x_new  = vSum(x, vProd(p, alpha, length), length);
  while(func(x_new, length) <= (func(x, length) + dotProd(vProd(gradCentralDiff(func, x, length), c*alpha, length), p, length))){
    x_new  = vSum(x, vProd(p, alpha, length), length);
    alpha  = alpha * rho;
  }
  return alpha;
}
