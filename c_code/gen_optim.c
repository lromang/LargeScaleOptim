/* #########################################
 * Luis Manuel Román García
 * luis.roangarci@gmail.com
 * #########################################
 *
 * -----------------------------------------
 * General purpose optimization functions.
 * -----------------------------------------
 *
 */

#include "line_alg.c"

/* -------------------------------------
 * Gradient approximation.
 * Approximates the gradient of func via
 * central differences.
 * IN
 * fun:    Funciton from which the gradient
 *         is to be obtained.
 * x:      Point where the gradient is to be
 *         obtained.
 * length: x's length.
 * OUT
 * res: A pointer to the gradient of func
 *      evaluated on x.
 * -------------------------------------
 */
double* gradCentralDiff(double (*func)(double*, int), double* x, int length){
  // Variable declaration.
  double *res, *laux, *raux;
  double epsilon, u;
  int i;
  // Variable initialization and space allocation.
  u       = 1.1e-16;
  epsilon = sqrt(u);
  res     = (double*)malloc(length * sizeof(double));
  laux    = (double*)malloc(length * sizeof(double));
  raux    = (double*)malloc(length * sizeof(double));
  // Auxiliar vectors fill in.
  for(i = 0; i < length; i++){
    laux[i] = x[i];
    raux[i] = x[i];
  }
  // Entry-wise gradient approximation.
  for(i = 0; i < length; i++){
    laux[i] = laux[i] + epsilon;
    raux[i] = raux[i] - epsilon;
    res[i]  = (func(laux, length) - func(raux, length)) / (2 * epsilon);
    raux[i] = x[i];
    laux[i] = x[i];
  }
  // Memory release.
  free(laux);
  free(raux);
  // Return result.
  return res;
}

/* -------------------------------------
 * Hessian vector product approximation
 * IN
 * fun:    Funciton from which the hessian
 *         is to be obtained.
 * x:      Point where the gradient is to be
 *         obtained.
 * p:      Vector that multiplies the
 *         Hessian.
 * length: x's length.
 * OUT
 * Pointer to result of product between
 * the Hessian and p evaluated on x.
 * -------------------------------------
 */
double* hessCentralDiff(double (*func)(double*, int), double* x, double* p, int length){
  // Variable declaration.
  double *grad, *aux_grad, *hess;
  double epsilon;
  // Epsilon initialization.
  epsilon = sqrt(1.1e-6);
  // Gradient and auxiliar gradient.
  grad     = gradCentralDiff(func, x, length);
  aux_grad = gradCentralDiff(func, vSum(x, vProd(p, epsilon, length), length), length);
  // Hessian approximation.
  hess = vProd(vSum(aux_grad, vProd(grad, -1, length), length), 1 / epsilon, length);
  // Memory release.
  free(grad);
  free(aux_grad);
  // Return result.
  return hess;
}

/* -------------------------------------
 * Backtrack
 * Calculates a step length via backtracking
 * IN
 * fun:    Funciton from which the step length
 *         is to be obtained.
 * x:      Point where the gradient is to be
 *         obtained.
 * p:      Vector that multiplies the
 *         Hessian.
 * length: x's length.
 * OUT
 * alpha: Step length.
 * -------------------------------------
 */
double backTrack(double (*func)(double*, int), double* x, double* p, int length, int verbose){
  // Variable declaration.
  double alpha, rho, c;
  int wolf_iters;
  // Variable initialization
  alpha = 1; // Return alpha = 1 whenever possible.
  rho   = 0.5; //((double) rand()/INT_MAX) + 1;
  c     = 1e-4;
  // Iterate
  // alpha > 1e-10 just to avoid numerical 0.
  wolf_iters = 0;
  while(func(vSum(x, vProd(p, alpha, length), length), length) > (func(x, length) +
                               dotProd(vProd(gradCentralDiff(func, x, length),
                                             c * alpha, length), p, length))){
    alpha      = alpha * rho;
    wolf_iters = wolf_iters + 1;
  }
  if(verbose){
    printf("\n backtrack iters = %d\n", wolf_iters);
  }
  // Return result.
  return alpha;
}
