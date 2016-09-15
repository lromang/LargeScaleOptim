#include <math.h>
#include "lbfgs.c"

double test_func(double*, int);
double dummy_func(double*, int);

int main(){
  // Variable declaration.
  double *optim_point_N, *optim_point_lbfgs;
  int length;

  // Ask for size of point.
  printf("Enter size of point:\n");
  scanf("%d", &length);

  /*
   * ###############################################################
   * Test Truncated Newton
   * ###############################################################
   */

  // Print result
  optim_point_N = NGC(test_func, length, 10, 1e-2);
  imprimeTit("Function minimum (NCG):");
  imprimeMatriz(optim_point_N, 1, length);

  /*
   * ###############################################################
   * Test LBFGS
   * ###############################################################
   */

  // Print result
  // optim_point_lbfgs = LBFGS(test_func, 10, length, 1e-2);
  imprimeTit("Function minimum (LBFGS):");
  // imprimeMatriz(optim_point_lbfgs, 1, length);


   return 0;
}

double test_func(double* x, int length){
  double res;
  int i;
  for(res = i = 0; i < length; i++){
    res = res + (3 - x[i])*(4 - x[i]) + x[i]*x[i]*x[i];
  }
  res = res;
  return res;
};
