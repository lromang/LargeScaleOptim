#include <math.h>
#include "utileries.c"

double test_func(double*, int);
double dummy_func(double*, int);

int main(){
  // Variable declaration.
  double *optim_point_N, *point;
  int length, i;




  // Ask for size of point.
  printf("Enter size of point:\n");
  scanf("%d", &length);

  /*
   * ###############################################################
   * Test Truncated Newton
   * ###############################################################
   */

  // Print results easy.
  optim_point_N = NGC(test_func, length, 10, 1e-2);
  imprimeTit("Function minimum (NCG):");
  imprimeMatriz(optim_point_N, 1, length);


  // Print results hard.
  optim_point_N = NGC(testFunc, length, 10, 1e-2);
  imprimeTit("Function minimum (NCG):");
  imprimeMatriz(optim_point_N, 1, length);

  // Test logistic.
  // Generate point
  optim_point_N = NGC(logistic_regression, min(length, 4), 10, 1e-2);
  imprimeTit("Function minimum (NCG):");
  imprimeMatriz(optim_point_N, 1, length);

  /*
   * ###############################################################
   * Test LBFGS
   * ###############################################################
   */

  // Print result
  //optim_point_lbfgs = LBFGS(test_func, length, 10, 1e-2);
  //imprimeTit("Function minimum (LBFGS):");
  //imprimeMatriz(optim_point_lbfgs, 1, length);

   return 0;
}
