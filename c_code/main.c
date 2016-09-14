#include <math.h>
#include "trunc_newton.c"


double test_func(double*, int);
double dummy_func(double*, int);

int main(){
  // Variable declaration.
  double *point, *p, *x;
  int length, i;

  // Ask for size of point.
  printf("Enter size of point:\n");
  scanf("%d", &length);

  // Space allocation.
  point = (double*) malloc(length * sizeof(double));
  p     = (double*) malloc(length * sizeof(double));

  // Point construction.
  for(i = 0; i < length; i++){
    point[i] = 1;
    p[i] = rand() % 10;
  }

  /*
   * ###############################################################
   * Tests
   * ###############################################################
   */
  x = point;

  imprimeTit("Initial point:");
  imprimeMatriz(x, 1, length);

  // Print result
  imprimeTit("Function minimum:");
  imprimeMatriz(NGC(test_func, x, length, 100, 1e-3), 1, length);
   return 0;
}

double test_func(double* x, int length){
  double res;
  int i;
  for(res = i = 0; i < length; i++){
    res = res + (3 - x[i])*(4 - x[i]);
  }
  res = res;
  return res;
};
