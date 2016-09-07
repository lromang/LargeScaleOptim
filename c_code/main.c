#include <math.h>
#include "trunc_newton.c"
#include "utileries.c"

double test_func(double*, int);

int main(){
  // Declaración de variables
  double* point;
  int length, i;
  /*
   *********************************************
   * Resolver sistema de ecuaciones Ax = b
   *********************************************
   */
  printf("Enter size of point:\n");
  scanf("%d", &length);

  point = (double*) malloc(length * sizeof(double));
  for(i = 0; i < length; i++){
    point[i] = rand() %10;
  }

  // Print point
  imprimeTit("El punto donde se evalualará la derivada es");
  imprimeMatriz(point, 1, length);

  // Print result
  imprimeTit("El resultado es:");
  imprimeMatriz(centralDiff(test_func, point, length), 1, length);



   return 0;
}

double test_func(double* x, int length){
  double res;
  int i;
  for(i = 0; i < length; i++){
    res = res + x[i]*x[i];
  }
  return res;
};
