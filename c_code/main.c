#include <math.h>
#include "trunc_newton.c"


double test_func(double*, int);
double dummy_func(double*, int);

int main(){
  // Declaración de variables
  double *point, *p, *x;
  int length, i;
  /*
   *********************************************
   * Resolver sistema de ecuaciones Ax = b
   *********************************************
   */
  printf("Enter size of point:\n");
  scanf("%d", &length);

  // Alocar espacio
  point = (double*) malloc(length * sizeof(double));
  p     = (double*) malloc(length * sizeof(double));

  for(i = 0; i < length; i++){
    point[i] = 1;
    p[i] = rand() % 10;
  }

  /* ###############################################################
   * Prueba NGC
   * ###############################################################
   */
  x = point;

  imprimeTit("El punto inicial donde se evalualará la función es:");
  imprimeMatriz(x, 1, length);

  // Print result
  imprimeTit("El mínimo de la función es:");
  imprimeMatriz(NGC(test_func, x, length, 100, 1e-3), 1, length);
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
