#include <math.h>
#include "trunc_newton.c"


double test_func(double*, int);
double dummy_func(double*, int);

int main(){
  // Declaración de variables
  double* point, *p;
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
    point[i] = rand() %10;
    p[i] = rand() % 10;
  }

  // Print point
  imprimeTit("El punto donde se evalualará la hessiana es");
  imprimeMatriz(point, 1, length);

  // Print vector
  imprimeTit("El vector por el que se multiplicará la hessiana es");
  imprimeMatriz(p, 1, length);

  // Print result
  imprimeTit("El resultado es:");
  imprimeMatriz(hessCentralDiff(test_func, point, p, length), 1, length);

   return 0;
}

double test_func(double* x, int length){
  double res;
  int i;
  for(res = i = 0; i < length; i++){
    res = res + x[i]*x[i];
  }
  res = res + x[0]*x[1];
  return res;
};
