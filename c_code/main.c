#include <stdio.h>
#include <stdlib.h>
#include "line_alg.c"
#include "utileries.c"

int main(){
  // Declaración de variables
  double *matrix, *vector;
  int    nrow, alpha;

  /*
   *********************************************
   * Solicitar a usuario tamaño de matriz y
   * vectores a utilizar en pruebas.
   *********************************************
   */
  printf("Escriba el tamaño de la matriz: \n");
  scanf("%d", &nrow);
  matrix = (double*) malloc((nrow * nrow) * sizeof(double));

  // Llenado de matriz y vector
  matrix = creaMatriz(nrow, nrow);
  vector = creaMatriz(1, nrow);

  // Imprimir matriz y vector.
   imprimeTit("Matriz del sistema");
   imprimeMatriz(matrix, nrow, nrow);
   imprimeTit("Vector del sistema");
   imprimeMatriz(vector, 1, nrow);

  /*
   *********************************************
   * Probar funciones.
   *********************************************
   */
   imprimeTit("El producto entre matriz y vector es");
   imprimeMatriz(mProd(matrix, vector, nrow, nrow), 1, nrow);
}
