#include <stdio.h>
#include <stdlib.h>

/* -------------------------------------
 * Imprime titulo
 * IN
 * title: Titulo que se quiere imprimir
 * -------------------------------------
 */
void imprimeTit(char * title){
  int k;
  printf("\n------------\n");
  for(k = 0; title[k] != '\0'; k ++){
    printf("%c", title[k]);
  }
  printf("\n------------\n");
}


/* -------------------------------------
 * Función de prueba.
 * -------------------------------------
 */
double testFunc(double* x, int length){
  int i;
  double sum;
  sum = 0;
  for(i = 0; i < length; i++){
    sum = sum  + (100 - i + 1) * (x[i] * x[i]);
  }
  return sum;
}
