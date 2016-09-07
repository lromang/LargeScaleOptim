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
