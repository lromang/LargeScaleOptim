#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* -------------------------------------
 * Print title
 * IN
 * title: Title to be printed.
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
 * Test function:
 * ## Characteristics ##
 * Ill conditioned.
 * -------------------------------------
 */
double testFunc(double* x, int length){
  int i;
  double sum;
  sum = 0;
  for(i = 0; i < length; i++){
    sum = sum  + (100 - i) * (x[i] * x[i]) + exp(x[i]);
  }
  // Return result.
  return sum;
}

/* -------------------------------------
 * Test function:
 * ## Characteristics ##
 * Easy
 * -------------------------------------
 */
double test_func(double* x, int length){
  double res;
  int i;
  for(res = i = 0; i < length; i++){
    res = res + (3 - x[i])*(4 - x[i]) + x[i]*x[i]*x[i]*x[i];
  }
  return res;
};

/* -------------------------------------
 * Read files
 * -------------------------------------
 */
const char* getfield(char* line, int num)
{
    const char* tok;
    for (tok = strtok(line, ",");
            tok && *tok;
            tok = strtok(NULL, ",\n"))
    {
        if (!--num)
            return tok;
    }
    // Return result.
    return NULL;
}
