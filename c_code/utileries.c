#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "trunc_newton.c"

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
 * Logistic regression.
 * ## Characteristics ##
 * Easy
 * -------------------------------------
 */
double logistic_regression(double* x, int length){
  // Variable declaration.
  int   *y;
  double res;
  int i, MAX_FILE_ROWS = 150;
  float w[MAX_FILE_ROWS][length];
  // Data file name.
  FILE *file = fopen("../data/iris", "r");
  // Allocate space
  y = (int*) malloc(MAX_FILE_ROWS * sizeof(int));
  // Read in values
  for(i = 0; i < MAX_FILE_ROWS; i++)
    {
      if (feof(file))
        break;
      // Read in values.
      fscanf(file, "%f %f %f %f %d", &(w[i][0]), &(w[i][1]), &(w[i][2]), &(w[i][3]), &(y[i]));
    }
  // Evaluate logistic.
  res = 0;
  for(i = 0; i < MAX_FILE_ROWS; i++){
    res = res + log(1 + exp(-dotProd(x, (double*) w[i], length)*y[i]));
  }
  return res;
};
