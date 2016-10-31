/* #########################################
 *
 * Luis Manuel Román García
 * luis.roangarci@gmail.com
 *
 * #########################################
 *
 * -----------------------------------------
 * General purpose utileries
 * -----------------------------------------
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Size of file.
int const MAX_FILE_ROWS = 2000;//7000000;
int const MAX_FILE_COLS = 27;
int const N_CLASS = 2;
// Value storage.
int    logistic_labels[2000];
double logistic_values[2000][27];
// Sample values
int    sample_logistic_labels[2000];
double sample_logistic_values[2000][27];
int run_logistic, SAMPLE, stocMode;

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
    res = res + x[i]*x[i]*x[i]*x[i] + (3 - x[i])*(4 - x[i]);
  }
  return res;
};

/* -------------------------------------
 * Create sample
 * -------------------------------------
 * This function receives a size and modifies
 * the structure of a global array of data
 */
void create_sample(int verbose){
  // Variable declaration.
  int* indexes;
  int i, j;
  // Modify global SAMPLE.
  if(verbose){
    printf("Tamaño muestra: %d \n", SAMPLE);
  }
  // Memory allocation.
  indexes = (int*) malloc(SAMPLE * sizeof(int));
  // Random indexes construction.
  for(i = 0; i < SAMPLE; i++){
    indexes[i] = rand() % SAMPLE;
  }
  // Fill in sample labels and sample values.
  for(i = 0; i < SAMPLE; i++){
    sample_logistic_labels[i] = logistic_labels[indexes[i]];
    for(j = 0; j < MAX_FILE_COLS; j++){
      sample_logistic_values[i][j] = logistic_values[indexes[i]][j];
    }
  }
}
