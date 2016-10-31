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

// ACTUAL DATA VARIABLES
int const MAX_FILE_ROWS = 2000;//7000000;
int const MAX_FILE_COLS = 27;
int    logistic_labels[2000];
double logistic_values[2000][27];
// SAMPLE VARIABLES
double sample_logistic_values[2000][27];
int    sample_logistic_labels[2000];
int SAMPLE;
// DEFAULT CONFIGURATION
double sampProp       = 1;
double regularization = .001;
int run_logistic      = 1;
int stocMode          = 1;
int seed              = 123454321;
int verbose           = 1;
int run_functions     = 0;
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

/* -------------------------------------
 * Print Configuration
 * -------------------------------------
 * This function receives a size and modifies
 * the structure of a global array of data
 */
void printConfig(){
  if(run_logistic){
  printf("\n\n#################################################");
  imprimeTit("LARGE SCALE OPTIMIZATION");
  printf("Actual configuration: \n");
  printf("1) EXECUTE LOGISTIC FUNCTION: %d\n", run_logistic);
  printf("2)       - REGULARIZATION: %lf\n", regularization);
  printf("3)       - STOCASTIC OPTIM: %d\n", stocMode);
  printf("4)       - RANDOM SEED: %d\n", seed);
  printf("5)       - PROPORTION OF SAMPLE: %lf\n", sampProp);
  printf("6) EXECUTE FUNCTIONS: %d\n", run_functions);
  printf("7) VERBOSE MODE: %d\n", verbose);
  printf("#################################################\n");
  }else{
  printf("\n\n#################################################");
  imprimeTit("LARGE SCALE OPTIMIZATION");
  printf("Actual configuration: \n");
  printf("1) EXECUTE LOGISTIC FUNCTION: %d\n", run_logistic);
  printf("2) EXECUTE FUNCTIONS: %d\n", run_functions);
  printf("3) VERBOSE MODE: %d\n", verbose);
  printf("#################################################\n");
  }
}
