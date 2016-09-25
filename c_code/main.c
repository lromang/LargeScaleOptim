/* #########################################
 *
 * Luis Manuel Román García
 * luis.roangarci@gmail.com
 *
 * #########################################
 *
 * -----------------------------------------
 * Main
 * -----------------------------------------
 *
 */


#include <math.h>
#include "utileries.c"

// Size of file.
int const MAX_FILE_ROWS = 150;
int const MAX_FILE_COLS = 4;

// Value storage.
int* logistic_labels;
float logistic_values[150][5];

double logistic_regression(double*, int);
double class_error(double*, int);


int main(){
  // Variable declaration.
  double *optim_point_N, *optim_point_lbfgs;
  double res;
  int i, length;
  // Logistic Variable declaration
  // Ask for size of point.
  printf("Enter size of point:\n");
  scanf("%d", &length);
  // Data file name.
  FILE *file = fopen("../data/iris", "r");
  // Allocate space variables
  logistic_labels = (int*) malloc(MAX_FILE_ROWS * sizeof(int));

  // Read in file
  for(i = 0; i < MAX_FILE_ROWS; i++){
      if (feof(file))
        break;
      logistic_values[i][0];
      fscanf(file, "%f %f %f %f %d",
             &(logistic_values[i][1]),
             &(logistic_values[i][2]),
             &(logistic_values[i][3]),
             &(logistic_values[i][4]),
             &(logistic_labels[i]));
      if(logistic_labels[i] == 0){logistic_labels[i] = -1;};
    }

  /*
   * ###############################################################
   * Test Truncated Newton
   * ###############################################################
   */

  // Print results easy.
  optim_point_N = NGC(test_func, length, 10, 1e-2);
  imprimeTit("Function minimum (NCG):");
  imprimeMatriz(optim_point_N, 1, length);


  // Print results hard.
   optim_point_N = NGC(testFunc, length, 10, 1e-2);
   imprimeTit("Function minimum (NCG):");
   imprimeMatriz(optim_point_N, 1, length);


  // Test logistic.
  optim_point_N = NGC(logistic_regression, 5, 10, 1e-2);
  imprimeTit("Function minimum (NCG):");
  imprimeMatriz(optim_point_N, 1, length);

  // Prediction error.
  imprimeTit("Class Error:");
  printf("%f \n", class_error(optim_point_N, length));

  /*
   * ###############################################################
   * Test LBFGS
   * ###############################################################
   */

  // Print result
  /*
   * optim_point_lbfgs = LBFGS(test_func, length, 20, 1e-2);
   * imprimeTit("Function minimum (LBFGS):");
   * imprimeMatriz(optim_point_lbfgs, 1, length);
  */
   return 0;
}


/* -------------------------------------
 * Logistic regression.
 * ## Characteristics ##
 * Easy
 * -------------------------------------
 */
double logistic_regression(double* x, int length){
  // Variable declaration.
  double res;
  int *y;
  int i;
  // Evaluate logistic.
  res = 0;
  for(i = 0; i < MAX_FILE_ROWS; i++){
    res = res + log(1 + exp(-dotProd(x, (double*) logistic_values[i], 5) * logistic_labels[i]));
  }
  return res;
};


/*
 * -------------------------------------
 * Eval function
 * -------------------------------------
 */
double class_error(double* coefs, int length){
  // Variable declaration.
  double class_error, entry_val;
  int *y;
  int i, pred;
  // Evaluate logistic.
  pred = 0;
  for(i = 0; i < MAX_FILE_ROWS; i++){
    entry_val = exp(- dotProd(coefs, (double*) logistic_values[i], length));
    entry_val = 1 / (1 + entry_val);
    // Classification threshold = .5
    pred < .5 ? pred = 1 : -1;
    class_error = class_error + (pred == logistic_labels[i] ? 1 : 0);
    //res = res + log(1 + exp(-dotProd(x, (double*) logistic_values[i], 5) * logistic_labels[i]));
  }
  return class_error / MAX_FILE_ROWS;
};
