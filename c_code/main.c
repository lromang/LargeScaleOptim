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
#include "slm.c"

// Size of file.
int const MAX_FILE_ROWS = 150;
int const MAX_FILE_COLS = 4;
int const N_CLASS = 4;

// Value storage.
int*   logistic_labels;
double  logistic_values[150][5];
double stochastic_logistic_regression(double*, int);
double softmax(double*, int);
double class_error(double*, int);


int main(){
  // Variable declaration.
  double *optim_point_N, *optim_point_lbfgs, *optim_point_slm_lbfgs;
  int i, length, seed, verbose;
  // Logistic Variable declaration
  // Ask for size of point.
  //printf("Enter size of point:\n");
  //scanf("%d", &length);
  verbose = 0;
  length  = 100;
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
   * Test Truncatd Newton
   * ###############################################################
   */

  // Print results easy.
  optim_point_N = NGC(test_func, length, 10, 1e-2, verbose);
  imprimeTit("Function minimum (NCG):");
  imprimeMatriz(optim_point_N, 1, length);


  // Print results hard.
  optim_point_N = NGC(testFunc, length, 10, 1e-2, verbose);
  imprimeTit("Function minimum (NCG):");
  imprimeMatriz(optim_point_N, 1, length);

  //printf("Input a random number seed: ");
  //scanf("%u", &seed);
  seed = 34234;
  srand(seed);

  // Test logistic.
  // optim_point_N = NGC(stochastic_logistic_regression, 5, 10, 1e-2);
  // imprimeTit("Function minimum (NCG):");
  // imprimeMatriz(optim_point_N, 1, length);


  // Prediction error.
  // imprimeTit("Class Error:");
  // printf(" %.5lf \n", class_error(optim_point_N, length));


  /*
   * ###############################################################
   * Test LBFGS
   * ###############################################################
   */

  // Print result easy.
  optim_point_lbfgs = LBFGS(test_func, length, 20, 1e-6, verbose);
  imprimeTit("Function minimum (LBFGS):");
  imprimeMatriz(optim_point_lbfgs, 1, length);

  // Print results hard.
  optim_point_lbfgs = LBFGS(testFunc, length, 20, 1e-6, verbose);
  imprimeTit("Function minimum (LBFGS):");
  imprimeMatriz(optim_point_lbfgs, 1, length);


  // Test logistic.
  // optim_point_N = LBFGS(stochastic_logistic_regression, 5, 10, 1e-2);
  // imprimeTit("Function minimum (LBFGS):");
  // imprimeMatriz(optim_point_N, 1, length);


  // Prediction error.
  // imprimeTit("Class Error:");
  // printf(" %.5lf \n", class_error(optim_point_N, length));


  /*
   * ###############################################################
   * Test Stochastically Initialized LFBGS
   * ###############################################################
   */
  // Print result easy.
  optim_point_slm_lbfgs = SLM_LBFGS(test_func, length, 20, 1e-4, 20, verbose);
  imprimeTit("Function minimum (SLM-LBFGS):");
  imprimeMatriz(optim_point_slm_lbfgs, 1, length);

  // Print results hard.
  optim_point_slm_lbfgs = SLM_LBFGS(testFunc, length, 20, 1e-4, 20, verbose);
  imprimeTit("Function minimum (SLM-LBFGS):");
  imprimeMatriz(optim_point_slm_lbfgs, 1, length);

  return 0;
}


/* -------------------------------------
 * Logistic regression.
 * ## Characterisics ##
 * Easy
 * -------------------------------------
 */
double stochastic_logistic_regression(double* x, int length){
  // Variable declaration.
  double res;
  int *y, *indexes;
  int i, samp_size, k;
  // Initialize samp_size and res.
  k = 10; // Size of sample ... proper subset of dataset.
  samp_size = rand() % (MAX_FILE_ROWS - k);
  res = 0;
  // Memory allocation.
  indexes = (int*) malloc(samp_size * sizeof(int));
  // Generate array of indexes.
  for(i = 0; i < samp_size; i++){
    indexes[i] = rand() % MAX_FILE_ROWS;
  }
  // Evaluate logistic Error.
  for(i = 0; i < samp_size; i++){
    res = res + log(1 +
                    exp(-dotProd(x, (double*) logistic_values[indexes[i]], 5) *
                        logistic_labels[indexes[i]]));
  }
  return res;
};

/* -------------------------------------
 * Multiclass Logistic Function
 * -------------------------------------
 * Theta = array of K x N. K = Númber of
 *         classes. N = Dimension of each
 *         observation.
 */
double softmax(double* theta, int length){
  double *theta_dot;
  double score, denom;
  int i, k, j;
  // Space allocation.
  theta_dot = (double*)malloc(length*sizeof(double));
  // Construct denom
  for(denom = i = 0; i < MAX_FILE_ROWS; i++){
    for(k = 0; k < N_CLASS; k++){
      // Theta dot construction.
      for(j = 0; j < length; j++){
        theta_dot[j] = theta[(length * k) + j];
      }
      denom = denom + exp(dotProd(theta_dot, logistic_values[i], length));
    }
  }

  // Construct score
  for(score = i = 0; i < MAX_FILE_ROWS; i++){
    for(k = 0; k < N_CLASS; k++){
      // Theta dot construction.
      for(j = 0; j < length; j++){
        theta_dot[j] = theta[length * k + j];
      }
      score = score + (double)(logistic_labels[i] == k) * log(exp(dotProd(theta_dot, logistic_values[i], length)) / denom);
    }
  }
  return -score;
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
