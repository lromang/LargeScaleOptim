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
int    logistic_labels[150];
double logistic_values[150][5];
double stochastic_softmax(double*, int);
double class_error(double*, int);


int main(){
  // Variable declaration.
  double *optim_point_N, *optim_point_lbfgs, *optim_point_slm_lbfgs;
  int i, length, seed, verbose, run_logistic;
  // Run logistic???
  printf("Run logistic?\n");
  scanf("%d", &run_logistic);
  printf("\nVerbose?\n");
  scanf("%d", &verbose);
  // Random seed for Softmax Regression.
  seed = 34234234;
  srand(seed);
  // Size of point.
  length  = 100;
  // Data file name.
  FILE *file = fopen("../data/iris", "r");
  // Read in file
  for(i = 0; i < MAX_FILE_ROWS; i++){
      if (feof(file))
        break;
      fscanf(file, "%lf %lf %lf %lf %d",
             &(logistic_values[i][0]),
             &(logistic_values[i][1]),
             &(logistic_values[i][2]),
             &(logistic_values[i][3]),
             &(logistic_labels[i]));
      if(verbose){
      printf("Entry: %d | col1 = %lf  col2 = %lf  col3 = %lf  col4 = %lf  col5 = %d \n",
             i,
             logistic_values[i][0],
             logistic_values[i][1],
             logistic_values[i][2],
             logistic_values[i][3],
             logistic_labels[i]);
      }
    }

  /*
   * ###############################################################
   * Test Truncatd Newton
   * ###############################################################
   */

  // Print results easy.
  optim_point_N = NGC(test_func, length, 10, 1e-2, verbose);
  imprimeTit("Problem 1 minimum (NCG):");
  imprimeMatriz(optim_point_N, 1, length);

  // Print results hard.
  optim_point_N = NGC(testFunc, length, 10, 1e-2, verbose);
  imprimeTit("Problem 2 minimum (NCG):");
  imprimeMatriz(optim_point_N, 1, length);

  /*
   * ###############################################################
   * Test LBFGS
   * ###############################################################
   */

  // Print result easy.
  optim_point_lbfgs = LBFGS(test_func, length, 20, 1e-6, verbose);
  imprimeTit("Problem 1 minimum (LBFGS):");
  imprimeMatriz(optim_point_lbfgs, 1, length);

  // Print results hard.
  optim_point_lbfgs = LBFGS(testFunc, length, 20, 1e-6, verbose);
  imprimeTit("Problem 2 minimum (LBFGS):");
  imprimeMatriz(optim_point_lbfgs, 1, length);

  /*
   * ###############################################################
   * Test Stochastically Initialized LFBGS
   * ###############################################################
   */
  // Print result easy.
  optim_point_slm_lbfgs = SLM_LBFGS(test_func, length, 20, 1e-4, 20, verbose);
  imprimeTit("Problem 1:  minimum (SLM-LBFGS):");
  imprimeMatriz(optim_point_slm_lbfgs, 1, length);

  // Print results hard.
  optim_point_slm_lbfgs = SLM_LBFGS(testFunc, length, 20, 1e-4, 20, verbose);
  imprimeTit("Problem 2: minimum (SLM-LBFGS):");
  imprimeMatriz(optim_point_slm_lbfgs, 1, length);

  /*
   * ###############################################################
   * Test Logistic
   * ###############################################################
   */

  if(run_logistic){
    // Test logistic.
    optim_point_N = NGC(stochastic_softmax, 5, 20, 1e-2, verbose);
    imprimeTit("Multinomial Logistic minimum (NCG):");
    imprimeMatriz(optim_point_N, 1, length);

    // Prediction error.
    imprimeTit("Class Error:");
    printf(" %.5lf \n", class_error(optim_point_N, length));

    // Test multinomial logistic.
    optim_point_N = LBFGS(stochastic_softmax, 5, 20, 1e-2, verbose);
    imprimeTit("Multinomial Logistic minimum (LBFGS):");
    imprimeMatriz(optim_point_N, 1, length);

    // Prediction error.
    imprimeTit("Class Error:");
    printf(" %.5lf \n", class_error(optim_point_N, length));

    // Test multinomial logistic.
    optim_point_N = SLM_LBFGS(stochastic_softmax, 5, 20, 1e-2, 20, verbose);
    imprimeTit("Multinomial Logistic minimum (SLM-LBFGS):");
    imprimeMatriz(optim_point_N, 1, length);

    // Prediction error.
    imprimeTit("Class Error:");
    printf(" %.5lf \n", class_error(optim_point_N, length));
  }

  return 0;
}

/* -------------------------------------
 * Multiclass Logistic Function
 * -------------------------------------
 * Theta = array of K x N. K = Númber of
 *         classes. N = Dimension of each
 *         observation.
 */
double stochastic_softmax(double* theta, int length){
  double *theta_dot;
  int *indexes;
  double score, denom;
  int i, k, j, samp_size, proper, deterministic;
  // Initialize samp_size and res.
  deterministic = 1;
  proper        = 10; // Upper bound in size for proper subset of dataset.
  samp_size     = rand() % (MAX_FILE_ROWS - proper);
  // Memory allocation.
  indexes = (int*) malloc(samp_size * sizeof(int));
  if(deterministic){
    // Generate array of indexes.
    for(i = 0; i < samp_size; i++){
      indexes[i] = rand() % MAX_FILE_ROWS;
    }
  }else{
    for(i = 0; i < MAX_FILE_ROWS; i++){
      indexes[i] = i;
    }
  }
  /* ------------------------------
   * Softmax error evaluation.
   * ------------------------------
   */

  // Space allocation.
  theta_dot = (double*)malloc(length*sizeof(double));
  // Construct denom
  for(denom = i = 0; i < samp_size; i++){
    for(k = 0; k < N_CLASS; k++){
      // Theta dot construction.
      for(j = 0; j < length; j++){
        theta_dot[j] = theta[(length * k) + j];
      }
      denom = denom + exp(dotProd(theta_dot, logistic_values[indexes[i]], length));
    }
  }

  // Construct score
  for(score = i = 0; i < samp_size; i++){
    for(k = 0; k < N_CLASS; k++){
      // Theta dot construction.
      for(j = 0; j < length; j++){
        theta_dot[j] = theta[length * k + j];
      }
      score = score + (double)(logistic_labels[indexes[i]] == k) *
        log(exp(dotProd(theta_dot, logistic_values[indexes[i]], length)) / denom);
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
