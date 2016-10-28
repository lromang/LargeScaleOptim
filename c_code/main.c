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
int const MAX_FILE_ROWS = 479;
int const MAX_FILE_COLS = 7;
int const N_CLASS = 3;

// Value storage.
int    logistic_labels[479];
double logistic_values[479][7];
double stochastic_softmax(double*, int);
double logActive(double*, int, int);
double logistic(double*, int);
double class_precision(double*, int, int);


int main(){
  // Variable declaration.
  double *optim_point_N, *optim_point_lbfgs, *optim_point_slm_lbfgs;
  double precision;
  int i, length, seed, verbose, run_logistic, run_functions;
  // Run logistic???
  printf("Run logistic?\n");
  scanf("%d", &run_logistic);
  printf("Run functions?\n");
  scanf("%d", &run_functions);
  printf("Verbose?\n");
  scanf("%d", &verbose);
  // Random seed for Softmax Regression.
  seed = 34234234;
  srand(seed);
  // Size of point.
  length  = 100;

  /*
   * ###############################################################
   * Test Functions
   * ###############################################################
   */

  if(run_functions){
    // Print results easy.
    optim_point_N = NGC(test_func, length, 10, 1e-2, verbose);
    imprimeTit("Problem 1 minimum (NCG):");
    imprimeMatriz(optim_point_N, 1, length);

    // Print results hard.
    optim_point_N = NGC(testFunc, length, 10, 1e-2, verbose);
    imprimeTit("Problem 2 minimum (NCG):");
    imprimeMatriz(optim_point_N, 1, length);

    // Print result easy.
    optim_point_lbfgs = LBFGS(test_func, length, 20, 1e-6, verbose);
    imprimeTit("Problem 1 minimum (LBFGS):");
    imprimeMatriz(optim_point_lbfgs, 1, length);

    // Print results hard.
    optim_point_lbfgs = LBFGS(testFunc, length, 20, 1e-6, verbose);
    imprimeTit("Problem 2 minimum (LBFGS):");
    imprimeMatriz(optim_point_lbfgs, 1, length);

    // Print result easy.
    optim_point_slm_lbfgs = SLM_LBFGS(test_func, length, 20, 1e-4, 20, verbose);
    imprimeTit("Problem 1:  minimum (SLM-LBFGS):");
    imprimeMatriz(optim_point_slm_lbfgs, 1, length);

    // Print results hard.
    optim_point_slm_lbfgs = SLM_LBFGS(testFunc, length, 20, 1e-4, 20, verbose);
    imprimeTit("Problem 2: minimum (SLM-LBFGS):");
    imprimeMatriz(optim_point_slm_lbfgs, 1, length);
  }

  /*
   * ###############################################################
   * Test Logistic
   * ###############################################################
   */

  if(run_logistic){
    length  = MAX_FILE_COLS;
    FILE *file = fopen("../data/ecoli", "r");
    // Read in file
    for(i = 0; i < MAX_FILE_ROWS; i++){
      if (feof(file))
        break;
      fscanf(file, "%lf %lf %lf %lf %lf %lf %lf %d",
             &(logistic_values[i][0]),
             &(logistic_values[i][1]),
             &(logistic_values[i][2]),
             &(logistic_values[i][3]),
             &(logistic_values[i][4]),
             &(logistic_values[i][5]),
             &(logistic_values[i][6]),
             &(logistic_labels[i]));
      if(verbose){
      printf("Entry: %d | col1 = %lf  col2 = %lf  col3 = %lf  col4 = %lf  col5 = %lf col6 = %lf  col7 = %lf col8 = %d \n",
             i,
             logistic_values[i][0],
             logistic_values[i][1],
             logistic_values[i][2],
             logistic_values[i][3],
             logistic_values[i][4],
             logistic_values[i][5],
             logistic_values[i][6],
             logistic_labels[i]);
      }
    }

    /*
    // Test logistic.
    optim_point_N = NGC(logistic, length, 100, 1e-2, verbose);
    imprimeTit("Multinomial Logistic minimum (NCG):");
    imprimeMatriz(optim_point_N, 1, length);

    // Prediction error.
    precision = class_precision(optim_point_N, length, verbose);
    printf("\n");
    imprimeTit("Classification Precision:");
    printf("%.5lf \n", precision);
    */

    // Test multinomial logistic.
    optim_point_N = LBFGS(logistic, length, 20, 1e-2, verbose);
    imprimeTit("Multinomial Logistic minimum (LBFGS):");
    imprimeMatriz(optim_point_N, 1, length);

    // Prediction error.
    precision = class_precision(optim_point_N, length, verbose);
    printf("\n");
    imprimeTit("Classification Precision (LBFGS):");
    printf("%.5lf \n", precision);
    /*
    // Test multinomial logistic.
    optim_point_N = SLM_LBFGS(logistic, length, 7, 1e-6, 100, verbose);
    imprimeTit("Multinomial Logistic minimum (SLM-LBFGS):");
    imprimeMatriz(optim_point_N, 1, length);

    // Prediction error.
    precision = class_precision(optim_point_N, length, verbose);
    printf("\n");
    imprimeTit("Classification Precision (SLM_LBFGS):");
    printf("%.5lf \n", precision);
  */
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
  double *theta_dot, *denom;
  double score;
  int i, k, j;
  // Space allocation.
  theta_dot = (double*) malloc(length * sizeof(double));
  denom     = (double*) malloc(length * sizeof(double));
  // Fill in denom empty values.
  for(i = 0; i < MAX_FILE_ROWS; i++){
    denom[i] = 0;
  }
  // Construct denom
  for(i = 0; i < MAX_FILE_ROWS; i++){
    for(k = 0; k < N_CLASS; k++){
      // Theta dot construction.
      for(j = 0; j < length; j++){
        theta_dot[j] = theta[(length * k) + j];
      }
      denom[i] = denom[i] + exp(dotProd(theta_dot, logistic_values[i], length));
    }
  }

  // Construct score
  for(score = i = 0; i < MAX_FILE_ROWS; i++){
    for(k = 0; k < N_CLASS; k++){
      // Theta dot construction.
      for(j = 0; j < length; j++){
        theta_dot[j] = theta[length * k + j];
      }
      score = score + (double)(logistic_labels[i] == k) *
        log(exp(dotProd(theta_dot, logistic_values[i], length)) / denom[i]);
    }
  }
  printf("Score: %lf\n", score);
  free(theta_dot);
  return -score;
};

/* -------------------------------------
 * Logistic Activation
 * -------------------------------------
 * Theta = array of K x N. K = Númber of
 *         classes. N = Dimension of each
 *         observation.
 */
double logActive(double* theta, int i, int length){
  return 1 / (1 + exp(-dotProd(logistic_values[i], theta, length)));
}

/* -------------------------------------
 * Binary Logistic
 * -------------------------------------
 * Theta = array of K x N. K = Númber of
 *         classes. N = Dimension of each
 *         observation.
 */
double logistic(double* theta, int length){
  double loss = 0;
  int *indexes;
  int i, stochastic, samp_size, proper;
  stochastic  = 0; // Sample size should be defined outside?
  proper      = 10; // N - proper = size of batch
  // Check if framework is stochastic
  if(stochastic){
    samp_size     = rand() % (MAX_FILE_ROWS - proper);
    printf("Tamaño muestra: %d \n", samp_size);
    // Memory allocation.
    indexes = (int*) malloc(samp_size * sizeof(int));
    // Indexes construction.
    for(i = 0; i < samp_size; i++){
      indexes[i] = rand() % samp_size;
    }
    // Logistic
    for(i = 0; i < samp_size; i++){
      loss = loss + logistic_labels[indexes[i]]*log(logActive(theta, indexes[i], length)) +
        (1 - logistic_labels[indexes[i]])*log(1 - logActive(theta, indexes[i], length));
    }
  }else{
    for(i = 0; i < MAX_FILE_ROWS; i++){
      loss = loss + logistic_labels[i]*log(logActive(theta, i, length)) +
        (1 - logistic_labels[i])*log(1 - logActive(theta, i, length));
    }
  }
  return -loss;
}


/*
 * -------------------------------------
 * Eval function
 * -------------------------------------
 */
double class_precision(double* coefs, int length, int verbose){
  // Variable declaration.
  double class_error, entry_val;
  int i, pred;
  // Evaluate logistic.
  pred = 0;
  for(i = 0; i < MAX_FILE_ROWS; i++){
    entry_val = exp(- dotProd(coefs, logistic_values[i], length));
    entry_val = 1 / (1 + entry_val);
    // Classification threshold = .5
    if(entry_val > .5){
      pred = 1;
    }else{
      pred = 0;
    }
    if(verbose){
      imprimeTit("PRED");
      printf("ACTIVATION VALUE: %lf | PREDICTION: %d | ACTUAL: %d", entry_val, pred, logistic_labels[i]);
    }
    class_error = class_error + (pred == logistic_labels[i] ? 1 : 0);
    //res = res + log(1 + exp(-dotProd(x, (double*) logistic_values[i], 5) * logistic_labels[i]));
  }
  return class_error / MAX_FILE_ROWS;
};
