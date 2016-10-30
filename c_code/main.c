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
int const MAX_FILE_ROWS = 20000;//7000000;
int const MAX_FILE_COLS = 27;
int const N_CLASS = 2;

// Value storage.
int    logistic_labels[20000];
double logistic_values[20000][27];
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
    FILE *file = fopen("../data/higgs/train_clean", "r");
    // Read in file
    for(i = 0; i < MAX_FILE_ROWS; i++){
      if (feof(file))
        break;
      fscanf(file, "%d %lf %lf %lf %lf %lf %lf %lf %lf  %lf"
             "%lf %lf %lf %lf %lf %lf %lf %lf  %lf"
             "%lf %lf %lf %lf %lf %lf %lf %lf  %lf",
             &(logistic_labels[i]),     &(logistic_values[i][0]),  &(logistic_values[i][1]),  &(logistic_values[i][2]),
             &(logistic_values[i][3]),  &(logistic_values[i][4]),  &(logistic_values[i][5]),  &(logistic_values[i][6]),
             &(logistic_values[i][7]),  &(logistic_values[i][8]),  &(logistic_values[i][9]),  &(logistic_values[i][10]),
             &(logistic_values[i][11]), &(logistic_values[i][12]), &(logistic_values[i][13]), &(logistic_values[i][14]),
             &(logistic_values[i][15]), &(logistic_values[i][16]), &(logistic_values[i][17]), &(logistic_values[i][18]),
             &(logistic_values[i][19]), &(logistic_values[i][20]), &(logistic_values[i][21]), &(logistic_values[i][22]),
             &(logistic_values[i][23]), &(logistic_values[i][24]), &(logistic_values[i][25]), &(logistic_values[i][26]));
      if(verbose && (i % 100) == 0){
      printf("Entry: %d | label = %d  col1 = %lf  col2 = %lf  col3 = %lf  col4 = %lf  col5 = %lf  col6 = %lf col7 = %lf col8 = %lf \n "
             "col9 = %lf col10 = %lf col11 = %lf col12 = %lf col13 = %lf col14 = %lf col15 = %lf col16 = %lf col17 = %lf col18 = %lf \n"
             "col19 = %lf col20 = %lf col21 = %lf col22 = %lf col23 = %lf col24 = %lf col25 = %lf col26 = %lf col27 = %lf \n",
             i,
             logistic_labels[i],     logistic_values[i][0],  logistic_values[i][1],  logistic_values[i][2],
             logistic_values[i][3],  logistic_values[i][4],  logistic_values[i][5],  logistic_values[i][6],
             logistic_values[i][7],  logistic_values[i][8],  logistic_values[i][9],  logistic_values[i][10],
             logistic_values[i][11], logistic_values[i][12], logistic_values[i][13], logistic_values[i][14],
             logistic_values[i][15], logistic_values[i][16], logistic_values[i][17], logistic_values[i][18],
             logistic_values[i][19], logistic_values[i][20], logistic_values[i][21], logistic_values[i][22],
             logistic_values[i][23], logistic_values[i][24], logistic_values[i][25], logistic_values[i][26]);
      }
    }

    imprimeTit("RUNNING NGC MODEL");

    // Test logistic.
    optim_point_N = NGC(logistic, length, 30, 6e-1, verbose);
    imprimeTit("Multinomial Logistic minimum (NCG):");
    imprimeMatriz(optim_point_N, 1, length);

    // Prediction error.
    precision = class_precision(optim_point_N, length, 0);
    printf("\n");
    imprimeTit("Classification Precision:");
    printf("%.5lf \n", precision);


    imprimeTit("RUNNING LBFGS MODEL");

    // Test multinomial logistic.
    optim_point_N = LBFGS(logistic, length, min(length, 30), 1e-1, verbose);
    imprimeTit("Multinomial Logistic minimum (LBFGS):");
    imprimeMatriz(optim_point_N, 1, length);

    // Prediction error.
    precision = class_precision(optim_point_N, length, 0);
    printf("\n");
    imprimeTit("Classification Precision (LBFGS):");
    printf("%.5lf \n", precision);

    /*
    imprimeTit("RUNNING SLM-LBFGS MODEL");

    // Test multinomial logistic.
    optim_point_N = SLM_LBFGS(logistic, length, 7, 1e-1, 100, verbose);
    imprimeTit("Multinomial Logistic minimum (SLM-LBFGS):");
    imprimeMatriz(optim_point_N, 1, length);

    // Prediction error.
    precision = class_precision(optim_point_N, length, 0);
    printf("\n");
    imprimeTit("Classification Precision (SLM_LBFGS):");
    printf("%.5lf \n", precision);
    */
  }
  return 0;
}

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
        (1 - logistic_labels[i])*log(1 - logActive(theta, i, length)) + .001*dotProd(theta, theta, length);
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
