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

// Function declaration.
double logActive(double*, int, int);
double logistic(double*, int);
double class_precision(double*, int, int);


int main(){
  // Variable declaration.
  double *optim_point_N, *optim_point_lbfgs, *optim_point_slm_lbfgs;
  double precision;
  int length;

  // Print options
  menu();

  /*
   * ###############################################################
   * Test Functions
   * ###############################################################
   */
  if(run_functions){
    // Size of point.
    length  = 100;
    // Print results easy.
    optim_point_N = NGC(test_func, length, 10, 1e-2, verbose, 100, 1);
    imprimeTit("Problem 1 minimum (NCG):");
    imprimeMatriz(optim_point_N, 1, length);

    // Print results hard.
    optim_point_N = NGC(testFunc, length, 10, 1e-2, verbose, 100, 1);
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
    // READ FILE
    readFile();
    // RUNNING NGC MODEL
    imprimeTit("RUNNING NGC MODEL");
    // Test logistic.
    optim_point_N = NGC(logistic, length, 30, 6e-1, verbose, 100, 1);
    imprimeTit("Logistic minimum (NCG):");
    imprimeMatriz(optim_point_N, 1, length);
    // Prediction error.
    precision = class_precision(optim_point_N, length, 0);
    printf("\n");
    imprimeTit("Classification Precision:");
    printf("%.5lf \n", precision);
    // RUNNING LBFGS MODEL
    /*
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
    */
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
  if(!stocMode){
    return 1 / (1 + exp(-dotProd(logistic_values[i], theta, length)));
  }
  return 1 / (1 + exp(-dotProd(sample_logistic_values[i], theta, length)));
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
  int i;
  if(!stocMode){
    SAMPLE = MAX_FILE_ROWS;
    for(i = 0; i < SAMPLE; i++){
      loss = loss + logistic_labels[i]*log(logActive(theta, i, length)) +
        (1 - logistic_labels[i])*log(1 - logActive(theta, i, length)) +
        regularization*dotProd(theta, theta, length);
    }
  }else{
    for(i = 0; i < SAMPLE; i++){
      loss = loss + sample_logistic_labels[i]*log(logActive(theta, i, length)) +
        (1 - sample_logistic_labels[i])*log(1 - logActive(theta, i, length)) +
        regularization*dotProd(theta, theta, length);
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
