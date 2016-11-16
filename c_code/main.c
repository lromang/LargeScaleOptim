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
double logSumExp(double*, int, int);

int main(){
  // Variable declaration.
  double *optim_point_N, *optim_point_lbfgs, *optim_point_slm_lbfgs;
  double precision;
  int length;
  // Print options
  // menu();
  /*
   * ###############################################################
   * Test Functions
   * ###############################################################
   */
  if(run_functions){
    // Size of point.
    length  = 100;
    // Print results easy.
    optim_point_N = NGC(test_func1, length, 10, 1e-2, verbose, 100, 1, .0001);
    imprimeTit("Problem 1 minimum (NCG):");
    imprimeMatriz(optim_point_N, 1, length);
    // Print results hard.
    optim_point_N = NGC(test_func2, length, 10, 1e-2, verbose, 100, 1, .0001);
    imprimeTit("Problem 2 minimum (NCG):");
    imprimeMatriz(optim_point_N, 1, length);
    // Print result easy.
    optim_point_lbfgs = LBFGS(test_func1, length, 20, 1e-6, verbose);
    imprimeTit("Problem 1 minimum (LBFGS):");
    imprimeMatriz(optim_point_lbfgs, 1, length);
    // Print results hard.
    optim_point_lbfgs = LBFGS(test_func2, length, 20, 1e-6, verbose);
    imprimeTit("Problem 2 minimum (LBFGS):");
    imprimeMatriz(optim_point_lbfgs, 1, length);
    // Print result easy.
    optim_point_slm_lbfgs = SLM_LBFGS(test_func1, length, 20, 1e-4, 10, verbose);
    imprimeTit("Problem 1:  minimum (SLM-LBFGS):");
    imprimeMatriz(optim_point_slm_lbfgs, 1, length);
    // Print results hard.
    optim_point_slm_lbfgs = SLM_LBFGS(test_func2, length, 20, 1e-4, 10, verbose);
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
    // Test logistic. // ADD THIS CONFIGURATIONS TO CODE
    optim_point_N = NGC(logistic, length, 2, 1e-2, verbose, 0, 100, .0001);
    imprimeTit("Logistic minimum (NCG):");
    imprimeMatriz(optim_point_N, 1, length);
    // Prediction error.
    precision = class_precision(optim_point_N, length, 0);
    printf("\n");
    imprimeTit("Classification Precision:");
    printf("%.5lf \n", precision);
    /**/
    /*
    // RUNNING LBFGS MODEL
    imprimeTit("RUNNING LBFGS MODEL");
    // Test multinomial logistic.
    optim_point_N = LBFGS(logistic, length, min(length, 15), 1e-5, verbose);
    imprimeTit("Multinomial Logistic minimum (LBFGS):");
    imprimeMatriz(optim_point_N, 1, length);
    // Prediction error.
    precision = class_precision(optim_point_N, length, 0);
    printf("\n");
    imprimeTit("Classification Precision (LBFGS):");
    printf("%.5lf \n", precision);
    */
    // RUNNING SLM-LBFGS MODEL
    /*
    imprimeTit("RUNNING SLM-LBFGS MODEL");
    // Test multinomial logistic.
    optim_point_N = SLM_LBFGS(logistic, length, 2, 1e-2, 5, verbose);
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
double logSumExp(double* theta, int i, int length){
  if(!stocMode){
    if(-log(1 + exp(pow(-1, logistic_labels[i])*
                    dotProd(logistic_values[i], theta, length))) <  -1e30){
      return -1e10; // Numerical stability...
    }
    return -log(1 + exp(pow(-1, logistic_labels[i])*
                        dotProd(logistic_values[i], theta, length)));
  }
  if(-log(1 + exp(pow(-1, sample_logistic_labels[i])*
                    dotProd(sample_logistic_values[i], theta, length))) <  -1e30){
      return -1e10; // Numerical stability...
  }
  return -log(1 + exp(pow(-1, sample_logistic_labels[i])*
                        dotProd(sample_logistic_values[i], theta, length)));
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
      // printf("logistic_label[%d] = %d | logSumExp[%d] = %lf\n", i, logistic_labels[i], i, logSumExp(theta,  i, length));
      loss = loss + logistic_labels[i]*logSumExp(theta, i, length) +
        (1 - logistic_labels[i])*logSumExp(theta, i, length) +
        regularization*dotProd(theta, theta, length);
    }
  }else{
    for(i = 0; i < SAMPLE; i++){
      loss = loss + sample_logistic_labels[i]*logSumExp(theta, i, length) +
        (1 - sample_logistic_labels[i])*logSumExp(theta, i, length) +
        regularization*dotProd(theta, theta, length);
    }
  }
  return -loss;
}
