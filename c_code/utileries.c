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
#define BIG_N 1000000

// ACTUAL DATA VARIABLES

int const MAX_FILE_ROWS = BIG_N;//7000000;
int const MAX_FILE_COLS = 27;
int    logistic_labels[BIG_N];
double logistic_values[BIG_N][27];
// SAMPLE VARIABLES GRADIENT
double sample_logistic_values[BIG_N][27];
int    sample_logistic_labels[BIG_N];

// DEFAULT CONFIGURATION
int SAMPLE = 0;
double sampProp       = 1; //.5
double regularization = .0001; //.0001
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
 * EASY
 * -------------------------------------
 */
double test_func1(double* x, int length){
  int i;
  double sum;
  sum = 0;
  for(i = 0; i < length; i++){
    sum = sum  + (100 - i) * (x[i] * x[i]);
  }
  // Return result.
  return sum;
};


/* -------------------------------------
 * Test function:
 * ## Characteristics ##
 * HARD
 * -------------------------------------
 */
double test_func2(double* x, int length){
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
    printf("Sample size: %d \n", SAMPLE);
  }
  // Memory allocation.
  indexes = (int*) malloc(SAMPLE * sizeof(int));
  // Random indexes construction.
  for(i = 0; i < SAMPLE; i++){
    indexes[i] = rand() % MAX_FILE_ROWS;
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
    if(stocMode){
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
      printf("2)       - REGULARIZATION: %lf\n", regularization);
      printf("3)       - STOCASTIC OPTIM: %d\n", stocMode);
      printf("4) EXECUTE FUNCTIONS: %d\n", run_functions);
      printf("5) VERBOSE MODE: %d\n", verbose);
      printf("#################################################\n");
    }
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

/* -------------------------------------
 * Print Configuration
 * -------------------------------------
 * Iteratively print menu until optimal
 * configuration.
 */
void menu(){
  int option;
  int correct = 0;
  char isCorrect;
  do{
  choose:
    printConfig();
    printf("Is this configuration correct (y/n)?\n");
    scanf(" %c", &isCorrect);
    if(isCorrect == 'y'){
      correct = 1;
    }else if(isCorrect == 'n'){
      correct = 0;
    }else{
      printf("\nPlease choose (y)es or (n)o:");
      goto choose;
    }
    if(!correct){
      printf("Choose the option you want to change:\n");
      scanf("%d", &option);
      switch(option){
      case 1:
        printf("Run logistic? (0/1)\n");
        scanf("%d", &run_logistic);
        break;
      case 2:
        if(!run_logistic){
          printf("Run functions? (0/1)\n");
          scanf("%d", &run_functions);
        }else{
          printf("Regularization Parameter?\n");
          scanf("%lf", &regularization);
        }
        break;
      case 3:
        if(!run_logistic){
          printf("Verbose? (0/1)\n");
        scanf("%d", &verbose);
        }else{
          printf("Stochastic Mode? (0/1)\n");
          scanf("%d", &stocMode);
        }
        break;
      case 4:
        if(!run_logistic){
          printf("Please choose a number between 1 and 3 \n");
        }else{
          if(stocMode){
            printf("Random seed?\n");
            scanf("%d", &seed);
          }else{
            printf("Run functions? (0/1)\n");
            scanf("%d", &run_functions);
          }
        }
        break;
      case 5:
        if(!run_logistic){
          printf("Please choose a number between 1 and 3 \n");
        }else{
          if(stocMode){
            printf("Sample proportion of dataset?\n");
            scanf("%lf", &sampProp);
          }else{
            printf("Verbose? (0/1)\n");
            scanf("%d", &verbose);
          }
        }
        break;
      case 6:
        if(!run_logistic){
          printf("Please choose a number between 1 and 3 \n");
        }else{
          if(stocMode){
            printf("Run functions? (0/1)\n");
            scanf("%d", &run_functions);
          }else{
            printf("Please choose a number between 1 and 5 \n");
          }
        }
        break;
      case 7:
        if(!run_logistic){
          printf("Please choose a number between 1 and 3 \n");
        }else{
          if(stocMode){
            printf("Verbose? (0/1)\n");
            scanf("%d", &verbose);
          }else{
            printf("Please choose a number between 1 and 5 \n");
          }
        }
        break;
      default:
        if(!run_logistic){
          printf("Please choose a number between 1 and 3 \n");
        }else{
          if(stocMode){
            printf(" Please choose a number between 1 and 7 \n");
          }else{
            printf("Please choose a number between 1 and 5 \n");
          }
        }
        break;
      }
    }
  }while(!correct);
  srand(seed);
}

/* -------------------------------------
 * READ FILE
 * -------------------------------------
 * Read file.
 * File that is to be read.
 */
void readFile(){
  int i;
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
    if(/*verbose*/0 && (i % (int)(MAX_FILE_ROWS*.1)) == 0){
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
}
