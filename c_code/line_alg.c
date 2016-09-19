/* #########################################
 * Luis Manuel Román García
 * luis.roangarci@gmail.com
 * #########################################
 *
 * -----------------------------------------
 * General purpose linear algebra funcitons.
 * -----------------------------------------
 *
 */
#include<limits.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

/* -------------------------------------
 * Matrix creation:
 * Generates a flat random  matrix
 * with the given dimensions.
 * IN
 * nRow:   Number of rows in matrix.
 * nCol:   Number of columns in matrix.
 * OUT
 * matrix: Pointer to array
 * -------------------------------------
 */
double * creaMatriz(int nRow, int nCol){
  // Variable declaration.
  double * matrix;
  int i, j, k;
  // Space allocation.
  matrix = (double*)(malloc((nRow * nCol) * sizeof(double)));
  k = 0;
  for(i = 0; i < nRow; i++){
    for(j = 0; j < nCol; j++){
      matrix[k] = rand() % 10;
      k = k + 1;
    }
  }
  // Return result.
  return matrix;
}

/* -------------------------------------
 * Matrix printing:
 * Prints a matrix with the given
 * dimensions
 * IN
 * A:    Matrix to be printed.
 * nRow: Number of rows of A.
 * nCol: Number of columns of A.
 * -------------------------------------
 */
void imprimeMatriz(double* A, int nCol, int nRow){
  // Variable declaration.
  int i, j, k;
  // Double loop to print each entry of A.
  k = 0;
  for(i = 0; i < nRow; i++){
    for(j = 0; j < nCol; j++){
      printf(" %.5lf ", A[k]);
      k = k + 1;
    }
    printf("\n");
  }
  printf("\n");
}


/* -------------------------------------
 * Product between scalar and vector:
 * Performs the product between a scalar
 * and a vector.
 * IN
 * v:      Vector to multiply.
 * alpha:  Scalar to multiply.
 * length: Length of vector.
 * OUT
 * prod_v: Vector with entries
 *         i = v_i*alpha.
 * -------------------------------------
 */
double* vProd(double* v, double alpha, int length){
  // Variable declaration.
  double* prod_v;
  int i;
  // Space allocation.
  prod_v = (double*) malloc(length * sizeof(double));
  // Entry-wise product.
  for(i = 0; i < length; i ++){
    prod_v[i] = v[i] * alpha;
  }
  // Return result.
  return prod_v;
};


/* -------------------------------------
 * Matrix transpose:
 * Generates the transpose matrix of A.
 * IN
 * A:    Matrix to be transposed.
 * nCol: Number of columns of A.
 * nRow: Number of rows of A.
 * OUT
 * A_trans: Pointer to A transpose.
 * -------------------------------------
 */
double* mTrans(double* A, int nCol, int nRow){
  // Variable declaration.
  double* A_trans;
  int i, j;
  // Space allocation.
  A_trans = (double*)malloc(nCol * nRow * sizeof(double));
  // Traspose Aij = Atji
  for(i = 0; i < nRow; i++){
    for(j = 0; j < nCol; j++){
      A_trans[i * nCol + j] = A[j * nCol + i];
    }
  }
  // Return result.
  return A_trans;
};

/* -------------------------------------
 * Vector-vector addition.
 * IN
 * v:      First vector to be added.
 * u:      Second vector to be added.
 * length: Vectors' length.
 * OUT
 * sum_v: Pointer to vector with entries
 *        i = u_i + v_i.
 * -------------------------------------
 */
double* vSum(double* v, double* u, int length){
  // Variable declaration.
  double* sum_v;
  int i;
  // Space allocation.
  sum_v = (double*)malloc(length * sizeof(double));
  // Entry-wise addition.
  for(i = 0; i < length; i ++){
    sum_v[i] = v[i] + u[i];
  }
  // Return result.
  return sum_v;
}


/* -------------------------------------
 * Vector comparison.
 * Determines if two vectors are equal.
 * IN
 * v:      Vector to be compared.
 * u:      Vector to be compared.
 * length: Vectors' length.
 * OUT
 * eq: 1 if the vectors are equal 0
 *     otherwise.
 * -------------------------------------
 */
int vEq(double* v, double* u, int length){
  // Variable declaration.
  int i;
  // Entry-wise comparison
  for(i = 0; (i < length) && (v[i] == u[i]); i ++);
  // Return result.
  return (i == length);
}

/* -------------------------------------
 * Dot product between two vectors.
 * IN
 * v:      Vector to be multiplied.
 * u:      Vector to be multiplied.
 * length: Vectors' length.
 * OUT
 * sum: Result of de dot product.
 * -------------------------------------
 */
double dotProd(double* v, double* u, int length){
  // Variable declaration.
  double sum;
  int i;
  // Addition.
  sum = 0;
  for(i = 0; i < length; i ++){
    sum = sum + v[i] * u[i];
  }
  // Return result.
  return sum;
}

/* -------------------------------------
 * Dot product between two vectors.
 * IN
 * v:      Vector to be multiplied.
 * u:      Vector to be multiplied.
 * length: Vectors' length.
 * OUT
 * a matrix with the dot prod of u and v
 * -------------------------------------
 */
double * dotProdMat(double* v, double* u, int length){
  // Variable declaration.
  double *dotProd, *aux_array;
  int i, j, k;
  // Space allocation.
  dotProd = (double * )malloc((length * length) * sizeof(double));
  for(k = i = 0; i < length; i++){
    aux_array = vProd(u, v[i], length);
    for(j = 0; j < length; j++){
      dotProd[k] = aux_array[j];
      k ++;
    }
  }
  // Return result.
  return dotProd;
}


/* -------------------------------------
 * Matrix vector product.
 * Carries out the product between a
 * matrix B and a vector v.
 * IN
 * B:    Matrix to be multiplied.
 * v:    Vector to be multiplied.
 * nRow: Number of rows of B.
 * nCol: Number of columns of B.
 * OUT
 * Bd: Pointer to product result.
 * -------------------------------------
 */
double* mProd(double* B, double* v, int nCol, int nRow){
  // Variable declaration.
  double* Bd;
  double sum;
  int i, j;
  // Space allocation.
  Bd = (double*)malloc(nRow * sizeof(double));
  // Double loop for product between rows of B and v.
  for(i = 0; i < nRow; i++){
    sum = 0;
    for(j = 0; j < nCol; j++){
      // Bj*v cumulative.
      sum = sum + B[i * nCol + j] * v[j];
    }
    Bd[i] = sum;
  }
  // Return result.
  return Bd;
}


/* -------------------------------------
 * Norm of vector.
 * Obtains the l2 norm of a vector.
 * IN
 * x:      Vector whose norm is to be
 *         obtained.
 * length: x's length.
 * OUT
 * l2 norm of x.
 * -------------------------------------
 */
double norm(double * x, int length){
  return sqrt(dotProd(x, x, length));
}


/* -------------------------------------
 * Identity matrix generation.
 * Generates an identity matrix with the
 * given direction.
 * IN
 * m: Number of columns and rows of the
 *    identity matrix.
 * OUT
 * id: Pointer to the identity matrix.
 * -------------------------------------
 */
double* identity(int m){
  double * id;
  int i, j, k;
  // Space allocation.
  id = (double*) malloc((m * m)* sizeof(double));
  for(k = i = 0; i < m; i ++){
    for(j = 0; j < m; j++){
      i == j ? id[k] = 1 : 0;
      k ++;
    }
  }
  // Return result.
  return id;
}


/* -------------------------------------
 * Minimum between two numbers.
 * Obtains the minimum between two
 * numbers.
 * IN
 * x: Number to be compared.
 * y: Number to be compared.
 * OUT
 * The minimum between x and y.
 * -------------------------------------
 */
double min (double x, double y){
  return x > y ? y : x;
}
