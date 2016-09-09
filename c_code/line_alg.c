/*
 * Luis Manuel Román García
 *
 * -------------------------------------
 * Rutinas para llevar a cabo
 * manipulación de vectores y matrices.
 * -------------------------------------
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

/* -------------------------------------
 * Crea Matriz
 * IN
 * nRow: número de filas que
 * debe tener la matriz
 * nCol: número de columnas que
 * debe tener la matriz
 * -------------------------------------
 */
double * creaMatriz(int nRow, int nCol){
  // Declarar variables.
  double * matriz;
  int i, j, k;
  // Alocar espacio.
  matriz = (double*)(malloc((nRow * nCol) * sizeof(double)));
  k = 0;
  for(i = 0; i < nRow; i++){
    for(j = 0; j < nCol; j++){
      matriz[k] = rand() % 10;
      k = k + 1;
    }
  }
  return matriz;
}

/* -------------------------------------
 * Producto escalar vector
 * IN
 * A: matriz a imprimir.
 * nCol: Número de columnas.
 * nRow: Número de filas
 * -------------------------------------
 */
void imprimeMatriz(double* A, int nCol, int nRow){
  // Declaración de variables.
  int i, j, k;

  // Imprimir matriz.
  k = 0;
  for(i = 0; i < nRow; i++){
    for(j = 0; j < nCol; j++){
      printf(" %.2lf ", A[k]);
      k = k + 1;
    }
    printf("\n");
  }
  printf("\n");
}


/* -------------------------------------
 * Producto escalar vector
 * IN
 * v: Vector a multiplicar.
 * alpha: Escalar a multiplicar
 * length: Longitud del vector.
 * OUT
 * prod_v: Arreglo que contiene el producto
 * de v y alpha entrada a entrada.
 * -------------------------------------
 */
double* vProd(double* v, double alpha, int length){
  // Declaración de variables.
  double* prod_v;
  int i;

  // Inicialización de arreglo
  prod_v = (double*) malloc(length * sizeof(double));

  // Multiplicación
  for(i = 0; i < length; i ++){
    prod_v[i] = v[i] * alpha;
  }

  return prod_v;
};


/* -------------------------------------
 * Trasponer matriz
 * IN
 * A: Matriz a transponer.
 * lengthRow: Número de columnas de A.
 * nRow: Número de renglones de A.
 * OUT
 * A_trans: A transpuesta.
 * -------------------------------------
 */
double* mTrans(double* A, int lengthRow, int nRow){
  // Declaración de variables.
  double* A_trans;
  int i, j;

  // Alocar espacio para matriz.
  A_trans = (double*)malloc(lengthRow * nRow * sizeof(double));

  // Trasponer
  for(i = 0; i < nRow; i++){
    for(j = 0; j < lengthRow; j++){
      A_trans[i*lengthRow + j] = A[j*lengthRow + i];
    }
  }

  return A_trans;
};

/* -------------------------------------
 * Suma vector-vector.
 * IN
 * v: Vector a sumar.
 * u: Vector a sumar.
 * length: Longitud de los vectores.
 * OUT
 * sum_v: Vector que contiene la suma
 * de u y v entrada a entrada.
 * -------------------------------------
 */
double* vSum(double* v, double* u, int length){
  // Declaración de variables.
  double* sum_v;
  int i;

  // Inicialización de arreglo
  sum_v = (double*)malloc(length * sizeof(double));

  // Suma
  for(i = 0; i < length; i ++){
    sum_v[i] = v[i] + u[i];
  }

  return sum_v;
}


/* -------------------------------------
 * Comparación entre dos vectores
 * IN
 * v: Vector a comparar.
 * u: Vector a comparar.
 * length: Longitud del vector.
 * OUT
 * eq: 1 si son iguales 0 en caso
 * contrario.
 * -------------------------------------
 */
int vEq(double* v, double* u, int length){
  // Declaración de variables.
  int eq, i;

  // Inicializar eq en 1.
  eq = 1;

  // Comparación
  for(i = 0; i < length; i ++){
    if(v[i] != u[i]){
      eq = 0;
      break;
    }
  }

  return eq;
}

/* -------------------------------------
 * Producto punto entre vectores
 * IN
 * v: Vector a multiplicar.
 * u: Vector a multiplicar.
 * length: Longitud de los vectores.
 * OUT
 * sum: Escalar que contiene el resultado
 * del producto punto entre ambos vectores.
 * -------------------------------------
 */
double dotProd(double* v, double* u, int length){
  // Declaración de variables.
  double sum;
  int i;

  // Suma
  sum = 0;
  for(i = 0; i < length; i ++){
    sum = sum + v[i] * u[i];
  }

  return sum;
}

/* -------------------------------------
 * Método para llevar a cabo producto de
 * matriz por vector.
 * IN
 * B: Matriz a multiplicar.
 * v: Vector a multiplicar.
 * lengthRow: Tamaño de cada fila.
 * nRow: Número de filas.
 * OUT
 * Bd: Resultado de multiplicar
 * matriz por vector.
 * -------------------------------------
 */
double* mProd(double* B, double* v, int lengthRow, int nRow){
  // Declaración de variables.
  int i, j;
  double* Bd;
  double sum;

  // Inicialización de variables.
  Bd = (double*)malloc(nRow * sizeof(double));

  for(i = 0; i < nRow; i++){
    sum = 0;
    for(j = 0; j < lengthRow; j++){
      sum = sum + B[i*lengthRow + j] * v[j];
    }
    Bd[i] = sum;
  }
  return Bd;
}


/* -------------------------------------
 * Función que calcula la norma de un
 * vector.
 * IN
 * x: vector al que se quiere calcular
 * la norma
 * norm: longitud del vecotr.
 * OUT
 * Norma del vector
 * -------------------------------------
 */
double norm(double * x, int length){
  return sqrt(dotProd(x, x, length));
}


/* -------------------------------------
 * Función que calcula el mínimo entre dos
 * números.
 * IN
 * x: double para comparar.
 * y: double para comparar.
 * OUT
 * El menor entre los dos números.
 * -------------------------------------
 */
double min (double x, double y){
  return x > y ? y : x;
}
