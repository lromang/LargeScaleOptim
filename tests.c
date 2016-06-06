/*
 * Luis Manuel Román García
 * 000117077
 * ----------------------------------
 * Rutinas de propósito general para
 * optimización numérica.
 * ----------------------------------
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

double* centralDiff(double (*func)(double*, int), double* x, int length);
double  testFunc(double* x, int length);
double* vProd(double* v, double alpha, int length);
double* vSum(double* v, double* u, int length);
double  dotProd(double* v, double* u, int length);
double  math(double x, double y);
double* mProd(double* B, double* v, int lengthRow, int nRow);
void    imprimeMatriz(double* A, int lengthRow, int nRow);
double* GC(double* A, double* b, int nrow);
double* NGC(double* A, double* b, int nrow, int N_max);
double* lSNCG(double* grad, double* B, int nrow);
double* mTrans(double* A, int lengthRow, int nRow);

int main(){
  // Declaración de variables.
  double *point,  *matrix, *sol, *matrix_symm, *dir_des;
  int    nrow, j, i, k;

  // Inicialización de variables.
  printf("--------------------------------\n");
  printf("Prueba Gradiente Conjugado\n");
  printf("--------------------------------\n\n");

  // Lectura de tamaño de matrix
  printf("Escriba el tamaño de la matriz: \n");
  scanf("%d", &nrow);
  point        = (double*) malloc(nrow * sizeof(double));
  matrix       = (double*) malloc((nrow * nrow) * sizeof(double));
  matrix_symm  = (double*) malloc((nrow * nrow) * sizeof(double));
  sol          = (double*) malloc(nrow * sizeof(double));
  dir_des      = (double*) malloc(nrow * sizeof(double));

  // Llenado de matriz y vector
  k = 0;
  for(j = 0; j < nrow; j++){
    for(i = 0; i < nrow; i++){
      matrix[k] = rand() % 10;
      k = k + 1;
    }
    point[j] = rand() % 10;
  }

  // Imprime Matriz.
  printf("\n------------\n");
  printf("Matriz del sistema");
  printf("\n------------\n");
  matrix_symm = vSum(mTrans(matrix, nrow, nrow),
                     matrix, (nrow*nrow));
  imprimeMatriz(matrix_symm, nrow, nrow);

  // Imprime Vector.
  printf("\n------------\n");
  printf("Vector del sistema");
  printf("\n------------\n");
  imprimeMatriz(point, 1, nrow);

  // Resolver con gradiente conjugado.
  sol = GC(matrix_symm, point, nrow);

  // Imprime Solución.
  printf("\n------------\n");
  printf("Solución del sistema");
  printf("\n------------\n");
  imprimeMatriz(sol, 1, nrow);

  // Prueba NGC
  printf("--------------------------------\n");
  printf("Prueba Newton Gradiente Conjugado\n");
  printf("--------------------------------\n\n");
  dir_des = NGC(matrix, point, nrow, 1e3);

  // Imprime Solución.
  printf("\n------------\n");
  printf("Dirección de descenso truncada");
  printf("\n------------\n");
  imprimeMatriz(dir_des, 1, nrow);

  return 0;
}

void imprimeMatriz(double* A, int lengthRow, int nRow){
  // Declaración de variables.
  int i, j, k;

  // Imprimir matriz.
  k = 0;
  for(i = 0; i < nRow; i++){
    for(j = 0; j < lengthRow; j++){
      printf(" %.2lf ", A[k]);
      k = k + 1;
    }
    printf("\n");
  }
  printf("\n");
}

/* -------------------------------------
 * Método de diferencias centrales para
 * aproximar el gradiente de una función
 * IN
 * fun: Apuntador a una función que recibe
 * un apuntador a un double y un entero.
 * x: Apuntador a un vector sobre el cual
 * se quiere aproximar el gradiente
 * length: Longitud del vector.
 * OUT
 * res: Gradiente de fun evaluado en x.
 * -------------------------------------
 */
double* centralDiff(double (*func)(double*, int), double* x, int length){
  // Declaración de variables.
  int i;
  double epsilon, u;
  double *res, *aux;

  // Inicialización de variables.
  u = 1.1e-16;
  epsilon = sqrt(u);
  res = (double*)malloc(length * sizeof(double));
  aux = (double*)malloc(length * sizeof(double));

  // Llenado de vector auxiliar.
  for(i = 0; i < length; i++){
    aux[i] = x[i];
  }

  // Evaluación de derivada.
  for(i = 0; i < length; i++){
    aux[i] = aux[i] + epsilon;
    res[i] = (func(aux, length) - func(x, length)) / epsilon;
    aux[i] = x[i];
  }

  // Liberar memoria
  free(aux);
  return res;
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
  prod_v = (double*)malloc(length * sizeof(double));

  // Multiplicación
  for(i = 0; i < length; i ++){
    prod_v[i] = v[i] * alpha;
  }

  return prod_v;
}


/* -------------------------------------
 * Transponer matriz
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
}


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

/* -------------------------------------
 * Método para buscar tamaño de paso
 * en la dirección de descenso.
 * IN
 * func: Apuntador a una función que recibe
 * cómo parámetros un apuntador a un double
 * y un entero.
 * x: Vector donde se lleva a cabo la
 * evaluación.
 * p: Dirección de descenso.
 * length: Longitud de los vectores.
 * OUT
 * alpha: Tamaño del paso en la dirección
 * de descenso.
 * -------------------------------------
 */
double backTrack(double (*func)(double*, int),
                 double* x,
                 double* p,
                 int length){

  // Decalaración de variables
  double alpha, rho, c, fx, fx_p;
  double *gradx, *aux_vprod, *aux_vsum;

  // Inicialización de variables
  alpha = rand() % 10;
  rho   = rand() % 1;
  c     = rand() % 1;
  aux_vprod = (double*)malloc(length * sizeof(double));
  aux_vsum  = (double*)malloc(length * sizeof(double));
  gradx     = (double*)malloc(length * sizeof(double));

  // Condiciones del Loop
  gradx = centralDiff(func, x, length);
  do{
    aux_vprod = vProd(p, alpha, length);
    aux_vsum  = vSum(x, aux_vprod, length);
    fx_p      = func(aux_vsum, length);
    fx        = func(x, length) + c * alpha * dotProd(gradx, p, length);
    alpha     = alpha * rho;
  } while(fx_p <= fx);

  // Liberar espacio
  free(aux_vprod);
  free(aux_vsum);
  free(gradx);

  return alpha;
}

/* -------------------------------------
 * Gradiente Conjugado
 * Resuelve equaciones de la forma
 * Ax = b
 * IN
 * A: Matríz simétrica positiva definida
 * b: Vector
 * length: longitud de b.
 * OUT
 * x
 * -------------------------------------
 */
double* GC(double* A, double* b, int nrow){
  double *r, *p, *x, *x_new, *r_new, *p_new;
  double alpha, beta, Tol;
  int k, i, N_max;

  // Número máximo de ciclos.
  N_max = 1e3;
  Tol   = 1e-3;

  // Alocar espacio para vectores.
  x     = (double*)malloc(sizeof(double)*nrow);
  x_new = (double*)malloc(sizeof(double)*nrow);
  r_new = (double*)malloc(sizeof(double)*nrow);
  p_new = (double*)malloc(sizeof(double)*nrow);

  // Punto inicial.
  for(i = 0; i < nrow; i++){
    x[i]    = rand() % 10;
  }

  // Inicializar r, p y k.
  // r = Ax0-b
  r = vSum(mProd(A, x, nrow, nrow), vProd(b, -1, nrow), nrow);
  // p = -r
  p = vProd(r, -1, nrow);
  k = 0;

  // Comenzar loop.
  // Mientras ||r|| > Tol && k < N_max
  while((sqrt(dotProd(r, r, nrow)) > Tol) && (k < N_max)){
    // alpha = r'r/p'Ap
    alpha = dotProd(r, r, nrow) / dotProd(p, mProd(A, p, nrow, nrow), nrow);
    // x_new = x + alpha*p
    x_new = vSum(x, vProd(p, alpha, nrow), nrow);
    x     = x_new;
    // r_new = r + alpha*A*p
    r_new = vSum(r, vProd(mProd(A, p, nrow, nrow), alpha, nrow), nrow);
    // beta  = r_new'r_new/r'r
    beta  = dotProd(r_new, r_new, nrow)/dotProd(r, r, nrow);
    r = r_new;
    // p_new = -r + beta*p
    p_new = vSum(vProd(r, -1, nrow), vProd(p, beta, nrow), nrow);
    p     = p_new;
    k++;
  }

  // Liberar Espacio.
  free(r);

  return x;
}

/* -------------------------------------
 * Newton Gradiente Conjugado -
 * IN
 * A: Hessian
 * b: -Gradient
 * OUT
 * p: Dirección de descenso que hace
 * que ||grad(f(x))|| < epsilon||f(x)||
 * -------------------------------------
 */
double* NGC(double* A, double* b, int nrow, int N_max){
  double *r, *p, *x, *z, *z_new, *r_new, *d, *d_new;
  double alpha, beta, Tol, epsilon, norm;
  int k, i;

  // Número máximo de ciclos.
  Tol   = 1e-3;

  // Alocar espacio para vectores.
  x     = (double*)malloc(sizeof(double)*nrow);
  z_new = (double*)malloc(sizeof(double)*nrow);
  r_new = (double*)malloc(sizeof(double)*nrow);
  d     = (double*)malloc(sizeof(double)*nrow);
  d_new = (double*)malloc(sizeof(double)*nrow);
  z     = (double*)malloc(sizeof(double)*nrow);

  // Inicializar x y z.
  for(i = 0; i < nrow; i++){
    x[i]    = rand() % 10;
    z[i]    = 0;
  }

  // Inicializar r y epsilon.
  // r = grad(fx) = Ax-b
  r       = vSum(mProd(A, x, nrow, nrow), vProd(b, -1, nrow), nrow);
  norm    = sqrt(dotProd(r, r, nrow));
  epsilon = min(.5, sqrt(norm)) * norm;

  // d = -r
  d = vProd(r, -1, nrow);
  k = 0;

  // Comenzar loop.
  // Mientras ||r|| > Tol && k < N_max
  while((sqrt(dotProd(r, r, nrow)) > Tol) && (k < N_max)){
    // Verificar si estamos en una dirección de descenso.
    if(dotProd(d, mProd(A, d, nrow, nrow), nrow) <= 0){
      if(k == 0){
        p = d;
        break; // Break sale de while.
      }else{
        p = z;
        break; // Break sale de while.
      }
    }
    // alpha = r'r/d'Ad
    alpha = dotProd(r, r, nrow) / dotProd(d, mProd(A, d, nrow, nrow), nrow);
    // x_new = x + alpha*d
    z_new = vSum(z, vProd(d, alpha, nrow), nrow);
    z     = z_new;
    // r_new = r + alpha*A*d
    r_new = vSum(r, vProd(mProd(A, d, nrow, nrow), alpha, nrow), nrow);
    // Truncamiento de GC.
    if(sqrt(dotProd(r_new, r_new, norm)) < epsilon){
      p = z;
      break; // Break sale de while.
    }
    // beta  = r_new'r_new/r'r
    beta  = dotProd(r_new, r_new, nrow) / dotProd(r, r, nrow);
    r = r_new;
    // p_new = -r + beta*p
    d_new = vSum(vProd(r, -1, nrow), vProd(d, beta, nrow), nrow);
    d     = d_new;
    k++;
  }

  // Liberar Espacio.
  free(r);

  return p;
}


/* -------------------------------------
 * Newton Gradiente Conjugado -
 * Búsqueda de Línea
 * IN
 *
 * OUT
 *
 * -------------------------------------
 */
double* lSNCG(double* grad, double* B, int nrow){
  // Declaración de variables.


}


/*
 * Función de prueba.
 */
double testFunc(double* x, int length){
  int i;
  double sum;
  sum = 0;
  for(i = 0; i < length; i++){
    sum = sum  + (100 - i + 1) * (x[i] * x[i]);
  }
  return sum;
}
