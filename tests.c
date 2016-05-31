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

int main(){
  // Declaración de variables.
  double *point, *derivative, *matriz;
  double evaluate;
  int    length, i, nrow;

  // Inicialización de variables.
  printf("Escriba el tamaño del vector: \n");
  scanf("%d", &length);
  point = (double*) malloc(length * sizeof(double));

  // Impresión de punto.
  printf("El punto sobre el que se hará la prueba es: \n");
  for(i = 0; i < length; i ++){
    point[i] = rand() % 10;
    printf("%.0lf\n", point[i]);
  }

  // Evaluación de función en el punto.
  evaluate = testFunc(point, length);

  // Impresión de valor de función evaluada en punto.
  printf("El valor de la función evaluada en ese punto es: %fl\n", evaluate);

  // Estimar gradiente de función.
  derivative = centralDiff(&testFunc, point, length);

  // Impresión de gradiente.
  printf("El valor de la derivada es: \n");
  for(i = 0; i < length; i ++){
    printf("%lf\n", derivative[i]);
  }

  // Prueba Matriz - Vector
  printf("Número de renglones: \n");
  scanf("%d", &nrow);
  printf("La matriz a multiplicar es: \n");
  matriz = (double*)malloc(length*nrow*sizeof(double));
  for(i = 0; i < nrow*length; i++){
    matriz[i] = rand() % 10;
  }

  // Imprimir matriz antes del producto
  imprimeMatriz(matriz, length, nrow);

  // Imprimir después del producto
  matriz = mProd(matriz, point, length, nrow);
  imprimeMatriz(matriz, length, 1);

  return 0;
}

void imprimeMatriz(double* A, int lengthRow, int nRow){
  // Declaración de variables.
  int i, j;

  // Imprimir matriz.
  for(i = 0; i < nRow; i++){
    for(j = 0; j < lengthRow; j++){
      printf(" %.2lf ", A[(i + 1) * j]);
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
      sum = sum + B[(i + 1) * j] * v[j];
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

  return alpha;
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
double* lSNCG(double* x, double* grad, double* B,
              int length, int TOL){
  // Declaración de variables.
  int i, j;
  double epsilon, normSquare, norm,
    dbd, alpha, normR, beta;
  double *r, *d, *z, *Bd, *oldR, *p;

  // Inicialización de variables.
  normSquare = dotProd(grad, grad, length);
  norm       = sqrt(normSquare);

  // Comienza loop externo
  for(i = 0; i < TOL; i++){
    epsilon = min(0.5, norm) * norm;
    z       = 0;
    r       = grad;
    d       = vProd(r, -1, length);
    // Comienza loop interno
    for(j = 0; j < TOL; j++){
      Bd      = mProd(B, d, length, length);
      dbd     = dotProd(Bd, d, length);
      if(dbd <= 0){
        if(j == 0){
          p = vProd(r, -1, length);
          break;
        }else{
          p = z;
          break;
        }
      }
      alpha = dotProd(r,r,length) / dbd;
      z     = vSum(z, vProd(d, alpha, length), length);
      oldR  = r;
      r     = vSum(r, vProd(Bd, alpha, length), length);
      normR = dotProd(r, r, length);
      if(normR < epsilon*epsilon){
        p = z;
        break;
      }
      beta = dotProd(r, r, length) / dotProd(oldR, oldR,  length);
      d    = vSum(vProd(r, -1, length), vProd(d, beta, length), length);
    }
    x = vSum(x, p, length);
  }
  return x;
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
