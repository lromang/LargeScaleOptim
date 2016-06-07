/*
 * Luis Manuel Román García
 *
 * ----------------------------------
 * Rutinas de propósito general para
 * optimización numérica. CUDA
 * ----------------------------------
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<cuda.h>

__global__ double* vSum(double* v, double* u, double* w, int length);

int main(){
  // Declaración de variables.
  double *d_v, *d_u, *d_w, *h_v, *h_u, *h_w;
  int i, size;

  // Lectura de longitud.
  printf("Escribir tamaño de vectores:\n");
  scanf("%d", &length);

  // Alocar espacio en host.
  h_v = (double*)malloc(length*sizeof(double));
  h_u = (double*)malloc(length*sizeof(double));

  // Alocar memoria en device.
  cudaError_t err0 = cudaMalloc((void**) &d_v, size);
  cudaError_t err1 = cudaMalloc((void**) &d_u, size);
  cudaError_t err2 = cudaMalloc((void**) &d_w, size);

  // Enviar argumentos a device.
  cudaMemcpy(d_v, h_v, size, cudaMemcpyHostToDevice);
  cudaMemcpy(d_u, h_u, size, cudaMemcpyHostToDevice);

  // Alocar espacio en device.
  cudaMalloc((void **) &v, length*sizeof(double));
  cudaMalloc((void **) &u, length*sizeof(double));
  cudaMalloc((void **) &w, length*sizeof(double));

  // Bloques y threads a levantar.
  dim3 DimGrid((n - 1)/256 + 1, 1, 1);
  dim3 DimBlock(256, 1, 1);

  //Ejecutar kernel.
  vSum<<<DimGrid, DimBlock>>>(d_u, d_v, d_w, n);

  // Regresar resultados.
  cudaMemcpy(h_w, d_w, size, cudaMemcpyDeviceToHost);

  // Liberar Memoria.
  cudaFree(d_v);
  cudaFree(d_u);
  cudaFree(d_w);

  // Verificar si hubo errores.
  if((err0 != cudaSuccess) || (err1 != cudaSuccess) || (err2 != cudaSuccess)){
    printf("%s en %s en línea %d \n",
           cudaGetErrorString(err),
           __FILE__,
           __LINE__);

    exit(EXIT_FAILURE);
  }


}

__global__ void vSum(double *v, double *u, double *w, int length){
  // Declaración de variables.
  int i, size;

  // Inicializar variables.
  size = length*sizeof(double);
  i    = blockIdx.x * blockDim.x + threadIdx.x;

  // Suma. Verificar que el índice es válido.
  if(i < length){
    w[i] = v[i] + u[i];
  }
}
