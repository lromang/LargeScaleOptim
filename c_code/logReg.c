
/*
 * logReg.c
 *
 *  Created on: April 28, 2010
 *      Author: YuchenWu
 */

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <cblas.h>
#include "logReg.h"

#define DEBUG_LOGREG


#ifdef USE_CUDA

#include "cublas.h"

static void logRegMallocCudaFloat(LogRegContext * lr, int size, LOGREG_REAL **var) {
    cublasStatus status;

    status = cublasAlloc(size, sizeof(LOGREG_REAL), (void **) var);
    if (status != CUBLAS_STATUS_SUCCESS) {
        fprintf(stderr,
            "!!! Device memory allocation error (feat by sample) \n");
        exit(-1);
    }
}

static void logRegSetVecCudaFloat(LogRegContext * lr, const int len,
    const LOGREG_REAL * const src, const int incs, LOGREG_REAL * const dest, const int incd) {

    cublasStatus status;

    status = cublasSetVector(len, sizeof(LOGREG_REAL), src, incs, dest, incd);
    if (status != CUBLAS_STATUS_SUCCESS) {
        fprintf(stderr, "!!! Device access error (write).\n");
        exit(-1);
    }
}

static void logRegGetVecCudaFloat(LogRegContext * lr, const int len,
    const LOGREG_REAL * const src, const int incs, LOGREG_REAL * const dest, const int incd) {
    cublasStatus status;

    status = cublasGetVector(len, sizeof(LOGREG_REAL), src, incs, dest, incd);
    if (status != CUBLAS_STATUS_SUCCESS) {
        fprintf(stderr, "!!! Device access error (write).\n");
        exit(-1);
    }
}

static void logRegSgemmCudaFloat(LogRegContext * lr, char transa, char transb,
    int m, int n, int k, LOGREG_REAL alpha, const LOGREG_REAL *A, int lda, const LOGREG_REAL *B,
    int ldb, LOGREG_REAL beta, LOGREG_REAL *C, int ldc) {
    cublasStatus status;

    cublasSgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    status = cublasGetError();
    if (status != CUBLAS_STATUS_SUCCESS) {
        fprintf(stderr, "!!!! Sgemm kernel execution error.\n");
        exit(-1);
    }
}

#else

#include <cblas.h>

#ifdef USE_DOUBLE

#define cblas_gemm cblas_dgemm
#define cblas_dot cblas_ddot
#define cblas_axpy cblas_daxpy

#else

#define cblas_gemm cblas_sgemm
#define cblas_dot cblas_sdot
#define cblas_axpy cblas_saxpy

#endif /* USE_DOUBLE */

#endif

static void logRegMallocReal(LogRegContext * lr, int size, LOGREG_REAL ** var) {
    *var = (LOGREG_REAL *) malloc(size * sizeof(LOGREG_REAL));
    if (*var == NULL) {
        fprintf(stderr, "!!! Host memory allocation error REAL. \n");
        exit(-1);
    }
}

static void logRegMallocInt(LogRegContext * lr, int size, int ** var) {
    *var = (int *) malloc(size * sizeof(int));
    if (*var == NULL) {
        fprintf(stderr, "!!! Host memory allocation error INT. \n");
        exit(-1);
    }
}

/**
 * Sample training points from daDataFeatures, daDataFeatures2
 * and naDataClassIdx into daFeatures, daFeatures2 and naClassIdx.
 */
static void logRegPrepareSample(LogRegContext * const lr) {
    int i;

    for (i = 0; i < lr->nSizeSample; i++) {
        if (lr->nCurrentSample >= lr->nNumTotal) {
        	//reset the sample index back to zero
            lr ->nCurrentSample = 0;
        }
        /* Assign the class label of */
        lr->naClassIdx[i] = lr->naDataClassIdx[lr->nCurrentSample];
        /* Find the starting point of the necessary feature vector and copy it: memcpy(dest,source,length)
         * 		Taking it from daDataFeatures: (pointer to the full data set), and placing them in daFeatures: (features used in the sample)*/
        memcpy(lr->daFeatures + i * lr->nNumFeat, lr->daDataFeatures + lr->nCurrentSample * lr->nNumFeat, lr->nNumFeat * sizeof(LOGREG_REAL));
        lr->nCurrentSample++;
    }
}

/**
 * Sample training points from daFeatures, daSamplePIJ, daExpIJ and daSumExpI
 * into daSubSampleFeatures, daSubSamplePIJ and daSubSampleExpIJ and daSubSampleExpJ.
 *
 * This is the sample for the subsampled values used in the Hessian .
 *
 */


static void logRegPrepareSubSample(LogRegContext * const lr) {
    int i, k;

    for (i = 0, k = 0; i < lr->nSizeSubSample; k++) {
        if (k >= lr->nSizeSample) {
            k = 0;
        }
        if (rand() < nearbyint(RAND_MAX * lr->dSubSampleRate)) {
            memcpy(lr->daSubSampleFeatures + i * lr->nNumFeat, lr->daFeatures + k * lr->nNumFeat, lr->nNumFeat * sizeof(LOGREG_REAL));
            memcpy(lr->daSubSamplePIJ + i * lr->nNumClass, lr->daPIJ + k * lr->nNumClass, lr->nNumClass * sizeof(LOGREG_REAL));
            memcpy(lr->daSubSampleExpIJ + i * lr->nNumClass, lr->daExpIJ + k * lr->nNumClass, lr->nNumClass * sizeof(LOGREG_REAL));
            lr->daSubSampleSumExpJ[i] = lr->daSumExpJ[k];
            i++;
        }
    }
}

/**
 *  Prepare UIJ, PIJ, QIJ, EXPIJ At The Point X
 */
static void logRegPrepareFunctionEval(LogRegContext * const lr,
    const int nNumFeat, const int nNumClass, const int nNumSample,
    const int * const naClassIdx, const LOGREG_REAL * const daFeatures,
    const LOGREG_REAL * const x, LOGREG_REAL * const daPIJ,
    LOGREG_REAL * const daQIJ, LOGREG_REAL * const daExpIJ, LOGREG_REAL * const daSumExpJ) {

    int i, j, k;

#if defined(USE_CUDA)
    // Using CUBLAS to compute the UIJ matrix.
    logRegSetVecCudaFloat(lr, nNumSample * nNumFeat, daFeatures,
        1, lr->cudaFeatBySample, 1);
    logRegSetVecCudaFloat(lr, lr->n, x, 1, lr->cudaFeatByClass, 1);
    logRegSgemmCudaFloat(lr, 'n', 'n', nNumClass, nNumSample,
        nNumFeat, 1.0, lr->cudaFeatByClass, nNumClass,
        lr->cudaFeatBySample, nNumFeat, 0.0, lr->cudaClassBySample,
        nNumClass);
    logRegGetVecCudaFloat(lr, nNumClass * nNumSample,
        lr->cudaClassBySample, 1, daUIJ, 1);
#else
    /*This computes w (or x)*f */
    cblas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nNumClass,
        nNumSample, nNumFeat, 1.0, x, nNumClass, daFeatures, nNumFeat, 0.0,
        daPIJ, nNumClass);
#endif

    // Transform the UIJ (stored in PIJ) matrix to ExpIJ and prepare sumExpIJ.
    memset(daSumExpJ, 0, nNumSample * sizeof(LOGREG_REAL));
    k = 0;
    for (j = 0; j < nNumSample; j++) {
        for (i = 0; i < nNumClass; i++) {
            daExpIJ[k] = exp(daPIJ[k]);
            daSumExpJ[j] += daExpIJ[k];
            k++;
        }
    }

    // Compute the PIJ matrix.
    // Compute the QIJ matrix to be
    //     QIJ = PIJ - delta(i, j)
    // where delta(i, j) is the indicator function
    // for sample j belonging to class i.
    k = 0;
    for (j = 0; j < nNumSample; j++) {
        for (i = 0; i < nNumClass; i++) {
            daPIJ[k] = daExpIJ[k] / daSumExpJ[j];
            if (i == naClassIdx[j]) {
                daQIJ[k] = daPIJ[k] - 1;
            } else {
                daQIJ[k] = daPIJ[k];
            }
            k++;
        }
    }
}

/* Given the PIJ matrix, compute objective function. */
static void logRegFormFunctionFromP(LogRegContext * const lr, const int nNumClass,
    const int nNumSample, const int * const naClassIdx,
    const LOGREG_REAL * const daPIJ, LOGREG_REAL * const f) {
    int i;

    *f = 0;
    for (i = 0; i < nNumSample; i++) {
        *f -= log2(daPIJ[i * lr->nNumClass + naClassIdx[i]]);
    }
    *f *= log(2);
    *f /= nNumSample;
}

#if 0

static void logRegFormFunction(LogRegContext * const lr, const int nNumClass,
    const int nNumSample, const int * const naClassIdx,
    const LOGREG_REAL * const daUIJ, const LOGREG_REAL * const daSumExpJ, LOGREG_REAL * const f) {
    int i;

    *f = 0;
    for (i = 0; i < nNumSample; i++) {
        *f += log2(daSumExpJ[i]);
    }
    *f = *f * log(2);
    for (i = 0; i < nNumSample; i++) {
        *f -= daUIJ[i * nNumClass + naClassIdx[i]];
    }
    *f /= (LOGREG_REAL) nNumSample;
}

#endif

static void logRegFormGradient(LogRegContext * const lr, const int nNumFeat,
    const int nNumClass, const int nNumSample, const LOGREG_REAL * const daFeatures,
    const LOGREG_REAL * const daQIJ, LOGREG_REAL * const g) {

	int i =0;
    if (nNumSample <= 0) {
        fprintf(stderr,
            "!!! nNumSample has to be positive in logRegFormGradient.\n");
        exit(-1);
    }

#if defined(USE_CUDA)
    // Compute the sample average gradient.
    // This essentially calls CUDA to compute
    //     QIJ * Features
    logRegSetVecCudaFloat(lr, nNumSample * nNumFeat, daFeatures,
        1, lr->cudaFeatBySample, 1);
    logRegSetVecCudaFloat(lr, nNumSample * nNumClass, daQIJ, 1,
        lr->cudaClassBySample, 1);
    cublasSgemm('n', 't', nNumClass, nNumFeat, nNumSample, 1.0,
        lr->cudaClassBySample, nNumClass, lr->cudaFeatBySample,
        nNumFeat, 0.0, lr->cudaFeatByClass, nNumClass);
    logRegGetVecCudaFloat(lr, lr->n, lr->cudaFeatByClass, 1, g, 1);
#else

    /* use CBLAS to compute: QIJ * Features
     * 		storing it in g*/
    cblas_gemm(CblasColMajor, CblasNoTrans, CblasTrans, nNumClass, nNumFeat,
        nNumSample, 1.0, daQIJ, nNumClass, daFeatures, nNumFeat, 0.0, g,
        nNumClass);

#endif

    /* Add up all the components and normalize by the number of sampled points*/
    for (i = 0; i < lr->n; i++) {
        g[i] /= (LOGREG_REAL) nNumSample;
    }

}



static void logRegFormGradientFromP(LogRegContext * const lr, const int nNumFeat,
    const int nNumClass, const int nNumSample, const LOGREG_REAL * const daFeatures,
    const int * const naClassIdx,
    const LOGREG_REAL * const daPIJ, LOGREG_REAL * const daQIJ, LOGREG_REAL * const g) {

    int i, j, k;

    k = 0;
    for (j = 0; j < nNumSample; j++) {
        for (i = 0; i < nNumClass; i++) {
            if (i == naClassIdx[j]) {
                daQIJ[k] = daPIJ[k] - 1;
            } else {
                daQIJ[k] = daPIJ[k];
            }
            k++;
        }
    }

    logRegFormGradient(lr, nNumFeat, nNumClass, nNumSample, daFeatures, daQIJ, g);
}

static void logRegFullFunctionGradientInternal(const LOGREG_REAL * const x, LOGREG_REAL * const g,
    LOGREG_REAL * const f, LOGREG_REAL * const daFeatures, int * const naClassIdx, int nNumTotal, LogRegContext * lr) {
    int i, j, k;
    LOGREG_REAL daSumExpJ;

#if defined(USE_CUDA)
#error "Function logRegFullFunctionGradientInternal is not implemented with USE_CUDA option."
#else
    cblas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, lr->nNumClass,
        nNumTotal, lr->nNumFeat, 1.0, x, lr->nNumClass, daFeatures, lr->nNumFeat, 0.0,
        lr->daPIJF, lr->nNumClass);
#endif

    for (j = 0; j < nNumTotal; j++) {
        daSumExpJ = 0;
        k = j * lr->nNumClass;
        for (i = 0; i < lr->nNumClass; i++) {
            daSumExpJ += exp(lr->daPIJF[k]);
            k++;
        }

        k = j * lr->nNumClass;
        for (i = 0; i < lr->nNumClass; i++) {
            lr->daPIJF[k] = exp(lr->daPIJF[k]) / daSumExpJ;
            k++;
        }
    }

    /* FUNCTION COMPUTATION STARTS HERE */
    logRegFormFunctionFromP(lr, lr->nNumClass, nNumTotal, naClassIdx, lr->daPIJF, f);

    /* GRADIENT COMPUTATION STARTS HERE */
    if (g != NULL) {
        logRegFormGradientFromP(lr, lr->nNumFeat, lr->nNumClass, nNumTotal,
            daFeatures, naClassIdx, lr->daPIJF, lr->daPIJF, g);
    }
}

static void logRegLoadFromFile(const char * const fileName, const int nNumFeat, const int nNumTotal, int * const naClassIdx, LOGREG_REAL * const daFeatures) {
    int i, j, k;
    FILE *fp;

    if (NULL == (fp = fopen(fileName, "r"))) {
        printf("!!! Failed to read data file %s.\n", fileName);
        exit(-1);
    }

    for (i = 0; i < nNumTotal; i++) {
        fscanf(fp, "%d", naClassIdx + i);
        for (j = 0; j < nNumFeat; j++) {
            k = i * nNumFeat + j;
#ifdef USE_DOUBLE
            fscanf(fp, "%lf", daFeatures + k);
#else
            fscanf(fp, "%f", daFeatures + k);
#endif
        }
        while (!feof(fp) && fgetc(fp) != '\n') {
        };
    }

    fclose(fp);
}

LogRegContext * newLogReg() {
    LogRegContext * lr;
    lr = (LogRegContext *) malloc(sizeof(*lr));

    if (lr == NULL) {
        printf("!!! Cannot allocate memory for logistic regression model.\n");
        exit(-1);
    }
    return lr;
}

void logRegInit(LogRegContext * lr, int nNumFeat, int nNumClass,
    int nSizeSubSample, int nSizeSample, int nNumTotal, int nNumTotalTest, int nNumSampleUnits){
	/* Changed this section to allow for dynamic sample sizes*/
    int i;
    int samplePerUnit;

    lr->nNumClass = nNumClass;
    lr->nNumFeat = nNumFeat;
    lr->nNumTotal = nNumTotal;

    lr->nNumTotalTest = nNumTotalTest;

    lr->nSizeSample = nSizeSample;
    lr->nSizeSubSample = nSizeSubSample;
    lr->n = nNumClass * nNumFeat;

    lr->nNumSampleUnits = nNumSampleUnits;
    lr->nCurrentSample = 0;

    lr->dSubSampleRate = (LOGREG_REAL) lr->nSizeSubSample / lr->nSizeSample;
    if (lr->dSubSampleRate <= 0 || lr->dSubSampleRate > 1) {
        fprintf(stderr, "!!! Unsupported sample rate: %g.\n",
            lr->dSubSampleRate);
        exit(-1);
    }

    if (lr->nSizeSample < lr->nNumSampleUnits || lr->nSizeSample
        % lr->nNumSampleUnits) {
        fprintf(stderr, "!!! Cannot split %d samples into %d units.\n",
            lr->nSizeSample, lr->nNumSampleUnits);
        exit(-1);
    } else {
        samplePerUnit = lr->nSizeSample / lr->nNumSampleUnits;
    }

    lr->daSampleInd = 0;

#if defined(USE_CUDA)
    /* Initialize CUDA */
    cublasStatus status;
    status = cublasInit();
    if (status != CUBLAS_STATUS_SUCCESS) {
        fprintf(stderr, "!!! CUBLAS initialization error.\n");
        exit(-1);
    }
    /* Allocate CUDA memory */
    logRegMallocCudaFloat(lr, lr->nNumFeat * lr->nNumTotal,
        &lr->cudaFeatBySample);
    logRegMallocCudaFloat(lr, lr->nNumClass * lr->nNumTotal,
        &lr->cudaClassBySample);
    logRegMallocCudaFloat(lr, lr->nNumFeat * lr->nNumClass,
        &lr->cudaFeatByClass);
    cublasGetError();
#endif

    /* Allocate Host memory */
    logRegMallocReal(lr, lr->nNumClass * lr->nSizeSample, &lr->daPIJ);
    logRegMallocReal(lr, lr->nNumClass * lr->nSizeSample, &lr->daQIJ);
    logRegMallocReal(lr, lr->nNumClass * lr->nSizeSample, &lr->daExpIJ);
    logRegMallocReal(lr, lr->nSizeSample, &lr->daSumExpJ);
    logRegMallocReal(lr, lr->nNumFeat * lr->nSizeSample, &lr->daFeatures);
    logRegMallocInt(lr, lr->nSizeSample, &lr->naClassIdx);

    logRegMallocReal(lr, lr->nNumClass * lr->nSizeSubSample,
        &lr->daSubSampleVIJ);
    logRegMallocReal(lr, lr->nNumClass * lr->nSizeSubSample,
        &lr->daSubSamplePIJ);
    logRegMallocReal(lr, lr->nNumClass * lr->nSizeSubSample,
        &lr->daSubSampleExpIJ);
    logRegMallocReal(lr, lr->nNumFeat * lr->nSizeSubSample,
        &lr->daSubSampleFeatures);
    logRegMallocReal(lr, lr->nSizeSubSample, &lr->daSubSampleSumExpJ);

    logRegMallocReal(lr, lr->nSizeSubSample, &lr->daTmpSubSample);
    logRegMallocReal(lr, lr->nNumClass * lr->nNumTotal, &lr->daPIJF);
    logRegMallocReal(lr, lr->nNumFeat * lr->nNumTotal, &lr->daDataFeatures);
    logRegMallocInt(lr, lr->nNumTotal, &lr->naDataClassIdx);

    lr->pdaGIJ = (LOGREG_REAL **) malloc(lr->nNumSampleUnits * sizeof(LOGREG_REAL *));
    logRegMallocInt(lr, lr->nNumSampleUnits, &lr->naGIJstartIdx);
    logRegMallocInt(lr, lr->nNumSampleUnits, &lr->naGIJsizes);

    for (i = 0; i < lr->nNumSampleUnits; i++) {
        logRegMallocReal(lr, lr->nNumClass * lr->nNumFeat, &lr->pdaGIJ[i]);
        lr->naGIJstartIdx[i] = i * samplePerUnit;
        if (i > 0) {
            lr->naGIJsizes[i-1] = lr->naGIJstartIdx[i] - lr->naGIJstartIdx[i - 1];
        }
    }
    lr->naGIJsizes[i - 1] = lr->nSizeSample - lr->naGIJstartIdx[i - 1];

    // Allocate space for testing data.
    if (lr->nNumTotalTest > 0) {
        logRegMallocReal(lr, lr->nNumTotalTest * lr->nNumFeat, &lr->daTestFeatures);
        logRegMallocInt(lr, lr->nNumTotalTest, &lr->naTestClassIdx);
    }

    lr->status = INIT;
}

#if 0
void logRegInit(LogRegContext * lr, int nNumFeat, int nNumClass,
    int nSizeSubSample, int nSizeSample, int nNumTotal, int nNumSampleUnits) {

    int i;
    int samplePerUnit;

    lr->nNumClass = nNumClass;
    lr->nNumFeat = nNumFeat;
    lr->nNumTotal = nNumTotal;
    lr->nSizeSample = nSizeSample;
    lr->nSizeSubSample = nSizeSubSample;
    lr->n = nNumClass * nNumFeat;

    lr->nNumSampleUnits = nNumSampleUnits;
    lr->nCurrentSample = 0;

    lr->dSubSampleRate = (LOGREG_REAL) lr->nSizeSubSample / lr->nSizeSample;
    if (lr->dSubSampleRate <= 0 || lr->dSubSampleRate > 1) {
        fprintf(stderr, "!!! Unsupported sample rate: %g.\n",
            lr->dSubSampleRate);
        exit(-1);
    }

    if (lr->nSizeSample < lr->nNumSampleUnits || lr->nSizeSample
        % lr->nNumSampleUnits) {
        fprintf(stderr, "!!! Cannot split %d samples into %d units.\n",
            lr->nSizeSample, lr->nNumSampleUnits);
        exit(-1);
    } else {
        samplePerUnit = lr->nSizeSample / lr->nNumSampleUnits;
    }

#if defined(USE_CUDA)
    /* Initialize CUDA */
    cublasStatus status;
    status = cublasInit();
    if (status != CUBLAS_STATUS_SUCCESS) {
        fprintf(stderr, "!!! CUBLAS initialization error.\n");
        exit(-1);
    }
    /* Allocate CUDA memory */
    logRegMallocCudaFloat(lr, lr->nNumFeat * lr->nNumTotal,
        &lr->cudaFeatBySample);
    logRegMallocCudaFloat(lr, lr->nNumClass * lr->nNumTotal,
        &lr->cudaClassBySample);
    logRegMallocCudaFloat(lr, lr->nNumFeat * lr->nNumClass,
        &lr->cudaFeatByClass);
    cublasGetError();
#endif

    /* Allocate Host memory */

    logRegMallocReal(lr, lr->nNumClass * lr->nNumTotal, &lr->daPIJF);
    logRegMallocReal(lr, lr->nNumClass * lr->nNumTotal, &lr->daPIJ);
    logRegMallocReal(lr, lr->nNumClass * lr->nNumTotal, &lr->daQIJ);
    logRegMallocReal(lr, lr->nNumClass * lr->nNumTotal, &lr->daExpIJ);
    logRegMallocReal(lr, lr->nNumTotal, &lr->daSumExpJ);
    logRegMallocReal(lr, lr->nNumFeat * lr->nNumTotal, &lr->daFeatures);
    logRegMallocInt(lr, lr->nNumTotal, &lr->naClassIdx);

    logRegMallocReal(lr, lr->nNumClass * lr->nSizeSubSample,
        &lr->daSubSampleVIJ);
    logRegMallocReal(lr, lr->nNumClass * lr->nSizeSubSample,
        &lr->daSubSamplePIJ);
    logRegMallocReal(lr, lr->nNumClass * lr->nSizeSubSample,
        &lr->daSubSampleExpIJ);
    logRegMallocReal(lr, lr->nNumFeat * lr->nSizeSubSample,
        &lr->daSubSampleFeatures);
    logRegMallocReal(lr, lr->nSizeSubSample, &lr->daSubSampleSumExpJ);
    logRegMallocReal(lr, lr->nSizeSubSample, &lr->daTmpSubSample);

    logRegMallocReal(lr, lr->nNumFeat * lr->nNumTotal, &lr->daDataFeatures);
    logRegMallocInt(lr, lr->nNumTotal, &lr->naDataClassIdx);

    lr->pdaGIJ = (LOGREG_REAL **) malloc(lr->nNumSampleUnits * sizeof(LOGREG_REAL *));
    logRegMallocInt(lr, lr->nNumSampleUnits, &lr->naGIJstartIdx);
    logRegMallocInt(lr, lr->nNumSampleUnits, &lr->naGIJsizes);

    for (i = 0; i < lr->nNumSampleUnits; i++) {
        logRegMallocReal(lr, lr->nNumClass * lr->nNumFeat, &lr->pdaGIJ[i]);
        lr->naGIJstartIdx[i] = i * samplePerUnit;
        if (i > 0) {
            lr->naGIJsizes[i-1] = lr->naGIJstartIdx[i] - lr->naGIJstartIdx[i - 1];
        }
    }
    lr->naGIJsizes[i - 1] = lr->nSizeSample - lr->naGIJstartIdx[i - 1];

    lr->status = INIT;
}
#endif

void freeLogReg(LogRegContext * lr) {
    int i;

#if defined(USE_CUDA)
    cublasFree(lr->cudaClassBySample);
    cublasFree(lr->cudaFeatByClass);
    cublasFree(lr->cudaFeatBySample);
    cublasShutdown();
#endif

    free(lr->daPIJF);
    free(lr->daPIJ);
    free(lr->daQIJ);
    free(lr->daExpIJ);
    free(lr->daSumExpJ);
    free(lr->naClassIdx);
    free(lr->daFeatures);

    free(lr->daSubSampleFeatures);
    free(lr->daSubSamplePIJ);
    free(lr->daSubSampleExpIJ);
    free(lr->daSubSampleVIJ);
    free(lr->daSubSampleSumExpJ);
    free(lr->daTmpSubSample);

    free(lr->daDataFeatures);
    free(lr->naDataClassIdx);

    for (i = 0; i < lr->nNumSampleUnits; i++) {
        free(lr->pdaGIJ[i]);
    }
    free(lr->pdaGIJ);
    free(lr->naGIJsizes);
    free(lr->naGIJstartIdx);

    if (lr->daTestFeatures) {
        free(lr->daTestFeatures);
    }

    if (lr->naTestClassIdx) {
        free(lr->naTestClassIdx);
    }

    free(lr);
}

/**
 * Load features from file to the LogRegContext.
 * Input:
 *      nNumTotal
 *      nNumFeat
 * Output:
 *      naDataClassIdx
 *      daDataFeatures
 *      daDataFeatures2
 */
void logRegLoadData(const LogRegContext * const lr, const char * const fileName) {
    logRegLoadFromFile(fileName, lr->nNumFeat, lr->nNumTotal, lr->naDataClassIdx, lr->daDataFeatures);
}

void logRegLoadTest(LogRegContext * const lr, const char * const fileName) {
    if (lr->nNumTotalTest == 0) {
        fprintf(stderr, "Total number of testing data is 0. Abort.\n");
        return;
    }

    logRegLoadFromFile(fileName, lr->nNumFeat, lr->nNumTotalTest, lr->naTestClassIdx, lr->daTestFeatures);
}

void logRegSampleGradientStd(const LOGREG_REAL * const g, LOGREG_REAL * const std,
    LogRegContext * lr) {
    if (lr->status != FUNCTION_GRADIENT_DONE) {
        fprintf(
            stderr,
            "!!! logRegSampleGradientVar called before logRegSampleFunctionGradient. Status = %d. \n",
            lr->status);
        exit(-1);
    }

    int i, j;

    memset(std, 0, lr->n * sizeof(LOGREG_REAL));
    for (j = 0; j < lr->n; j++) {
        for (i = 0; i < lr->nNumSampleUnits; i++) {
            std[j] += pow(lr->pdaGIJ[i][j], 2);
        }
        std[j] /= lr->nNumSampleUnits;
        std[j] -= pow(g[j], 2);

        if (std[j] < -1e-8) {
            fprintf(stderr, "!!! Observed negative variance: %g\n", std[j]);
            exit(-1);
        }

        std[j] /= lr->nNumSampleUnits;
        std[j] = sqrt(std[j] > 0 ? std[j] : 0);
    }
}

/**
 * Compute d = I * v.
 *
 * TESTING
 */
void logRegIdentitySampleCovarianceVector(const LOGREG_REAL * const x,
    const LOGREG_REAL * const gbar, const LOGREG_REAL * const v, LOGREG_REAL * const d,
    LogRegContext * lr) {
    memcpy(d, v, lr->n * sizeof(LOGREG_REAL));
}

/**
 * Compute d = Cov(g) * v, where x is
 * the current point, gbar is the current gradient estimate.
 */
void logRegSampleCovarianceVector(const LOGREG_REAL * const x,
    const LOGREG_REAL * const gbar, const LOGREG_REAL * const v, LOGREG_REAL * const d,
    LogRegContext * lr) {

    int i;
    LOGREG_REAL gbarv, gv, sumgv;

    if (!(lr->status == FUNCTION_GRADIENT_DONE || lr->status
        == HESSIAN_VEC_DONE)) {
        fprintf(
            stderr,
            "!!! Hessian vector products performed before function and gradient evaluations.\n");
        exit(-1);
    }

    memset(d, 0, lr->n * sizeof(LOGREG_REAL));
    gbarv = cblas_dot(lr->n, gbar, 1, v, 1);
    sumgv = 0;

#if 1
    for (i = 0; i < lr->nNumSampleUnits; i++) {
        gv = cblas_dot(lr->n, lr->pdaGIJ[i], 1, v, 1) - gbarv;
        sumgv += gv;
        cblas_axpy(lr->n, gv, lr->pdaGIJ[i], 1, d, 1);
    }
    cblas_axpy(lr->n, -sumgv, gbar, 1, d, 1);

    for (i = 0; i < lr->n; i++) {
        d[i] /= lr->nNumSampleUnits * lr->nNumSampleUnits;
    }
#else
    for (i = 0; i < lr->nSizeSample; i++) {
        logRegFormGradient(lr, lr->nNumFeat, lr->nNumClass, 1, lr->daFeatures
            + i * lr->nNumFeat, lr->daQIJ + i * lr->nNumClass, lr->daTmpN1);
        gv = cblas_dot(lr->n, lr->daTmpN1, 1, v, 1) - gbarv;
        sumgv += gv;
        cblas_axpy(lr->n, gv, lr->daTmpN1, 1, d, 1);
    }
    cblas_axpy(lr->n, -sumgv, gbar, 1, d, 1);

    cblas_sscal(lr->n, pow(lr->nSizeSample, -2), d, 1);
#endif

}

void logRegSampleHessianVector(const LOGREG_REAL * const x, const LOGREG_REAL * const v,
    LOGREG_REAL * const daHv, LogRegContext * lr) {
    int i, j, k;

    /*
     * IF NECESSARY, RANDOMLY SAMPLE THE TRAINNING DATA, PIJ, EXPIJ.
     */
    switch (lr->status) {
    case FUNCTION_GRADIENT_DONE:
        // Evaluate HV for the first time.
        // Sample Feature matrix and PIJ
        logRegPrepareSubSample(lr);
        break;
    case HESSIAN_VEC_DONE:
        break;
    default:
        fprintf(
            stderr,
            "!!! Hessian vector products performed before function and gradient evaluations.\n");
        double tmp;
        logRegSampleFunctionGradient(x, NULL, &tmp, 0, lr);
        logRegPrepareSubSample(lr);
        break;
    }

    // compute VIJ = v * Feat

#if defined(USE_CUDA)
    // Using CUBLAS to compute the first stage of VIJ matrix.
    logRegSetVecCudaFloat(lr, lr->nSizeSubSample * lr->nNumFeat,
        lr->daSubSampleFeatures, 1, lr->cudaFeatBySample, 1);
    logRegSetVecCudaFloat(lr, lr->n, v, 1, lr->cudaFeatByClass, 1);
    logRegSgemmCudaFloat(lr, 'n', 'n', lr->nNumClass, lr->nSizeSubSample,
        lr->nNumFeat, 1.0, lr->cudaFeatByClass, lr->nNumClass,
        lr->cudaFeatBySample, lr->nNumFeat, 0.0, lr->cudaClassBySample,
        lr->nNumClass);
    logRegGetVecCudaFloat(lr, lr->nNumClass * lr->nSizeSubSample,
        lr->cudaClassBySample, 1, lr->daSubSampleVIJ, 1);
#else
    cblas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, lr->nNumClass,
        lr->nSizeSubSample, lr->nNumFeat, 1.0, v, lr->nNumClass,
        lr->daSubSampleFeatures, lr->nNumFeat, 0.0, lr->daSubSampleVIJ,
        lr->nNumClass);
#endif

    // Finish preparing VIJ matrix
    memset(lr->daTmpSubSample, 0, lr->nSizeSubSample * sizeof(LOGREG_REAL));

    for (j = 0; j < lr->nSizeSubSample; j++) {
        for (i = 0; i < lr->nNumClass; i++) {
            k = j * lr->nNumClass + i;
            lr->daTmpSubSample[j] += lr->daSubSampleExpIJ[k]
                * lr->daSubSampleVIJ[k];
        }
        lr->daTmpSubSample[j] /= lr->daSubSampleSumExpJ[j];
        for (i = 0; i < lr->nNumClass; i++) {
            k = j * lr->nNumClass + i;
            lr->daSubSampleVIJ[k] -= lr->daTmpSubSample[j];
            lr->daSubSampleVIJ[k] *= lr->daSubSamplePIJ[k];
        }
    }

#if defined(USE_CUDA)
    // Using CUBLAS to compute FEATURE * VIJ
    logRegSetVecCudaFloat(lr, lr->nSizeSubSample * lr->nNumFeat,
        lr->daSubSampleFeatures, 1, lr->cudaFeatBySample, 1);
    logRegSetVecCudaFloat(lr, lr->nSizeSubSample * lr->nNumClass,
        lr->daSubSampleVIJ, 1, lr->cudaClassBySample, 1);
    cublasSgemm('n', 't', lr->nNumClass, lr->nNumFeat, lr->nSizeSubSample, 1.0,
        lr->cudaClassBySample, lr->nNumClass, lr->cudaFeatBySample,
        lr->nNumFeat, 0.0, lr->cudaFeatByClass, lr->nNumClass);
    logRegGetVecCudaFloat(lr, lr->n, lr->cudaFeatByClass, 1, daHv, 1);
#else
    cblas_gemm(CblasColMajor, CblasNoTrans, CblasTrans, lr->nNumClass,
        lr->nNumFeat, lr->nSizeSubSample, 1.0, lr->daSubSampleVIJ,
        lr->nNumClass, lr->daSubSampleFeatures, lr->nNumFeat, 0.0, daHv,
        lr->nNumClass);
#endif

    for (i = 0; i < lr->n; i++) {
        daHv[i] /= lr->nSizeSubSample;
    }

    lr->status = HESSIAN_VEC_DONE;
}

void logRegSampleFunctionGradient(const LOGREG_REAL * const x, LOGREG_REAL * const g,
    LOGREG_REAL * const f, int bResample, LogRegContext * lr) {
    int i;

    /* SAMPLE THE TRAINING DATA. */
    if (lr->nSizeSample == lr->nNumTotal) {
        // In case of full sampling, resample only once.
        if (lr->status == INIT) {
            logRegPrepareSample(lr);
        }
    } else if (bResample) {
        logRegPrepareSample(lr);
    }

    /* PREPARE UIJ, PIJ, QIJ, EXPIJ AT THE POINT X */
    logRegPrepareFunctionEval(lr, lr->nNumFeat, lr->nNumClass, lr->nSizeSample,
        lr->naClassIdx, lr->daFeatures, x, lr->daPIJ, lr->daQIJ,
        lr->daExpIJ, lr->daSumExpJ);

    /* FUNCTION COMPUTATION STARTS HERE */
    logRegFormFunctionFromP(lr, lr->nNumClass, lr->nSizeSample, lr->naClassIdx, lr->daPIJ, f);

    /* GRADIENT COMPUTATION STARTS HERE */
    if (g != NULL) {
        memset(g, 0, lr->n * sizeof(LOGREG_REAL));
        for (i = 0; i < lr->nNumSampleUnits; i++) {
            logRegFormGradient(lr, lr->nNumFeat, lr->nNumClass,
                lr->naGIJsizes[i], lr->daFeatures + lr->nNumFeat
                    * lr->naGIJstartIdx[i], lr->daQIJ + lr->nNumClass
                    * lr->naGIJstartIdx[i], lr->pdaGIJ[i]);
            cblas_axpy(lr->n, 1.0, lr->pdaGIJ[i], 1, g, 1);
        }
        for (i = 0; i < lr->n; i++) {
            g[i] /= lr->nNumSampleUnits;
        }
#if 0
        logRegFormGradient(lr, lr->nNumFeat, lr->nNumClass, lr->nSizeSample,
                    lr->daFeatures, lr->daQIJ, g);
#endif
    }

    lr->status = FUNCTION_GRADIENT_DONE;
}

void logRegValidate(const LOGREG_REAL * const x, LOGREG_REAL * const f, LogRegContext *lr) {
    if (lr->nNumTotalTest <= 0 ||
        NULL == lr->daTestFeatures ||
        NULL == lr->naTestClassIdx) {
        fprintf(stderr, "Testing data is not initialized. \n");
        *f = -999;
        return;
    }

    logRegFullFunctionGradientInternal(x, NULL, f, lr->daTestFeatures, lr->naTestClassIdx, lr->nNumTotalTest, lr);
}

void logRegFullFunctionGradient(const LOGREG_REAL * const x, LOGREG_REAL * const g,
    LOGREG_REAL * const f, LogRegContext * lr) {
    logRegFullFunctionGradientInternal(x, g, f, lr->daDataFeatures, lr->naDataClassIdx, lr->nNumTotal, lr);
}
