/*
 * main.c
 *
 *  Created on: May 28, 2010
 *      Author: YuchenWu
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <cblas.h>
#include <logRegDouble.h>
#include "gpcg.h"

static void evalFG(const double * const x,
    double * const g, double * const f, const int bResample,
    void * userParam) {
    LogRegContext * lr = (LogRegContext *) userParam;
    logRegSampleFunctionGradient(x, g, f, bResample, lr);
}

static void evalFullFG(const double * const x,
    double * const g, double * const f, void * const userParam) {
    LogRegContext * lr = (LogRegContext *) userParam;
    logRegFullFunctionGradient(x, g, f, lr);
}

static void evalHessMv(const double * const x, const double * const v,
    double * const d, void * userParam){
    LogRegContext * lr = (LogRegContext *) userParam;
    logRegSampleHessianVector(x, v, d, lr);
}

#define N 30315     /* Total number of features (problem dimension)
						 Number of classes 129
						 Number of features 235*/

#define T 191607    // Total number of data points
#define S 191607    // Sample size for gradients

int main(int argc, char **argv) {
    LogRegContext * lr;
    GpcgOptions * opt;
    GpcgStatus status;
    int nCgIter = 5;
    double p = 0.05;
    double xinit[N];
    int i;

    for (i=0; i<N; i++) {
        xinit[i] = 0;
    }

    if (argc > 1) {
        nCgIter = atoi(argv[1]);
        printf("Changed nCgIter to %d.\n", nCgIter);
    }

    // Initialize the logistic model.
    lr = newLogReg();
    logRegInit(lr, 235, 129, (int)(p*S), S, T, 0, 1);
    logRegLoadData(lr, "../train_features-500.squared.mle");

    // Set the GPCG options.
    opt = gpcgGetDefaultOptions();
    opt->dHessRho = 1.0;
    opt->dEyeRho = 0;
    opt->dL1 = 10.0 / T;
    opt->nCgMaxIter = nCgIter;
    opt->nMaxIteration = 50;
    opt->dCgTol = 1e-4;

    opt->bResample = (S < T);

    // Run the GPCG algorithm.
    status = gpcg(lr->n, xinit, evalFG, evalHessMv, evalFullFG, opt, lr); 

    printf("Termination status: %d\n", status);

    freeLogReg(lr);
    return 0;
}
