# Large Scale Optimization

This repository contains a set of scripts to perform:

* Limited Memory BFGS
* Stochastic Conjugate Gradient Truncated Newton
* Stochastic Limited Memory BFGS

Both this methods are meant to work on a large scale optimization environment.

## Compilation

make

## Execution

During the execution, the program is going to ask the user for the dimention of the problem. The three functions that the code is going to optimize are the following:


* Characteristics: Convex.

```
double test_func(double* x, int length){
  double res;
  int i;
  for(res = i = 0; i < length; i++){
    res = res + (3 - x[i])*(4 - x[i]) + x[i]*x[i]*x[i]*x[i];
  }
  return res;
};
```

* Characteristics: Ill conditioned.

```
double testFunc(double* x, int length){
  int i;
  double sum;
  sum = 0;
  for(i = 0; i < length; i++){
    sum = sum  + (100 - i) * (x[i] * x[i]) + exp(x[i]);
  }
  // Return result.
  return sum;
};
```

* Logistic regression on iris dataset.

