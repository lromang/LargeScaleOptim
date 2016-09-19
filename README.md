# Large Scale Optimization

This repository contains a set of scripts to perform:

* Stochastic Conjugate Gradient Truncated Newton
* Stochastic Limited Memory BFGS

Both this methods are meant to work on a large scale optimization environment.

## Compilation

gcc -o code main.c -L gen_optim.c -L line_alg.c -L utileries.c -lm

## Execution

After execution, the program is going to ask you for the size of the the dimention of the problem. The two functions that it is going to optimize are the following:


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


