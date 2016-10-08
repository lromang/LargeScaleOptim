CC=gcc

all: gen_optim

gen_optim:
	$(CC) -o code main.c -L gen_optim.c -L line_alg.c -L utileries.c -L lbfgs.c -L trunc_newton.c -L slm.c -lm
clean:
	rm code
