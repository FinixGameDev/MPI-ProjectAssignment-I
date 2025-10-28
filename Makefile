a.out: main.c
	mpicc main.c -lm && mpirun -np 4 --oversubscribe a.out input6