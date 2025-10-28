#Poject Assignement 1

## Objective
Solve the All-Pairs Shortest Path problem using MPI and Fox's Algorithm.

## Implementation

### Phase 1 - Environment Setup
We start by setting up a base structure for our D1 Matrices like so:
```c
typedef struct
{
int size;
int* data;
} Matrix;
```
And setting up some read, write and saving functions:
```c
Matrix read_file(char *filename) {

	Matrix matrix;
	FILE *file = fopen(filename, "r");

	if (file == NULL) {

		printf("Error: Could not open file %s\n", filename);
		Matrix empty_matrix = {0, NULL};
		return empty_matrix;
	}

	fscanf(file, "%d", &matrix.size);
	matrix.data = (int*)malloc(matrix.size * matrix.size * sizeof(int));

	for (int i = 0; i < matrix.size * matrix.size; i++) {
		fscanf(file, "%d", &matrix.data[i]);
	} 

	fclose(file);
	return matrix;
}
```
```c
void print_matrix(Matrix matrix) {
	printf("%d\n", matrix.size);
	for (int i = 0; i < matrix.size; i++) {
		for (int j = 0; j < matrix.size; j++) {
			printf("%d ", matrix.data[i * matrix.size + j]);
		}
	printf("\n");
	}
}
```
```c
void save_matrix(char *filename, Matrix matrix) {
	FILE *file = fopen(filename, "w");

	fprintf(file, "%d\n", matrix.size);
	for (int i = 0; i < matrix.size; i++) {
		for (int j = 0; j < matrix.size; j++) {
			fprintf(file, "%d ", matrix.data[i * matrix.size + j]);
		}
	fprintf(file, "\n");
	}
}
```
And then setup a basic makefile for debugging:
```makefile
a.out: main.c

mpicc main.c && mpirun -np 4 --oversubscribe a.out input6
```
For the main function we setup the basic MPI boiler plate and test our functions in Rank 0 to avoid duplication.
```c
int main(int argc, char **argv) {

	int my_rank;
	int num_procs;

	Matrix matrix;«
  
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	if (my_rank == 0)
	{
		matrix = read_file(argv[1]);
		if (matrix.size == 0) {
			return -1;
		}
		
		print_matrix(matrix);
		char output_filename[20];
		sprintf(output_filename, "output%d", matrix.size);
		save_matrix(output_filename, matrix);
	}
	MPI_Finalize();
	return 0;
}
```
After compiling and running we get the output both in the terminal and in a new file called "output6":
```
6
0 2 0 5 0 0 
0 0 0 0 0 0 
0 2 0 0 0 5 
0 0 0 0 1 0 
3 9 3 0 0 0 
0 0 0 0 1 0 
```

### Phase 2 Process Distribution
First we put some safeguards, like checking if we have an integer Q:
```c
grid_size = (int)sqrt(num_procs);
if (my_rank == 0)
{
	if (num_procs != grid_size * grid_size || matrix.size % grid_size != 0) {
		printf("Error: Number of processes must be a perfect square and matrix size must be divisible by sqrt(number of processes).\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
}
```
We star getting some compiler errors here: like using MPI_Abort causing the program to prematurely kill itself just by existing ? Not even being called, and the math.h (that contains sqrt()) library requiring to add -lm in the end of the compilation process.
