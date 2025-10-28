#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>
#include <stddef.h>
#include <math.h>
#include <stdbool.h>


typedef struct 
{
    int size;
    int* data;
} Matrix;

Matrix read_file(char *filename);
void print_matrix(Matrix matrix);
void save_matrix(char *filename, Matrix matrix);
void FoxAlgorithm(Matrix *A, Matrix *B, Matrix *C, int rank, int size);
bool checkresults(Matrix *result, char* expected_result_filename);

    
int main(int argc, char **argv) {

    //1: Initialize MPI environment and get the input file name from command line arguments
    int my_rank;
    int num_procs;
    int grid_size;

    Matrix matrix;

    Matrix result;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    //1.5 Check if matrix size is compatible with number of processes
    grid_size = (int)sqrt(num_procs);
    if (my_rank == 0)
    {
        if (num_procs != grid_size * grid_size || matrix.size % grid_size != 0) {
            printf("Error: Number of processes must be a perfect square and matrix size must be divisible by sqrt(number of processes).\n");
            //MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }
    

    if (my_rank == 0)
    {
        matrix = read_file(argv[1]);
        if (matrix.size == 0) {
            printf("Error: Matrix size is zero.\n");
        }
        else
        {
            print_matrix(matrix);
            char output_filename[20];
            sprintf(output_filename, "output%d", matrix.size);
            save_matrix(output_filename, matrix);
            printf("Grid size: %d\n", grid_size);

            if (checkresults(&matrix, "output6")) {
                printf("Results are correct.\n");
            } else {
                printf("Results are incorrect.\n");
            }
        }
    
    }

    //2 Distribute Matrix among processes

    //Apply Fox's Algorithm

    //Gather results

    //Print Matrix and output to file

    
    MPI_Finalize();

    return 0;
}

void FoxAlgorithm(Matrix *A, Matrix *B, Matrix *C, int rank, int size) {
    // Implement Fox's algorithm for matrix multiplication here
}

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

void print_matrix(Matrix matrix) {
    printf("%d\n", matrix.size);
    for (int i = 0; i < matrix.size; i++) {
        for (int j = 0; j < matrix.size; j++) {
            printf("%d ", matrix.data[i * matrix.size + j]);
        }
        printf("\n");
    }
}

void save_matrix(char *filename, Matrix matrix) {
    FILE *file = fopen(filename, "w");

    fprintf(file, "%d\n", matrix.size);
    for (int i = 0; i < matrix.size; i++) {
        for (int j = 0; j < matrix.size - 1; j++) {
            fprintf(file, "%d ", matrix.data[i * matrix.size + j]);
        }
        fprintf(file, "%d\n", matrix.data[i * matrix.size + matrix.size - 1]);
    }
}

bool checkresults(Matrix *result, char* expected_result_filename){
    Matrix expected_result = read_file(expected_result_filename);
    if (expected_result.size != result->size) {
        return false;
    }
    for (int i = 0; i < result->size * result->size; i++) {
        if (expected_result.data[i] != result->data[i]) {
            return false;
        }
    }
    return true;
}