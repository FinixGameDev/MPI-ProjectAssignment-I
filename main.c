#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <limits.h>
#include <string.h>

#define MATRIX unsigned long
#define MPI_MATRIX MPI_UNSIGNED_LONG

// Function to read the input file, allocate memory for the blocks and matrices and init the data
void proc_init (MATRIX** A, MATRIX** B, MATRIX** C, MATRIX** local_A, MATRIX** local_B, MATRIX** local_C,MATRIX** inital_A, int *size, int *block_size, char* input );
// Function to initialize the matrix with the data inside the input file
void data_init(MATRIX* A, MATRIX* B, MATRIX* C, int size, char * input);
// Free all allocated memory
void proc_end (MATRIX* A, MATRIX* B, MATRIX* C, MATRIX* local_A, MATRIX* local_B, MATRIX* local_C,MATRIX* inital_A);
// Function for creating the two-dimensional grid communicator and communicators for each row and each column of the grid
void setup_grid();
// Self-explanatory function to print the matrix
void print_matrix(MATRIX *matrix, int size);
// Distribution of data among the processes
void distribute_data(MATRIX* A, MATRIX* B, MATRIX* C, MATRIX* inital_A, MATRIX* local_B, MATRIX* local_C, int size, int block_size);
// Function for checkerboard matrix decomposition
void matrix_scatter(MATRIX* matrix, MATRIX* block, int size, int block_size);
// Function for parallel execution of the Fox method
void fox_algorithmn(MATRIX* local_A, MATRIX* inital_A, MATRIX* local_B, MATRIX* local_C, int block_size);
// Join final results in Matrix C
void gather_results (MATRIX* C, MATRIX* local_C, int size, int block_size);
// Set INT_MAX values back to 0 in the result matrix
void adjust_values(MATRIX * C, int size);
// Function to compare the result matrix with an expected output file
bool check_results(MATRIX *result, char *expected_result_filename, int size);
// Save Matrix to output
void save_matrix(char *filename, MATRIX* matrix, int size);

//Global variables

// Genral MPI variables
int num_procs = 0;
int my_rank = 0;

// Grid settings and comunicators
int grid_size;
MPI_Comm comm;
int coords[2];
MPI_Comm col_comm;
MPI_Comm row_comm;

int main(int argc, char *argv[]) {
    
    // General vars
    MATRIX* A;
    MATRIX* B;
    MATRIX* C;
    int size;
    
    // Sub matrices vars
    int block_size;
    MATRIX *inital_A;
    MATRIX *local_A;
    MATRIX *local_B;
    MATRIX *local_C;

    // iteration until Df = DgxDg
    int iteration = 1;

    //I/O vars
    char* input;               

    // Profiling vars
    double start, end, time; 

    // Inititialize MPI environment
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    // Validating the grid size to be the perfect sqrt of the number of processes
    grid_size = sqrt((double)num_procs);
    if (num_procs != grid_size*grid_size) {
        if (my_rank == 0) {
            printf ("ERROR: Invalid configuration! \n");
        }

        MPI_Finalize();
        return -1;
    }

    if (my_rank == 0)
    {
        if (argc < 2){
            printf("ERROR: You must run the program as 'mpirun -np {numprocess} program {filename}");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        else{
            input = argv[1];
        }

    }
        
    MPI_Barrier(MPI_COMM_WORLD);

    //Initialize matrices and sub-matrices    
    proc_init(&A, &B, &C, &local_A, &local_B, &local_C, &inital_A, &size, &block_size, input);
    setup_grid();
        
    MPI_Barrier(MPI_COMM_WORLD);
        
    start = MPI_Wtime();
    
    
    do {
        distribute_data(A, B, C, inital_A, local_B, local_C, size, block_size);
        fox_algorithmn(local_A, inital_A, local_B, local_C, block_size);
        gather_results(C, local_C, size, block_size);
            
        if (my_rank == 0){
            memcpy(A, C, sizeof(MATRIX) * size * size);
            memcpy(A, C, sizeof(MATRIX) * size * size);
            iteration = iteration * 2;
        }
        MPI_Bcast(&iteration, 1, MPI_INT, 0, MPI_COMM_WORLD);
            
    }while (iteration < size);
        
    end = MPI_Wtime();

    // Print final results
    if (my_rank == 0) {
        adjust_values(C, size);
        printf("Final Result Matrix:\n");
        print_matrix(C, size);
    }

                
    time = end-start;
    if (my_rank == 0) {
        //Pofile and check results to a cheat file
        printf("Time of execution = %f\n", time);

        //Result verification
        printf("Would you like to compare the results with an expected output file? (y/n): \n");
        char response;
        scanf("%c", &response);

        if (response == 'y' || response == 'Y') {
            char expected_result_filename[256];
            printf("Enter the expected result filename: \n");
            scanf("%s", expected_result_filename);

            if (check_results(C, expected_result_filename, size)) 
                printf("Results are correct.\n");
            else 
                printf("Results are incorrect.\n");
        }
        else
           printf("Skipping result comparison.\n");

        printf("Saving Matrix to file 'output':");
        save_matrix("output", C, size);
        printf("Matrix saved to file 'output'.");
    }
        
    // Free memory
    proc_end(A, B, C, local_A, local_B, local_C, inital_A);
    
    MPI_Finalize();
    return 0;
}

void adjust_values(MATRIX * C, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if(C[i * size + j] == INT_MAX)
                C[i * size + j] = 0;
        }
    }
}

void fox_algorithmn(MATRIX* local_A, MATRIX* inital_A, MATRIX* local_B, MATRIX* local_C, int block_size) {
    for (int count = 0; count < grid_size; count ++) {
        
        //Send the block of matrix A to the process that will use it
        int pivot = (coords[0] + count) % grid_size;
        // Copy the block of matrix A to be broadcasted
        if (coords[1] == pivot) {
            for (int i = 0; i < block_size*block_size; i++)
            local_A[i] = inital_A[i];
        }
        // Block broadcasting
        MPI_Bcast(local_A, block_size * block_size, MPI_MATRIX, pivot, row_comm);

        // Block multiplication
        for (int i=0; i<block_size; i++) {
            for (int j=0; j<block_size; j++) {
                for (int k=0; k<block_size; k++) {
                    int new_min = local_A[i*block_size+k] + local_B[k*block_size+j];

                    if (new_min < local_C[i*block_size+j])
                        local_C[i*block_size+j] = new_min;
                }
            }
        }

        
        // Shift block of matrix B upwards by one
         MPI_Status Status;
    
        int next_proc = coords[0] + 1;
        if ( coords[0] == grid_size-1 )
            next_proc = 0;
    
        int prev_proc = coords[0] - 1;
        if ( coords[0] == 0 )
            prev_proc = grid_size-1;

        MPI_Sendrecv_replace( local_B, block_size*block_size, MPI_MATRIX, next_proc, 0, prev_proc, 0, col_comm, &Status);
    }
}

void gather_results (MATRIX* C, MATRIX* local_C, int size, int block_size) {
    MATRIX * result_row = malloc(sizeof(MATRIX) * size * block_size);
    for (int i=0; i < block_size; i++) {
        MPI_Gather( &local_C[i*block_size], block_size, MPI_MATRIX, &result_row[i*size], block_size, MPI_MATRIX, 0, row_comm);
    }
    if (coords[1] == 0) {
        MPI_Gather(result_row, block_size*size, MPI_MATRIX, C, block_size*size, MPI_MATRIX, 0, col_comm);
    }
    free(result_row);
}

void distribute_data(MATRIX* A, MATRIX* B, MATRIX* C, MATRIX* inital_A, MATRIX* local_B, MATRIX* local_C, int size, int block_size) {
    matrix_scatter(A, inital_A, size, block_size);
    matrix_scatter(B, local_B, size, block_size);
    matrix_scatter(C, local_C, size, block_size);
}

void matrix_scatter(MATRIX* matrix, MATRIX* block, int size, int block_size)
{
    MATRIX * matrix_row = malloc(sizeof(MATRIX) * block_size * size);
    if (coords[1] == 0) {
        MPI_Scatter(matrix, block_size * size, MPI_MATRIX, matrix_row, block_size * size, MPI_MATRIX, 0, col_comm);
    }
    for (int i = 0; i < block_size; i++) {
        MPI_Scatter(&matrix_row[i*size], block_size, MPI_MATRIX, &(block[ i * block_size]), block_size, MPI_MATRIX, 0, row_comm);
    }
    free(matrix_row);
}

void setup_grid() {
    int dims[2];  // Number of processes per dimension of the grid
    int periodic[2]; // =1, if the grid dimension should be periodic
    int sub_dims[2]; // =1, if the grid dimension should be fixed
    
    dims[0] = grid_size;
    dims[1] = grid_size;
    periodic[0] = 1;
    periodic[1] = 1;
    
    // Creation of the Cartesian communicator grid_size x grid_size with "paradox"
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periodic, 1, &comm);
    
    // Determination of the cartesian coordinates for every process
    MPI_Cart_coords(comm, my_rank, 2, coords);
    
    // Createing communicators for rows
    sub_dims[0] = 0;
    sub_dims[1] = 1;
    MPI_Cart_sub(comm, sub_dims, &row_comm);
    
    // Creating communicators for columns
    sub_dims[0] = 1;
    sub_dims[1] = 0; 
    MPI_Cart_sub(comm, sub_dims, &col_comm);
}


void print_matrix(MATRIX *matrix, int size){
    int i, j;  // Loop variables
    for (i=0; i < size; i++) {
        for (j = 0; j < size; j++) {
            printf("%ld ", matrix[i*size+j]);
        }
        printf("\n");
    }
}

// Function for memory allocation and initialization of matrix elements
void proc_init (MATRIX** A, MATRIX** B, MATRIX** C, MATRIX** local_A, MATRIX** local_B, MATRIX** local_C, MATRIX** inital_A, int *size, int *block_size, char* input ){
    // Setting the size of matrices
    
    if (my_rank == 0){
        FILE * file;
        file = fopen(input, "r");
        
        if (file == NULL)
        {
            printf("Error! Could not open file\n");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }

        char str[128];
        fscanf(file, "%s", str);
        *size = atoi(str);

        fclose(file);
        
        // Validation
        if ((*size) % grid_size != 0){
            printf("\nERROR:size of matrices must be divisible by the grid size!\n");
            printf("ERROR: Invalid configuration!\n");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Bcast(size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Allocating the memmory
    *block_size = (*size)/grid_size;
    (*local_A) = malloc(sizeof(MATRIX) * (*block_size) * (*block_size));
    (*local_B) = malloc(sizeof(MATRIX) * (*block_size) * (*block_size));
    (*local_C) = malloc(sizeof(MATRIX) * (*block_size) * (*block_size));
    (*inital_A) = malloc(sizeof(MATRIX) * (*block_size) * (*block_size));
    
    if (my_rank == 0) {
        (*A) = malloc(sizeof(MATRIX) * (*size) * (*size));
        (*B) = malloc(sizeof(MATRIX) * (*size) * (*size));
        (*C) = malloc(sizeof(MATRIX) * (*size) * (*size));
        
        // Initialize the matrices
        data_init(*A, *B, *C, *size, input);
    }
    
    // Init the result block of each process (c)
    for (int i = 0; i < (*block_size) * (*block_size); i++) {
        (*local_C)[i] = 0;
    }
}

void data_init(MATRIX* A, MATRIX* B, MATRIX* C, int size, char * input){
    int i, j;  // Loop variables
    
    FILE * file;
    file = fopen(input, "r");
    
    if (file == NULL)
    {
        printf("Error! Could not open file\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    
    // reading the size (again) but we don't need it now
    char str[128];
    fscanf(file, "%s", str);
    
    for (i=0; i < size; i++){
        for (j=0; j<size; j++) {
            fscanf(file, "%s", str);
            
            MATRIX number = atoi(str);
            if(number == 0 && i != j){
                A[i*size+j] = INT_MAX;
                B[i*size+j] = INT_MAX;
                C[i*size+j] = INT_MAX;
            }
            else {
                A[i*size+j] = number;
                B[i*size+j] = number;
                C[i*size+j] = number;
            }
        }
    }
    fclose(file);
}


// Function for computational process termination
void proc_end (MATRIX* A, MATRIX* B, MATRIX* C, MATRIX* local_A, MATRIX* local_B, MATRIX* local_C,MATRIX* inital_A){
    if (my_rank == 0) {
        free(A);
        free(B);
        free(C);
    }
    free(local_A);
    free(local_B);
    free(local_C);
    free(inital_A);
}

bool check_results(MATRIX *result, char *expected_result_filename, int size) {
    
    FILE* file;
    file = fopen(expected_result_filename, "r");
    if (file == NULL) {
        printf("Error: Could not open file %s\n", expected_result_filename);
        return false;
    }

    char str[128];

    MATRIX* expected_result = malloc(sizeof(MATRIX) * size * size);   
    
    for (int i=0; i < size; i++){
        for (int j=0; j<size; j++) {
            fscanf(file, "%s", str);
            
            MATRIX number = atoi(str);
            expected_result[i*size+j] = number;
        }
    }
    
    fclose(file);
    for (int i = 0; i < size * size; i++) {
        if (result[i] != expected_result[i]) {
            printf("Mismatch at index %d: expected %ld, got %ld\n", i, expected_result[i], result[i]);
            return false;
        }
    }
    return true;
}

void save_matrix(char *filename, MATRIX* matrix, int size) {
    FILE *file = fopen(filename, "w");

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size - 1; j++) {
            fprintf(file, "%ld ", matrix[i * size + j]);
        }
        fprintf(file, "%ld\n", matrix[i * size + size - 1]);
    }

    fclose(file);
}