#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>

// fox - MPI implementation of min-plus matrix multiplication using Fox's algorithm
// followed by repeated squaring to compute all-pairs shortest paths.

#define INF (INT_MAX/4)

static void print_matrix(int *A, int N) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (A[i*N + j] >= INF/2) printf("0");
            else printf("%d", A[i*N + j]);
            if (j+1 < N) printf(" ");
        }
        printf("\n");
    }
}

// Read matrix from stdin on rank 0; broadcast N then matrix
static int *read_matrix(int *N_out, int world_rank) {
    int N = 0;
    int *A = NULL;
    if (world_rank == 0) {
        if (scanf("%d", &N) != 1) {
            fprintf(stderr, "Failed to read N\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        A = (int*)malloc(N*N*sizeof(int));
        if (!A) { fprintf(stderr, "malloc failed\n"); MPI_Abort(MPI_COMM_WORLD, 1); }
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                int x; if (scanf("%d", &x) != 1) { fprintf(stderr, "bad matrix input\n"); MPI_Abort(MPI_COMM_WORLD,1); }
                if (x == 0 && i != j) A[i*N + j] = INF; else A[i*N + j] = x;
            }
        }
    }
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (world_rank != 0) A = (int*)malloc(N*N*sizeof(int));
    MPI_Bcast(A, N*N, MPI_INT, 0, MPI_COMM_WORLD);
    *N_out = N; return A;
}

// Min-plus multiply C = min_over_k (A(i,k) + B(k,j)) for blocks
static void minplus_mul_block(int *A, int *B, int *C, int bs) {
    int elems = bs*bs;
    for (int i = 0; i < elems; ++i) C[i] = INF;
    for (int k = 0; k < bs; ++k) {
        for (int i = 0; i < bs; ++i) {
            int a = A[i*bs + k];
            if (a >= INF/2) continue;
            for (int j = 0; j < bs; ++j) {
                int b = B[k*bs + j];
                if (b >= INF/2) continue;
                int s = a + b;
                if (s < C[i*bs + j]) C[i*bs + j] = s;
            }
        }
    }
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int world_rank, world_size; MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int N; int *D = read_matrix(&N, world_rank);

    // Check P is a perfect square Q*Q and N % Q == 0
    int Q = -1;
    for (int q = 1; q*q <= world_size; ++q) if (q*q == world_size) { Q = q; break; }
    if (Q == -1) {
        if (world_rank == 0) fprintf(stderr, "Number of processes must be a perfect square (P=Q*Q)\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (N % Q != 0) {
        if (world_rank == 0) fprintf(stderr, "Matrix dimension N must be divisible by Q (N %% Q == 0)\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int bs = N / Q; // block size

    // Create cartesian grid QxQ
    int dims[2] = {Q, Q}; int periods[2] = {0, 0};
    MPI_Comm grid; MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &grid);
    int grid_rank; MPI_Comm_rank(grid, &grid_rank);
    int coords[2]; MPI_Cart_coords(grid, grid_rank, 2, coords);
    int prow = coords[0], pcol = coords[1];

    // allocate local blocks
    int block_elems = bs*bs;
    int *Ablock = malloc(block_elems * sizeof(int));
    int *Bblock = malloc(block_elems * sizeof(int));
    int *Cblock = malloc(block_elems * sizeof(int));
    int *tempA = malloc(block_elems * sizeof(int));
    int *tempC = malloc(block_elems * sizeof(int));
    if (!Ablock||!Bblock||!Cblock||!tempA||!tempC) { if (world_rank==0) fprintf(stderr,"alloc fail\n"); MPI_Abort(MPI_COMM_WORLD,1); }

    // distribute global D into local Ablock (initial matrix)
    if (grid_rank == 0) {
        // root prepares and sends
        for (int bi = 0; bi < Q; ++bi) for (int bj = 0; bj < Q; ++bj) {
            int dest_rank; int dest_coords[2] = {bi, bj}; MPI_Cart_rank(grid, dest_coords, &dest_rank);
            int *buf = malloc(block_elems*sizeof(int));
            for (int i = 0; i < bs; ++i) for (int j = 0; j < bs; ++j) {
                int gi = bi*bs + i, gj = bj*bs + j; buf[i*bs + j] = D[gi*N + gj];
            }
            if (dest_rank == 0) memcpy(Ablock, buf, block_elems*sizeof(int)); else MPI_Send(buf, block_elems, MPI_INT, dest_rank, 0, MPI_COMM_WORLD);
            free(buf);
        }
    } else {
        MPI_Recv(Ablock, block_elems, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // initial Bblock := Ablock (we will square D)
    memcpy(Bblock, Ablock, block_elems*sizeof(int));

    // create row and column communicators
    MPI_Comm row_comm, col_comm;
    int remain_row[2] = {0,1}; MPI_Cart_sub(grid, remain_row, &row_comm); // keep columns -> communicator contains processes in same row
    int remain_col[2] = {1,0}; MPI_Cart_sub(grid, remain_col, &col_comm); // keep rows -> communicator contains processes in same column
    int row_rank, col_rank; MPI_Comm_rank(row_comm, &row_rank); MPI_Comm_rank(col_comm, &col_rank);

    int exponent = 1;
    while (exponent < N) {
        // compute Cblock = Ablock (distributed) * Bblock (distributed) using Fox
        for (int i = 0; i < block_elems; ++i) Cblock[i] = INF;

        for (int stage = 0; stage < Q; ++stage) {
            int k = (prow + stage) % Q;
            // root in this row that holds the A block to broadcast has column index k -> its row_rank == k
            if (row_rank == k) {
                // broadcast Ablock from this process to row
                MPI_Bcast(Ablock, block_elems, MPI_INT, k, row_comm);
                // use Ablock as Aused
                // multiply Ablock (as Aused) with local Bblock
                minplus_mul_block(Ablock, Bblock, tempC, bs);
            } else {
                MPI_Bcast(tempA, block_elems, MPI_INT, k, row_comm);
                minplus_mul_block(tempA, Bblock, tempC, bs);
            }
            // accumulate
            for (int t = 0; t < block_elems; ++t) if (tempC[t] < Cblock[t]) Cblock[t] = tempC[t];

            // rotate Bblock upward by 1 in column communicator
            int src = (col_rank + 1) % Q;
            int dst = (col_rank - 1 + Q) % Q;
            MPI_Sendrecv_replace(Bblock, block_elems, MPI_INT, dst, 0, src, 0, col_comm, MPI_STATUS_IGNORE);
        }

        // after local multiply, set Ablock := Cblock; Bblock := Cblock for next squaring
        memcpy(Ablock, Cblock, block_elems*sizeof(int));
        memcpy(Bblock, Cblock, block_elems*sizeof(int));
        exponent *= 2;
    }

    // gather blocks back to grid root (grid_rank==0)
    int *gatherbuf = NULL;
    if (grid_rank == 0) gatherbuf = malloc(world_size * block_elems * sizeof(int));
    MPI_Gather(Ablock, block_elems, MPI_INT, gatherbuf, block_elems, MPI_INT, 0, grid);

    if (world_rank == 0) {
        int *out = malloc(N*N*sizeof(int));
        // gatherbuf is ordered by cartesian ranks 0..P-1
        for (int r = 0; r < world_size; ++r) {
            int rc[2]; MPI_Cart_coords(grid, r, 2, rc);
            int bi = rc[0], bj = rc[1];
            int *blk = &gatherbuf[r*block_elems];
            for (int i = 0; i < bs; ++i) for (int j = 0; j < bs; ++j) {
                int gi = bi*bs + i, gj = bj*bs + j; out[gi*N + gj] = blk[i*bs + j];
            }
        }
        print_matrix(out, N);
        free(out); free(gatherbuf);
    }

    free(Ablock); free(Bblock); free(Cblock); free(tempA); free(tempC); free(D);
    MPI_Finalize(); return 0;
}