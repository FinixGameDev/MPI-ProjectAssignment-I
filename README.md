fox - MPI min-plus matrix multiplication (Fox) + repeated squaring

Build:

    make

Run (example using 4 processes):

    mpirun -np 4 ./fox < example.in

Input format: first line N, then NxN integers (0 means no edge except diagonal)

This program requires the number of processes P to be a perfect square Q*Q and N %% Q == 0.
