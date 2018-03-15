# NMFSolver
This C++ library proposes a Non-negative Matrix Factorization solver. 

## Installation
First, this library requires Armadillo.

### Method 1
    Build the lib.
    Add the libNMFSolver.a and the headers files to your project.
    Add the lib '-lNMFSolver' as library for the compilation
    Include "NMFSolver.hpp"

### Method 2

    Copy the files to your source directory
    Include "NMFSolver.hpp"

## Use

The principal class of this library is the NMFSolver class.
The constructor takes as arguments:

    A The matrix to be factorize
    W The W matrix such that *A = WH*
    H The H matrix such that *A = WH*

Note that the dimensions of *W* and *H* must be already defined and guide the solving.

The other arguments are the following:

    gradient_method: 0 for KL multiplicative update, 2 for L2 additive, 21 for L2 without coordinate descent W/H (Default 0)
    init_method: 0 initialize W with random columns of A, 1 initialize W with random values [0,1], both initialize H with random values [0,1] (Default 1)
    sparsity_coefficient: Sparsity coefficient for the gradient update (Default 0.001)
    time_out_in_second: Stop criteria in seconds (Default 3)
    number_of_iteration_step: Stop criteria in number of iterations (Default 200)
    convergence_stop: Stop criteria of convergence (Default 1e-6)
    W_fix: If true, the algorithm does not modify W (Default false)
    H_fix: If true, the algorithm does not modify H (Default false)
    verbose: If true, print the current cost at each iteration (Default false)

To run the solver, call the method solve

    Solver s(A,W,S);
    s.solve();


To extract the cost of the solution

    double cost = s.lossFunction()




