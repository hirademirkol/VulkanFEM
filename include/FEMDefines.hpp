// Uncomment to set a max number of iterations
// Defaults are twice the number of columns on a matrix (Eigen), or 30000 (Kompute)
// #define MAX_ITER 5000

// Uncomment to set a tolerance for the solver
// Default is 1e-16
#define TOLERANCE 1e-6

// Comment to use explicitly built Matrix (Very high memory usage!)
// Only implemented on CPU
#define MATRIX_FREE

#ifdef MATRIX_FREE
    #include "MatrixFreeSparse.hpp"

    // Uncommment to use multigrid preconditiner
    #define MULTIGRID
#endif