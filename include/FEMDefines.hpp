// #define MAX_ITER 5000
#define TOLERANCE 1e-6

#define MATRIX_FREE

#ifdef MATRIX_FREE
    #include "MatrixFreeSparse.hpp"

    #define MULTIGRID

    #ifdef MULTIGRID
        #define NUM_LEVELS 3
    #endif

#endif