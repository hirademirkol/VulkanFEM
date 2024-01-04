#include "FEMDefines.hpp"
#include "MatrixFreeSparse.hpp"

#include <vector>

/// @brief Solve the system with kompute package on the GPU
/// @tparam scalar Only double implemented
/// @param systemMatrix (in) System matrix to solve
/// @param f (in) RHS vector
/// @param u (out) Solution vector
template< typename scalar>
void solveWithKompute(const MatrixFreeSparse<scalar>& systemMatrix, const std::vector<scalar>& f, std::vector<scalar>& u);