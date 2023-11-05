#include "FEMDefines.hpp"
#include "MatrixFreeSparse.hpp"

#include <vector>

template< typename scalar>
void solveWithKompute(const MatrixFreeSparse<scalar>& systemMatrix, const std::vector<scalar>& f, std::vector<scalar>& u);