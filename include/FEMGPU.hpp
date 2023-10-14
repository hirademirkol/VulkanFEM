#include "FEMDefines.hpp"

#include <vector>

#include <kompute/Kompute.hpp>

#include "MatrixFreeSparse.hpp"

// Compiled shaders
#include "MatxVec.hpp"
#include "FixedNodes.hpp"
#include "VecDotVec.hpp"
#include "GaussSeidel.hpp"
#include "ConjugateGradient_1.hpp"
#include "ConjugateGradient_2.hpp"

#ifdef MULTIGRID
#include "CWiseMult.hpp"
#include "CWiseAdd.hpp"
#include "Smooth.hpp"
#include "Restrict.hpp"
#include "Interpolate.hpp"
#endif

typedef std::shared_ptr<kp::Tensor> Tensor;
typedef std::shared_ptr<kp::Sequence> Sequence;
typedef std::shared_ptr<kp::Algorithm> Algorithm;
using TensorDataTypes = kp::Tensor::TensorDataTypes;
using TensorTypes = kp::Tensor::TensorTypes;

template< typename scalar>
void solveWithKompute(const MatrixFreeSparse<scalar>& systemMatrix, const std::vector<scalar>& f, std::vector<scalar>& u);