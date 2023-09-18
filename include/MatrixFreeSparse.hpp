#ifndef __MATRIX_FREE_SPARSE_HPP__
#define __MATRIX_FREE_SPARSE_HPP__

#include <Eigen/Sparse>

class MatrixFreeSparse;
using Eigen::SparseMatrix;
 
namespace Eigen {
namespace internal {
  // MatrixFreeSparse looks-like a SparseMatrix, so let's inherits its traits:
  template<>
  struct traits<MatrixFreeSparse> :  public Eigen::internal::traits<Eigen::SparseMatrix<double> >
  {};
}
}
 
// Example of a matrix-free wrapper from a user type to Eigen's compatible type
// For the sake of simplicity, this example simply wrap a Eigen::SparseMatrix.
class MatrixFreeSparse : public Eigen::EigenBase<MatrixFreeSparse> {
public:
  // Required typedefs, constants, and method:
  typedef double Scalar;
  typedef double RealScalar;
  typedef int StorageIndex;
  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic,
    IsRowMajor = false,
  };
 
  EIGEN_CONSTEXPR Index rows() const EIGEN_NOEXCEPT { return numElements; }
  EIGEN_CONSTEXPR Index cols() const EIGEN_NOEXCEPT { return numElements; }
 
  template<typename Rhs>
  Eigen::Product<MatrixFreeSparse,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
    return Eigen::Product<MatrixFreeSparse,Rhs,Eigen::AliasFreeProduct>(*this, x.derived());
  }

  MatrixFreeSparse(){}

  MatrixFreeSparse(StorageIndex _numElements, Eigen::Matrix<double, 24, 24> _elementStiffnessMat, Eigen::Array<int, Eigen::Dynamic, 8> _elementToNode, Eigen::ArrayXi _fixedNodes) 
                  : numElements(_numElements),
                    elementStiffnessMat(_elementStiffnessMat),
                    elementToNode(_elementToNode),
                    fixedNodes(_fixedNodes) {}


  const Eigen::Matrix<double, 24, 24> elementStiffnessMat;
	Eigen::Array<int, Eigen::Dynamic, 8> elementToNode;
  Eigen::ArrayXi fixedNodes;

  //Necessary structs for Multigrid 
  void PrepareMultigrid(int _numLevels,
                        std::vector<Eigen::Array<int, Eigen::Dynamic, 8>> _elementToNodeMatrices,
                        std::vector<Eigen::Array<int, Eigen::Dynamic, 27>> _restrictionMappings,
                        std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> _restrictionCoefficients,
                        std::vector<Eigen::SparseMatrix<double>> _interpolationMatrices,
                        std::vector<Eigen::VectorXd> _invDiagKOnLevels,
                        Eigen::SparseMatrix<double> _Kc,
                        Eigen::ArrayXi _coarseFreeDoFs)
  {
    numLevels = _numLevels;
    elementToNodeMatrices = _elementToNodeMatrices;
    restrictionMappings = _restrictionMappings;
    restrictionCoefficients = _restrictionCoefficients;
    interpolationMatrices = _interpolationMatrices;
    invDiagKOnLevels = _invDiagKOnLevels;
    Kc = _Kc;
    coarseFreeDoFs = _coarseFreeDoFs;
  }

  int numLevels;
  std::vector<Eigen::Array<int, Eigen::Dynamic, 8>> elementToNodeMatrices;
  std::vector<Eigen::Array<int, Eigen::Dynamic, 27>> restrictionMappings;
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> restrictionCoefficients;
  std::vector<Eigen::SparseMatrix<double>> interpolationMatrices;
  std::vector<Eigen::VectorXd> invDiagKOnLevels;
  Eigen::SparseMatrix<double> Kc;
  Eigen::ArrayXi coarseFreeDoFs;
  
  StorageIndex numElements;
};


namespace Eigen {
namespace internal {

  template<typename Rhs>
  struct generic_product_impl<MatrixFreeSparse, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
  : generic_product_impl_base<MatrixFreeSparse,Rhs,generic_product_impl<MatrixFreeSparse,Rhs> >
  {
    typedef typename Product<MatrixFreeSparse,Rhs>::Scalar Scalar;
 
    template<typename Dest>
    static void scaleAndAddTo(Dest& dst, const MatrixFreeSparse& lhs, const Rhs& rhs, const Scalar& alpha)
    {
      const Array<int, 1, 24> c   {0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2};
      const Array<int, 1, 24> xInd{0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7};

      for(auto line : lhs.elementToNode.rowwise())
      {
        Array<int, 1, 24> xs = 3 * line(xInd) + c;
        dst(xs) += alpha * lhs.elementStiffnessMat * rhs(xs);
      }

      dst(3 * lhs.fixedNodes    ) = Eigen::ArrayXd::Zero(lhs.fixedNodes.size());
      dst(3 * lhs.fixedNodes + 1) = Eigen::ArrayXd::Zero(lhs.fixedNodes.size());
      dst(3 * lhs.fixedNodes + 2) = Eigen::ArrayXd::Zero(lhs.fixedNodes.size());
    };
  };
}
}

#endif // __MATRIX_FREE_SPARSE_HPP__