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
 
  Index rows() const { return numElements; }
  Index cols() const { return numElements; }
 
  template<typename Rhs>
  Eigen::Product<MatrixFreeSparse,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
    return Eigen::Product<MatrixFreeSparse,Rhs,Eigen::AliasFreeProduct>(*this, x.derived());
  }

  MatrixFreeSparse(){}

  MatrixFreeSparse(StorageIndex _numElements, Eigen::Matrix<double, 24, 24> _elementStiffnessMat, Eigen::Array<int, Eigen::Dynamic, 8> _elementToNode) 
                  : numElements(_numElements),
                    elementStiffnessMat(_elementStiffnessMat),
                    elementToNode(_elementToNode) {}

  const Eigen::Matrix<double, 24, 24> elementStiffnessMat;
	Eigen::Array<int, Eigen::Dynamic, 8> elementToNode;
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
      // This method should implement "dst += alpha * lhs * rhs" inplace,
      // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
      assert(alpha==Scalar(1) && "scaling is not implemented");
      EIGEN_ONLY_USED_FOR_DEBUG(alpha);
 
      const Array<int, 1, 24> c{0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2};
      const Array<int, 1, 24> xInd{0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7};

      for(auto line : lhs.elementToNode.rowwise())
      {
        Array<int, 1, 24> xs = 3 * line(xInd) + c;
        dst(xs) += lhs.elementStiffnessMat * rhs(xs);       
      }
    };
  };
}
}

#endif // __MATRIX_FREE_SPARSE_HPP__