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

  MatrixFreeSparse() : elementStiffnessMat(0), elementToNode(0) {}

  MatrixFreeSparse(StorageIndex _numElements, Eigen::Matrix<double, 300, 1> _elementStiffnessMat, std::vector<std::array<uint64_t, 8>> _elementToNode) 
                  : numElements(_numElements),
                    elementStiffnessMat(_elementStiffnessMat),
                    elementToNode(_elementToNode) {}

  const Eigen::Matrix<double, 300, 1> elementStiffnessMat;
	std::vector<std::array<uint64_t, 8>> elementToNode;
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
 
      // Here we could simply call dst.noalias() += lhs.my_matrix() * rhs,
      // but let's do something fancier (and less efficient):
      for(auto line : lhs.elementToNode)
      {
        for(int i = 0; i < 8; i++)
          for(int j = 0; j < 8; j++)
          {
            uint64_t node1 = line[i];
            uint64_t node2 = line[j];
            
            if(node1 == -1 || node2 == -1)
					  continue;

            for(int c1 = 0; c1 < 3; c1++)
            {
              for(int c2 = 0; c2 < 3; c2++)
              {
                int iMatrix = node1*3 + c1;
                int jMatrix = node2*3 + c2;
                int elementIndex = i <= j ? (i*(47-i))/2 + j + 1 : (j*(47-j))/2 + i + 1;

                dst(i) += lhs.elementStiffnessMat(elementIndex) * rhs(j);
              }
            }
          }
      }
    };
  };
}
}

#endif // __MATRIX_FREE_SPARSE_HPP__