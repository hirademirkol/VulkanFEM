#ifndef __MATRIX_FREE_SPARSE_HPP__
#define __MATRIX_FREE_SPARSE_HPP__

#include <Eigen/Sparse>

// template<typename scalar, int options, typename storageIndex>
// class MatrixFreeSparse;

// using Eigen::SparseMatrix;

// namespace Eigen {
// namespace internal {
//   // MatrixFreeSparse looks-like a SparseMatrix, so let's inherits its traits:
//   template<typename scalar, int options, typename storageIndex>
//   struct traits<MatrixFreeSparse<scalar, options, storageIndex>> :  public traits<Eigen::SparseMatrix<scalar, options, storageIndex> >
//   {};

//   template<>
//   struct traits<MatrixFreeSparse<double, 0, int>> :  public traits<Eigen::SparseMatrix<double, 0, int> >
//   {};
//   template<>
//   struct traits<MatrixFreeSparse<float, 0, int>> :  public traits<Eigen::SparseMatrix<double, 0, int> >
//   {};
// }
// }

// template<typename scalar, int options = 0, typename storageIndex = int>
// class MatrixFreeSparse : public Eigen::EigenBase<MatrixFreeSparse<scalar, options, storageIndex>> {
// public:
//   // Required typedefs, constants, and method:
//   typedef scalar Scalar;
//   typedef scalar RealScalar;
//   typedef storageIndex StorageIndex;
//   enum {
//     ColsAtCompileTime = Eigen::Dynamic,
//     MaxColsAtCompileTime = Eigen::Dynamic,
//     IsRowMajor = false
//   };

//   Eigen::Index rows() const { return numElements; }
//   Eigen::Index cols() const { return numElements; }
 
//   template<typename Rhs>
//   Eigen::Product<MatrixFreeSparse<Scalar, options, storageIndex>,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
//     return Eigen::Product<MatrixFreeSparse,Rhs,Eigen::AliasFreeProduct>(*this, x.derived());
//   }
 
//   // Custom API:
//   MatrixFreeSparse() : elementStiffnessMat(0), elementToNode(0) {}

//   MatrixFreeSparse(StorageIndex _numElements, Eigen::Matrix<scalar, 300, 1> _elementStiffnessMat, std::vector<std::array<uint64_t, 8>> _elementToNode) 
//                   : numElements(_numElements),
//                     elementStiffnessMat(_elementStiffnessMat),
//                     elementToNode(_elementToNode) {}

// private:
//   const Eigen::Matrix<scalar, 300, 1> elementStiffnessMat;
// 	std::vector<std::array<uint64_t, 8>> elementToNode;
//   StorageIndex numElements;

// };

// template class MatrixFreeSparse<double, 0, int>;
// template class MatrixFreeSparse<float, 0, int>;

// template Eigen::Product<MatrixFreeSparse<double, 0, int>,
//                         Eigen::Matrix<double,-1,1,0,-1,1>,
//                         Eigen::AliasFreeProduct> MatrixFreeSparse<double, 0, int>::operator*
//                         (const Eigen::MatrixBase<Eigen::Matrix<double,-1,1,0,-1,1>>& x) const;

// template Eigen::Product<MatrixFreeSparse<float, 0, int>,
//                         Eigen::Matrix<float,-1,1,0,-1,1>,
//                         Eigen::AliasFreeProduct> MatrixFreeSparse<float, 0, int>::operator*
//                         (const Eigen::MatrixBase<Eigen::Matrix<float,-1,1,0,-1,1>>& x) const;

class MatrixFreeSparse;
using Eigen::SparseMatrix;
 
namespace Eigen {
namespace internal {
  // MatrixFreeSparse looks-like a SparseMatrix, so let's inherits its traits:
  template<>
  struct traits<MatrixFreeSparse> :  public Eigen::internal::traits<Eigen::SparseMatrix<double> >
  {};

  template<>
  struct evaluator_traits<MatrixFreeSparse> : public Eigen::internal::evaluator_traits<Eigen::SparseMatrix<double> >
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
 
   // Custom API:
  MatrixFreeSparse() : elementStiffnessMat(0), elementToNode(0) {}

  MatrixFreeSparse(StorageIndex _numElements, Eigen::Matrix<double, 300, 1> _elementStiffnessMat, std::vector<std::array<uint64_t, 8>> _elementToNode) 
                  : numElements(_numElements),
                    elementStiffnessMat(_elementStiffnessMat),
                    elementToNode(_elementToNode) {}

private:
  const Eigen::Matrix<double, 300, 1> elementStiffnessMat;
	std::vector<std::array<uint64_t, 8>> elementToNode;
  StorageIndex numElements;

};

#endif // __MATRIX_FREE_SPARSE_HPP__