#ifndef __MATRIX_FREE_SPARSE_HPP__
#define __MATRIX_FREE_SPARSE_HPP__

#include <Eigen/Sparse>

// Declaration of Matrix-Free form matrix
template <typename scalar>
class MatrixFreeSparse;

template<>
class MatrixFreeSparse<double>;

// template<>
// class MatrixFreeSparse<float>;

using Eigen::SparseMatrix;

// Inherit the traits of Sparse Matrix from Eigen for usage in solvers
namespace Eigen {
namespace internal {
  template<>
  struct traits<MatrixFreeSparse<double>> :  public Eigen::internal::traits<Eigen::SparseMatrix<double>> {};

  template<>
  struct traits<MatrixFreeSparse<float>> :  public Eigen::internal::traits<Eigen::SparseMatrix<float>> {};
  }
}

// Definition of Matrix-Free form matrix
// Based on the example by Eigen: https://eigen.tuxfamily.org/dox/group__MatrixfreeSolverExample.html
template<>
class MatrixFreeSparse<double> : public Eigen::EigenBase<MatrixFreeSparse<double>> {
public:
  typedef double Scalar;
  typedef double RealScalar;
  typedef int StorageIndex;
  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic,
    IsRowMajor = false,
  };
 
  template<typename Rhs>
  Eigen::Product<MatrixFreeSparse<double>,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
    return Eigen::Product<MatrixFreeSparse<double>,Rhs,Eigen::AliasFreeProduct>(*this, x.derived());
  }

  MatrixFreeSparse(){}

  MatrixFreeSparse(StorageIndex _numElements, Eigen::Matrix<double, 24, 24> _elementStiffnessMat, Eigen::Array<int, Eigen::Dynamic, 4> _elementToNode, Eigen::ArrayXi _fixedNodes) 
                  : numElements(_numElements),
                    elementStiffnessMat(_elementStiffnessMat),
                    elementToNode(_elementToNode),
                    fixedNodes(_fixedNodes) {}


  EIGEN_CONSTEXPR Index rows() const EIGEN_NOEXCEPT { return numElements; }
  EIGEN_CONSTEXPR Index cols() const EIGEN_NOEXCEPT { return numElements; }

  //Necessary data for Matrix-free Solver
  StorageIndex numElements;
  const Eigen::Matrix<double, 24, 24> elementStiffnessMat;
	Eigen::Array<int, Eigen::Dynamic, 4> elementToNode;
  Eigen::ArrayXi fixedNodes;

  //Necessary data for Multigrid Preconditioner
  void PrepareMultigrid(int _numLevels,
                        int _skipLevels,
                        std::vector<Eigen::Array<int, Eigen::Dynamic, 4>> _elementToNodeMatrices,
                        std::vector<Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic>> _restrictionMappings,
                        std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> _restrictionCoefficients,
                        std::vector<Eigen::VectorXd> _invDiagKOnLevels,
                        Eigen::SparseMatrix<double> _Kc,
                        Eigen::ArrayXi _coarseFreeDoFs)
  {
    numLevels = _numLevels;
    skipLevels = _skipLevels;
    elementToNodeMatrices = _elementToNodeMatrices;
    restrictionMappings = _restrictionMappings;
    restrictionCoefficients = _restrictionCoefficients;
    invDiagKOnLevels = _invDiagKOnLevels;
    Kc = _Kc;
    coarseFreeDoFs = _coarseFreeDoFs;
  }

  int numLevels;
  int skipLevels;
  std::vector<Eigen::Array<int, Eigen::Dynamic, 4>> elementToNodeMatrices;
  std::vector<Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic>> restrictionMappings;
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> restrictionCoefficients;
  std::vector<Eigen::VectorXd> invDiagKOnLevels;
  Eigen::SparseMatrix<double> Kc;
  Eigen::ArrayXi coarseFreeDoFs;
};

// template<>
// class MatrixFreeSparse<float> : public Eigen::EigenBase<MatrixFreeSparse<float>> {
// public:
//   typedef float Scalar;
//   typedef float RealScalar;
//   typedef int StorageIndex;
//   enum {
//     ColsAtCompileTime = Eigen::Dynamic,
//     MaxColsAtCompileTime = Eigen::Dynamic,
//     IsRowMajor = false,
//   };
 
//   template<typename Rhs>
//   Eigen::Product<MatrixFreeSparse<float>,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
//     return Eigen::Product<MatrixFreeSparse<float>,Rhs,Eigen::AliasFreeProduct>(*this, x.derived());
//   }

//   MatrixFreeSparse(){}

//   MatrixFreeSparse(StorageIndex _numElements, Eigen::Matrix<float, 24, 24> _elementStiffnessMat, Eigen::Array<int, Eigen::Dynamic, 8> _elementToNode, Eigen::ArrayXi _fixedNodes) 
//                   : numElements(_numElements),
//                     elementStiffnessMat(_elementStiffnessMat),
//                     elementToNode(_elementToNode),
//                     fixedNodes(_fixedNodes) {}


//   EIGEN_CONSTEXPR Index rows() const EIGEN_NOEXCEPT { return numElements; }
//   EIGEN_CONSTEXPR Index cols() const EIGEN_NOEXCEPT { return numElements; }

//   const Eigen::Matrix<float, 24, 24> elementStiffnessMat;
// 	Eigen::Array<int, Eigen::Dynamic, 8> elementToNode;
//   Eigen::ArrayXi fixedNodes;

//   //Necessary structs for Multigrid 
//   void PrepareMultigrid(int _numLevels,
//                         std::vector<Eigen::Array<int, Eigen::Dynamic, 8>> _elementToNodeMatrices,
//                         std::vector<Eigen::Array<int, Eigen::Dynamic, 27>> _restrictionMappings,
//                         std::vector<Eigen::Matrix<float, Eigen::Dynamic, 1>> _restrictionCoefficients,
//                         std::vector<Eigen::VectorXd> _invDiagKOnLevels,
//                         Eigen::SparseMatrix<float> _Kc,
//                         Eigen::ArrayXi _coarseFreeDoFs)
//   {
//     numLevels = _numLevels;
//     elementToNodeMatrices = _elementToNodeMatrices;
//     restrictionMappings = _restrictionMappings;
//     restrictionCoefficients = _restrictionCoefficients;
//     invDiagKOnLevels = _invDiagKOnLevels;
//     Kc = _Kc;
//     coarseFreeDoFs = _coarseFreeDoFs;
//   }

//   int numLevels;
//   std::vector<Eigen::Array<int, Eigen::Dynamic, 8>> elementToNodeMatrices;
//   std::vector<Eigen::Array<int, Eigen::Dynamic, 27>> restrictionMappings;
//   std::vector<Eigen::Matrix<float, Eigen::Dynamic, 1>> restrictionCoefficients;
//   std::vector<Eigen::VectorXd> invDiagKOnLevels;
//   Eigen::SparseMatrix<float> Kc;
//   Eigen::ArrayXi coarseFreeDoFs;
  
//   StorageIndex numElements;
// };


namespace Eigen {
namespace internal {

  // template<typename Rhs>
  // struct generic_product_impl<MatrixFreeSparse<float>, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
  // : generic_product_impl_base<MatrixFreeSparse<float>,Rhs,generic_product_impl<MatrixFreeSparse<float>,Rhs> >
  // {
  //   typedef typename Product<MatrixFreeSparse<float>,Rhs>::Scalar Scalar;
 
  //   template<typename Dest>
  //   static void scaleAndAddTo(Dest& dst, const MatrixFreeSparse<float>& lhs, const Rhs& rhs, const Scalar& alpha)
  //   {
  //     const Array<int, 1, 24> c   {0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2};
  //     const Array<int, 1, 24> xInd{0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7};

  //     for(auto line : lhs.elementToNode.rowwise())
  //     {
  //       Array<int, 1, 24> xs = 3 * line(xInd) + c;
  //       dst(xs) += alpha * lhs.elementStiffnessMat * rhs(xs);
  //     }

  //     dst(3 * lhs.fixedNodes    ) = Eigen::ArrayXf::Zero(lhs.fixedNodes.size());
  //     dst(3 * lhs.fixedNodes + 1) = Eigen::ArrayXf::Zero(lhs.fixedNodes.size());
  //     dst(3 * lhs.fixedNodes + 2) = Eigen::ArrayXf::Zero(lhs.fixedNodes.size());
  //   };
  // };

  // Implementation of a*A*b for iterative solvers,
  // Where  a -> alpha (scalar),
  //        A -> lhs (matrix),
  // and    b -> rhs (vector).
  // Based on the example by Eigen: https://eigen.tuxfamily.org/dox/group__MatrixfreeSolverExample.html
  template<typename Rhs>
  struct generic_product_impl<MatrixFreeSparse<double>, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
  : generic_product_impl_base<MatrixFreeSparse<double>,Rhs,generic_product_impl<MatrixFreeSparse<double>,Rhs> >
  {
    typedef typename Product<MatrixFreeSparse<double>,Rhs>::Scalar Scalar;
 
    template<typename Dest>
    static void scaleAndAddTo(Dest& dst, const MatrixFreeSparse<double>& lhs, const Rhs& rhs, const Scalar& alpha)
    {
      const Array<int, 1, 24> c   {0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5};
      const Array<int, 1, 24> xInd{0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3};

      // Loop over each element to add its contribution to the multiplication
      for(auto line : lhs.elementToNode.rowwise())
      {
        Array<int, 1, 24> xs = 3 * line(xInd) + c;
        dst(xs) += alpha * lhs.elementStiffnessMat * rhs(xs);
      }

      // Zero the fixed nodes back
      dst(3 * lhs.fixedNodes    ) = Eigen::ArrayXd::Zero(lhs.fixedNodes.size());
      dst(3 * lhs.fixedNodes + 1) = Eigen::ArrayXd::Zero(lhs.fixedNodes.size());
      dst(3 * lhs.fixedNodes + 2) = Eigen::ArrayXd::Zero(lhs.fixedNodes.size());
    };
  };
}
}

#endif // __MATRIX_FREE_SPARSE_HPP__