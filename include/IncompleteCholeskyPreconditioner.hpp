#ifndef __ICHOL_HPP__
#define __ICHOL_HPP__

#include "FEM.hpp"

namespace Eigen 
{
  // Incomplete Cholesky Preconditioner implementation, based on the basic preconditioner implementations from Eigen
  template <typename scalar>
  class IncompleteCholeskyPreconditioner
  {
      typedef Matrix<scalar, Dynamic,1> Vector;
    public:
      typedef typename Vector::StorageIndex StorageIndex;
      enum {
        ColsAtCompileTime = Dynamic,
        MaxColsAtCompileTime = Dynamic
      };

      IncompleteCholeskyPreconditioner() : m_isInitialized(false) {}

      template<typename MatType>
      explicit IncompleteCholeskyPreconditioner(const MatType& mat)
      {
        compute(mat);
      }

      EIGEN_CONSTEXPR Index rows() const EIGEN_NOEXCEPT { return m_LL.rows(); }
      EIGEN_CONSTEXPR Index cols() const EIGEN_NOEXCEPT { return m_LL.cols(); }

      template<typename MatType>
      IncompleteCholeskyPreconditioner& analyzePattern(const MatType& )
      {
        return *this;
      }

      template<typename MatType>
      IncompleteCholeskyPreconditioner& factorize(const MatType& mat)
      {
        // Initialize the IncompleteCholesky solver from the matrix and write down the factorization matrix
        IncompleteCholesky<scalar> iChol(mat);
        m_LL = iChol.matrixL() * iChol.matrixL().transpose();
        m_LL.cwiseSqrt();

        m_isInitialized = true;
        return *this;
      }

      template<typename MatType>
      IncompleteCholeskyPreconditioner& compute(const MatType& mat)
      {
        return factorize(mat);
      }

      /** \internal */
      template<typename Rhs, typename Dest>
      void _solve_impl(const Rhs& b, Dest& x) const
      {
        // Return L*L'*b for preconditioning
        x = m_LL * b;
      }

      template<typename Rhs> inline const Solve<IncompleteCholeskyPreconditioner, Rhs>
      solve(const MatrixBase<Rhs>& b) const
      {
        eigen_assert(m_isInitialized && "IncompleteCholeskyPreconditioner is not initialized.");
        eigen_assert(m_LL.rows()==b.rows()
                  && "IncompleteCholeskyPreconditioner::solve(): invalid number of rows of the right hand side matrix b");
        return Solve<IncompleteCholeskyPreconditioner, Rhs>(*this, b.derived());
      }

      ComputationInfo info() { return Success; }

    private:
      SparseMatrix<scalar> m_LL;
      bool m_isInitialized;
  };
}

#endif // __ICHOL_HPP__