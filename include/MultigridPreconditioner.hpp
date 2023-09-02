#ifndef __MULTIGRID_HPP__
#define __MULTIGRID_HPP__

#include "FEM.hpp"
#include <Eigen/SparseCholesky>

namespace Eigen 
{
  class MultigridPreconditioner
  {
      typedef Matrix<double, Dynamic, 1> Vector;
    public:
      typedef typename Vector::StorageIndex StorageIndex;
      enum {
        ColsAtCompileTime = Dynamic,
        MaxColsAtCompileTime = Dynamic
      };

      MultigridPreconditioner() : m_isInitialized(false) {}

      explicit MultigridPreconditioner(const MatrixFreeSparse& mat) : ldltSolver(mat.Kc)
      {
        compute(mat);
      }

      EIGEN_CONSTEXPR Index rows() const EIGEN_NOEXCEPT { return matrix->rows(); }
      EIGEN_CONSTEXPR Index cols() const EIGEN_NOEXCEPT { return matrix->cols(); }

      MultigridPreconditioner& analyzePattern(const MatrixFreeSparse& )
      {
        return *this;
      }

      MultigridPreconditioner& factorize(const MatrixFreeSparse& mat)
      {
        matrix = &mat;
        numLevels = mat.numLevels;
        ldltSolver.compute(mat.Kc);

        m_isInitialized = true;
        return *this;
      }

      MultigridPreconditioner& compute(const MatrixFreeSparse& mat)
      {
        return factorize(mat);
      }

      /** \internal */
      void _solve_impl(const Vector& b, Vector& x) const
      {
        std::vector<Vector> rs(numLevels);
        std::vector<Vector> dxs(numLevels);

        rs[0] = b;

        for(int i = 0; i < numLevels - 1; i++)
        {
          rs[i+1] = matrix->restrictionMatrices[i] * rs[i];
        }

        dxs[numLevels - 1] = ldltSolver.solve(rs[numLevels - 1]);

        for(int i = numLevels - 1; i > 0; i--)
        {
          dxs[i-1] = matrix->interpolationMatrices[i-1] * dxs[i];
        }

        x = b - *matrix * dxs[0];
      }

      inline const Solve<MultigridPreconditioner, Vector>
      solve(const MatrixBase<Vector>& b) const
      {
        eigen_assert(m_isInitialized && "MultigridPreconditioner is not initialized.");
        eigen_assert(matrix->rows()==b.rows()
                  && "MultigridPreconditioner::solve(): invalid number of rows of the right hand side matrix b");
        return Solve<MultigridPreconditioner, Vector>(*this, b.derived());
      }

      ComputationInfo info() { return Success; }

    private:
      const MatrixFreeSparse* matrix;
      SimplicialLDLT<Eigen::SparseMatrix<double>> ldltSolver;
      int numLevels;

      bool m_isInitialized;
  };
}

#endif // __MULTIGRID_HPP__