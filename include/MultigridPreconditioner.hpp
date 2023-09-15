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
		enum
		{
			ColsAtCompileTime = Dynamic,
			MaxColsAtCompileTime = Dynamic
		};

		MultigridPreconditioner() : m_isInitialized(false) {}

		explicit MultigridPreconditioner(const MatrixFreeSparse &mat) : ldltSolver(mat.Kc)
		{
			compute(mat);
		}

		EIGEN_CONSTEXPR Index rows() const EIGEN_NOEXCEPT { return matrix->rows(); }
		EIGEN_CONSTEXPR Index cols() const EIGEN_NOEXCEPT { return matrix->cols(); }

		MultigridPreconditioner &analyzePattern(const MatrixFreeSparse &)
		{
			return *this;
		}

		MultigridPreconditioner &factorize(const MatrixFreeSparse &mat)
		{
			matrix = &mat;
			numLevels = mat.numLevels;
			ldltSolver.compute(mat.Kc);

			m_isInitialized = true;
			return *this;
		}

		MultigridPreconditioner &compute(const MatrixFreeSparse &mat)
		{
			return factorize(mat);
		}

		inline void smooth(Vector &r, Vector &x, const Vector &invDiag) const
		{
			Vector dx = r.cwiseProduct(invDiag);
			x += dx;
			// r -= *matrix*dx;
		}

		/** \internal */
		void _solve_impl(const Vector &b, Vector &x) const
		{
			std::vector<Vector> xs(numLevels);
			std::vector<Vector> rs(numLevels);

			rs[0] = b;
			xs[0] = Vector::Zero(rs[0].rows());

			for (int i = 0; i < numLevels - 1; i++)
			{
				smooth(rs[i], xs[i], matrix->invDiagKOnLevels[i]);
				rs[i + 1] = matrix->restrictionMatrices[i] * rs[i];
				xs[i + 1] = Vector::Zero(rs[i + 1].rows());
			}

			Vector solution = ldltSolver.solve(rs[numLevels - 1](matrix->coarseFreeDoFs));
			xs[numLevels - 1](matrix->coarseFreeDoFs) = solution;

			for (int i = numLevels - 1; i > 0; i--)
			{
				xs[i - 1] += matrix->interpolationMatrices[i - 1] * xs[i];
				smooth(rs[i - 1], xs[i - 1], matrix->invDiagKOnLevels[i - 1]);
			}

			x = xs[0];

			x(3 * matrix->fixedNodes) = Eigen::ArrayXd::Zero(matrix->fixedNodes.size());
			x(3 * matrix->fixedNodes + 1) = Eigen::ArrayXd::Zero(matrix->fixedNodes.size());
			x(3 * matrix->fixedNodes + 2) = Eigen::ArrayXd::Zero(matrix->fixedNodes.size());
		}

		inline const Solve<MultigridPreconditioner, Vector>
		solve(const MatrixBase<Vector> &b) const
		{
			eigen_assert(m_isInitialized && "MultigridPreconditioner is not initialized.");
			eigen_assert(matrix->rows() == b.rows() && "MultigridPreconditioner::solve(): invalid number of rows of the right hand side matrix b");
			return Solve<MultigridPreconditioner, Vector>(*this, b.derived());
		}

		ComputationInfo info() { return Success; }

	private:
		const MatrixFreeSparse *matrix;
		SimplicialLDLT<Eigen::SparseMatrix<double>> ldltSolver;
		int numLevels;

		bool m_isInitialized;
	};
}

#endif // __MULTIGRID_HPP__