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

			restrictionOperator << 	0, 0.5000, 1.0000, 0, 0.2500, 0.5000, 0, 0, 0, 0, 0.2500, 0.5000, 0, 0.1250, 0.2500, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
									0, 0, 0, 0, 0.2500, 0.5000, 0, 0.5000, 1.0000, 0, 0, 0, 0, 0.1250, 0.2500, 0, 0.2500, 0.5000, 0, 0, 0, 0, 0, 0, 0, 0, 0,
									0, 0, 0,0.5000, 0.2500, 0, 1.0000, 0.5000, 0, 0, 0, 0, 0.2500, 0.1250, 0, 0.5000, 0.2500, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
									1.0000, 0.5000, 0, 0.5000, 0.2500, 0, 0, 0, 0, 0.5000, 0.2500, 0, 0.2500, 0.1250, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
									0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2500, 0.5000, 0, 0.1250, 0.2500, 0, 0, 0, 0, 0.5000, 1.0000, 0, 0.2500, 0.5000, 0, 0, 0,
									0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1250, 0.2500, 0, 0.2500, 0.5000, 0, 0, 0, 0, 0.2500, 0.5000, 0, 0.5000, 1.0000,
									0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2500, 0.1250, 0, 0.5000, 0.2500, 0, 0, 0, 0, 0.5000, 0.2500, 0, 1.0000, 0.5000, 0,
									0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5000, 0.2500, 0, 0.2500, 0.1250, 0, 0, 0, 0, 1.0000, 0.5000, 0, 0.5000, 0.2500, 0, 0, 0, 0;



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
		}

		inline Vector restrict(const Vector &x, int level) const
		{
			// return matrix->restrictionMatrices[level] * x;

			Vector result = Vector::Zero(matrix->invDiagKOnLevels[level + 1].rows());
			Vector tempX = x.cwiseProduct(matrix->restrictionCoefficients[level]);
      		// const Array<int, 1, 24> c   {0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2};
      		// const Array<int, 1, 24> xInd{0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7};

			int num = 0;
			for(auto line : matrix->elementToNodeMatrices[level].rowwise())
			{
				// Array<int, 1, 24> xs = 3 * line(xInd) + c;
				Array<int, 27, 1> mapping = matrix->restrictionMappings[level].row(num);
				// Matrix<double, 8, 1> res = restrictionOperator * x(mapping);
				Matrix<double, 27, 1> mask = Matrix<double, 27, 1>::Ones(mapping.rows());
				
				for(int i = 0; i < 27; i++)
					if(mapping(i) == -1)
					{
						mask(i) = 0;
						mapping(i) = 0;
					}

				for(int c = 0; c < 3; c++)
				{
					// Matrix<double, 27, 1> temp = tempX(3 * mapping + c);
					// temp = temp.cwiseProduct(mask);
					result(3 * line + c) += restrictionOperator *  tempX(3 * mapping + c).cwiseProduct(mask);
				}

				num++;
			}

			return result;
		}

		inline Vector interpolate(const Vector &x, int level) const
		{
			return matrix->interpolationMatrices[level] * x;
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
				rs[i + 1] = restrict(rs[i], i);
				xs[i + 1] = Vector::Zero(rs[i + 1].rows());
			}

			Vector solution = ldltSolver.solve(rs[numLevels - 1](matrix->coarseFreeDoFs));
			xs[numLevels - 1](matrix->coarseFreeDoFs) = solution;

			for (int i = numLevels - 1; i > 0; i--)
			{
				xs[i - 1] += interpolate(xs[i], i - 1);
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

		Matrix<double, 8, 27> restrictionOperator;
	};
}

#endif // __MULTIGRID_HPP__