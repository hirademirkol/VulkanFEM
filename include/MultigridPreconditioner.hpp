#ifndef __MULTIGRID_HPP__
#define __MULTIGRID_HPP__

#include "FEM.hpp"
#include "RestrictionOperator.hpp"

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

		explicit MultigridPreconditioner(const MatrixFreeSparse<double> &mat) : ldltSolver(mat.Kc)
		{
			compute(mat);
		}

		EIGEN_CONSTEXPR Index rows() const EIGEN_NOEXCEPT { return matrix->rows(); }
		EIGEN_CONSTEXPR Index cols() const EIGEN_NOEXCEPT { return matrix->cols(); }

		MultigridPreconditioner &analyzePattern(const MatrixFreeSparse<double> &)
		{
			return *this;
		}

		MultigridPreconditioner &factorize(const MatrixFreeSparse<double> &mat)
		{
			matrix = &mat;
			numLevels = mat.numLevels;
			ldltSolver.compute(mat.Kc);

			restrictionOperator = Eigen::Map<Eigen::Matrix<double, 8, 27>>(const_cast<double*>(restrictionOperatorValues.data()), 8, 27).transpose();
			restrictionOperator4 = Eigen::Map<Eigen::Matrix<double, 8, 125>>(const_cast<double*>(restrictionOperator4Values.data()), 8, 125).transpose();

			m_isInitialized = true;
			return *this;
		}

		MultigridPreconditioner &compute(const MatrixFreeSparse<double> &mat)
		{
			return factorize(mat);
		}

		inline void smooth(Vector &r, Vector &x, const Vector &invDiag) const
		{
			Vector dx = r.cwiseProduct(invDiag);
			x += dx;
		}

		inline Vector restrict(const Vector &r, int level) const
		{
			Vector result = Vector::Zero(matrix->invDiagKOnLevels[level + 1].rows());
			Vector tempR = r.cwiseProduct(matrix->restrictionCoefficients[level]);

			int num = 0;
			for(auto line : matrix->elementToNodeMatrices[level].rowwise())
			{
				Array<int, Eigen::Dynamic, 1> mapping = matrix->restrictionMappings[level].row(num);
				Matrix<double, Eigen::Dynamic, 1> mask = Matrix<double, Eigen::Dynamic, 1>::Ones(mapping.rows());

				for(int i = 0; i < mapping.rows(); i++)
				{
					if(mapping(i) == -1)
					{
						mask(i) = 0;
						mapping(i) = 0;
					}
				}

				const Array<int, 1, 8> xInd{0, 0, 1, 1, 2, 2, 3, 3};
				const Array<int, 1, 8> offset{0, 1, 0, 1, 0, 1, 0, 1};

				for(int c = 0; c < 3; c++)
				{
					if(level == 0 && matrix->skipLevels == 1)
						result(3 * (line(xInd) + offset) + c) += tempR(3 * mapping + c).cwiseProduct(mask).transpose() * restrictionOperator4;
					else
						result(3 * (line(xInd) + offset) + c) += tempR(3 * mapping + c).cwiseProduct(mask).transpose() * restrictionOperator;
				}

				num++;
			}
			return result;
		}

		inline Vector interpolate(const Vector &x, int level) const
		{
			Vector temp = Vector::Zero(matrix->invDiagKOnLevels[level].rows() + 3);
			int overRow = matrix->invDiagKOnLevels[level].rows() / 3;

			int num = 0;
			for(auto line : matrix->elementToNodeMatrices[level].rowwise())
			{
				Array<int, Eigen::Dynamic, 1> mapping = matrix->restrictionMappings[level].row(num);

				for(int i = 0; i < mapping.rows(); i++)
				{
					if(mapping(i) == -1)
					{
						mapping(i) = overRow;
					}
				}

				const Array<int, 1, 8> xInd{0, 0, 1, 1, 2, 2, 3, 3};
				const Array<int, 1, 8> offset{0, 1, 0, 1, 0, 1, 0, 1};

				for(int c = 0; c < 3; c++)
				{
					if(level == 0 && matrix->skipLevels == 1)
						temp(3 * mapping + c) += restrictionOperator4 *  x(3 * (line(xInd) + offset) + c);
					else
						temp(3 * mapping + c) += restrictionOperator *  x(3 * (line(xInd) + offset) + c);
				}

				num++;
			}

			Vector result(temp.rows() - 3);
			std::copy(temp.begin(), temp.end() - 3, result.begin());

			return result.cwiseProduct(matrix->restrictionCoefficients[level]);
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
		const MatrixFreeSparse<double> *matrix;
		SimplicialLDLT<Eigen::SparseMatrix<double>> ldltSolver;
		int numLevels;

		bool m_isInitialized;

		Matrix<double, 27, 8> restrictionOperator;
		Matrix<double, 125, 8> restrictionOperator4;
	};
}

#endif // __MULTIGRID_HPP__