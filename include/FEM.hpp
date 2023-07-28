#include "utils.hpp"

#include <utility>
#include <map>
#include <unordered_map>
#include <array>

#include <math.h>

#include <Eigen/Sparse>

#define MAX_ITER 200
#define index(i,j,n) (i)*(n)+(j)
#define get_symmetric(A,i,j) (i) <= (j) ? A[(i)][(j)] : A[(j)][(i)]

//TODO: pass the used vertices out for applyBoundaryConditions to map global nodes to matrix nodes 
template<typename scalar>
Eigen::SparseMatrix<scalar> assembleSystemMatrix(int* voxelModel, Vec3i voxelGridDimensions, scalar elementStiffness[24][24], const std::set<uint64_t>& fixedNodes)
{
	Vec3i vertexGridDimensions = voxelGridDimensions + Vec3i(1);

	// int levels = 3;
	std::vector<Vec3i> usedElements;
	// std::vector<Vec3i> usedElementsInLevel[levels];


	Vec3i element;
	FOR3(element,  Vec3i(0), voxelGridDimensions)
	{
		if(voxelModel[Linearize(element, voxelGridDimensions)] == 1)
		{
			usedElements.push_back(element);
		}
	}

	Vec3i::MAX = std::max(vertexGridDimensions.x,vertexGridDimensions.y);
	Vec3i::MAX = std::max(vertexGridDimensions.z, Vec3i::MAX);

	std::map<uint64_t, uint64_t> usedNodes, freeNodes;
	Vec3i vert;

	for(auto element : usedElements)
	{
		FOR3(vert, Vec3i(0), Vec3i(2))
		{
			auto val = Linearize(element + vert, vertexGridDimensions);
			if(!usedNodes.contains(val))
			{
				usedNodes[val] = val;
			}
		}
	}

	uint64_t index = 0, expectedIndex = 0;
	for (auto line : usedNodes)
	{
		if(fixedNodes.find(3*expectedIndex    ) == fixedNodes.end() ||
		   fixedNodes.find(3*expectedIndex + 1) == fixedNodes.end() ||
		   fixedNodes.find(3*expectedIndex + 2) == fixedNodes.end())
			freeNodes[line.first] = index++;
		else
		{ expectedIndex++; }
	}

	std::vector<std::array<uint64_t, 8>> elementToGlobal;
	for (auto element : usedElements)
	{
		std::array<uint64_t, 8> verts;
		FOR3(vert, Vec3i(0), Vec3i(2))
		{
			verts[Linearize(vert, Vec3i(2))] = freeNodes[Linearize(element + vert, vertexGridDimensions)];
		}
		elementToGlobal.push_back(verts);
	}

	// saveMatrix<uint64_t, 8>(elementToGlobal);

	int size = freeNodes.size() * 3;

	Eigen::SparseMatrix<scalar> systemMatrix(size, size);
	std::vector<Eigen::Triplet<scalar>> triplets;

	for(auto line : elementToGlobal)
	{
		FOR3(vert, Vec3i(0), Vec3i(2))
		{
			Vec3i vert2; 
			FOR3(vert2, Vec3i(0), Vec3i(2))
			{
				auto i = Linearize(vert,  Vec3i(2));
				auto j = Linearize(vert2, Vec3i(2));

				auto iGlobal = line[i];
				auto jGlobal = line[j];

				for(int c1 = 0; c1 < 3; c1++)
				{
					for(int c2 = 0; c2 < 3; c2++)
					{
						triplets.push_back(Eigen::Triplet<scalar>(iGlobal*3 + c1, jGlobal*3 + c2, get_symmetric(elementStiffness,i*3 + c1,j*3 + c2)));
					}
				}
			}
		}
	}

	systemMatrix.setFromTriplets(triplets.begin(), triplets.end());
	systemMatrix.makeCompressed();
	std::cout << "System Matrix : Size:" << systemMatrix.rows() << "x" << systemMatrix.cols() << ", Non-Zero:" << systemMatrix.nonZeros() << ", " << ((float)systemMatrix.nonZeros()/systemMatrix.rows()) << " full per row" << std::endl;

	return systemMatrix;
}

template <typename scalar>
void applyBoundaryConditions(std::vector<scalar>& f, std::map<uint64_t, Vec3<scalar>>& loadedNodes)
{
	for(auto node : loadedNodes)
	{
		f[3*node.first] = node.second.x;
		f[3*node.first + 1] = node.second.y;
		f[3*node.first + 2] = node.second.z;
	}
}

#undef index
namespace Eigen {

	template <typename scalar>
	class MultigridPreconditioner
	{
		typedef Matrix<scalar, Dynamic, 1> Vector;

		public:
		typedef typename Vector::StorageIndex StorageIndex;
		
		enum {
		ColsAtCompileTime = Dynamic,
		MaxColsAtCompileTime = Dynamic
		};

		MultigridPreconditioner() : m_isInitialized(false) {}

		explicit MultigridPreconditioner(const SparseMatrix<scalar>& mat) : m_invdiag(mat.cols())
		{
			compute(mat);
		}

		static void SetMatrices(const SparseMatrix<scalar> mat[2],
										 const SparseMatrix<scalar> restriction[],
										 const SparseMatrix<scalar> interpolation[],
										 int numLevels)
		{
			s_matrices = mat;
			s_restrictionMatrices = restriction;
			s_interpolationMatrices = interpolation;
			s_numLevels = numLevels;
		}


		EIGEN_CONSTEXPR Index rows() const EIGEN_NOEXCEPT { return m_invdiag.size(); }
		EIGEN_CONSTEXPR Index cols() const EIGEN_NOEXCEPT { return m_invdiag.size(); }

		MultigridPreconditioner& analyzePattern(const SparseMatrix<scalar>& )
		{
			return *this;
		}

		MultigridPreconditioner& factorize(const SparseMatrix<scalar>& mat)
		{
			m_invdiag.resize(mat.cols());

			for(int j=0; j<mat.outerSize(); ++j)
			{
				typename SparseMatrix<scalar>::InnerIterator it(mat,j);
				while(it && it.index()!=j) ++it;
				if(it && it.index()==j && it.value()!=scalar(0))
				m_invdiag(j) = scalar(1)/it.value();
				else
				m_invdiag(j) = scalar(1);
			}



			m_isInitialized = true;
			return *this;
		}

		MultigridPreconditioner& compute(const SparseMatrix<scalar>& mat)
		{
			return factorize(mat);
		}

		/** \internal */
		void _solve_impl(const Vector& b, Vector& x) const
		{
			// Vector r = b - s_matrices[0] * x;

			// r *= m_invdiag;

			// for(int i = 0; i < m_numLevels)
			// {
			// 	r *= s_restrictionMatrices[i+1];
			// }

			// SimplicialLLT<SparseMatrix<scalar>> subSolver(m_matrix[2]);

			// subSolver.solve(r);

			// for(int i = 0; i < m_numLevels)
			// {
			// 	r *= s_interpolationMatrices[i+1];
			// }

			// x += r;
		}

		template<typename Rhs> inline const Solve<MultigridPreconditioner, Rhs>
		solve(const MatrixBase<Rhs>& b) const
		{
			eigen_assert(m_isInitialized && "MultigridPreconditioner is not initialized.");
			eigen_assert(m_invdiag.size()==b.rows()
						&& "MultigridPreconditioner::solve(): invalid number of rows of the right hand side matrix b");
			return Solve<MultigridPreconditioner, Rhs>(*this, b.derived());
		}

		ComputationInfo info() { return Success; }

		protected:
		Vector m_invdiag;
		bool m_isInitialized;

		private:
		static SparseMatrix<scalar>* s_matrices;
		static SparseMatrix<scalar>* s_restrictionMatrices;
		static SparseMatrix<scalar>* s_interpolationMatrices;
		static int s_numLevels;
	};
}

#define index(i,j,n) (i)*(n)+(j)

template <typename scalar>
void solveWithCG(const Eigen::SparseMatrix<scalar>& A, const std::vector<scalar>& b, std::vector<scalar>& x)
{
	Eigen::Matrix<scalar, Eigen::Dynamic, 1> b_eig(b.size());
	Eigen::Matrix<scalar, Eigen::Dynamic, 1> x_eig(x.size());

	memcpy(b_eig.data(), b.data(), b.size()*sizeof(scalar)); 

	// int levels = 3; 
	// Eigen::SparseMatrix<scalar> matrices[2];
	// Eigen::SparseMatrix<scalar> restrictionMatrices[levels];
	// Eigen::SparseMatrix<scalar> interpolationMatrices[levels];

	// matrices[0] = A;
	// restrictionMatrices


	// for (int i = 0; i < levels; i++)
	// {
		
	// }

	// Eigen::MultigridPreconditioner<scalar>::SetMatrices(matrices, restrictionMatrices, interpolationMatrices, levels);
	Eigen::ConjugateGradient<Eigen::SparseMatrix<scalar>> solver(A);

	solver.setMaxIterations((int)(MAX_ITER/10));

	x_eig = solver.solve(b_eig);

	std::cout << "#iterations:     " << (int)(MAX_ITER/10) << std::endl;
	std::cout << "estimated error: " << solver.error()      << std::endl;
	
	for (int i = (int)(MAX_ITER/10); i < MAX_ITER; i+=(int)(MAX_ITER/10))
	{
		x_eig = solver.solveWithGuess(b_eig, x_eig);
		std::cout << "#iterations:     " << i + (int)(MAX_ITER/10) << std::endl;
		std::cout << "estimated error: " << solver.error()      << std::endl;
	}

	std::cout << "Final #iterations:     " << MAX_ITER << std::endl;
	std::cout << "Final estimated error: " << solver.error()      << std::endl;

	std::memcpy(x.data(), x_eig.data(), x.size()*sizeof(scalar)); 
}