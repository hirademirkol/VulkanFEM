#include "utils.hpp"

#include <utility>
#include <map>
#include <unordered_map>
#include <array>

#include <math.h>

// #define CSPARSE
#ifdef CSPARSE
#include <cs.h>
#endif

#define EIGEN
#ifdef EIGEN
#include <Eigen/Sparse>
#endif

#define MAX_ITER 200
#define index(i,j,n) (i)*(n)+(j)
#define get_symmetric(A,i,j) (i) <= (j) ? A[(i)][(j)] : A[(j)][(i)]

template <typename scalar> 
struct Mat{
	std::map<uint64_t, scalar> data;
	uint64_t rows;
	const scalar zero = (scalar)0.0;

	Mat<scalar>(uint64_t size)
	{
		rows = size;
	}

	const scalar& operator[](uint64_t index) const
	{
		if(auto a = data.find(index); a != data.end())
			return a->second;
		else
			return zero;
	}

	scalar& operator[](uint64_t index)
	{
		if(auto a = data.find(index); a != data.end())
			return a->second;
		else
			data[index] = 0; return data[index];
	}

	std::vector<uint64_t> indices()
	{
		std::vector<uint64_t> keys;

		for(auto it = data.begin(); it != data.end(); it++)
			keys.push_back(it->first);

		return keys;
	}

	std::vector<scalar> values()
	{
		std::vector<scalar> values;

		for(auto it = data.begin(); it != data.end(); it++)
			values.push_back(it->second);

		return values;
	}

	void clearZeros()
	{
		std::erase_if(data, [](const auto& item){
			auto const& [key, value] = item;
			return value == 0;
		});
	}
};

typedef Mat<double> Matd;
typedef Mat<float> Matf;

inline bool stationary(Vec3i element)
{
	//Return true if at staionary boundary
	return element.z == 0;
}

template<typename scalar>
// #ifdef CSPARSE
// cs* assembleSystemMatrix(int* voxelModel, Vec3i voxelGridDimensions, scalar elementStiffness[24][24])
// #else
#ifdef EIGEN
Eigen::SparseMatrix<scalar> assembleSystemMatrix(int* voxelModel, Vec3i voxelGridDimensions, scalar elementStiffness[24][24], const std::set<uint64_t>& fixedNodes)
//TODO: pass the used vertices out for applyBoundaryConditions to map global nodes to matrix nodes 
#else
Mat<scalar> assembleSystemMatrix(int* voxelModel, Vec3i voxelGridDimensions, scalar elementStiffness[24][24])
#endif
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
	// std::vector<int> usedNodes;
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


// #ifdef CSPARSE
// 	cs* systemMatrix = new cs();
// 	systemMatrix->m = size;
// 	systemMatrix->n = size;
#ifdef EIGEN
	Eigen::SparseMatrix<scalar> systemMatrix(size, size);
	std::vector<Eigen::Triplet<scalar>> triplets;
	// systemMatrix.reserve(75*size);
#else
	Mat<scalar> systemMatrix(size);
#endif
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
// #ifdef CSPARSE
// 						cs_entry(systemMatrix, iGlobal, jGlobal, get_symmetric(elementStiffness,i*3 + c1,j*3 + c2));
#ifdef EIGEN
						// systemMatrix.coeffRef(iGlobal*3 + c1, jGlobal*3 + c2) += get_symmetric(elementStiffness,i*3 + c1,j*3 + c2);
						triplets.push_back(Eigen::Triplet<scalar>(iGlobal*3 + c1, jGlobal*3 + c2, get_symmetric(elementStiffness,i*3 + c1,j*3 + c2)));
#else
						systemMatrix[index(iGlobal*3 + c1, jGlobal*3 + c2, size)] += get_symmetric(elementStiffness,i*3 + c1,j*3 + c2);
#endif
					}
				}
			}
		}
	}

// #ifdef CSPARSE
// 	cs_dupl(systemMatrix);
// 	std::cout << "System Matrix : Size:" << size << "x" << size << ", Non-Zero:" << systemMatrix->nz << ", " << ((float)systemMatrix->nz/size) << " full per row" << std::endl;
#ifdef EIGEN
	systemMatrix.setFromTriplets(triplets.begin(), triplets.end());
	systemMatrix.makeCompressed();
	std::cout << "System Matrix : Size:" << systemMatrix.rows() << "x" << systemMatrix.cols() << ", Non-Zero:" << systemMatrix.nonZeros() << ", " << ((float)systemMatrix.nonZeros()/systemMatrix.rows()) << " full per row" << std::endl;
#else
	std::cout << "System Matrix : Size:" << size << "x" << size << ", Non-Zero:" << systemMatrix.data.size() << ", " << ((float)systemMatrix.data.size()/size) << " full per row" << std::endl;
#endif
	return systemMatrix;
}

template <typename scalar>
#ifdef EIGEN
void applyBoundaryConditions(std::vector<scalar>& f, std::map<uint64_t, Vec3<scalar>>& loadedNodes)
{
	for(auto node : loadedNodes)
	{
		f[3*node.first] = node.second.x;
		f[3*node.first + 1] = node.second.y;
		f[3*node.first + 2] = node.second.z;
	}
}
#else
void applyBoundaryConditions(Mat<scalar>& systemMatrix, std::vector<scalar> f, const std::vector<uint64_t>& fixedNodes)
{
	for(auto iter = fixedNodes.cend()-1; iter > fixedNodes.cbegin(); iter--)
	{
		for(auto node : fixedNodes)
		{
			systemMatrix.data.erase(index(*iter, node, systemMatrix.rows));
		}
	}
#endif

// template void applyBoundaryConditions<double>(Matd& systemMatrix, std::vector<double> f, const std::vector<uint64_t>& fixedNodes);
// template void applyBoundaryConditions<float>(Matf& systemMatrix, std::vector<float> f, const std::vector<uint64_t>& fixedNodes);

template <typename scalar>
std::vector<scalar> multiply(const Mat<scalar>& A, const std::vector<scalar>& b)
{
	uint64_t size = b.size();
	std::vector<scalar> f;
	f.resize(size, 0);

	//#pragma omp parallel
	for(auto it = A.data.begin(); it != A.data.end(); it++)
	{
		uint64_t i = it->first / size;
		uint64_t j = it->first % size;
		
		f[i] += it->second * b[j];
	}

	return f;
}

template <typename scalar>
std::vector<scalar> multiply(const std::vector<scalar>& a, const std::vector<scalar>& b)
{
	uint64_t size = a.size();
	std::vector<scalar> result(a);

	//#pragma omp parallel for
	for(uint64_t i = 0; i < size; i++)
	{
		result[i] *= b[i];
	}

	return result;
}

template <typename scalar>
void add(std::vector<scalar>& addee, const std::vector<scalar>& addend)
{
	uint64_t size = addee.size();

	//#pragma omp parallel for
	for(uint64_t i = 0; i < size; i++)
	{
		addee[i] += addend[i];
	}
}

template <typename scalar>
void subtract(std::vector<scalar>& subtractee, const std::vector<scalar>& subtracted)
{
	uint64_t size = subtractee.size();

	//#pragma omp parallel for
	for(uint64_t i = 0; i < size; i++)
	{
		subtractee[i] -= subtracted[i];
	}
}

template <typename scalar>
scalar multiplyTranspose(const std::vector<scalar>& a, const std::vector<scalar>& b)
{
	uint64_t size = a.size();
	scalar result = 0;
	
	//#pragma omp parallel for
	for(int64_t i = 0; i < size; i++)
	{
		result += a[i]*b[i];
	}

	return result;
}

template <typename scalar>
std::vector<scalar> multiply(const scalar& a, const std::vector<scalar>& b)
{
	uint64_t size = b.size();
	std::vector<scalar> result(b);
	
	//#pragma omp parallel for
	for(int64_t i = 0; i < size; i++)
	{
		result[i] *= a;
	}

	return result;
}

template <typename scalar>
scalar two_norm(const std::vector<scalar>& b)
{
	scalar norm = (scalar)0;

	for(auto val : b)
	{
		norm += val*val;
	}

	return sqrt(norm);
}

#ifdef EIGEN
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
#endif

template <typename scalar>
#ifdef EIGEN
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
	// solver.setMaxIterations(MAX_ITER);
	x_eig = solver.solve(b_eig);

	std::cout << "#iterations:     " << solver.iterations() << std::endl;
	std::cout << "estimated error: " << solver.error()      << std::endl;

	std::memcpy(x.data(), x_eig.data(), x.size()*sizeof(scalar)); 
}
#else
void solveWithCG(const Mat<scalar>& A, const std::vector<scalar>& b, std::vector<scalar>& x, const std::vector<uint64_t> &fixedNodes)
{
	uint64_t size = x.size();
	std::vector<scalar> r(b);
	std::vector<scalar> r2(r);

	scalar lambda = (scalar)0;
	subtract(r, multiply(A, x));

	// smoothWithJacobi(A, r, x);
	std::vector<scalar> p(r);

	// subtract(r, multiply(A, x));

	for(int i = 0; i < MAX_ITER; i++)
	{
		scalar n1 = multiplyTranspose(r,r);
		lambda = n1 / multiplyTranspose(p, multiply(A, p));

		// lambda *= (scalar)0.01;
		add(x, multiply(lambda, p));

		subtract(r, multiply(lambda, multiply(A, p)));
		r2 = r;
		
		for(int i = 0; i < r.size(); i++)
		{
			auto val = A[index(i,i,A.rows)];
			if(val != 0.0)
				r[i] += r[i] / val;
		}

		scalar n2 = multiplyTranspose(r2, r2);

		add(r2, multiply(n2/n1, p));
		p = r2;

		r.swap(r2);
		// smoothWithJacobi(A, r, x);
		// subtract(r, multiply(A, x));

		std::cout << "iteration: " << i << ", " << "Norm: " << two_norm(r2) << std::endl;
	}
}
#endif

// template void solveWithCG<double>(const Matd& A, const std::vector<double>& b, std::vector<double>& x, const std::vector<uint64_t> &fixedNodes);
// template void solveWithCG<float>(const Matf& A, const std::vector<float>& b, std::vector<float>& x, const std::vector<uint64_t> &fixedNodes);