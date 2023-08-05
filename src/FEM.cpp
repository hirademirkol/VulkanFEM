#include "out.hpp"
#include "FEM.hpp"
#include "IncompleteCholeskyPreconditioner.hpp"

//TODO: pass the used vertices out for applyBoundaryConditions to map global nodes to matrix nodes
#ifdef MATRIX_FREE
template<typename scalar>
MatrixFreeSparse assembleSystemMatrix(int* voxelModel, Vec3i voxelGridDimensions, double elementStiffness[24][24], const std::set<uint64_t>& fixedNodes)
#else
template<typename scalar>
Eigen::SparseMatrix<scalar> assembleSystemMatrix(int* voxelModel, Vec3i voxelGridDimensions, scalar elementStiffness[24][24], const std::set<uint64_t>& fixedNodes)
#endif
{

	std::cout << "Assembling the system matrix" << std::endl;
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

	uint64_t index = 1, expectedIndex = 0;
	for (auto line : usedNodes)
	{
		if(fixedNodes.find(expectedIndex) == fixedNodes.end())
			freeNodes[line.first] = index++;
		
		expectedIndex++;
	}

	int size = freeNodes.size() * 3;

	std::vector<std::array<uint64_t, 8>> elementToGlobal;
	for (auto element : usedElements)
	{
		std::array<uint64_t, 8> verts;
		FOR3(vert, Vec3i(0), Vec3i(2))
		{
			verts[Linearize(vert, Vec3i(2))] = freeNodes[Linearize(element + vert, vertexGridDimensions)] - 1;
		}
		elementToGlobal.push_back(verts);
	}
	// saveMatrix<uint64_t, 8>(elementToGlobal);

#ifdef MATRIX_FREE
	Eigen::Matrix<scalar, 300, 1> elementStiffnessMatrix;
	for(int i = 0; i < 24; i++)
	for(int j = 0; j < i + 1; i++)
	{
		elementStiffnessMatrix((j*(47-j))/2 + i + 1) = elementStiffness[i][j];
	}


	MatrixFreeSparse systemMatrix(size, elementStiffnessMatrix, elementToGlobal);
#else

	Eigen::SparseMatrix<scalar> systemMatrix(size, size);
	std::vector<Eigen::Triplet<scalar>> triplets;
	// std::vector<scalar> ind,jnd,v;
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

				if(iGlobal == -1 || jGlobal == -1)
					continue;

				for(int c1 = 0; c1 < 3; c1++)
				{
					for(int c2 = 0; c2 < 3; c2++)
					{
						triplets.push_back(Eigen::Triplet<scalar>(iGlobal*3 + c1, jGlobal*3 + c2, get_symmetric(elementStiffness,i*3 + c1,j*3 + c2)));
						// ind.push_back(iGlobal*3 + c1);
						// jnd.push_back(jGlobal*3 + c2);
						// v.push_back(get_symmetric(elementStiffness,i*3 + c1,j*3 + c2));
					}
				}
			}
		}
	}

	// saveVector(ind, "i");
	// saveVector(jnd, "j");
	// saveVector(v, "v");

	systemMatrix.setFromTriplets(triplets.begin(), triplets.end());
	systemMatrix.makeCompressed();
	std::cout << "System Matrix : Size:" << systemMatrix.rows() << "x" << systemMatrix.cols() << ", Non-Zero:" << systemMatrix.nonZeros() << ", " << ((float)systemMatrix.nonZeros()/systemMatrix.rows()) << " full per row" << std::endl;
#endif
	return systemMatrix;
}

#ifdef MATRIX_FREE
template MatrixFreeSparse assembleSystemMatrix<double>(int* voxelModel, Vec3i voxelGridDimensions, double elementStiffness[24][24], const std::set<uint64_t>& fixedNodes);
#else
template Eigen::SparseMatrix<double> assembleSystemMatrix<double>(int* voxelModel, Vec3i voxelGridDimensions, double elementStiffness[24][24], const std::set<uint64_t>& fixedNodes);
template Eigen::SparseMatrix<float> assembleSystemMatrix<float>(int* voxelModel, Vec3i voxelGridDimensions, float elementStiffness[24][24], const std::set<uint64_t>& fixedNodes);
#endif

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

template void applyBoundaryConditions<double>(std::vector<double>& f, std::map<uint64_t, Vec3<double>>& loadedNodes);
template void applyBoundaryConditions<float>(std::vector<float>& f, std::map<uint64_t, Vec3<float>>& loadedNodes);

#ifdef MATRIX_FREE
template <typename scalar>
void solveWithCG(const MatrixFreeSparse& A, const std::vector<scalar>& b, std::vector<scalar>& x)
#else
template <typename scalar>
void solveWithCG(const Eigen::SparseMatrix<scalar>& A, const std::vector<scalar>& b, std::vector<scalar>& x)
#endif
{
	Eigen::Matrix<scalar, Eigen::Dynamic, 1> b_eig(b.size());
	Eigen::Matrix<scalar, Eigen::Dynamic, 1> x_eig(x.size());

	memcpy(b_eig.data(), b.data(), b.size()*sizeof(scalar)); 

	std::cout << "Setting up the CG solver with Incomplete Cholesky Preconditioner" << std::endl;

#ifdef MATRIX_FREE
	Eigen::ConjugateGradient<MatrixFreeSparse, Eigen::Lower | Eigen::Upper, Eigen::IdentityPreconditioner> solver(A);
#else
	Eigen::ConjugateGradient<Eigen::SparseMatrix<scalar>, Eigen::Lower, Eigen::IncompleteCholeskyPreconditioner<scalar>> solver(A);
#endif
	
#ifdef MAX_ITER
	solver.setMaxIterations(MAX_ITER);
#endif

	std::cout << "Starting the solver" << std::endl;
	x_eig = solver.solve(b_eig);

	std::cout << "#iterations:     " << solver.iterations() << std::endl;
	std::cout << "Estimated error: " << solver.error()      << std::endl;


	std::memcpy(x.data(), x_eig.data(), x.size()*sizeof(scalar)); 
}

#ifdef MATRIX_FREE
template void solveWithCG<double>(const MatrixFreeSparse& A, const std::vector<double>& b, std::vector<double>& x);
#else
template void solveWithCG<double>(const Eigen::SparseMatrix<double>& A, const std::vector<double>& b, std::vector<double>& x);
template void solveWithCG<float>(const Eigen::SparseMatrix<float>& A, const std::vector<float>& b, std::vector<float>& x);
#endif