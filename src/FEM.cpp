#include "out.hpp"
#include "FEM.hpp"
#include "IncompleteCholeskyPreconditioner.hpp"


#ifdef MATRIX_FREE
#ifdef MULTIGRID
	#include "MultigridPreconditioner.hpp"
#endif
#endif

#include <chrono>

template<typename scalar>
inline Eigen::SparseMatrix<scalar> assembleK(int size, std::vector<std::array<uint64_t, 8>>& elementToGlobal, double elementStiffness[24][24], scalar multiplier = 1.0)
{
	Eigen::SparseMatrix<scalar> systemMatrix(size, size);
	Vec3i node, node2;

	std::vector<Eigen::Triplet<scalar>> triplets;
	for(auto line : elementToGlobal)
	{
		FOR3(node, Vec3i(0), Vec3i(2))
		{
			FOR3(node2, Vec3i(0), Vec3i(2))
			{
				auto i = Linearize(node,  Vec3i(2));
				auto j = Linearize(node2, Vec3i(2));

				auto iGlobal = line[i];
				auto jGlobal = line[j];

				if(iGlobal == -1 || jGlobal == -1)
					continue;

				for(int c1 = 0; c1 < 3; c1++)
				{
					for(int c2 = 0; c2 < 3; c2++)
					{
						triplets.push_back(Eigen::Triplet<scalar>(iGlobal*3 + c1, jGlobal*3 + c2, multiplier * (get_symmetric(elementStiffness,i*3 + c1,j*3 + c2))));
					}
				}
			}
		}
	}

	// saveMatrix(triplets, "K3");

	systemMatrix.setFromTriplets(triplets.begin(), triplets.end());
	systemMatrix.makeCompressed();

	return systemMatrix;
}

inline void EnlistUsedElements(int* voxelModel, Vec3i voxelGridDimensions, Vec3i nodeGridDimensions,
							   std::vector<Vec3i>& usedElements, std::map<uint64_t, uint64_t>& usedNodes)
{
	Vec3i element;
	FOR3(element,  Vec3i(0), voxelGridDimensions)
	{
		if(voxelModel[Linearize(element, voxelGridDimensions)] == 1)
		{
			usedElements.push_back(element);
		}
	}

	Vec3i::MAX = std::max(nodeGridDimensions.x, nodeGridDimensions.y);
	Vec3i::MAX = std::max(nodeGridDimensions.z, Vec3i::MAX);

	Vec3i node;

	for(auto element : usedElements)
	{
		FOR3(node, Vec3i(0), Vec3i(2))
		{
			auto val = Linearize(element + node, nodeGridDimensions);
			if(!usedNodes.contains(val))
			{
				usedNodes[val] = val;
			}
		}
	}
}

#ifdef MULTIGRID

template <typename scalar>
inline Eigen::VectorXd GetInverseDiagonal(int size, Eigen::Matrix<scalar, 24, 24> elementStiffnessMat, Eigen::Array<int, Eigen::Dynamic, 8> elementToNode)
{
	Eigen::VectorXd diag(size);

	const Eigen::Array<int, 1, 24> c{0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2};
	const Eigen::Array<int, 1, 24> xInd{0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7};

	for (auto line : elementToNode.rowwise())
	{	
		Eigen::Array<int, 1, 24> xs = 3 * line(xInd) + c;
		diag(xs) += elementStiffnessMat.diagonal();
	}
	Eigen::VectorXd invDiag = diag.cwiseInverse();

	return invDiag;
}

#endif

//TODO: pass the used nodes out for applyBoundaryConditions to map global nodes to matrix nodes
#ifdef MATRIX_FREE
template<typename scalar>
MatrixFreeSparse<scalar> assembleSystemMatrix(int* voxelModel, Vec3i voxelGridDimensions, double elementStiffness[24][24], const std::set<uint64_t>& fixedNodes)
#else
template<typename scalar>
Eigen::SparseMatrix<scalar> assembleSystemMatrix(int* voxelModel, Vec3i voxelGridDimensions, scalar elementStiffness[24][24], const std::set<uint64_t>& fixedNodes)
#endif
{
	std::cout << "Assembling the system matrix" << std::endl;
	Vec3i nodeGridDimensions = voxelGridDimensions + Vec3i(1);

	std::vector<Vec3i> usedElements;
	std::map<uint64_t, uint64_t> usedNodes;
	EnlistUsedElements(voxelModel, voxelGridDimensions, nodeGridDimensions, usedElements, usedNodes);

#ifdef MATRIX_FREE

	uint64_t index = 0;
	for (auto line : usedNodes)
	{
		usedNodes[line.first] = index++;
	}

	int size = usedNodes.size() * 3;

	int line = 0; Vec3i node;
	Eigen::Array<int, Eigen::Dynamic, 8> elementToGlobal(usedElements.size(), 8);
	for (auto element : usedElements)
	{
		FOR3(node, Vec3i(0), Vec3i(2))
		{
			elementToGlobal(line, Linearize(node, Vec3i(2))) = (int)usedNodes[Linearize(element + node, nodeGridDimensions)];
		}
		line++;
	}

	Eigen::Matrix<scalar, 24, 24> elementStiffnessMatrix;
	for(int i = 0; i < 24; i++)
	{
		elementStiffnessMatrix(i,i) = (scalar)elementStiffness[i][i];
		for(int j = 0; j < i ; j++)
		{
			elementStiffnessMatrix(i,j) = elementStiffness[j][i];
			elementStiffnessMatrix(j,i) = elementStiffness[j][i];
		}
	}

	Eigen::ArrayXi fixed(fixedNodes.size());

	index = 0;
	for(auto val : fixedNodes)
	{
		fixed[index++] = (int)val;
	}
	MatrixFreeSparse<scalar> systemMatrix(size, elementStiffnessMatrix, elementToGlobal, fixed);

	#ifdef MULTIGRID

	std::cout << "Preparing the Multigrid structure" << std::endl;

	//Prepare multigrid related matrices
	int numLevels = NUM_LEVELS;
	std::vector<Vec3i> levelDims;
	std::vector<std::vector<Vec3i>> levelElements;
	std::vector<std::map<uint64_t, uint64_t>> usedNodesInLevels;
	std::vector<std::set<uint64_t>> fixedNodesInLevels;
	std::vector<Eigen::VectorXd> invDiagKOnLevels;
  	std::vector<Eigen::Array<int, Eigen::Dynamic, 8>> elementToNodeMatrices;
  	std::vector<Eigen::Array<int, Eigen::Dynamic, 27>> restrictionMappings;
	std::vector<Eigen::Matrix<scalar, Eigen::Dynamic, 1>> restrictionCoefficients;

	levelDims.push_back(voxelGridDimensions);
	levelElements.push_back(usedElements);
	usedNodesInLevels.push_back(usedNodes);
	fixedNodesInLevels.push_back(fixedNodes);
	invDiagKOnLevels.push_back(GetInverseDiagonal(size, elementStiffnessMatrix, elementToGlobal));

	for(int i = 1; i < numLevels; i++)
	{
		Vec3i newLevelDims;
		newLevelDims.x = levelDims[i-1].x / 2 + (levelDims[i-1].x % 2);
		newLevelDims.y = levelDims[i-1].y / 2 + (levelDims[i-1].y % 2);
		newLevelDims.z = levelDims[i-1].z / 2 + (levelDims[i-1].z % 2);
		levelDims.push_back(newLevelDims);

		std::set<Vec3i> elements, nodes;
		for(auto finerElement: levelElements[i-1])
		{
			int x = finerElement.x / 2;
			int y = finerElement.y / 2;
			int z = finerElement.z / 2;

			elements.insert(Vec3i(x,y,z));
		}

		Vec3i nodeDims = levelDims[i] + Vec3i(1);
		Vec3i finerNodeDims = levelDims[i-1] + Vec3i(1);
		std::map<uint64_t, uint64_t> usedNodesInLevel;
		std::map<Vec3i, Vec3i> coarseToFinerNodeMap;

		for(auto element : elements)
		{
			uint64_t ind;
			Vec3i globalNode;

			FOR3(node, Vec3i(0), Vec3i(2))
			{
				globalNode = element + node;
				nodes.insert(globalNode);

				ind = Linearize(globalNode, nodeDims);
				if(!usedNodesInLevel.contains(ind))
				{
					usedNodesInLevel[ind] = ind;

					if(!coarseToFinerNodeMap.contains(globalNode))
						coarseToFinerNodeMap[globalNode] = Vec3i(globalNode.x * 2, globalNode.y * 2, globalNode.z * 2);
				}
			}
		}

		index = 0;
		for (auto line : usedNodesInLevel)
		{
			usedNodesInLevel[line.first] = index++;
		}

		std::vector<Eigen::Triplet<scalar>> restrtictionTriplets;

		for(auto globalNode : nodes)
		{
			auto finerNode = coarseToFinerNodeMap[globalNode];

			Vec3i diff;
			FOR3(diff, Vec3i(-1), Vec3i(2))
			{
				if(usedNodesInLevels[i-1].contains(Linearize(finerNode + diff, finerNodeDims)))
				{
					auto dist = abs(diff.x) + abs(diff.y) + abs(diff.z);
					scalar val = pow(0.5, 3 + dist);

					for(int c = 0; c < 3; c++)
					{
						restrtictionTriplets.push_back(Eigen::Triplet<scalar>(3 * usedNodesInLevel[Linearize(globalNode, nodeDims)] + c,
																			  3 * usedNodesInLevels[i-1][Linearize(finerNode + diff, finerNodeDims)] + c,
																			  val));
					}
				}
			}
		}

		std::vector<Vec3i> elementsVector(elements.size());
		std::copy(elements.begin(), elements.end(), elementsVector.begin());
 
		levelElements.push_back(elementsVector);
		usedNodesInLevels.push_back(usedNodesInLevel);

		auto restriction = Eigen::SparseMatrix<scalar>(usedNodesInLevels[i].size() * 3,usedNodesInLevels[i-1].size() * 3);
		restriction.setFromTriplets(restrtictionTriplets.begin(), restrtictionTriplets.end());
		auto interpolation = restriction.transpose() * 8;

		//Find fixed nodes on the level
		Eigen::ArrayXi fixedNodesAbove(fixedNodesInLevels[i-1].size());
		index = 0;
		for(auto val : fixedNodesInLevels[i-1])
		{
			fixedNodesAbove[index++] = (int)val;
		}

		Eigen::VectorXd fixingForcesAbove = Eigen::VectorXd::Zero(usedNodesInLevels[i-1].size() * 3);

		fixingForcesAbove(3 * fixedNodesAbove    ) = Eigen::ArrayXd::Ones(fixedNodesAbove.size());
		fixingForcesAbove(3 * fixedNodesAbove + 1) = Eigen::ArrayXd::Ones(fixedNodesAbove.size());
		fixingForcesAbove(3 * fixedNodesAbove + 2) = Eigen::ArrayXd::Ones(fixedNodesAbove.size());

		Eigen::VectorXd fixingForcesOnLevel = restriction * fixingForcesAbove;
		std::set<uint64_t> fixedNodesOnLevel;

		for(index = 0; index < fixingForcesOnLevel.rows(); index++)
		{
			int val = index % 3;

			if(fixingForcesOnLevel(index) != 0.0 && (fixedNodesOnLevel.find(val) == fixedNodesOnLevel.end())) 
				fixedNodesOnLevel.insert(val);
		}
		fixedNodesInLevels.push_back(fixedNodesOnLevel);

		int levelSize = usedNodesInLevel.size() * 3;
		Eigen::Array<int, Eigen::Dynamic, 8> elementToGlobalOnLevel(levelElements[i].size(), 8);
		Eigen::Array<int, Eigen::Dynamic, 27> restrictionMapping(levelElements[i].size(), 27);
		Eigen::VectorXd restrictionCoeffVector = Eigen::VectorXd::Zero(3 * usedNodesInLevels[i - 1].size());
		line = 0;
		for (auto element : levelElements[i])
		{
			int restrictionCoeff = 0;
			FOR3(node, Vec3i(0), Vec3i(2))
			{
				elementToGlobalOnLevel(line, Linearize(node, Vec3i(2))) = (int)usedNodesInLevel[Linearize(element + node, nodeDims)];
			}

			FOR3(node, Vec3i(0), Vec3i(3))
			{
				auto val = Linearize(coarseToFinerNodeMap[element] + node, finerNodeDims);
				if(usedNodesInLevels[i - 1].contains(val))
				{
					restrictionMapping(line, Linearize(node, Vec3i(3))) = (int)usedNodesInLevels[i - 1][val];
					restrictionCoeffVector(3 * usedNodesInLevels[i - 1][val]    )++;
					restrictionCoeffVector(3 * usedNodesInLevels[i - 1][val] + 1)++;
					restrictionCoeffVector(3 * usedNodesInLevels[i - 1][val] + 2)++;
				}
				else
					restrictionMapping(line, Linearize(node, Vec3i(3))) = -1;
			
			}

			line++;
		}
		elementToNodeMatrices.push_back(elementToGlobalOnLevel);
		restrictionMappings.push_back(restrictionMapping);
		restrictionCoefficients.push_back(restrictionCoeffVector.cwiseInverse());
		invDiagKOnLevels.push_back(GetInverseDiagonal<scalar>(levelSize, pow(0.5, 2*i) * elementStiffnessMatrix, elementToGlobalOnLevel));

		if(newLevelDims.x <= 2 || newLevelDims.y <= 2 || newLevelDims.z <= 2)
		{
			numLevels = i + 1;
			std::cout << "Enough coarseness reached, using number of levels: " << numLevels << std::endl;
			break;
		}
	}

	std::map<uint64_t, uint64_t> freeNodes;
	std::vector<uint64_t> freeIndices;
	index = 1; 
	uint64_t expectedIndex = 0;
	for (auto line : usedNodesInLevels[numLevels - 1])
	{
		if(fixedNodesInLevels[numLevels - 1].find(expectedIndex) == fixedNodesInLevels[numLevels - 1].end())
		{
			freeNodes[line.first] = index++;
			freeIndices.push_back(expectedIndex);
		}

		expectedIndex++;
	}
	
	Eigen::ArrayXi freeDoFsOnCoarsest(freeNodes.size() * 3);
	index = 0;
	for(auto val : freeIndices)
	{
		freeDoFsOnCoarsest[index++] = 3 * (int)val;
		freeDoFsOnCoarsest[index++] = 3 * (int)val + 1;
		freeDoFsOnCoarsest[index++] = 3 * (int)val + 2;
	}

	int Ksize = freeNodes.size() * 3;
	std::vector<std::array<uint64_t, 8>> elementToGlobalCoarsest;
	for (auto element : levelElements[numLevels - 1])
	{
		std::array<uint64_t, 8> nodes;
		FOR3(node, Vec3i(0), Vec3i(2))
		{
			nodes[Linearize(node, Vec3i(2))] = freeNodes[Linearize(element + node, levelDims[numLevels - 1] + Vec3i(1))] - 1;
		}
		elementToGlobalCoarsest.push_back(nodes);
	}

	auto Kc = assembleK<scalar>(Ksize, elementToGlobalCoarsest, elementStiffness, pow(0.5, 2*(numLevels - 1)));

	systemMatrix.PrepareMultigrid(numLevels, elementToNodeMatrices, restrictionMappings, restrictionCoefficients, invDiagKOnLevels, Kc, freeDoFsOnCoarsest);
	#endif // MULTIGRID

#else // MATRIX_FREE

	std::map<uint64_t, uint64_t> freeNodes;
	uint64_t index = 1, expectedIndex = 0;
	for (auto line : usedNodes)
	{
		if(fixedNodes.find(expectedIndex) == fixedNodes.end())
			freeNodes[line.first] = index++;
		
		expectedIndex++;
	}

	int size = freeNodes.size() * 3;

	Vec3i node;
	std::vector<std::array<uint64_t, 8>> elementToGlobal;
	for (auto element : usedElements)
	{
		std::array<uint64_t, 8> nodes;
		FOR3(node, Vec3i(0), Vec3i(2))
		{
			nodes[Linearize(node, Vec3i(2))] = freeNodes[Linearize(element + node, nodeGridDimensions)] - 1;
		}
		elementToGlobal.push_back(nodes);
	}

	auto systemMatrix = assembleK<scalar>(size, elementToGlobal, elementStiffness);

	std::cout << "System Matrix : Size:" << systemMatrix.rows() << "x" << systemMatrix.cols() << ", Non-Zero:" << systemMatrix.nonZeros() << ", " << ((float)systemMatrix.nonZeros()/systemMatrix.rows()) << " full per row" << std::endl;
#endif // MATRIX_FREE

	return systemMatrix;
}

#ifdef MATRIX_FREE
template MatrixFreeSparse<double> assembleSystemMatrix<double>(int* voxelModel, Vec3i voxelGridDimensions, double elementStiffness[24][24], const std::set<uint64_t>& fixedNodes);
// template MatrixFreeSparse<float> assembleSystemMatrix<float>(int* voxelModel, Vec3i voxelGridDimensions, double elementStiffness[24][24], const std::set<uint64_t>& fixedNodes);
#else
template Eigen::SparseMatrix<double> assembleSystemMatrix<double>(int* voxelModel, Vec3i voxelGridDimensions, double elementStiffness[24][24], const std::set<uint64_t>& fixedNodes);
template Eigen::SparseMatrix<float> assembleSystemMatrix<float>(int* voxelModel, Vec3i voxelGridDimensions, float elementStiffness[24][24], const std::set<uint64_t>& fixedNodes);
#endif

template <typename scalar>
void applyBoundaryConditions(std::vector<scalar>& f, std::map<uint64_t, Vec3<scalar>>& loadedNodes, const std::set<uint64_t>& fixedNodes)
{
	for(auto node : loadedNodes)
	{
#ifdef MATRIX_FREE
        uint64_t ind = node.first;
#else
        uint64_t ind = node.first - fixedNodes.size();
#endif
		f[3*ind    ] = node.second.x;
		f[3*ind + 1] = node.second.y;
		f[3*ind + 2] = node.second.z;
	}
}

template void applyBoundaryConditions<double>(std::vector<double>& f, std::map<uint64_t, Vec3<double>>& loadedNodes, const std::set<uint64_t>& fixedNodes);
template void applyBoundaryConditions<float>(std::vector<float>& f, std::map<uint64_t, Vec3<float>>& loadedNodes, const std::set<uint64_t>& fixedNodes);

#ifdef MATRIX_FREE
template <typename scalar>
void solveWithEigen(const MatrixFreeSparse<scalar>& A, const std::vector<scalar>& b, std::vector<scalar>& x)
#else
template <typename scalar>
void solveWithEigen(const Eigen::SparseMatrix<scalar>& A, const std::vector<scalar>& b, std::vector<scalar>& x)
#endif
{
	Eigen::Matrix<scalar, Eigen::Dynamic, 1> b_eig(b.size());
	Eigen::Matrix<scalar, Eigen::Dynamic, 1> x_eig(x.size());

	memcpy(b_eig.data(), b.data(), b.size()*sizeof(scalar)); 

#ifdef MATRIX_FREE
	#ifdef MULTIGRID
	std::cout << "Setting up the CG solver with Matrix Free Formulation and Multigrid Preconditioner" << std::endl;
	Eigen::ConjugateGradient<MatrixFreeSparse<scalar>, Eigen::Lower | Eigen::Upper, Eigen::MultigridPreconditioner> solver(A);
	#else
	std::cout << "Setting up the CG solver with Matrix Free Formulation" << std::endl;
	Eigen::ConjugateGradient<MatrixFreeSparse<scalar>, Eigen::Lower | Eigen::Upper, Eigen::IdentityPreconditioner> solver(A);
	#endif
#else
	std::cout << "Setting up the CG solver with Incomplete Cholesky Preconditioner" << std::endl;
	Eigen::ConjugateGradient<Eigen::SparseMatrix<scalar>, Eigen::Lower, Eigen::IncompleteCholeskyPreconditioner<scalar>> solver(A);
#endif
	
#ifdef MAX_ITER
	solver.setMaxIterations(MAX_ITER);
#endif

#ifdef TOLERANCE
	solver.setTolerance(TOLERANCE);
#endif

	std::cout << "Starting the solver" << std::endl;

	auto start = std::chrono::system_clock::now();
	x_eig = solver.solve(b_eig);
	auto end = std::chrono::system_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);

	std::cout << "Solving took:    " << duration.count() << " s" << std::endl;
	std::cout << "#iterations:     " << solver.iterations()      << std::endl;
	std::cout << "Estimated error: " << solver.error()           << std::endl;

	std::memcpy(x.data(), x_eig.data(), x.size()*sizeof(scalar)); 
}

#ifdef MATRIX_FREE
template void solveWithEigen<double>(const MatrixFreeSparse<double>& A, const std::vector<double>& b, std::vector<double>& x);
// template void solveWithEigen<float>(const MatrixFreeSparse<float>& A, const std::vector<float>& b, std::vector<float>& x);
#else
template void solveWithEigen<double>(const Eigen::SparseMatrix<double>& A, const std::vector<double>& b, std::vector<double>& x);
template void solveWithEigen<float>(const Eigen::SparseMatrix<float>& A, const std::vector<float>& b, std::vector<float>& x);
#endif