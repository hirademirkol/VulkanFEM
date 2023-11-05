#include "out.hpp"
#include "FEM.hpp"
#include "IncompleteCholeskyPreconditioner.hpp"


#ifdef MATRIX_FREE
#ifdef MULTIGRID
	#include "MultigridPreconditioner.hpp"
	#include "RestrictionOperator.hpp"
#endif
#endif

#include <chrono>
#include <algorithm>
#include <execution>
#include <numeric>

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
							   std::vector<Vec3i>& usedElements, std::vector<int>& elementIndices, std::map<uint64_t, uint64_t>& usedNodes)
{
	Vec3i element; int index = 0;
	FOR3(element,  Vec3i(0), voxelGridDimensions)
	{
		if(voxelModel[Linearize(element, voxelGridDimensions)] == 1)
		{
			usedElements.push_back(element);
			elementIndices.push_back(index++);
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
inline Eigen::VectorXd GetInverseDiagonal(int size, Eigen::Matrix<scalar, 24, 24> elementStiffnessMat, Eigen::Array<int, Eigen::Dynamic, 4> elementToNode)
{
	Eigen::VectorXd diag = Eigen::VectorXd::Zero(size);

	const Eigen::Array<int, 1, 24> c   {0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5};
	const Eigen::Array<int, 1, 24> xInd{0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3};

	for (auto line : elementToNode.rowwise())
	{	
		Eigen::Array<int, 1, 24> xs = 3 * line(xInd) + c;
		diag(xs) += elementStiffnessMat.diagonal();
	}
	Eigen::VectorXd invDiag = diag.cwiseInverse();

	return invDiag;
}

template <typename scalar>
inline Eigen::VectorXd GetInverseDiagonal(int size, Eigen::Matrix<scalar, 24, 24> elementStiffnessMat, Eigen::Array<int, Eigen::Dynamic, 4> elementToNode, std::map<int, Eigen::Matrix<double, 24, 24>> uniqueElementStiffnessMatrices)
{
	Eigen::VectorXd diag = Eigen::VectorXd::Zero(size);

	const Eigen::Array<int, 1, 24> c   {0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5};
	const Eigen::Array<int, 1, 24> xInd{0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3};

	int index = 0;
	for (auto line : elementToNode.rowwise())
	{	
		Eigen::Array<int, 1, 24> xs = 3 * line(xInd) + c;
		
		if(uniqueElementStiffnessMatrices.contains(index))
			diag(xs) += uniqueElementStiffnessMatrices[index].diagonal();
		else
			diag(xs) += elementStiffnessMat.diagonal();

		index++;
	}

	Eigen::VectorXd invDiag = diag.cwiseInverse();

	return invDiag;
}

#endif

#ifdef MATRIX_FREE
	#ifdef MULTIGRID
template<typename scalar>
MatrixFreeSparse<scalar> assembleSystemMatrix(int* voxelModel, Vec3i voxelGridDimensions, double elementStiffness[24][24], const std::set<uint64_t>& fixedNodes, int numLevels)
    #else
template<typename scalar>
MatrixFreeSparse<scalar> assembleSystemMatrix(int* voxelModel, Vec3i voxelGridDimensions, double elementStiffness[24][24], const std::set<uint64_t>& fixedNodes)
	#endif
#else
//TODO: pass the used nodes out for applyBoundaryConditions to map global nodes to matrix nodes
template<typename scalar>
Eigen::SparseMatrix<scalar> assembleSystemMatrix(int* voxelModel, Vec3i voxelGridDimensions, scalar elementStiffness[24][24], const std::set<uint64_t>& fixedNodes)
#endif
{
	std::cout << "Assembling the system matrix" << std::endl;
	Vec3i nodeGridDimensions = voxelGridDimensions + Vec3i(1);

	std::vector<Vec3i> usedElements;
	std::vector<int> elementIndices;
	std::map<uint64_t, uint64_t> usedNodes;
	EnlistUsedElements(voxelModel, voxelGridDimensions, nodeGridDimensions, usedElements, elementIndices, usedNodes);

#ifdef MATRIX_FREE

	uint64_t index = 0;
	for (auto line : usedNodes)
	{
		usedNodes[line.first] = index++;
	}

	int size = usedNodes.size() * 3;

	Eigen::Array<int, Eigen::Dynamic, 4> elementToGlobal(usedElements.size(), 4);

	std::for_each(std::execution::par_unseq, elementIndices.begin(), elementIndices.end(), [&](int index)
	{
		Vec3i node;
		FOR3(node, Vec3i(0), Vec3i(1,2,2))
		{
			elementToGlobal(index, Linearize(node, Vec3i(1,2,2))) = (int)usedNodes[Linearize(usedElements[index] + node, nodeGridDimensions)];
		}
	});

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

	std::array<int, 64> restrictionOnElementData =
	{
		0, 1, 3, 4, 9, 10, 12, 13,
		1, 2, 4, 5, 10, 11, 13, 14,
		3, 4, 6, 7, 12, 13, 15, 16,
		4, 5, 7, 8, 13, 14, 16, 17,
		9, 10, 12, 13, 18, 19, 21, 22,
		10, 11, 13, 14, 19, 20, 22, 23,
		12, 13, 15, 16, 21, 22, 24, 25,
		13, 14, 16, 17, 22, 23, 25, 26
	};
	Eigen::Array<int, 8, 8> restrictionOnElement =  Eigen::Map<Eigen::Array<int, 8, 8>>(restrictionOnElementData.data());
	Eigen::Matrix<double, 27, 8> restrictionOperator = Eigen::Map<Eigen::Matrix<double, 8, 27>>(const_cast<double*>(restrictionOperatorValues.data()), 8, 27).transpose();
	Eigen::Matrix<double, 81, 24> restrictionOperatorDoF = Eigen::Map<Eigen::Matrix<double, 24, 81>>(const_cast<double*>(restrictionOperatorDoFValues.data()), 24, 81).transpose();

	//Prepare multigrid related matrices
	std::vector<Vec3i> levelDims;
	std::vector<std::vector<Vec3i>> levelElements;
	std::vector<std::map<uint64_t, uint64_t>> usedNodesInLevels;
	std::vector<std::set<uint64_t>> fixedNodesInLevels;
	std::vector<Eigen::VectorXd> invDiagKOnLevels;
  	std::vector<Eigen::Array<int, Eigen::Dynamic, 4>> elementToNodeMatrices;
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
		std::map<Vec3i, std::vector<Vec3i>> coarseToFinerElementMap;
		for(auto finerElement: levelElements[i-1])
		{
			int x = finerElement.x / 2;
			int y = finerElement.y / 2;
			int z = finerElement.z / 2;

			elements.insert(Vec3i(x,y,z));
			coarseToFinerElementMap[Vec3i(x,y,z)].push_back(finerElement);
		}

		Vec3i nodeDims = levelDims[i] + Vec3i(1);
		Vec3i finerNodeDims = levelDims[i-1] + Vec3i(1);
		std::map<uint64_t, uint64_t> usedNodesInLevel;
		std::map<Vec3i, Vec3i> coarseToFinerNodeMap;
		std::map<int, Eigen::Matrix<double, 24, 24>> uniqueKsOnLevel;

		for(auto element : elements)
		{
			uint64_t ind;
			Vec3i globalNode;
			Vec3i node;

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
		
		std::vector<Vec3i> elementsVector(elements.size());
		std::copy(elements.begin(), elements.end(), elementsVector.begin());
		std::vector<int> levelElementIndices(elementsVector.size());
		std::iota(levelElementIndices.begin(), levelElementIndices.end(), 0);

		levelElements.push_back(elementsVector);
		usedNodesInLevels.push_back(usedNodesInLevel);

		int levelSize = usedNodesInLevel.size() * 3;
		Eigen::Array<int, Eigen::Dynamic, 4> elementToGlobalOnLevel(levelElements[i].size(), 4);
		Eigen::Array<int, Eigen::Dynamic, 27> restrictionMapping(levelElements[i].size(), 27);
		Eigen::VectorXd restrictionCoeffVector = Eigen::VectorXd::Zero(3 * usedNodesInLevels[i - 1].size());

		std::mutex mapMutex;

		std::for_each(std::execution::par_unseq, levelElementIndices.begin(), levelElementIndices.end(), [&](int index)
		{
			Vec3i node; int restrictionCoeff = 0;
			FOR3(node, Vec3i(0), Vec3i(1,2,2))
			{
				elementToGlobalOnLevel(index, Linearize(node, Vec3i(1,2,2))) = (int)usedNodesInLevel[Linearize(elementsVector[index] + node, nodeDims)];
			}

			FOR3(node, Vec3i(0), Vec3i(3))
			{
				auto val = Linearize(coarseToFinerNodeMap[elementsVector[index]] + node, finerNodeDims);
				if(usedNodesInLevels[i - 1].contains(val))
				{
					restrictionMapping(index, Linearize(node, Vec3i(3))) = (int)usedNodesInLevels[i - 1][val];
					restrictionCoeffVector(3 * usedNodesInLevels[i - 1][val]    )++;
					restrictionCoeffVector(3 * usedNodesInLevels[i - 1][val] + 1)++;
					restrictionCoeffVector(3 * usedNodesInLevels[i - 1][val] + 2)++;
				}
				else
				{
					restrictionMapping(index, Linearize(node, Vec3i(3))) = -1;
				}
			}

			if(coarseToFinerElementMap[elementsVector[index]].size() < 8)
			{
				Eigen::Matrix<scalar, 81, 81> tmpKfine = Eigen::Matrix<scalar, 81, 81>::Zero(81, 81);

				const Eigen::Array<int, 1, 24> c   {0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2};
				const Eigen::Array<int, 1, 24> xInd{0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7};

				Vec3i firstFinerElement = Vec3i(elementsVector[index].x * 2, elementsVector[index].y * 2, elementsVector[index].z * 2);
				for(auto finerElement : coarseToFinerElementMap[elementsVector[index]])
				{
					Vec3i element = finerElement - firstFinerElement;
					int elementIndex = Linearize(element, Vec3i(2));

					tmpKfine(3 * restrictionOnElement(elementIndex, xInd) + c, 3 * restrictionOnElement(elementIndex, xInd) + c) += pow(2, i - 1) * elementStiffnessMatrix;
				}

				std::lock_guard<std::mutex> lock(mapMutex);
				uniqueKsOnLevel[index] = restrictionOperatorDoF.transpose() * tmpKfine * restrictionOperatorDoF;
			}
		});

		elementToNodeMatrices.push_back(elementToGlobalOnLevel);
		restrictionMappings.push_back(restrictionMapping);
		restrictionCoefficients.push_back(restrictionCoeffVector.cwiseInverse());

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

		//Restrict fixing forces
		Eigen::Matrix<double, Eigen::Dynamic, 1> fixingForcesOnLevel = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(usedNodesInLevel.size() * 3);
		Eigen::Matrix<double, Eigen::Dynamic, 1> tempR = fixingForcesAbove.cwiseProduct(restrictionCoefficients[i - 1]);
		
		int num = 0;
		for(auto line : elementToGlobalOnLevel.rowwise())
		{
			Eigen::Array<int, 27, 1> mapping = restrictionMapping.row(num);
			Eigen::Matrix<double, 27, 1> mask = Eigen::Matrix<double, 27, 1>::Ones(mapping.rows());

			for(int i = 0; i < 27; i++)
			{
				if(mapping(i) == -1)
				{
					mask(i) = 0;
					mapping(i) = 0;
				}
			}
			
			const Eigen::Array<int, 1, 8> xInd{0, 0, 1, 1, 2, 2, 3, 3};
			const Eigen::Array<int, 1, 8> offset{0, 1, 0, 1, 0, 1, 0, 1};

			for(int c = 0; c < 3; c++)
			{
				fixingForcesOnLevel(3 * (line(xInd) + offset) + c) += tempR(3 * mapping + c).cwiseProduct(mask).transpose() * restrictionOperator;
			}

			num++;
		}

		//Enlist fixed nodes on level
		std::set<uint64_t> fixedNodesOnLevel;

		for(index = 0; index < fixingForcesOnLevel.rows(); index++)
		{
			int val = index / 3;

			if(fixingForcesOnLevel(index) != 0.0 && (fixedNodesOnLevel.find(val) == fixedNodesOnLevel.end())) 
				fixedNodesOnLevel.insert(val);
		}
		fixedNodesInLevels.push_back(fixedNodesOnLevel);

		invDiagKOnLevels.push_back(GetInverseDiagonal<scalar>(levelSize, pow(2, i) * elementStiffnessMatrix, elementToGlobalOnLevel, uniqueKsOnLevel));

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
		Vec3i node;
		FOR3(node, Vec3i(0), Vec3i(2))
		{
			nodes[Linearize(node, Vec3i(2))] = freeNodes[Linearize(element + node, levelDims[numLevels - 1] + Vec3i(1))] - 1;
		}
		elementToGlobalCoarsest.push_back(nodes);
	}

	auto Kc = assembleK<scalar>(Ksize, elementToGlobalCoarsest, elementStiffness, pow(2, numLevels - 1));

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
	#ifdef MULTIGRID
template MatrixFreeSparse<double> assembleSystemMatrix<double>(int* voxelModel, Vec3i voxelGridDimensions, double elementStiffness[24][24], const std::set<uint64_t>& fixedNodes, int numLevels);
	#else
template MatrixFreeSparse<double> assembleSystemMatrix<double>(int* voxelModel, Vec3i voxelGridDimensions, double elementStiffness[24][24], const std::set<uint64_t>& fixedNodes);
// template MatrixFreeSparse<float> assembleSystemMatrix<float>(int* voxelModel, Vec3i voxelGridDimensions, double elementStiffness[24][24], const std::set<uint64_t>& fixedNodes);
	#endif
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
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

	std::cout << "Solving took:    " << (float)duration.count() / 1000 << " s" << std::endl;
	std::cout << "#iterations:     " << solver.iterations()      		<< std::endl;
	std::cout << "Estimated error: " << solver.error()           		<< std::endl;

	std::memcpy(x.data(), x_eig.data(), x.size()*sizeof(scalar)); 
}

#ifdef MATRIX_FREE
template void solveWithEigen<double>(const MatrixFreeSparse<double>& A, const std::vector<double>& b, std::vector<double>& x);
// template void solveWithEigen<float>(const MatrixFreeSparse<float>& A, const std::vector<float>& b, std::vector<float>& x);
#else
template void solveWithEigen<double>(const Eigen::SparseMatrix<double>& A, const std::vector<double>& b, std::vector<double>& x);
template void solveWithEigen<float>(const Eigen::SparseMatrix<float>& A, const std::vector<float>& b, std::vector<float>& x);
#endif