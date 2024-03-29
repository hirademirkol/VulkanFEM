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

// Assemble the system matrix from the elementToGlobal matrix and one element stiffness matrix
template<typename scalar>
inline Eigen::SparseMatrix<scalar> assembleK(int size, std::vector<std::array<uint64_t, 8>>& elementToGlobal, scalar elementStiffness[24][24], scalar multiplier = 1.0)
{
	Eigen::SparseMatrix<scalar> systemMatrix(size, size);
	Vec3i node, node2;

	std::vector<Eigen::Triplet<scalar>> triplets;

	// Loop over each element
	for(auto line : elementToGlobal)
	{
		// Loop over each set of two nodes on the element
		FOR3(node, Vec3i(0), Vec3i(2))
		{
			FOR3(node2, Vec3i(0), Vec3i(2))
			{
				auto i = Linearize(node,  Vec3i(2));
				auto j = Linearize(node2, Vec3i(2));

				// Get the global node index
				auto iGlobal = line[i];
				auto jGlobal = line[j];

				// Skip if the node is not used
				if(iGlobal == -1 || jGlobal == -1)
					continue;

				// Add the row, column and the coefficient to the triplet list
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

	// Create the sparse matrix from the triplets
	systemMatrix.setFromTriplets(triplets.begin(), triplets.end());
	systemMatrix.makeCompressed();

	return systemMatrix;
}

// Create a list of used elements and nodes from the voxel model
inline void EnlistUsedElements(int* voxelModel, Vec3i voxelGridDimensions, Vec3i nodeGridDimensions,
							   std::vector<Vec3i>& usedElements, std::vector<int>& elementIndices, std::map<uint64_t, uint64_t>& usedNodes)
{
	Vec3i element; int index = 0;

	// Loop over each element in the voxel model
	FOR3(element,  Vec3i(0), voxelGridDimensions)
	{
		// Pick filled ones
		if(voxelModel[Linearize(element, voxelGridDimensions)] == 1)
		{
			usedElements.push_back(element);
			elementIndices.push_back(index++);
		}
	}

	// MAX is set for Vec3's to use when being hashed, largest coordinate size is picked for MAX for unique hashes
	Vec3i::MAX = std::max(nodeGridDimensions.x, nodeGridDimensions.y);
	Vec3i::MAX = std::max(nodeGridDimensions.z, Vec3i::MAX);

	Vec3i node;

	// Loop over each corner of each element and create a unique set of used nodes' global indices
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

// Calculate the inverse diagonal into a vector from element stiffness matrix and elementToNode array
template <typename scalar>
inline Eigen::VectorXd GetInverseDiagonal(int size, Eigen::Matrix<scalar, 24, 24> elementStiffnessMat, Eigen::Array<int, Eigen::Dynamic, 4> elementToNode)
{
	Eigen::VectorXd diag = Eigen::VectorXd::Zero(size);

	const Eigen::Array<int, 1, 24> c   {0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5};
	const Eigen::Array<int, 1, 24> xInd{0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3};

	// Add the contribution of each element to the diagonal
	for (auto line : elementToNode.rowwise())
	{	
		Eigen::Array<int, 1, 24> xs = 3 * line(xInd) + c;
		diag(xs) += elementStiffnessMat.diagonal();
	}
	Eigen::VectorXd invDiag = diag.cwiseInverse();

	return invDiag;
}

// Calculate the inverse diagonal into a vector from unique element stiffness matrices and elementToNode array
template <typename scalar>
inline Eigen::VectorXd GetInverseDiagonal(int size, Eigen::Matrix<scalar, 24, 24> elementStiffnessMat, Eigen::Array<int, Eigen::Dynamic, 4> elementToNode, std::map<int, Eigen::Matrix<double, 24, 24>> uniqueElementStiffnessMatrices)
{
	Eigen::VectorXd diag = Eigen::VectorXd::Zero(size);

	const Eigen::Array<int, 1, 24> c   {0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5};
	const Eigen::Array<int, 1, 24> xInd{0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3};

	int index = 0;
	
	// Add the contribution of each element to the diagonal
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
MatrixFreeSparse<scalar> assembleSystemMatrix(int* voxelModel, Vec3i voxelGridDimensions, double elementStiffness[24][24], const std::set<uint64_t>& fixedNodes, int numLevels, int skipLevels)
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

	// Each node is given an ordered index 1 to numIndices
	uint64_t index = 0;
	for (auto line : usedNodes)
	{
		usedNodes[line.first] = index++;
	}

	int size = usedNodes.size() * 3;

	Eigen::Array<int, Eigen::Dynamic, 4> elementToGlobal(usedElements.size(), 4);

	// elementToGlobal is filled with global index of each node in each element
	std::for_each(std::execution::par_unseq, elementIndices.begin(), elementIndices.end(), [&](int index)
	{
		Vec3i node;
		FOR3(node, Vec3i(0), Vec3i(1,2,2))
		{
			elementToGlobal(index, Linearize(node, Vec3i(1,2,2))) = (int)usedNodes[Linearize(usedElements[index] + node, nodeGridDimensions)];
		}
	});

	// Element stiffness matrix is copied into an Eigen matrix
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

	// Fixed nodes are copied into an Eigen array
	Eigen::ArrayXi fixed(fixedNodes.size());
	index = 0;
	for(auto val : fixedNodes)
	{
		fixed[index++] = (int)val;
	}

	// System matrix is created in matrix-free form
	MatrixFreeSparse<scalar> systemMatrix(size, elementStiffnessMatrix, elementToGlobal, fixed);

	#ifdef MULTIGRID

	std::cout << "Preparing the Multigrid structure" << std::endl;

	// Restriction operators are created from data
	Eigen::Array<int, 8, 8> restrictionOnElement =  Eigen::Map<Eigen::Array<int, 8, 8>>(const_cast<int*>(restrictionOnElementValues.data()));
	Eigen::Array<int, 64, 8> restrictionOnElement4 =  Eigen::Map<Eigen::Array<int, 8, 64>>(const_cast<int*>(restrictionOnElement4Values.data()), 8, 64).transpose();
	Eigen::Matrix<double, 27, 8> restrictionOperator = Eigen::Map<Eigen::Matrix<double, 8, 27>>(const_cast<double*>(restrictionOperatorValues.data()), 8, 27).transpose();
	Eigen::Matrix<double, 125, 8> restrictionOperator4 = Eigen::Map<Eigen::Matrix<double, 8, 125>>(const_cast<double*>(restrictionOperator4Values.data()), 8, 125).transpose();
	Eigen::Matrix<double, 81, 24> restrictionOperatorDoF = Eigen::Map<Eigen::Matrix<double, 24, 81>>(const_cast<double*>(restrictionOperatorDoFValues.data()), 24, 81).transpose();
	Eigen::Matrix<double, 375, 24> restrictionOperatorDoF4 = Eigen::Map<Eigen::Matrix<double, 24, 375>>(const_cast<double*>(restrictionOperatorDoF4Values.data()), 24, 375).transpose();

	// Multigrid related matrix and arrays are prepared
	std::vector<Vec3i> levelDims;
	std::vector<std::vector<Vec3i>> levelElements;
	std::vector<std::map<uint64_t, uint64_t>> usedNodesInLevels;
	std::vector<std::set<uint64_t>> fixedNodesInLevels;
	std::vector<Eigen::VectorXd> invDiagKOnLevels;
  	std::vector<Eigen::Array<int, Eigen::Dynamic, 4>> elementToNodeMatrices;
  	std::vector<Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic>> restrictionMappings;
	std::vector<Eigen::Matrix<scalar, Eigen::Dynamic, 1>> restrictionCoefficients;

	// First level matrices are already created
	levelDims.push_back(voxelGridDimensions);
	levelElements.push_back(usedElements);
	usedNodesInLevels.push_back(usedNodes);
	fixedNodesInLevels.push_back(fixedNodes);
	invDiagKOnLevels.push_back(GetInverseDiagonal(size, elementStiffnessMatrix, elementToGlobal));

	// Loop over each Multigrid level
	for(int i = 1; i < numLevels; i++)
	{
		// Calculation of new dimensions
		Vec3i newLevelDims;
		newLevelDims.x = levelDims[i-1].x / 2 + (levelDims[i-1].x % 2);
		newLevelDims.y = levelDims[i-1].y / 2 + (levelDims[i-1].y % 2);
		newLevelDims.z = levelDims[i-1].z / 2 + (levelDims[i-1].z % 2);
		levelDims.push_back(newLevelDims);

		std::set<Vec3i> elements, nodes;
		std::map<Vec3i, std::vector<Vec3i>> coarseToFinerElementMap;

		// Calculation of divisions, can be 4 in the case of level skipping
		int divider = 2;
		int numberOfChildren = 8;
		if(i == 1 && skipLevels > 0)
		{
			divider = pow(2, skipLevels + 1);
			numberOfChildren = pow(divider, 3);
		}

		// Each fine element is mapped to coarse elements
		for(auto finerElement: levelElements[i-1])
		{
			int x = finerElement.x / divider;
			int y = finerElement.y / divider;
			int z = finerElement.z / divider;

			elements.insert(Vec3i(x,y,z));
			coarseToFinerElementMap[Vec3i(x,y,z)].push_back(finerElement);
		}

		Vec3i nodeDims = levelDims[i] + Vec3i(1);
		Vec3i finerNodeDims = levelDims[i-1] + Vec3i(1);
		std::map<uint64_t, uint64_t> usedNodesInLevel;
		std::map<Vec3i, Vec3i> coarseToFinerNodeMap;
		std::map<int, Eigen::Matrix<double, 24, 24>> uniqueKsOnLevel;

		// Each element is looped over and used nodes are listed, these nodes are also mapped to their corresponing finer nodes
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
						coarseToFinerNodeMap[globalNode] = Vec3i(globalNode.x * divider, globalNode.y * divider, globalNode.z * divider);
				}
			}
		}

		// Each node is given an ordered index 1 to numIndices on level
		index = 0;
		for (auto line : usedNodesInLevel)
		{
			usedNodesInLevel[line.first] = index++;
		}
		
		// Elements set is copied into a vector and a corresponding element index vector is created
		std::vector<Vec3i> elementsVector(elements.size());
		std::copy(elements.begin(), elements.end(), elementsVector.begin());
		std::vector<int> levelElementIndices(elementsVector.size());
		std::iota(levelElementIndices.begin(), levelElementIndices.end(), 0);

		levelElements.push_back(elementsVector);
		usedNodesInLevels.push_back(usedNodesInLevel);

		int levelSize = usedNodesInLevel.size() * 3;
		int numBoundingNodesOnAxis = divider + 1;
		Eigen::Array<int, Eigen::Dynamic, 4> elementToGlobalOnLevel(levelElements[i].size(), 4);
		Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic> restrictionMapping(levelElements[i].size(), (int)pow(numBoundingNodesOnAxis, 3));
		Eigen::VectorXd restrictionCoeffVector = Eigen::VectorXd::Zero(3 * usedNodesInLevels[i - 1].size());

		std::mutex mapMutex;

		// Parallelly, each element is used to be mapped to their nodes' global indices,
		// Restriction coefficient and mapping is calculated depending on the weights of each finer node inside the element
		std::for_each(std::execution::par_unseq, levelElementIndices.begin(), levelElementIndices.end(), [&](int index)
		{
			Vec3i node; int restrictionCoeff = 0;
			FOR3(node, Vec3i(0), Vec3i(1,2,2))
			{
				elementToGlobalOnLevel(index, Linearize(node, Vec3i(1,2,2))) = (int)usedNodesInLevel[Linearize(elementsVector[index] + node, nodeDims)];
			}

			FOR3(node, Vec3i(0), Vec3i(numBoundingNodesOnAxis))
			{
				auto val = Linearize(coarseToFinerNodeMap[elementsVector[index]] + node, finerNodeDims);
				if(usedNodesInLevels[i - 1].contains(val))
				{
					restrictionMapping(index, Linearize(node, Vec3i(numBoundingNodesOnAxis))) = (int)usedNodesInLevels[i - 1][val];
					restrictionCoeffVector(3 * usedNodesInLevels[i - 1][val]    )++;
					restrictionCoeffVector(3 * usedNodesInLevels[i - 1][val] + 1)++;
					restrictionCoeffVector(3 * usedNodesInLevels[i - 1][val] + 2)++;
				}
				else
				{
					restrictionMapping(index, Linearize(node, Vec3i(numBoundingNodesOnAxis))) = -1;
				}
			}

			// If a coarse element has non filled child (finer) elements, its unique element stiffness matrix is calculated
			if(coarseToFinerElementMap[elementsVector[index]].size() < numberOfChildren)
			{
				Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> tmpKfine = Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(3*pow(numBoundingNodesOnAxis, 3), 3*pow(numBoundingNodesOnAxis, 3));
				
				const Eigen::Array<int, 1, 24> c   {0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2};
				const Eigen::Array<int, 1, 24> xInd{0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7};

				Vec3i firstFinerElement = Vec3i(elementsVector[index].x * divider, elementsVector[index].y * divider, elementsVector[index].z * divider);

				// Each element's contribution is added into a temporary K, as values from the global stiffness matrix on level
				for(auto finerElement : coarseToFinerElementMap[elementsVector[index]])
				{
					Vec3i element = finerElement - firstFinerElement;
					int elementIndex = Linearize(element, Vec3i(divider));

					if(i == 1 && skipLevels == 1)
						tmpKfine(3 * restrictionOnElement4(elementIndex, xInd) + c, 3 * restrictionOnElement4(elementIndex, xInd) + c) += pow(2, i - 1) * elementStiffnessMatrix;
					else
						tmpKfine(3 * restrictionOnElement(elementIndex, xInd) + c, 3 * restrictionOnElement(elementIndex, xInd) + c) += pow(2, i - 1) * elementStiffnessMatrix;
				}

				// Temporary K is restricted to create the unique element stiffness matrix
				Eigen::Matrix<scalar, 24, 24> K;
				if(i == 1 && skipLevels == 1)
					K = restrictionOperatorDoF4.transpose() * tmpKfine * restrictionOperatorDoF4;
				else
					K = restrictionOperatorDoF.transpose() * tmpKfine * restrictionOperatorDoF;

				std::lock_guard<std::mutex> lock(mapMutex);
				uniqueKsOnLevel[index] = K;
			}
		});

		elementToNodeMatrices.push_back(elementToGlobalOnLevel);
		restrictionMappings.push_back(restrictionMapping);
		restrictionCoefficients.push_back(restrictionCoeffVector.cwiseInverse());

		// The fixed nodes on the level above is written in an array
		Eigen::ArrayXi fixedNodesAbove(fixedNodesInLevels[i-1].size());
		index = 0;
		for(auto val : fixedNodesInLevels[i-1])
		{
			fixedNodesAbove[index++] = (int)val;
		}

		Eigen::VectorXd fixingForcesAbove = Eigen::VectorXd::Zero(usedNodesInLevels[i-1].size() * 3);

		// Fixing forces on level above are created as a vector of ones
		fixingForcesAbove(3 * fixedNodesAbove    ) = Eigen::ArrayXd::Ones(fixedNodesAbove.size());
		fixingForcesAbove(3 * fixedNodesAbove + 1) = Eigen::ArrayXd::Ones(fixedNodesAbove.size());
		fixingForcesAbove(3 * fixedNodesAbove + 2) = Eigen::ArrayXd::Ones(fixedNodesAbove.size());

		// Fixing forces are restricted down to this level to pick the fixed nodes on the level
		Eigen::Matrix<double, Eigen::Dynamic, 1> fixingForcesOnLevel = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(usedNodesInLevel.size() * 3);
		Eigen::Matrix<double, Eigen::Dynamic, 1> tempR = fixingForcesAbove.cwiseProduct(restrictionCoefficients[i - 1]);
		
		int num = 0;
		for(auto line : elementToGlobalOnLevel.rowwise())
		{
			Eigen::Array<int, Eigen::Dynamic, 1> mapping = restrictionMapping.row(num);
			Eigen::Matrix<double, Eigen::Dynamic, 1> mask = Eigen::Matrix<double, Eigen::Dynamic, 1>::Ones(mapping.rows());

			for(int i = 0; i < pow(numBoundingNodesOnAxis, 3); i++)
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
				if(i == 1 && skipLevels == 1)
					fixingForcesOnLevel(3 * (line(xInd) + offset) + c) += tempR(3 * mapping + c).cwiseProduct(mask).transpose() * restrictionOperator4;
				else
					fixingForcesOnLevel(3 * (line(xInd) + offset) + c) += tempR(3 * mapping + c).cwiseProduct(mask).transpose() * restrictionOperator;
			}

			num++;
		}

		// The fixed nodes on this level are picked
		std::set<uint64_t> fixedNodesOnLevel;

		for(index = 0; index < fixingForcesOnLevel.rows(); index++)
		{
			int val = index / 3;

			if(fixingForcesOnLevel(index) != 0.0 && (fixedNodesOnLevel.find(val) == fixedNodesOnLevel.end())) 
				fixedNodesOnLevel.insert(val);
		}
		fixedNodesInLevels.push_back(fixedNodesOnLevel);

		// Inverse diagonal of the system matrix on the level is calculated
		invDiagKOnLevels.push_back(GetInverseDiagonal<scalar>(levelSize, pow(2, i + skipLevels) * elementStiffnessMatrix, elementToGlobalOnLevel, uniqueKsOnLevel));

		if(newLevelDims.x <= 2 || newLevelDims.y <= 2 || newLevelDims.z <= 2)
		{
			numLevels = i + 1;
			std::cout << "Enough coarseness reached, using number of levels: " << numLevels << std::endl;
			break;
		}
	}

	// Free nodes and DoFs on the coarsest level are listed
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

	// ElementToGlobal matrix for the coarsest level is created
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

	// Explicit system matrix for the coarsest level is assembled
	auto Kc = assembleK<scalar>(Ksize, elementToGlobalCoarsest, elementStiffness, pow(2, numLevels + skipLevels - 1));

	// Multigrid related data is saved into the system
	systemMatrix.PrepareMultigrid(numLevels, skipLevels, elementToNodeMatrices, restrictionMappings, restrictionCoefficients, invDiagKOnLevels, Kc, freeDoFsOnCoarsest);
	#endif // MULTIGRID

#else // MATRIX_FREE

	// Free Nodes are enlisted for usage in the matrix solver
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

	// Each element is looped over and their free nodes' global indices are listed in the elementToGlobal matrix
	for (auto element : usedElements)
	{
		std::array<uint64_t, 8> nodes;
		FOR3(node, Vec3i(0), Vec3i(2))
		{
			nodes[Linearize(node, Vec3i(2))] = freeNodes[Linearize(element + node, nodeGridDimensions)] - 1;
		}
		elementToGlobal.push_back(nodes);
	}

	Eigen::SparseMatrix<scalar> systemMatrix = assembleK<scalar>(size, elementToGlobal, elementStiffness);

	std::cout << "System Matrix : Size:" << systemMatrix.rows() << "x" << systemMatrix.cols() << ", Non-Zero:" << systemMatrix.nonZeros() << ", " << ((float)systemMatrix.nonZeros()/systemMatrix.rows()) << " full per row" << std::endl;
#endif // MATRIX_FREE

	return systemMatrix;
}

#ifdef MATRIX_FREE
	#ifdef MULTIGRID
template MatrixFreeSparse<double> assembleSystemMatrix<double>(int* voxelModel, Vec3i voxelGridDimensions, double elementStiffness[24][24], const std::set<uint64_t>& fixedNodes, int numLevels, int skipLevels);
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
	// Add each loading condition to the RHS vector
	for(auto node : loadedNodes)
	{
#ifdef MATRIX_FREE
        uint64_t ind = node.first;
#else
		// For built matrix, fixed nodes are excluded. Here it is assumed that all fixed nodes removed are come before the loaded nodes
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
	// Data is copied into Eigen matrices
	Eigen::Matrix<scalar, Eigen::Dynamic, 1> b_eig(b.size());
	Eigen::Matrix<scalar, Eigen::Dynamic, 1> x_eig(x.size());

	memcpy(b_eig.data(), b.data(), b.size()*sizeof(scalar)); 

	// CG Solver is created with corresponding preconditioner and matrix form
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