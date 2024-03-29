#include <memory>
#include <vector>
#include <chrono>

#include "utils.hpp"
#include "in.hpp"
#include "out.hpp"

#include "FEM.hpp"
#include "FEMDefines.hpp"

// Uncomment to use the CPU solver from Eigen
// #define CPU

#ifndef CPU

#include "FEMGPU.hpp"

// Only double precision is implemented for the GPU solver
#define real double

#else

#define real double //or float

#endif

int main(int argc, char* argv[])
{
	Vec3i voxelModelSize; // Number of elements on each axis

	if(argc < 2)
	{
		std::cout << "An input file is necessary" << std::endl;
		return EXIT_FAILURE;
	}
	std::string fileName(argv[1]);
	std::string modelName = fileName.substr(0, fileName.find('.', 1)); // Name of the file to be used for further operations, in path format

	auto start = std::chrono::system_clock::now();
	int* voxelModel = loadVoxel(modelName, voxelModelSize.x, voxelModelSize.y, voxelModelSize.z);

#ifdef MULTIGRID
	// Parse the necessary variables for Multigrid Preconditioner
	int numLevels = 3;
	int skipLevels = 0;
	if(argc > 2)
	{
		if(std::string(argv[2]) == "-n") numLevels = std::stoi(argv[3]);
	}
	if(argc > 4)
	{
		if(std::string(argv[4]) == "--skip") skipLevels = std::stoi(argv[5]);
	}
	std::cout << "Number of levels for Multigrid: " << numLevels;
	if(skipLevels != 0) std::cout << " (+ " << skipLevels << " jump(s))";
	std::cout << std::endl;
#endif

	// Get the element stiffness matrix and the boundary conditions from the files
	double Ke[24][24];
	getKe(modelName, Ke);

	std::set<uint64_t> fixedNodes;
	std::map<uint64_t, Vec3<real>> loadedNodes;
	getBoundaryConditions(modelName, fixedNodes, loadedNodes);

	// Assemble the system matrix according to the given configurations
#ifndef MATRIX_FREE
	Eigen::SparseMatrix<real> systemMatrix = assembleSystemMatrix<real>(voxelModel, voxelModelSize, Ke, fixedNodes);
#else
	#ifndef MULTIGRID
	MatrixFreeSparse<real> systemMatrix = assembleSystemMatrix<real>(voxelModel, voxelModelSize, Ke, fixedNodes);
	#else
	MatrixFreeSparse<real> systemMatrix = assembleSystemMatrix<real>(voxelModel, voxelModelSize, Ke, fixedNodes, numLevels, skipLevels);
	#endif
#endif

	uint64_t numDoF = systemMatrix.rows();

	std::cout << "System statistics:" << std::endl;
	std::cout << "\tVoxel model resolution: " << voxelModelSize.x << " x " << voxelModelSize.y << " x " << voxelModelSize.z << std::endl;
#ifdef MATRIX_FREE
	std::cout << "\tElements: " << systemMatrix.elementToNode.rows() << std::endl;
#endif // MATRIX_FREE
	std::cout << "\tDegrees of freedom: " << numDoF << std::endl;

	std::vector<real> f, u;
	f.resize(numDoF, 0.0);
	u.resize(numDoF, 0.0);

	applyBoundaryConditions(f, loadedNodes, fixedNodes);

	auto end = std::chrono::system_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	std::cout << "Preparation took:    " << (float)duration.count() / 1000 << " s" << std::endl;

	// Solve the system with the defined solver
#ifdef CPU
	solveWithEigen<real>(systemMatrix, f, u);
#else

#ifndef MATRIX_FREE
	std::cout << "Sparse matrix is not implemented for kompute." << std::endl;
	return EXIT_FAILURE;
#else
	solveWithKompute<real>(systemMatrix, f, u);
#endif

#endif
	
	// Calculate and print the compliance
	double compliance = 0;
	for(int i = 0; i < numDoF; i++)
		compliance += f[i] * u[i];

	std::cout << "Compliance:      " << compliance << std::endl;

	// Delete the voxel data
	delete[] voxelModel;

	// Save the solution
	std::cout << "Saving the solution" << std::endl;
	saveVector(u, modelName);

	return EXIT_SUCCESS;
}