#include <memory>
#include <vector>

#define CPU

#ifndef CPU
#include <kompute/Kompute.hpp>

// Compiled shaders
#include "MatxVec.hpp"
#include "GaussSeidel.hpp"
#include "CG.hpp"

typedef std::shared_ptr<kp::Tensor> Tensor;
using TensorDataTypes = kp::Tensor::TensorDataTypes;
using TensorTypes = kp::Tensor::TensorTypes;
#endif

#include "utils.hpp"
#include "in.hpp"
#include "out.hpp"

#include "Ke.h"
#include "FEM.hpp"

int main()
{
	Vec3i voxelModelSize; // Number of elements on each axis
	Vec3d elementSize = Vec3d(1e-3); //m

	int* voxelModel = loadVoxel(voxelModelSize.x, voxelModelSize.y, voxelModelSize.z);

	double Ke[24][24];
	// ComputeKe(elementSize.x, elementSize.y, elementSize.z, 2e9, 0.394, Ke);
	getKe(Ke);
/*
#pragma region // TODO: Check this and do the decoloring on the shader accordingly
	// Preparation of grids
	// int nActiveElements = elementGridSize.multiplyComponents();
	// int nAcviteVertices = vertexGridSize.multiplyComponents();

	// std::vector<int> elementGrid;
	// std::vector<int> vertexGrid;
	// std::vector<int8_t> vertexFixed; //boolean
	// std::vector<int> vertexPosition;
	// std::vector<float> vertexEnergy;
	// std::vector<float> vertexMass;
	// std::vector<float> vertexF;
	// std::vector<float> vertexU;

	// Prepare elementGrid and vertexGrid
	Vec3i elementGridSize, vertexGridSize;

	std::vector<int> coloredVertices;

	int vert = 0;
	Vec3i vi;
	for (int color = 0; color < 8; color++)
	{
		coloredVertices.push_back(vert);
		Vec3i vi0 = Vectorize(color, 2);
		Vec3i vj0, vj1;

		for (vi.z = vi0.z; vi.z < vertexGridSize.z; vi.z += 2)
		{
			vj0.z = std::max(-1, -vi.z);
			vj1.z = std::min(1, elementGridSize.z - vi.z);

			for (vi.y = vi0.y; vi.y < vertexGridSize.y; vi.y += 2)
			{
				vj0.y = std::max(-1, -vi.y);
				vj1.y = std::min(1, elementGridSize.y - vi.y);

				for (vi.x = vi0.x; vi.x < vertexGridSize.x; vi.x += 2)
				{
					vj0.x = std::max(-1, -vi.x);
					vj1.x = std::min(1, elementGridSize.x - vi.x);
					bool isVertexRequired = false;
					Vec3i vj;
					FOR3 (vj, vj0, vj1)
					{
						if (elementGrid[Linearize(vi + vj, elementGridSize)] != -1)
						{
							isVertexRequired = true;
						}
					}
					if (isVertexRequired)
					{
						vertexGrid[Linearize(vi, vertexGridSize)] = vert;
						vertexPosition[vert] = PackPosition(vi);
						vert++;
					}
				}
			}
		}
	}
#pragma endregion
*/

	std::set<uint64_t> fixedNodes;
	std::map<uint64_t, Vec3d> loadedNodes;
	getBoundaryConditions(fixedNodes, loadedNodes);


	auto systemMatrix = assembleSystemMatrix<double>(voxelModel, voxelModelSize, Ke, fixedNodes);

	uint64_t numDoF = systemMatrix.rows();

	std::vector<double> f, u;
	f.resize(numDoF, 0.0);
	u.resize(numDoF, 0.0);

	applyBoundaryConditions(f, loadedNodes, fixedNodes);

#ifdef CPU
	solveWithCG<double>(systemMatrix, f, u);
#else

	// float dt = 0.05f;
	// int oneOverDt = (int)1/dt;

	std::vector<double> ap;
	ap.resize(numDoF, 0.0);

	std::vector<double> r(f);
	std::vector<double> p(f);

	// Setting Kompute  

	// Data for Matrix-free algorithm
	/*
	Tensor elementGridTensor = mgr.tensor((void*)elementGrid.data(), (uint64_t)nElements, sizeof(int), eInt, eStorage);
	Tensor vertexGridTensor = mgr.tensor((void*)vertexGrid.data(), (uint64_t)nVertices, sizeof(int), eInt, eStorage);
	Tensor vertexFixedTensor = mgr.tensor((void*)vertexFixed.data(), (uint64_t)nVertices, sizeof(bool), eBool, eStorage);
	Tensor positionTensor = mgr.tensor((void*)vertexPosition.data(), (uint64_t)nVertices, sizeof(int), eInt, eStorage);
	Tensor energyTensor = mgr.tensor((void*)vertexEnergy.data(), (uint64_t)nVertices, sizeof(float), eFloat, eStorage);
	Tensor massTensor = mgr.tensor((void*)vertexMass.data(), (uint64_t)nVertices, sizeof(float), eFloat, eStorage);
	Tensor KeTensor = mgr.tensor((void*)Ke, (uint64_t)(24*24), sizeof(double), eDouble, eStorage);
	Tensor fTensor = mgr.tensor((void*)vertexF.data(), (uint64_t)nVertices, sizeof(float), eFloat, eDevice);
	Tensor uTensor = mgr.tensor((void*)vertexU.data(), (uint64_t)nVertices, sizeof(float), eFloat, eDevice);

	std::vector<int> constants{elementGridSize.x, 
							elementGridSize.y,
							elementGridSize.z,
							vertexGridSize.x,
							vertexGridSize.y, 
							vertexGridSize.z,
							nVertices,
							oneOverDt};

	
	const std::vector<std::shared_ptr<kp::Tensor>> params = {elementGridTensor,
															 vertexGridTensor,
															 vertexFixedTensor,
															 positionTensor,
															 energyTensor,
															 massTensor,
															 KeTensor,
															 fTensor,
															 uTensor};

	const std::vector<uint64_t> gaussSeidelShader = std::vector<uint64_t>(
		shaders::GAUSSSEIDEL_COMP_SPV.begin(),	shaders::GAUSSSEIDEL_COMP_SPV.end());


	std::vector<int> pushConstants{0};

	int blockSize = 128;
	// int numGroups
	const kp::Workgroup workgroup({128, 128, 1});

	std::shared_ptr<kp::Algorithm> algoGaussSeidel = mgr.algorithm<int, int>(params, gaussSeidelShader, workgroup, constants, pushConstants);

	std::shared_ptr<kp::Sequence> seq = mgr.sequence();
	
	seq->record<kp::OpTensorSyncDevice>({uTensor, fTensor});
	for(int i = 0; i < 8; i++)
	{
		kp::Constants color{(float)i};
		seq->record<kp::OpAlgoDispatch>(algoGaussSeidel, color);
		seq->record<kp::OpTensorSyncLocal>({uTensor, fTensor});
	}

	for(int i = 0; i < 20; i++)
		seq->eval();
	*/

	kp::Manager mgr(0, {}, { "VK_EXT_shader_atomic_float" });

	// Data for Matrix based algorithm
	auto indices = systemMatrix.indices();
	auto values = systemMatrix.values();
	uint64_t numIndices = indices.size();

	double *d1, *n1, *n2;

	d1 = new double(0.0);
	n1 = new double(0.0);
	n2 = new double(0.0);

	for(int i = 0; i < numDoF; i++)
		*n2 += r[i] * r[i];


	Tensor systemMatrixIndices = mgr.tensor((void*)indices.data(), numIndices, sizeof(uint64_t), TensorDataTypes::eUnsignedInt, TensorTypes::eStorage);
	Tensor systemMatrixValues = mgr.tensor((void*)values.data(), numIndices, sizeof(double), TensorDataTypes::eDouble, TensorTypes::eStorage);
	Tensor pTensor = mgr.tensor((void*)p.data(), (uint64_t)numDoF, sizeof(double), TensorDataTypes::eDouble, TensorTypes::eDevice);
	Tensor apTensor = mgr.tensor((void*)ap.data(), (uint64_t)numDoF, sizeof(double), TensorDataTypes::eDouble, TensorTypes::eDevice);
	Tensor rTensor = mgr.tensor((void*)r.data(), (uint64_t)numDoF, sizeof(double), TensorDataTypes::eDouble, TensorTypes::eDevice);
	Tensor fTensor = mgr.tensor((void*)f.data(), (uint64_t)numDoF, sizeof(double), TensorDataTypes::eDouble, TensorTypes::eDevice);
	Tensor uTensor = mgr.tensor((void*)u.data(), (uint64_t)numDoF, sizeof(double), TensorDataTypes::eDouble, TensorTypes::eDevice);
	Tensor d1Tensor = mgr.tensor((void*)d1, 1, sizeof(double), TensorDataTypes::eDouble, TensorTypes::eDevice);
	Tensor n1Tensor = mgr.tensor((void*)n1, 1, sizeof(double), TensorDataTypes::eDouble, TensorTypes::eDevice);
	Tensor n2Tensor = mgr.tensor((void*)n2, 1, sizeof(double), TensorDataTypes::eDouble, TensorTypes::eDevice);
	

	const std::vector<Tensor> paramsMatxVec = {systemMatrixIndices,
											   systemMatrixValues,
											   pTensor,
											   apTensor};

	const std::vector<Tensor> paramsCG = {pTensor,
										  apTensor,
										  rTensor,
										  fTensor,
										  uTensor,
										  d1Tensor,
										  n1Tensor,
										  n2Tensor};

	const kp::Workgroup perNonZeroWorkgroup({numIndices, 1, 1});
	const kp::Workgroup perDoFWorkgroup({numDoF, 1, 1});

	const std::vector<uint32_t> MatxVecShader = std::vector<uint32_t>(
		shaders::MATXVEC_COMP_SPV.begin(),	shaders::MATXVEC_COMP_SPV.end());

	const std::vector<uint32_t> CGShader = std::vector<uint32_t>(
		shaders::CONJUGATEGRADIENT_COMP_SPV.begin(),	shaders::CONJUGATEGRADIENT_COMP_SPV.end());

	std::shared_ptr<kp::Algorithm> algoMatxVec = mgr.algorithm<uint64_t, uint64_t>(paramsMatxVec, MatxVecShader, {}, {numDoF, numIndices}, {});
	std::shared_ptr<kp::Algorithm> algoCG = mgr.algorithm(paramsCG, CGShader, perDoFWorkgroup);

	std::shared_ptr<kp::Sequence> seq = mgr.sequence();

	seq->record<kp::OpTensorSyncDevice>({uTensor, fTensor, pTensor, rTensor, n2Tensor});
	seq->record<kp::OpAlgoDispatch>(algoMatxVec);
	seq->record<kp::OpTensorSyncLocal>({pTensor, apTensor});
	// seq->record<kp::OpAlgoDispatch>(algoCG);
	// seq->record<kp::OpTensorSyncLocal>({uTensor, fTensor, pTensor, rTensor, n2Tensor, apTensor});

	// for(int i = 0; i < MAX_ITER; i++)//while(n2 > tolerance)
	// {
		seq->eval();
	// 	KP_LOG_INFO("Iteration: {}, Norm: {}", i, *n2);
	// }

	ap.assign(apTensor->data<double>(), apTensor->data<double>() + numDoF);
	// saveMatrix(systemMatrix.data, numDoF, "A");
	// saveMatrix(std::vector<std::vector<double>>({p}), "p");
	saveMatrix(std::vector<std::vector<double>>({ap}), "ap");

	delete d1, n1, n2;
#endif

	delete voxelModel;
	
	saveVector(u, "u");
	// printModel(system.data(), vertexGridSize.multiplyComponents()*vertexGridSize.multiplyComponents());

	return EXIT_SUCCESS;
}