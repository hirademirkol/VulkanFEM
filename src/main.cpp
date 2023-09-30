#include <memory>
#include <vector>
#include <chrono>

// #define CPU

#ifndef CPU
#include <kompute/Kompute.hpp>

// Compiled shaders
#include "MatxVec.hpp"
#include "FixedNodes.hpp"
#include "VecDotVec.hpp"
#include "GaussSeidel.hpp"
#include "ConjugateGradient_1.hpp"
#include "ConjugateGradient_2.hpp"

typedef std::shared_ptr<kp::Tensor> Tensor;
using TensorDataTypes = kp::Tensor::TensorDataTypes;
using TensorTypes = kp::Tensor::TensorTypes;
#endif

#include "utils.hpp"
#include "in.hpp"
#include "out.hpp"

#include "Ke.h"
#include "FEM.hpp"

int main(int argc, char* argv[])
{
	Vec3i voxelModelSize; // Number of elements on each axis
	Vec3d elementSize = Vec3d(1e-3); //m

	if(argc < 2)
	{
		std::cout << "An input file is necessary" << std::endl;
		return EXIT_FAILURE;
	}
	std::string fileName(argv[1]);
	std::string modelName = fileName.substr(0, fileName.find('.', 1));
	int* voxelModel = loadVoxel(modelName, voxelModelSize.x, voxelModelSize.y, voxelModelSize.z);

	double Ke[24][24];
	// ComputeKe(elementSize.x, elementSize.y, elementSize.z, 2e9, 0.394, Ke);
	getKe(Ke);

	std::set<uint64_t> fixedNodes;
	std::map<uint64_t, Vec3d> loadedNodes;
	getBoundaryConditions(modelName, fixedNodes, loadedNodes);

	auto systemMatrix = assembleSystemMatrix<double>(voxelModel, voxelModelSize, Ke, fixedNodes);

	uint64_t numDoF = systemMatrix.rows();

	std::vector<double> f, u;
	f.resize(numDoF, 0.0);
	u.resize(numDoF, 0.0);

	applyBoundaryConditions(f, loadedNodes, fixedNodes);

#ifdef CPU
	solveWithCG<double>(systemMatrix, f, u);
#else

	std::vector<double> r(f);
	std::vector<double> p(f);
	std::vector<double> Atp;
	Atp.resize(numDoF, 0.0);

	// Setting Kompute
	kp::Manager mgr(0, {}, { "VK_EXT_shader_atomic_float" });

	// double *norm1, *norm2, *dotP, *zero;

	double *norm1 = new double(0.0);
	double *norm2 = new double(0.0);
	double *dotP = new double(0.0);
	double *zero = new double(0.0);

	for(int i = 0; i < numDoF; i++)
		*norm2 += r[i] * r[i];

	Eigen::Array<int, 8, Eigen::Dynamic> elementToNodeArray = systemMatrix.elementToNode.transpose();

	Tensor elementStiffnessTensor = mgr.tensor((void*)systemMatrix.elementStiffnessMat.data(), 24*24, sizeof(double), TensorDataTypes::eDouble);
	Tensor elementToGlobalTensor = mgr.tensor((void*)elementToNodeArray.data(), systemMatrix.elementToNode.rows()*8, sizeof(int), TensorDataTypes::eInt);
	Tensor fixedNodesTensor = mgr.tensor((void*)systemMatrix.fixedNodes.data(), systemMatrix.fixedNodes.rows(), sizeof(int), TensorDataTypes::eInt);
	Tensor pTensor = mgr.tensor((void*)p.data(), (uint64_t)numDoF, sizeof(double), TensorDataTypes::eDouble, TensorTypes::eDevice);
	Tensor AtpTensor = mgr.tensor((void*)Atp.data(), (uint64_t)numDoF, sizeof(double), TensorDataTypes::eDouble, TensorTypes::eDevice);
	Tensor rTensor = mgr.tensor((void*)r.data(), (uint64_t)numDoF, sizeof(double), TensorDataTypes::eDouble, TensorTypes::eDevice);
	Tensor uTensor = mgr.tensor((void*)u.data(), (uint64_t)numDoF, sizeof(double), TensorDataTypes::eDouble, TensorTypes::eDevice);
	Tensor norm1Tensor = mgr.tensor((void*)norm1, 1, sizeof(double), TensorDataTypes::eDouble, TensorTypes::eDevice);
	Tensor norm2Tensor = mgr.tensor((void*)norm2, 1, sizeof(double), TensorDataTypes::eDouble, TensorTypes::eDevice);
	Tensor dotPTensor = mgr.tensor((void*)dotP, 1, sizeof(double), TensorDataTypes::eDouble, TensorTypes::eDevice);
	Tensor zeroTensor = mgr.tensor((void*)zero, 1, sizeof(double), TensorDataTypes::eDouble, TensorTypes::eDevice);


	const std::vector<Tensor> paramsMatxVec = {elementStiffnessTensor,
											   elementToGlobalTensor,
											   pTensor,
											   AtpTensor};

	const std::vector<Tensor> paramsFixedNodes = {AtpTensor,
											   fixedNodesTensor};

	const std::vector<Tensor> paramsVecDotVec = {pTensor,
										 		 AtpTensor,
												 dotPTensor};

	const std::vector<Tensor> paramsVecNorm = {rTensor,
										 		 rTensor,
												 norm2Tensor};

	const std::vector<Tensor> paramsCG_1 = {pTensor,
										  AtpTensor,
										  rTensor,
										  uTensor,
										  norm2Tensor,
										  dotPTensor};

	const std::vector<Tensor> paramsCG_2 = {pTensor,
										  rTensor,
										  norm1Tensor,
										  norm2Tensor};


	// const kp::Workgroup perElementDoFWorkgroup({(uint32_t)systemMatrix.elementToNode.rows(), 24, 24});
	const kp::Workgroup perElementDoFWorkgroup({(uint32_t)systemMatrix.elementToNode.rows(), 8, 8});
	const kp::Workgroup perFixedNodeWorkgroup({(uint32_t)systemMatrix.fixedNodes.rows(), 1, 1});
	const kp::Workgroup perDoFWorkgroup({(uint32_t)numDoF, 1, 1});

	const std::vector<uint32_t> MatxVecShader = std::vector<uint32_t>(
		shaders::MATXVEC_COMP_SPV.begin(),	shaders::MATXVEC_COMP_SPV.end());

	const std::vector<uint32_t> FixedNodesShader = std::vector<uint32_t>(
		shaders::FIXEDNODES_COMP_SPV.begin(),	shaders::FIXEDNODES_COMP_SPV.end());

	const std::vector<uint32_t> VecDotVecShader = std::vector<uint32_t>(
		shaders::VECDOTVEC_COMP_SPV.begin(),	shaders::VECDOTVEC_COMP_SPV.end());

	const std::vector<uint32_t> CGShader_1 = std::vector<uint32_t>(
		shaders::CONJUGATEGRADIENT_1_COMP_SPV.begin(),	shaders::CONJUGATEGRADIENT_1_COMP_SPV.end());

	const std::vector<uint32_t> CGShader_2 = std::vector<uint32_t>(
		shaders::CONJUGATEGRADIENT_2_COMP_SPV.begin(),	shaders::CONJUGATEGRADIENT_2_COMP_SPV.end());

	std::shared_ptr<kp::Algorithm> algoMatxVec = mgr.algorithm(paramsMatxVec, MatxVecShader, perElementDoFWorkgroup);
	std::shared_ptr<kp::Algorithm> algoFixedNodes = mgr.algorithm(paramsFixedNodes, FixedNodesShader, perFixedNodeWorkgroup);
	std::shared_ptr<kp::Algorithm> algoVecDotVec = mgr.algorithm(paramsVecDotVec, VecDotVecShader, perDoFWorkgroup);
	std::shared_ptr<kp::Algorithm> algoVecNorm = mgr.algorithm(paramsVecNorm, VecDotVecShader, perDoFWorkgroup);
	std::shared_ptr<kp::Algorithm> algoCG_1 = mgr.algorithm(paramsCG_1, CGShader_1, perDoFWorkgroup);
	std::shared_ptr<kp::Algorithm> algoCG_2 = mgr.algorithm(paramsCG_2, CGShader_2, perDoFWorkgroup);

	mgr.sequence()->eval<kp::OpTensorSyncDevice>(paramsMatxVec);
	mgr.sequence()->eval<kp::OpTensorSyncDevice>(paramsCG_1);
	mgr.sequence()->eval<kp::OpTensorSyncDevice>(paramsCG_2);

	std::shared_ptr<kp::Sequence> seq1 =
		mgr.sequence()->record<kp::OpTensorSyncDevice>({AtpTensor, dotPTensor})
					  ->record<kp::OpAlgoDispatch>(algoMatxVec)
					  ->record<kp::OpAlgoDispatch>(algoFixedNodes)
					  ->record<kp::OpAlgoDispatch>(algoVecDotVec)
					  ->record<kp::OpAlgoDispatch>(algoCG_1);

	std::shared_ptr<kp::Sequence> seq2 =
		mgr.sequence()->record<kp::OpTensorSyncDevice>({norm1Tensor, norm2Tensor})
					  ->record<kp::OpAlgoDispatch>(algoVecNorm)
				  	  ->record<kp::OpTensorSyncLocal>({norm2Tensor});

	std::shared_ptr<kp::Sequence> seq3 =
		mgr.sequence()->record<kp::OpAlgoDispatch>(algoCG_2);


#ifdef MAX_ITER
	int maxIterations = MAX_ITER;
#else
	int maxIterations = 10000;
#endif

#ifdef TOLERANCE
	double tolerance = TOLERANCE;
#else
	double tolerance = 1e-16;
#endif

	double threshold = tolerance * tolerance * *norm2;

	std::cout << "Starting the GPU Solver" << std::endl;

	int iter;
	auto start = std::chrono::system_clock::now();
	for(iter = 0; iter < maxIterations; iter++)
	{
		memset((void*)AtpTensor->data<double>(), 0, Atp.size()*sizeof(double));
		*dotPTensor->data<double>() = 0.0;

		seq1->eval();

		*norm1Tensor->data<double>() = *norm2Tensor->data<double>();
		*norm2Tensor->data<double>() = 0.0;

		seq2->eval();

		if(iter%100 == 0)
			std::cout << "Iteration: "<< iter <<", Norm: " << (*norm2Tensor->data<double>()) << std::endl;

		if(*norm2Tensor->data<double>() < threshold)
		{
			break;
		}

		seq3->eval();
	}
	auto end = std::chrono::system_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
	
	std::cout << "Solving took:    " << duration.count() << " s" 					<< std::endl;
	std::cout << "#iterations:     " << iter	    								<< std::endl;
	std::cout << "Estimated error: " << sqrt(*norm2Tensor->data<double>() / *norm2)	<< std::endl;

	mgr.sequence()->eval<kp::OpTensorSyncLocal>({uTensor});
	u = uTensor->vector<double>();
	delete norm2, dotP;

#endif

	delete voxelModel;

	std::cout << "Saving the solution" << std::endl;
	saveVector(u, modelName);

	return EXIT_SUCCESS;
}