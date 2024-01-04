#include "FEMGPU.hpp"
#include "RestrictionOperator.hpp"

#include <iostream>
#include <chrono>

#include <kompute/Kompute.hpp>

// Compiled shaders
#include "MatxVec.hpp"
#include "FixedNodes.hpp"
#include "VecDotVec.hpp"
#include "ConjugateGradient_1.hpp"
#include "ConjugateGradient_2.hpp"

#ifdef MULTIGRID
#include "CWiseMult.hpp"
#include "PreSmoothAndAXPY.hpp"
#include "Smooth.hpp"
#include "Restrict.hpp"
#include "Interpolate.hpp"
#include "RestrictDouble.hpp"
#include "InterpolateDouble.hpp"
#endif

typedef std::shared_ptr<kp::Tensor> Tensor;
typedef std::shared_ptr<kp::Sequence> Sequence;
typedef std::shared_ptr<kp::Algorithm> Algorithm;
using TensorDataTypes = kp::Tensor::TensorDataTypes;
using TensorTypes = kp::Tensor::TensorTypes;

template< typename scalar>
void solveWithKompute(const MatrixFreeSparse<scalar>& systemMatrix, const std::vector<scalar>& f, std::vector<scalar>& u)
{
	std::cout << "Initializing GPU Solver" << std::endl;
    uint64_t numDoF = f.size();

	uint64_t memoryUsage = 0;

    std::vector<scalar> r(f);
	std::vector<scalar> p(f);
	std::vector<scalar> Atp;
	Atp.resize(numDoF, 0.0);

	// Setting up Kompute
	kp::Manager mgr(0, {}, { "VK_EXT_shader_atomic_float" });

	scalar *norm1 = new scalar(0.0);
	scalar *norm2 = new scalar(0.0);
	scalar *dotP = new scalar(0.0);

	for(int i = 0; i < numDoF; i++)
		*norm2 += r[i] * r[i];

	// A transpose of the non-symmetric matrices is required, since Kompute assigns the data in column-major format
	Eigen::Array<int, 4, Eigen::Dynamic> elementToNodeArray = systemMatrix.elementToNode.transpose();

	// Creating the data tensors(buffers) with initial data from the prepared Eigen matrices
	Tensor elementStiffnessTensor = mgr.tensor((void*)systemMatrix.elementStiffnessMat.data(), 24*24, sizeof(scalar), TensorDataTypes::eDouble);
	Tensor elementToGlobalTensor = mgr.tensor((void*)elementToNodeArray.data(), systemMatrix.elementToNode.rows()*4, sizeof(int), TensorDataTypes::eInt);
	Tensor fixedNodesTensor = mgr.tensor((void*)systemMatrix.fixedNodes.data(), systemMatrix.fixedNodes.rows(), sizeof(int), TensorDataTypes::eInt);
	Tensor pTensor = mgr.tensor((void*)p.data(), (uint64_t)numDoF, sizeof(scalar), TensorDataTypes::eDouble, TensorTypes::eDevice);
	Tensor tempTensor = mgr.tensor((void*)Atp.data(), (uint64_t)numDoF, sizeof(scalar), TensorDataTypes::eDouble, TensorTypes::eDevice);
	Tensor rTensor = mgr.tensor((void*)r.data(), (uint64_t)numDoF, sizeof(scalar), TensorDataTypes::eDouble, TensorTypes::eDevice);
	Tensor uTensor = mgr.tensor((void*)u.data(), (uint64_t)numDoF, sizeof(scalar), TensorDataTypes::eDouble, TensorTypes::eDevice);
	Tensor norm1Tensor = mgr.tensor((void*)norm1, 1, sizeof(scalar), TensorDataTypes::eDouble, TensorTypes::eDevice);
	Tensor norm2Tensor = mgr.tensor((void*)norm2, 1, sizeof(scalar), TensorDataTypes::eDouble, TensorTypes::eDevice);
	Tensor dotPTensor = mgr.tensor((void*)dotP, 1, sizeof(scalar), TensorDataTypes::eDouble, TensorTypes::eDevice);

	// Memory usage of the tensors are noted down
	memoryUsage += elementStiffnessTensor->size() * sizeof(scalar);
	memoryUsage += elementToGlobalTensor->size() * sizeof(int);
	memoryUsage += fixedNodesTensor->size() * sizeof(int);
	memoryUsage += pTensor->size() * sizeof(scalar);
	memoryUsage += tempTensor->size() * sizeof(scalar);
	memoryUsage += rTensor->size() * sizeof(scalar);
	memoryUsage += uTensor->size() * sizeof(scalar);
	memoryUsage += norm1Tensor->size() * sizeof(scalar);
	memoryUsage += norm2Tensor->size() * sizeof(scalar);
	memoryUsage += dotPTensor->size() * sizeof(scalar);

	// Parameter sets for each algorithm, must be in the order written in the shader
	const std::vector<Tensor> paramsMatxVec = {elementStiffnessTensor,
											   elementToGlobalTensor,
											   pTensor,
											   tempTensor};

	const std::vector<Tensor> paramsFixedNodes = {tempTensor,
											   fixedNodesTensor};

	const std::vector<Tensor> paramsVecDotVec = {pTensor,
										 		 tempTensor,
												 dotPTensor};

	const std::vector<Tensor> paramsVecNorm = {rTensor,
										 		 rTensor,
												 norm2Tensor};

	const std::vector<Tensor> paramsCG_1 = {pTensor,
										  tempTensor,
										  rTensor,
										  uTensor,
										  norm2Tensor,
										  dotPTensor};

	// Workgroup sizes for each different type dispatches
	const kp::Workgroup perElementWorkgroup({(uint32_t)systemMatrix.elementToNode.rows(), 1, 1});
	const kp::Workgroup perFixedNodeWorkgroup({(uint32_t)systemMatrix.fixedNodes.rows()/64 + 1, 1, 1});
	const kp::Workgroup perDoFWorkgroup({(uint32_t)numDoF/64 + 1, 1, 1});

	// Precompiled shaders are read from the header files included above
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

	// Algorithms are created with corresponding parameters, shaders and workgroup sizes
	Algorithm algoMatxVec = mgr.algorithm(paramsMatxVec, MatxVecShader, perElementWorkgroup);
	Algorithm algoFixedNodes = mgr.algorithm(paramsFixedNodes, FixedNodesShader, perFixedNodeWorkgroup);
	Algorithm algoVecDotVec = mgr.algorithm(paramsVecDotVec, VecDotVecShader, perDoFWorkgroup);
	Algorithm algoVecNorm = mgr.algorithm(paramsVecNorm, VecDotVecShader, perDoFWorkgroup);
	Algorithm algoCG_1 = mgr.algorithm(paramsCG_1, CGShader_1, perDoFWorkgroup);

	// All parameter tensors' data needs to be synced to the device
	mgr.sequence()->eval<kp::OpTensorSyncDevice>(paramsMatxVec);
	mgr.sequence()->eval<kp::OpTensorSyncDevice>(paramsFixedNodes);
	mgr.sequence()->eval<kp::OpTensorSyncDevice>(paramsCG_1);

	// Sequences(sets of ordered algorithms) are created
	Sequence seqXUpdate =
		mgr.sequence()->record<kp::OpTensorSyncDevice>({tempTensor, dotPTensor})
					  ->record<kp::OpAlgoDispatch>(algoMatxVec)
					  ->record<kp::OpAlgoDispatch>(algoFixedNodes)
					  ->record<kp::OpAlgoDispatch>(algoVecDotVec)
					  ->record<kp::OpAlgoDispatch>(algoCG_1);

	Sequence seqNormCalc =
		mgr.sequence()->record<kp::OpTensorSyncDevice>({norm1Tensor, norm2Tensor})
					  ->record<kp::OpAlgoDispatch>(algoVecNorm)
				  	  ->record<kp::OpTensorSyncLocal>({norm2Tensor});

#ifndef MULTIGRID

	// Second part of the CG for non-preconditied version
	const std::vector<Tensor> paramsCG_2 = {pTensor,
											rTensor,
											norm1Tensor,
											norm2Tensor};

	mgr.sequence()->eval<kp::OpTensorSyncDevice>(paramsCG_2);

	Algorithm algoCG_2 = mgr.algorithm(paramsCG_2, CGShader_2, perDoFWorkgroup);
	
	Sequence seqPUpdate =
		mgr.sequence()->record<kp::OpAlgoDispatch>(algoCG_2);


#else

	// Shaders for the Multigrid Preconditioning
	const std::vector<uint32_t> CWiseMultShader = std::vector<uint32_t>(
		shaders::CWISEMULT_COMP_SPV.begin(),	shaders::CWISEMULT_COMP_SPV.end());
		
	const std::vector<uint32_t> PreSmoothAndAxpyShader = std::vector<uint32_t>(
		shaders::PRESMOOTHANDAXPY_COMP_SPV.begin(),	shaders::PRESMOOTHANDAXPY_COMP_SPV.end());

	const std::vector<uint32_t> SmoothShader = std::vector<uint32_t>(
		shaders::SMOOTH_COMP_SPV.begin(),	shaders::SMOOTH_COMP_SPV.end());
		
	const std::vector<uint32_t> RestrictShader = std::vector<uint32_t>(
		shaders::RESTRICT_COMP_SPV.begin(),	shaders::RESTRICT_COMP_SPV.end());
		
	const std::vector<uint32_t> InterpolateShader = std::vector<uint32_t>(
		shaders::INTERPOLATE_COMP_SPV.begin(),	shaders::INTERPOLATE_COMP_SPV.end());

	const std::vector<uint32_t> RestrictDoubleShader = std::vector<uint32_t>(
		shaders::RESTRICTDOUBLE_COMP_SPV.begin(),	shaders::RESTRICTDOUBLE_COMP_SPV.end());
		
	const std::vector<uint32_t> InterpolateDoubleShader = std::vector<uint32_t>(
		shaders::INTERPOLATEDOUBLE_COMP_SPV.begin(),	shaders::INTERPOLATEDOUBLE_COMP_SPV.end());

	// Vectors for keeping tensors and sequences on each Multigrid level
	std::vector<Tensor> elementToNodeTensors;
	std::vector<Tensor> restrictionMappingTensors;
	std::vector<Tensor> restrictionCoefficientTensors;
	std::vector<Tensor> invDiagKOnLevelTensors;
	std::vector<Tensor> rTensors;
	std::vector<Tensor> tempTensors;
	std::vector<std::vector<double>> rs;
	std::vector<std::vector<double>> temps;

	std::vector<Sequence> seqRestricts;
	std::vector<Sequence> seqInterpolates;

	// Tensors for keeping the restriction operators
	auto restrictionOperatorTensor = mgr.tensor((void*)restrictionOperatorValues.data(), 27*8, sizeof(double), TensorDataTypes::eDouble);
	auto restrictionOperator4Tensor = mgr.tensor((void*)restrictionOperator4Values.data(), 125*8, sizeof(double), TensorDataTypes::eDouble);

	// These tensors are recycled during the multigrid section, hence also added here
	rTensors.push_back(rTensor);
	tempTensors.push_back(tempTensor);

	int numLevels = systemMatrix.numLevels;

	// For each level, necessary tensors and sequences are created
	for(int i = 1; i < numLevels; i++)
	{
		Eigen::Array<int, 4, Eigen::Dynamic> matrix = systemMatrix.elementToNodeMatrices[i - 1].transpose();
		auto elementToNodeTensor = mgr.tensor((void*)matrix.data(), matrix.cols()*4, sizeof(int), TensorDataTypes::eInt);
		elementToNodeTensors.push_back(elementToNodeTensor);

		Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic> matrix2 = systemMatrix.restrictionMappings[i - 1].transpose();
		auto restrictionMappingTensor = mgr.tensor((void*)matrix2.data(), matrix2.cols()*matrix2.rows(), sizeof(int), TensorDataTypes::eInt);
		restrictionMappingTensors.push_back(restrictionMappingTensor);

		// Restriction coefficients are converted to float for memory optimization
		std::vector<float> dataVector(systemMatrix.restrictionCoefficients[i - 1].rows());
		for(int j = 0; j < systemMatrix.restrictionCoefficients[i - 1].rows(); j++)
		{
			dataVector[j] = (float)systemMatrix.restrictionCoefficients[i - 1](j);
		}
		auto restrictionCoefficientTensor = mgr.tensor((void*)dataVector.data(), dataVector.size(), sizeof(float), TensorDataTypes::eFloat);
		restrictionCoefficientTensors.push_back(restrictionCoefficientTensor);

		// Inverse diagonal coefficients are converted to float for memory optimization
		std::vector<float> dataVector2(systemMatrix.invDiagKOnLevels[i - 1].rows());
		for(int j = 0; j < systemMatrix.invDiagKOnLevels[i - 1].rows(); j++)
		{
			dataVector2[j] = (float)systemMatrix.invDiagKOnLevels[i - 1](j);
		}
		auto invDiagKOnLevelTensor = mgr.tensor((void*)dataVector2.data(), dataVector2.size(), sizeof(float), TensorDataTypes::eFloat);
		invDiagKOnLevelTensors.push_back(invDiagKOnLevelTensor);

		std::vector<double> r(systemMatrix.invDiagKOnLevels[i].rows());
		rs.push_back(r);
		auto rTensorOnLevelBelow = mgr.tensor((void*)r.data(), systemMatrix.invDiagKOnLevels[i].rows(), sizeof(double), TensorDataTypes::eDouble);
		rTensors.push_back(rTensorOnLevelBelow);

		std::vector<double> temp(systemMatrix.invDiagKOnLevels[i].rows());
		temps.push_back(temp);
		auto tempTensorOnLevelBelow = mgr.tensor((void*)temp.data(), systemMatrix.invDiagKOnLevels[i].rows(), sizeof(double), TensorDataTypes::eDouble);
		tempTensors.push_back(tempTensorOnLevelBelow);

		// Parameter sets for each algorithm
		std::vector<Tensor> paramsRestrict;
		std::vector<Tensor> paramsInterpolate;
		if(i == 1 && systemMatrix.skipLevels == 1)
		{
			paramsRestrict = { restrictionOperator4Tensor,
							   elementToNodeTensor,
							   restrictionMappingTensor,
							   tempTensors[i - 1],
							   rTensorOnLevelBelow};

			paramsInterpolate = { restrictionOperator4Tensor,
								  elementToNodeTensor,
								  restrictionMappingTensor,
								  tempTensorOnLevelBelow,
								  tempTensors[i - 1]};
		}
		else
		{
			paramsRestrict = { restrictionOperatorTensor,
							   elementToNodeTensor,
							   restrictionMappingTensor,
							   tempTensors[i - 1],
							   rTensorOnLevelBelow};

			paramsInterpolate = { restrictionOperatorTensor,
								  elementToNodeTensor,
								  restrictionMappingTensor,
								  tempTensorOnLevelBelow,
								  tempTensors[i - 1]};
		}

		const std::vector<Tensor> paramsSmooth = { tempTensors[i - 1],
													rTensors[i - 1],
													invDiagKOnLevelTensors[i - 1]};

		const std::vector<Tensor> paramsCWiseMult = { rTensors[i - 1],
													restrictionCoefficientTensor,
													tempTensors[i - 1]};

		const std::vector<Tensor> paramsPreSmoothAndAxpy = { tempTensors[i - 1],
															 rTensors[i - 1],
															 invDiagKOnLevelTensors[i - 1],
															 restrictionCoefficientTensors[i - 1]};

		// Workgroup sizes for each different types of algorithms
		const kp::Workgroup perElementOnLevelWorkgroup({(uint32_t)matrix.cols(), 1, 1});
		const kp::Workgroup perDoFOnLevelWorkgroup({(uint32_t)systemMatrix.invDiagKOnLevels[i - 1].rows()/64 + 1, 1, 1});

		// Creating algorithms and sequences for restriction and interpolation on the level
		Algorithm algoRestrict;
		Algorithm algoInterpolate;

		if(i == 1 && systemMatrix.skipLevels == 1)
		{
			algoRestrict = mgr.algorithm(paramsRestrict, RestrictDoubleShader, perElementOnLevelWorkgroup);
			algoInterpolate = mgr.algorithm(paramsInterpolate, InterpolateDoubleShader, perElementOnLevelWorkgroup);
		}
		else
		{
			algoRestrict = mgr.algorithm(paramsRestrict, RestrictShader, perElementOnLevelWorkgroup);
			algoInterpolate = mgr.algorithm(paramsInterpolate, InterpolateShader, perElementOnLevelWorkgroup);
		}

		Algorithm algoSmooth = mgr.algorithm(paramsSmooth, SmoothShader, perDoFOnLevelWorkgroup);
		Algorithm algoCWiseMult = mgr.algorithm(paramsCWiseMult, CWiseMultShader, perDoFOnLevelWorkgroup);
		Algorithm algoPreSmoothAndAxpy = mgr.algorithm(paramsPreSmoothAndAxpy, PreSmoothAndAxpyShader, perDoFOnLevelWorkgroup);

		std::vector<Tensor> syncTensors = {rTensorOnLevelBelow};
		auto seqRestrict =
			mgr.sequence()->record<kp::OpTensorSyncDevice>(syncTensors)
						  ->record<kp::OpAlgoDispatch>(algoCWiseMult)
						  ->record<kp::OpAlgoDispatch>(algoRestrict);

		auto seqInterpolate = 
			mgr.sequence()->record<kp::OpTensorSyncDevice>({tempTensors[i - 1]})
						  ->record<kp::OpAlgoDispatch>(algoInterpolate)
						  ->record<kp::OpAlgoDispatch>(algoPreSmoothAndAxpy)
						  ->record<kp::OpAlgoDispatch>(algoSmooth);

		seqRestricts.push_back(seqRestrict);
		seqInterpolates.push_back(seqInterpolate);
	}

	// Syncing the data of each tensor related to Multigrid
	mgr.sequence()->eval<kp::OpTensorSyncDevice>(elementToNodeTensors);
	mgr.sequence()->eval<kp::OpTensorSyncDevice>(restrictionMappingTensors);
	mgr.sequence()->eval<kp::OpTensorSyncDevice>(restrictionCoefficientTensors);
	mgr.sequence()->eval<kp::OpTensorSyncDevice>(invDiagKOnLevelTensors);
	mgr.sequence()->eval<kp::OpTensorSyncDevice>(rTensors);
	mgr.sequence()->eval<kp::OpTensorSyncDevice>(tempTensors);
	mgr.sequence()->eval<kp::OpTensorSyncDevice>({restrictionOperatorTensor, restrictionOperator4Tensor});

	// Writing down the memory usage for each tensor
	for( auto tensor : elementToNodeTensors)
		memoryUsage += tensor->size()*sizeof(int);

	for( auto tensor : restrictionMappingTensors)
		memoryUsage += tensor->size()*sizeof(int);

	for( auto tensor : restrictionCoefficientTensors)
		memoryUsage += tensor->size()*sizeof(float);

	for( auto tensor : invDiagKOnLevelTensors)
		memoryUsage += tensor->size()*sizeof(float);

	for( auto tensor : rTensors)
		memoryUsage += tensor->size()*sizeof(scalar);

	memoryUsage -= rTensor->size()*sizeof(scalar); // Double counted above

	for( auto tensor : tempTensors)
		memoryUsage += tensor->size()*sizeof(scalar);

	memoryUsage -= tempTensor->size()*sizeof(scalar); // Double counted above

	memoryUsage += restrictionOperatorTensor->size()*sizeof(scalar);
	memoryUsage += restrictionOperator4Tensor->size()*sizeof(scalar);

	const std::vector<Tensor> paramsFixedNodesZ = { tempTensors[0],
											  fixedNodesTensor};

	const std::vector<Tensor> paramsRdotZ = {rTensor,
											 tempTensors[0],
											 norm2Tensor};

	const std::vector<Tensor> paramsCG_2Z = {pTensor,
										  tempTensors[0],
										  norm1Tensor,
										  norm2Tensor};


	Algorithm algoFixedNodesZ = mgr.algorithm(paramsFixedNodesZ, FixedNodesShader, perFixedNodeWorkgroup);
	Algorithm algoRdotZ = mgr.algorithm(paramsRdotZ, VecDotVecShader, perDoFWorkgroup);
	Algorithm algoCG_2Z = mgr.algorithm(paramsCG_2Z, CGShader_2, perDoFWorkgroup);

	// Second part of the CG with preconditioning
	Sequence seqPUpdateWithZ =
		mgr.sequence()->record<kp::OpTensorSyncDevice>({norm2Tensor})
					  ->record<kp::OpAlgoDispatch>(algoFixedNodesZ)
					  ->record<kp::OpAlgoDispatch>(algoRdotZ)
					  ->record<kp::OpAlgoDispatch>(algoCG_2Z);

	// CPU Cholesky solver for the coarsest level
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldltSolver(systemMatrix.Kc);

#endif

	std::cout << "\tData assigned GPU memory: " << (double)memoryUsage/1024/1024 << " MB" << std::endl;

#ifdef MAX_ITER
	int maxIterations = MAX_ITER;
#else
	int maxIterations = 30000;
#endif

#ifdef TOLERANCE
	scalar tolerance = TOLERANCE;
#else
	scalar tolerance = 1e-16;
#endif

	scalar threshold = tolerance * tolerance * *norm2;

	int iter;
	std::cout << "Starting the GPU Solver" << std::endl;
	auto start = std::chrono::system_clock::now();

#ifdef MULTIGRID
	// Preconditioning for first update direction
	{
		// Go down the V-Cycle
		for(int i = 0; i < numLevels - 1; i++)
		{
			memset((void*)rTensors[i + 1]->data<scalar>(), 0, systemMatrix.invDiagKOnLevels[i + 1].rows()*sizeof(scalar));
			seqRestricts[i]->eval();
		}

		// Sync coarsest level residual to the CPU
		mgr.sequence()->eval<kp::OpTensorSyncLocal>({rTensors[numLevels - 1]});

		// Create Eigen Vector from residual data and empty result data
		Eigen::Matrix<scalar, Eigen::Dynamic, 1> vector1 = Eigen::Map<Eigen::Matrix<scalar, Eigen::Dynamic, 1>> (rTensors[numLevels - 1]->data<scalar>(), rTensors[numLevels - 1]->vector<scalar>().size());
		Eigen::Matrix<scalar, Eigen::Dynamic, 1> vector2 = Eigen::Matrix<scalar, Eigen::Dynamic, 1>::Zero(vector1.rows());

		// Solve the system for free nodes, put the solution into result vector
		Eigen::Matrix<scalar, Eigen::Dynamic, 1> solution = ldltSolver.solve(vector1(systemMatrix.coarseFreeDoFs));
		vector2(systemMatrix.coarseFreeDoFs) = solution;

		// Copy the result vector into the corresponding tensor and sync it to the GPU
		memcpy((void*)tempTensors[numLevels - 1]->data<scalar>(), (void*)vector2.data(), vector2.rows()*sizeof(scalar));
		mgr.sequence()->eval<kp::OpTensorSyncDevice>({tempTensors[numLevels - 1]});

		// Go up the V-Cycle
		for (int i = numLevels - 1; i > 0; i--)
		{
			// Zero the memory of tempTensors since they are used before, it will be synced to GPU during the sequence
			memset((void*)tempTensors[i - 1]->data<scalar>(), 0, systemMatrix.invDiagKOnLevels[i - 1].rows()*sizeof(scalar));
			seqInterpolates[i - 1]->eval();
		}

		// Calculation of the norm for CG
		*norm2Tensor->data<scalar>() = 0.0;
		mgr.sequence()->record<kp::OpTensorSyncDevice>({norm2Tensor})
					  ->record<kp::OpAlgoDispatch>(algoFixedNodesZ)
					  ->record<kp::OpAlgoDispatch>(algoRdotZ)
					  ->record<kp::OpTensorSyncLocal>({norm2Tensor})
					  ->eval();

		// Resulting vector is copied into pTensor to be used as initial search direction
		mgr.sequence()->eval<kp::OpTensorSyncLocal>({tempTensors[0]});
		memcpy((void*)pTensor->data<scalar>(), (void*)tempTensors[0]->data<scalar>(), tempTensors[0]->size()*sizeof(scalar));
		mgr.sequence()->eval<kp::OpTensorSyncDevice>({pTensor});		
	}
#endif

	// Solving cycle
	for(iter = 0; iter < maxIterations; iter++)
	{
		// Zero the memory of tempTensor and dotPTensor since it is be used before, it will be synced to GPU during the sequence
		memset((void*)tempTensor->data<scalar>(), 0, Atp.size()*sizeof(scalar));
		*dotPTensor->data<scalar>() = 0.0;

		// Dispatch for the first part of CG
		seqXUpdate->eval();

		// Norms are exchanged
		*norm1Tensor->data<scalar>() = *norm2Tensor->data<scalar>();
		*norm2Tensor->data<scalar>() = 0.0;

		// New norm calculation
		seqNormCalc->eval();

		// Print the status every 100th iteration
		if(iter%100 == 0)
			std::cout << "Iteration: "<< iter <<", Norm: " << (*norm2Tensor->data<scalar>()) << std::endl;

		// Conclude the solving if the error is small enough
		if(*norm2Tensor->data<scalar>() < threshold)
		{
			break;
		}

#ifdef MULTIGRID
		//Preconditioning

		// Go down the V-Cycle
		for(int i = 0; i < numLevels - 1; i++)
		{
			// Zero the memory of rTensors since they are used before, it will be synced to GPU during the sequence
			memset((void*)rTensors[i + 1]->data<scalar>(), 0, systemMatrix.invDiagKOnLevels[i + 1].rows()*sizeof(scalar));
			seqRestricts[i]->eval();
		}

		// Sync coarsest level residual to the CPU
		mgr.sequence()->eval<kp::OpTensorSyncLocal>({rTensors[numLevels - 1]});

		// Create Eigen Vector from residual data and empty result data
		Eigen::Matrix<scalar, Eigen::Dynamic, 1> vector1 = Eigen::Map<Eigen::Matrix<scalar, Eigen::Dynamic, 1>> (rTensors[numLevels - 1]->data<scalar>(), rTensors[numLevels - 1]->vector<scalar>().size());
		Eigen::Matrix<scalar, Eigen::Dynamic, 1> vector2 = Eigen::Matrix<scalar, Eigen::Dynamic, 1>::Zero(vector1.rows());

		// Solve the system for free nodes, put the solution into result vector
		Eigen::Matrix<scalar, Eigen::Dynamic, 1> solution = ldltSolver.solve(vector1(systemMatrix.coarseFreeDoFs));
		vector2(systemMatrix.coarseFreeDoFs) = solution;

		// Copy the result vector into the corresponding tensor and sync it to the GPU
		memcpy((void*)tempTensors[numLevels - 1]->data<scalar>(), (void*)vector2.data(), vector2.rows()*sizeof(scalar));
		mgr.sequence()->eval<kp::OpTensorSyncDevice>({tempTensors[numLevels - 1]});

		*norm2Tensor->data<scalar>() = 0.0;

		// Go up the V-Cycle
		for (int i = numLevels - 1; i > 0; i--)
		{
			// Zero the memory of tempTensors since they are used before, it will be synced to GPU during the sequence
			memset((void*)tempTensors[i - 1]->data<scalar>(), 0, systemMatrix.invDiagKOnLevels[i - 1].rows()*sizeof(scalar));
			seqInterpolates[i - 1]->eval();
		}

		// Dispatch for the second part of CG
		seqPUpdateWithZ->eval();
		mgr.sequence()->eval<kp::OpTensorSyncLocal>({norm2Tensor});

#else
		// Dispatch for the second part of CG without preconditioning
		seqPUpdate->eval();
#endif
	}
	auto end = std::chrono::system_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	
	std::cout << "Solving took:    " << (float)duration.count() / 1000 << " s" 		<< std::endl;
	std::cout << "#iterations:     " << iter	    								<< std::endl;
	std::cout << "Estimated error: " << sqrt(*norm2Tensor->data<scalar>() / *norm2)	<< std::endl;

	// Sync the result to the CPU
	mgr.sequence()->eval<kp::OpTensorSyncLocal>({uTensor});
	u = uTensor->vector<scalar>();

	delete norm1, norm2, dotP;
}

template void solveWithKompute<double>(const MatrixFreeSparse<double>& A, const std::vector<double>& b, std::vector<double>& x);