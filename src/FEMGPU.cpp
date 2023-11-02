#include "FEMGPU.hpp"
#include "RestrictionOperator.hpp"

#include <iostream>
#include <chrono>

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

	// Setting Kompute
	kp::Manager mgr(0, {}, { "VK_EXT_shader_atomic_float" });

	scalar *norm1 = new scalar(0.0);
	scalar *norm2 = new scalar(0.0);
	scalar *dotP = new scalar(0.0);

	for(int i = 0; i < numDoF; i++)
		*norm2 += r[i] * r[i];

	Eigen::Array<int, 4, Eigen::Dynamic> elementToNodeArray = systemMatrix.elementToNode.transpose();

	Tensor elementStiffnessTensor = mgr.tensor((void*)systemMatrix.elementStiffnessMat.data(), 24*24, sizeof(scalar), TensorDataTypes::eDouble);
	Tensor elementToGlobalTensor = mgr.tensor((void*)elementToNodeArray.data(), systemMatrix.elementToNode.rows()*4, sizeof(int), TensorDataTypes::eInt);
	Tensor fixedNodesTensor = mgr.tensor((void*)systemMatrix.fixedNodes.data(), systemMatrix.fixedNodes.rows(), sizeof(int), TensorDataTypes::eInt);
	Tensor pTensor = mgr.tensor((void*)p.data(), (uint64_t)numDoF, sizeof(scalar), TensorDataTypes::eDouble, TensorTypes::eDevice);
	Tensor AtpTensor = mgr.tensor((void*)Atp.data(), (uint64_t)numDoF, sizeof(scalar), TensorDataTypes::eDouble, TensorTypes::eDevice);
	Tensor rTensor = mgr.tensor((void*)r.data(), (uint64_t)numDoF, sizeof(scalar), TensorDataTypes::eDouble, TensorTypes::eDevice);
	Tensor uTensor = mgr.tensor((void*)u.data(), (uint64_t)numDoF, sizeof(scalar), TensorDataTypes::eDouble, TensorTypes::eDevice);
	Tensor norm1Tensor = mgr.tensor((void*)norm1, 1, sizeof(scalar), TensorDataTypes::eDouble, TensorTypes::eDevice);
	Tensor norm2Tensor = mgr.tensor((void*)norm2, 1, sizeof(scalar), TensorDataTypes::eDouble, TensorTypes::eDevice);
	Tensor dotPTensor = mgr.tensor((void*)dotP, 1, sizeof(scalar), TensorDataTypes::eDouble, TensorTypes::eDevice);

	memoryUsage += elementStiffnessTensor->size() * sizeof(scalar);
	memoryUsage += elementToGlobalTensor->size() * sizeof(int);
	memoryUsage += fixedNodesTensor->size() * sizeof(int);
	memoryUsage += pTensor->size() * sizeof(scalar);
	memoryUsage += AtpTensor->size() * sizeof(scalar);
	memoryUsage += rTensor->size() * sizeof(scalar);
	memoryUsage += uTensor->size() * sizeof(scalar);
	memoryUsage += norm1Tensor->size() * sizeof(scalar);
	memoryUsage += norm2Tensor->size() * sizeof(scalar);
	memoryUsage += dotPTensor->size() * sizeof(scalar);

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

	const kp::Workgroup perElementWorkgroup({(uint32_t)systemMatrix.elementToNode.rows(), 1, 1});
	const kp::Workgroup perFixedNodeWorkgroup({(uint32_t)systemMatrix.fixedNodes.rows()/64 + 1, 1, 1});
	const kp::Workgroup perDoFWorkgroup({(uint32_t)numDoF/64 + 1, 1, 1});

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

	Algorithm algoMatxVec = mgr.algorithm(paramsMatxVec, MatxVecShader, perElementWorkgroup);
	Algorithm algoFixedNodes = mgr.algorithm(paramsFixedNodes, FixedNodesShader, perFixedNodeWorkgroup);
	Algorithm algoVecDotVec = mgr.algorithm(paramsVecDotVec, VecDotVecShader, perDoFWorkgroup);
	Algorithm algoVecNorm = mgr.algorithm(paramsVecNorm, VecDotVecShader, perDoFWorkgroup);
	Algorithm algoCG_1 = mgr.algorithm(paramsCG_1, CGShader_1, perDoFWorkgroup);
	Algorithm algoCG_2 = mgr.algorithm(paramsCG_2, CGShader_2, perDoFWorkgroup);

	mgr.sequence()->eval<kp::OpTensorSyncDevice>(paramsMatxVec);
	mgr.sequence()->eval<kp::OpTensorSyncDevice>(paramsFixedNodes);
	mgr.sequence()->eval<kp::OpTensorSyncDevice>(paramsCG_1);
	mgr.sequence()->eval<kp::OpTensorSyncDevice>(paramsCG_2);

	Sequence seqXUpdate =
		mgr.sequence()->record<kp::OpTensorSyncDevice>({AtpTensor, dotPTensor})
					  ->record<kp::OpAlgoDispatch>(algoMatxVec)
					  ->record<kp::OpAlgoDispatch>(algoFixedNodes)
					  ->record<kp::OpAlgoDispatch>(algoVecDotVec)
					  ->record<kp::OpAlgoDispatch>(algoCG_1);

	Sequence seqNormCalc =
		mgr.sequence()->record<kp::OpTensorSyncDevice>({norm1Tensor, norm2Tensor})
					  ->record<kp::OpAlgoDispatch>(algoVecNorm)
				  	  ->record<kp::OpTensorSyncLocal>({norm2Tensor});

	Sequence seqPUpdate =
		mgr.sequence()->record<kp::OpAlgoDispatch>(algoCG_2);


#ifdef MULTIGRID

	const std::vector<uint32_t> CWiseMultShader = std::vector<uint32_t>(
		shaders::CWISEMULT_COMP_SPV.begin(),	shaders::CWISEMULT_COMP_SPV.end());
		
	const std::vector<uint32_t> CWiseAddShader = std::vector<uint32_t>(
		shaders::CWISEADD_COMP_SPV.begin(),	shaders::CWISEADD_COMP_SPV.end());

	const std::vector<uint32_t> SmoothShader = std::vector<uint32_t>(
		shaders::SMOOTH_COMP_SPV.begin(),	shaders::SMOOTH_COMP_SPV.end());
		
	const std::vector<uint32_t> RestrictShader = std::vector<uint32_t>(
		shaders::RESTRICT_COMP_SPV.begin(),	shaders::RESTRICT_COMP_SPV.end());
		
	const std::vector<uint32_t> InterpolateShader = std::vector<uint32_t>(
		shaders::INTERPOLATE_COMP_SPV.begin(),	shaders::INTERPOLATE_COMP_SPV.end());


	std::vector<Tensor> elementToNodeTensors;
	std::vector<Tensor> restrictionMappingTensors;
	std::vector<Tensor> restrictionCoefficientTensors;
	std::vector<Tensor> invDiagKOnLevelTensors;
	std::vector<Tensor> rTensors;
	std::vector<Tensor> tempTensors;
	std::vector<Tensor> resultTensors;
	std::vector<std::vector<double>> rs;
	std::vector<std::vector<double>> temps;
	std::vector<std::vector<double>> results;

	std::vector<Sequence> seqRestricts;
	std::vector<Sequence> seqInterpolates;

	auto restrictionOperatorTensor = mgr.tensor((void*)restrictionOperatorValues.data(), 27*8, sizeof(double), TensorDataTypes::eDouble);

	rTensors.push_back(rTensor);

	std::vector<double> temp(systemMatrix.invDiagKOnLevels[0].rows());
	temps.push_back(temp);
	auto tempTensor = mgr.tensor((void*)temp.data(), systemMatrix.invDiagKOnLevels[0].rows(), sizeof(double), TensorDataTypes::eDouble);
	tempTensors.push_back(tempTensor);

	std::vector<double> result(systemMatrix.invDiagKOnLevels[0].rows());
	results.push_back(result);
	auto resultTensor = mgr.tensor((void*)result.data(), systemMatrix.invDiagKOnLevels[0].rows(), sizeof(double), TensorDataTypes::eDouble);
	resultTensors.push_back(resultTensor);

	int numLevels = systemMatrix.numLevels;

	for(int i = 1; i < numLevels; i++)
	{
		Eigen::Array<int, 4, Eigen::Dynamic> matrix = systemMatrix.elementToNodeMatrices[i - 1].transpose();
		auto elementToNodeTensor = mgr.tensor((void*)matrix.data(), matrix.cols()*4, sizeof(int), TensorDataTypes::eInt);
		elementToNodeTensors.push_back(elementToNodeTensor);

		Eigen::Array<int, 27, Eigen::Dynamic> matrix2 = systemMatrix.restrictionMappings[i - 1].transpose();
		auto restrictionMappingTensor = mgr.tensor((void*)matrix2.data(), matrix2.cols()*27, sizeof(int), TensorDataTypes::eInt);
		restrictionMappingTensors.push_back(restrictionMappingTensor);

		auto restrictionCoefficientTensor = mgr.tensor((void*)systemMatrix.restrictionCoefficients[i - 1].data(), systemMatrix.restrictionCoefficients[i - 1].rows(), sizeof(double), TensorDataTypes::eDouble);
		restrictionCoefficientTensors.push_back(restrictionCoefficientTensor);

		auto invDiagKOnLevelTensor = mgr.tensor((void*)systemMatrix.invDiagKOnLevels[i - 1].data(), systemMatrix.invDiagKOnLevels[i - 1].rows(), sizeof(double), TensorDataTypes::eDouble);
		invDiagKOnLevelTensors.push_back(invDiagKOnLevelTensor);

		std::vector<double> r(systemMatrix.invDiagKOnLevels[i].rows());
		rs.push_back(r);
		auto rTensorOnLevelBelow = mgr.tensor((void*)r.data(), systemMatrix.invDiagKOnLevels[i].rows(), sizeof(double), TensorDataTypes::eDouble);
		rTensors.push_back(rTensorOnLevelBelow);

		std::vector<double> temp(systemMatrix.invDiagKOnLevels[i].rows());
		temps.push_back(temp);
		auto tempTensorOnLevelBelow = mgr.tensor((void*)temp.data(), systemMatrix.invDiagKOnLevels[i].rows(), sizeof(double), TensorDataTypes::eDouble);
		tempTensors.push_back(tempTensorOnLevelBelow);

		std::vector<double> result(systemMatrix.invDiagKOnLevels[i].rows());
		results.push_back(result);
		auto resultTensorOnLevelBelow = mgr.tensor((void*)result.data(), systemMatrix.invDiagKOnLevels[i].rows(), sizeof(double), TensorDataTypes::eDouble);
		resultTensors.push_back(resultTensorOnLevelBelow);

		const std::vector<Tensor> paramsRestrict = { restrictionOperatorTensor,
													elementToNodeTensor,
													restrictionMappingTensor,
													tempTensors[i - 1],
													rTensorOnLevelBelow};

		const std::vector<Tensor> paramsInterpolate = { restrictionOperatorTensor,
														elementToNodeTensor,
														restrictionMappingTensor,
														resultTensorOnLevelBelow,
														tempTensors[i - 1]};

		const std::vector<Tensor> paramsSmooth = { resultTensors[i - 1],
													rTensors[i - 1],
													invDiagKOnLevelTensors[i - 1]};

		const std::vector<Tensor> paramsCWiseMult = { rTensors[i - 1],
													restrictionCoefficientTensor,
													tempTensors[i - 1]};

		const std::vector<Tensor> paramsCWiseMult2 = { tempTensors[i - 1],
													restrictionCoefficientTensor,
													tempTensors[i - 1]};

		const std::vector<Tensor> paramsCWiseAdd = { tempTensors[i - 1],
													resultTensors[i - 1],
													resultTensors[i - 1]};
												  
		const kp::Workgroup perElementOnLevelWorkgroup({(uint32_t)matrix.cols(), 1, 1});
		const kp::Workgroup perDoFOnLevelWorkgroup({(uint32_t)systemMatrix.invDiagKOnLevels[i - 1].rows()/64 + 1, 1, 1});

		Algorithm algoRestrict = mgr.algorithm(paramsRestrict, RestrictShader, perElementOnLevelWorkgroup);
		Algorithm algoInterpolate = mgr.algorithm(paramsInterpolate, InterpolateShader, perElementOnLevelWorkgroup);
		Algorithm algoSmooth = mgr.algorithm(paramsSmooth, SmoothShader, perDoFOnLevelWorkgroup);
		Algorithm algoCWiseMult = mgr.algorithm(paramsCWiseMult, CWiseMultShader, perDoFOnLevelWorkgroup);
		Algorithm algoCWiseMult2 = mgr.algorithm(paramsCWiseMult2, CWiseMultShader, perDoFOnLevelWorkgroup);
		Algorithm algoCWiseAdd = mgr.algorithm(paramsCWiseAdd, CWiseAddShader, perDoFOnLevelWorkgroup);

		std::vector<Tensor> syncTensors = {resultTensors[i - 1], rTensorOnLevelBelow};
		auto seqRestrict = 
			mgr.sequence()->record<kp::OpTensorSyncDevice>(syncTensors)
						  ->record<kp::OpAlgoDispatch>(algoSmooth)
						  ->record<kp::OpAlgoDispatch>(algoCWiseMult)
						  ->record<kp::OpAlgoDispatch>(algoRestrict);

						  
		auto seqInterpolate = 
			mgr.sequence()->record<kp::OpTensorSyncDevice>({tempTensors[i - 1]})
						  ->record<kp::OpAlgoDispatch>(algoInterpolate)
						  ->record<kp::OpAlgoDispatch>(algoCWiseMult2)
						  ->record<kp::OpAlgoDispatch>(algoCWiseAdd)
						  ->record<kp::OpAlgoDispatch>(algoSmooth);

		seqRestricts.push_back(seqRestrict);
		seqInterpolates.push_back(seqInterpolate);
	}
	
	mgr.sequence()->eval<kp::OpTensorSyncDevice>(elementToNodeTensors);
	mgr.sequence()->eval<kp::OpTensorSyncDevice>(restrictionMappingTensors);
	mgr.sequence()->eval<kp::OpTensorSyncDevice>(restrictionCoefficientTensors);
	mgr.sequence()->eval<kp::OpTensorSyncDevice>(invDiagKOnLevelTensors);
	mgr.sequence()->eval<kp::OpTensorSyncDevice>(rTensors);
	mgr.sequence()->eval<kp::OpTensorSyncDevice>(tempTensors);
	mgr.sequence()->eval<kp::OpTensorSyncDevice>(resultTensors);
	mgr.sequence()->eval<kp::OpTensorSyncDevice>({restrictionOperatorTensor});

	for( auto tensor : elementToNodeTensors)
		memoryUsage += tensor->size()*sizeof(int);

	for( auto tensor : restrictionMappingTensors)
		memoryUsage += tensor->size()*sizeof(int);

	for( auto tensor : restrictionCoefficientTensors)
		memoryUsage += tensor->size()*sizeof(scalar);

	for( auto tensor : invDiagKOnLevelTensors)
		memoryUsage += tensor->size()*sizeof(scalar);

	for( auto tensor : rTensors)
		memoryUsage += tensor->size()*sizeof(scalar);

	for( auto tensor : tempTensors)
		memoryUsage += tensor->size()*sizeof(scalar);

	for( auto tensor : resultTensors)
		memoryUsage += tensor->size()*sizeof(scalar);

	memoryUsage += restrictionOperatorTensor->size()*sizeof(scalar);
	memoryUsage -= rTensor->size()*sizeof(scalar); // Double counted above

	const std::vector<Tensor> paramsFixedNodesZ = { resultTensors[0],
											  fixedNodesTensor};

	const std::vector<Tensor> paramsRdotZ = {rTensor,
											 resultTensors[0],
											 norm2Tensor};

	const std::vector<Tensor> paramsCG_2Z = {pTensor,
										  resultTensors[0],
										  norm1Tensor,
										  norm2Tensor};


	Algorithm algoFixedNodesZ = mgr.algorithm(paramsFixedNodesZ, FixedNodesShader, perFixedNodeWorkgroup);
	Algorithm algoRdotZ = mgr.algorithm(paramsRdotZ, VecDotVecShader, perDoFWorkgroup);
	Algorithm algoCG_2Z = mgr.algorithm(paramsCG_2Z, CGShader_2, perDoFWorkgroup);

	Sequence seqPUpdateWithZ =
		mgr.sequence()->record<kp::OpTensorSyncDevice>({norm2Tensor})
					  ->record<kp::OpAlgoDispatch>(algoFixedNodesZ)
					  ->record<kp::OpAlgoDispatch>(algoRdotZ)
					  ->record<kp::OpAlgoDispatch>(algoCG_2Z);

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
	//Preconditioning for first update direction
	{
		for(int i = 0; i < numLevels - 1; i++)
		{
			memset((void*)resultTensors[i]->data<scalar>(), 0, systemMatrix.invDiagKOnLevels[i].rows()*sizeof(scalar));
			memset((void*)rTensors[i + 1]->data<scalar>(), 0, systemMatrix.invDiagKOnLevels[i + 1].rows()*sizeof(scalar));
			seqRestricts[i]->eval();
		}

		mgr.sequence()->eval<kp::OpTensorSyncLocal>({rTensors[numLevels - 1]});

		Eigen::Matrix<scalar, Eigen::Dynamic, 1> vector1 = Eigen::Map<Eigen::Matrix<scalar, Eigen::Dynamic, 1>> (rTensors[numLevels - 1]->data<scalar>(), rTensors[numLevels - 1]->vector<scalar>().size());
		Eigen::Matrix<scalar, Eigen::Dynamic, 1> vector2 = Eigen::Matrix<scalar, Eigen::Dynamic, 1>::Zero(vector1.rows());

		Eigen::Matrix<scalar, Eigen::Dynamic, 1> solution = ldltSolver.solve(vector1(systemMatrix.coarseFreeDoFs));
		vector2(systemMatrix.coarseFreeDoFs) = solution;

		memcpy((void*)resultTensors[numLevels - 1]->data<scalar>(), (void*)vector2.data(), vector2.rows()*sizeof(scalar));

		mgr.sequence()->eval<kp::OpTensorSyncDevice>({resultTensors[numLevels - 1]});

		for (int i = numLevels - 1; i > 0; i--)
		{
			memset((void*)tempTensors[i - 1]->data<scalar>(), 0, systemMatrix.invDiagKOnLevels[i - 1].rows()*sizeof(scalar));
			seqInterpolates[i - 1]->eval();
		}

		*norm2Tensor->data<scalar>() = 0.0;
		mgr.sequence()->record<kp::OpTensorSyncDevice>({norm2Tensor})
					  ->record<kp::OpAlgoDispatch>(algoFixedNodesZ)
					  ->record<kp::OpAlgoDispatch>(algoRdotZ)
					  ->record<kp::OpTensorSyncLocal>({norm2Tensor})
					  ->eval();

		mgr.sequence()->eval<kp::OpTensorSyncLocal>({resultTensors[0]});
		memcpy((void*)pTensor->data<scalar>(), (void*)resultTensors[0]->data<scalar>(), resultTensors[0]->size()*sizeof(scalar));
		mgr.sequence()->eval<kp::OpTensorSyncDevice>({pTensor});		
	}
#endif

	for(iter = 0; iter < maxIterations; iter++)
	{
		memset((void*)AtpTensor->data<scalar>(), 0, Atp.size()*sizeof(scalar));
		*dotPTensor->data<scalar>() = 0.0;

		seqXUpdate->eval();

		*norm1Tensor->data<scalar>() = *norm2Tensor->data<scalar>();
		*norm2Tensor->data<scalar>() = 0.0;

		seqNormCalc->eval();

		if(iter%100 == 0)
			std::cout << "Iteration: "<< iter <<", Norm: " << (*norm2Tensor->data<scalar>()) << std::endl;

		if(*norm2Tensor->data<scalar>() < threshold)
		{
			break;
		}
#ifdef MULTIGRID

		for(int i = 0; i < numLevels - 1; i++)
		{
			memset((void*)resultTensors[i]->data<scalar>(), 0, systemMatrix.invDiagKOnLevels[i].rows()*sizeof(scalar));
			memset((void*)rTensors[i + 1]->data<scalar>(), 0, systemMatrix.invDiagKOnLevels[i + 1].rows()*sizeof(scalar));
			seqRestricts[i]->eval();
		}

		mgr.sequence()->eval<kp::OpTensorSyncLocal>({rTensors[numLevels - 1]});

		Eigen::Matrix<scalar, Eigen::Dynamic, 1> vector1 = Eigen::Map<Eigen::Matrix<scalar, Eigen::Dynamic, 1>> (rTensors[numLevels - 1]->data<scalar>(), rTensors[numLevels - 1]->vector<scalar>().size());
		Eigen::Matrix<scalar, Eigen::Dynamic, 1> vector2 = Eigen::Matrix<scalar, Eigen::Dynamic, 1>::Zero(vector1.rows());

		Eigen::Matrix<scalar, Eigen::Dynamic, 1> solution = ldltSolver.solve(vector1(systemMatrix.coarseFreeDoFs));
		vector2(systemMatrix.coarseFreeDoFs) = solution;

		memcpy((void*)resultTensors[numLevels - 1]->data<scalar>(), (void*)vector2.data(), vector2.rows()*sizeof(scalar));

		mgr.sequence()->eval<kp::OpTensorSyncDevice>({resultTensors[numLevels - 1]});

		*norm2Tensor->data<scalar>() = 0.0;

		for (int i = numLevels - 1; i > 0; i--)
		{
			memset((void*)tempTensors[i - 1]->data<scalar>(), 0, systemMatrix.invDiagKOnLevels[i - 1].rows()*sizeof(scalar));
			seqInterpolates[i - 1]->eval();
		}

		seqPUpdateWithZ->eval();
		mgr.sequence()->eval<kp::OpTensorSyncLocal>({norm2Tensor});

#else
		seqPUpdate->eval();
#endif
	}
	auto end = std::chrono::system_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	
	std::cout << "Solving took:    " << (float)duration.count() / 1000 << " s" 		<< std::endl;
	std::cout << "#iterations:     " << iter	    								<< std::endl;
	std::cout << "Estimated error: " << sqrt(*norm2Tensor->data<scalar>() / *norm2)	<< std::endl;

	mgr.sequence()->eval<kp::OpTensorSyncLocal>({uTensor});
	u = uTensor->vector<scalar>();

	delete norm1, norm2, dotP;
}

template void solveWithKompute<double>(const MatrixFreeSparse<double>& A, const std::vector<double>& b, std::vector<double>& x);