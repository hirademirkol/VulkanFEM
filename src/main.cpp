#include <memory>
#include <vector>

#include "out.hpp"

#include "MatxVec.hpp"
#include "VecxVec.hpp"
#include <kompute/Kompute.hpp>

int main()
{
	kp::Manager mgr(0, {}, {"VK_EXT_shader_atomic_float"});

	int32_t n = 2;
	int32_t m = 3;

	std::shared_ptr<kp::TensorT<float>> A = mgr.tensor({2.0, 4.0, 6.0, 2.0, 4.0, 1.0});
	std::shared_ptr<kp::TensorT<float>> b = mgr.tensor({0.0, 1.0, 2.0});
	std::shared_ptr<kp::TensorT<float>> tensorOut = mgr.tensor({0.0, 0.0});

	const std::vector<std::shared_ptr<kp::Tensor>> params = {A,
															 b,
															 tensorOut};

	const std::vector<uint32_t> matxVecShader = std::vector<uint32_t>(
		shaders::MATRIXOPS_MATXVEC_SPV.begin(), shaders::MATRIXOPS_MATXVEC_SPV.end());

		
	const std::vector<uint32_t> vecxVecShader = std::vector<uint32_t>(
		shaders::VECXVEC_SPV.begin(), shaders::VECXVEC_SPV.end());

	std::shared_ptr<kp::Algorithm> algoMatxVec = mgr.algorithm(params, matxVecShader, {}, std::vector<int32_t>({n, m}), {});
	std::shared_ptr<kp::Algorithm> algoVecxVec = mgr.algorithm(params, vecxVecShader);

	mgr.sequence()
		->record<kp::OpTensorSyncDevice>(params)
		->record<kp::OpAlgoDispatch>(algoMatxVec)
		->record<kp::OpTensorSyncLocal>(params)
		->record<kp::OpAlgoDispatch>(algoVecxVec)
		->record<kp::OpTensorSyncLocal>(params)
		->eval();

	// prints A
	// prints "A:
	// 		   { 2 4 6 }
	// 		   { 2 4 1 }"
	std::cout << "A:" << std::endl;
	printMatrix(A->vector(), n, m);

	// prints b
	// prints "b:
	// 		   { 0 }
	// 		   { 1 }
	// 		   { 2 }"
	std::cout << "b:" << std::endl;
	printMatrix(b->vector(), m, 1);

	// prints A.b
	// prints "Output:
	// 		   { 256 } 
	// 		   { 36 }"
	std::cout << "Output:" << std::endl;
	printMatrix(tensorOut->vector(), n, 1);
}