#include "utils.hpp"

#include <utility>
#include <map>
#include <set>
#include <unordered_map>
#include <array>
#include <iostream>

#include <math.h>

#include <Eigen/Sparse>

#define MAX_ITER 200
#define index(i,j,n) (i)*(n)+(j)
#define get_symmetric(A,i,j) (i) <= (j) ? A[(i)][(j)] : A[(j)][(i)]

template<typename scalar>
Eigen::SparseMatrix<scalar> assembleSystemMatrix(int* voxelModel, Vec3i voxelGridDimensions, scalar elementStiffness[24][24], const std::set<uint64_t>& fixedNodes);

template <typename scalar>
void applyBoundaryConditions(std::vector<scalar>& f, std::map<uint64_t, Vec3<scalar>>& loadedNodes);

template <typename scalar>
void solveWithCG(const Eigen::SparseMatrix<scalar>& A, const std::vector<scalar>& b, std::vector<scalar>& x);

template Eigen::SparseMatrix<double> assembleSystemMatrix<double>(int* voxelModel, Vec3i voxelGridDimensions, double elementStiffness[24][24], const std::set<uint64_t>& fixedNodes);
template Eigen::SparseMatrix<float> assembleSystemMatrix<float>(int* voxelModel, Vec3i voxelGridDimensions, float elementStiffness[24][24], const std::set<uint64_t>& fixedNodes);

template void applyBoundaryConditions<double>(std::vector<double>& f, std::map<uint64_t, Vec3<double>>& loadedNodes);
template void applyBoundaryConditions<float>(std::vector<float>& f, std::map<uint64_t, Vec3<float>>& loadedNodes);

template void solveWithCG<double>(const Eigen::SparseMatrix<double>& A, const std::vector<double>& b, std::vector<double>& x);
template void solveWithCG<float>(const Eigen::SparseMatrix<float>& A, const std::vector<float>& b, std::vector<float>& x);