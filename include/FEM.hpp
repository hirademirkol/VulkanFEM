#ifndef __FEM_HPP__
#define __FEM_HPP__

#include "utils.hpp"

#include <utility>
#include <map>
#include <set>
#include <unordered_map>
#include <array>
#include <iostream>

#include <math.h>

#include <Eigen/Sparse>

// #define MAX_ITER 100
#define index(i,j,n) (i)*(n)+(j)
#define get_symmetric(A,i,j) (i) <= (j) ? A[(i)][(j)] : A[(j)][(i)]

template<typename scalar>
Eigen::SparseMatrix<scalar> assembleSystemMatrix(int* voxelModel, Vec3i voxelGridDimensions, scalar elementStiffness[24][24], const std::set<uint64_t>& fixedNodes);

template <typename scalar>
void applyBoundaryConditions(std::vector<scalar>& f, std::map<uint64_t, Vec3<scalar>>& loadedNodes);

template <typename scalar>
void solveWithCG(const Eigen::SparseMatrix<scalar>& A, const std::vector<scalar>& b, std::vector<scalar>& x);

#endif // __FEM_HPP__
