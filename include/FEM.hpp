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

#include "FEMDefines.hpp"

#define index(i,j,n) (i)*(n)+(j)
#define get_symmetric(A,i,j) (i) <= (j) ? A[(i)][(j)] : A[(j)][(i)]

#ifdef MATRIX_FREE
    #ifdef MULTIGRID
    template<typename scalar>
    MatrixFreeSparse<scalar> assembleSystemMatrix(int* voxelModel, Vec3i voxelGridDimensions, double elementStiffness[24][24], const std::set<uint64_t>& fixedNodes, int numLevels, int skipLevels);
    #else
    template<typename scalar>
    MatrixFreeSparse<scalar> assembleSystemMatrix(int* voxelModel, Vec3i voxelGridDimensions, double elementStiffness[24][24], const std::set<uint64_t>& fixedNodes);
    #endif
#else
template<typename scalar>
Eigen::SparseMatrix<scalar> assembleSystemMatrix(int* voxelModel, Vec3i voxelGridDimensions, scalar elementStiffness[24][24], const std::set<uint64_t>& fixedNodes);
#endif

template <typename scalar>
void applyBoundaryConditions(std::vector<scalar>& f, std::map<uint64_t, Vec3<scalar>>& loadedNodes, const std::set<uint64_t>& fixedNodes);


#ifdef MATRIX_FREE
template<typename scalar>
void solveWithEigen(const MatrixFreeSparse<scalar>& A, const std::vector<scalar>& b, std::vector<scalar>& x);
#else
template <typename scalar>
void solveWithEigen(const Eigen::SparseMatrix<scalar>& A, const std::vector<scalar>& b, std::vector<scalar>& x);
#endif

#endif // __FEM_HPP__
