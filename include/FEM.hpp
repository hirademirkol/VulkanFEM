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

// #define MAX_ITER 1000
#define TOLERANCE 1e-6

#define MATRIX_FREE

#ifdef MATRIX_FREE
    #include "MatrixFreeSparse.hpp"

    #define MULTIGRID

    #ifdef MULTIGRID
        #define NUM_LEVELS 3
    #endif

#endif

#define index(i,j,n) (i)*(n)+(j)
#define get_symmetric(A,i,j) (i) <= (j) ? A[(i)][(j)] : A[(j)][(i)]

#ifdef MATRIX_FREE
template<typename scalar>
MatrixFreeSparse assembleSystemMatrix(int* voxelModel, Vec3i voxelGridDimensions, double elementStiffness[24][24], const std::set<uint64_t>& fixedNodes);
#else
template<typename scalar>
Eigen::SparseMatrix<scalar> assembleSystemMatrix(int* voxelModel, Vec3i voxelGridDimensions, scalar elementStiffness[24][24], const std::set<uint64_t>& fixedNodes);
#endif

template <typename scalar>
void applyBoundaryConditions(std::vector<scalar>& f, std::map<uint64_t, Vec3<scalar>>& loadedNodes, const std::set<uint64_t>& fixedNodes);


#ifdef MATRIX_FREE
template<typename scalar>
void solveWithCG(const MatrixFreeSparse& A, const std::vector<scalar>& b, std::vector<scalar>& x);
#else
template <typename scalar>
void solveWithCG(const Eigen::SparseMatrix<scalar>& A, const std::vector<scalar>& b, std::vector<scalar>& x);
#endif

#endif // __FEM_HPP__
