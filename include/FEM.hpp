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
    /// @brief Assemble the system matrix in matrix-free form and prepare necessary data for Multigrid Preconditioner
    /// @tparam scalar Only double is implemented
    /// @param voxelModel Raw voxel data
    /// @param voxelGridDimensions Dimensions of the raw voxel data
    /// @param elementStiffness Element stiffness matrix
    /// @param fixedNodes Fixed boundary condition nodes
    /// @param numLevels Number of levels to use for Multigrid Preconditioner
    /// @param skipLevels Number of levels to skip after first level for Multigrid Preconditioner
    /// @return The assembled system matrix
    template<typename scalar>
    MatrixFreeSparse<scalar> assembleSystemMatrix(int* voxelModel, Vec3i voxelGridDimensions, double elementStiffness[24][24], const std::set<uint64_t>& fixedNodes, int numLevels, int skipLevels);
    #else
    /// @brief Assemble the system matrix in matrix-free form
    /// @tparam scalar Only double is implemented
    /// @param voxelModel Raw voxel data
    /// @param voxelGridDimensions Dimensions of the raw voxel data
    /// @param elementStiffness Element stiffness matrix
    /// @param fixedNodes Fixed boundary condition nodes
    /// @return The assembled system matrix
    template<typename scalar>
    MatrixFreeSparse<scalar> assembleSystemMatrix(int* voxelModel, Vec3i voxelGridDimensions, double elementStiffness[24][24], const std::set<uint64_t>& fixedNodes);
    #endif
#else
/// @brief Assemble the system matrix in Sparse Matrix form
/// @tparam scalar double or float
/// @param voxelModel Raw voxel data
/// @param voxelGridDimensions Dimensions of the raw voxel data
/// @param elementStiffness Element stiffness matrix
/// @param fixedNodes Fixed boundary condition nodes
/// @return The assembled system matrix
template<typename scalar>
Eigen::SparseMatrix<scalar> assembleSystemMatrix(int* voxelModel, Vec3i voxelGridDimensions, scalar elementStiffness[24][24], const std::set<uint64_t>& fixedNodes);
#endif

/// @brief Apply the boundary conditions to the RHS vector
/// @tparam scalar double or float
/// @param f RHS vector
/// @param loadedNodes Loaded nodes to apply
/// @param fixedNodes Fixed nodes, only used with explicitly built matrix where fixed nodes are removed from the matrix. Assumes all fixed nodes come before the loaded nodes.
template <typename scalar>
void applyBoundaryConditions(std::vector<scalar>& f, std::map<uint64_t, Vec3<scalar>>& loadedNodes, const std::set<uint64_t>& fixedNodes);

#ifdef MATRIX_FREE
/// @brief Solve the system matrix with CG Solver from Eigen on the CPU
/// @tparam scalar Only double implemented
/// @param A (in) System matrix
/// @param b (in) RHS vector
/// @param x (out) Solution vector
template<typename scalar>
void solveWithEigen(const MatrixFreeSparse<scalar>& A, const std::vector<scalar>& b, std::vector<scalar>& x);
#else
/// @brief Solve the system matrix with CG Solver from Eigen on the CPU
/// @tparam scalar double or float
/// @param A (in) System matrix
/// @param b (in) RHS vector
/// @param x (out) Solution vector
template <typename scalar>
void solveWithEigen(const Eigen::SparseMatrix<scalar>& A, const std::vector<scalar>& b, std::vector<scalar>& x);
#endif

#endif // __FEM_HPP__
