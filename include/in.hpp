#ifndef __IN_HPP__
#define __IN_HPP__

#include <fstream>
#include <string>
#include <filesystem>

#include <set>
#include <map>

#include "utils.hpp"

/// @brief Load the voxel model from the data file
/// @param model_name (in) Name of the file
/// @param sizeX (out) Size of the model (X)
/// @param sizeY (out) Size of the model (Y)
/// @param sizeZ (out) Size of the model (Z)
/// @return Array pointer to the raw data, must be deleted manually
int* loadVoxel(std::string& model_name, int &sizeX, int &sizeY, int &sizeZ);

/// @brief Get the element stiffness matrix from the data file
/// @tparam scalar double or float
/// @param model_name (in) Name of the model element stiffness matrix belongs to, if exists. If not, data/Ke0.dat will be used.
/// @param Ke (out) Element stiffness matrix as an array
template <typename scalar>
void getKe(std::string& model_name, scalar Ke[24][24]);

/// @brief Get the boundary condiditions from the data files
/// @tparam scalar double or float
/// @param model_name (in) Name of the model boundary conditions belong to
/// @param fixedNodes (out) Set of fixed nodes in linearized form
/// @param loadedNodes (out) Map of loaded nodes in linearized form and their loads
template <typename scalar>
void getBoundaryConditions(std::string& model_name, std::set<uint64_t>& fixedNodes, std::map<uint64_t, Vec3<scalar>>& loadedNodes);

#endif // __IN_HPP__
