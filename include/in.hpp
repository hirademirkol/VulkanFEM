#ifndef __IN_HPP__
#define __IN_HPP__

#include <fstream>
#include <string>
#include <filesystem>

#include <set>
#include <map>

#include "utils.hpp"

int* loadVoxel(int &sizeX, int &sizeY, int &sizeZ);

template <typename scalar>
void getKe(scalar Ke[24][24]);

template <typename scalar>
void getBoundaryConditions(std::set<uint64_t>& fixedNodes, std::map<uint64_t, Vec3<scalar>>& loadedNodes);


template void getKe<double>(double Ke[24][24]);
template void getKe<float>(float Ke[24][24]);

template void getBoundaryConditions<double>(std::set<uint64_t>& fixedNodes, std::map<uint64_t, Vec3<double>>& loadedNodes);
template void getBoundaryConditions<float>(std::set<uint64_t>& fixedNodes, std::map<uint64_t, Vec3<float>>& loadedNodes);

#endif // __IN_HPP__
