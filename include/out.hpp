#ifndef __OUT_HPP__
#define __OUT_HPP__

#include <iostream>
#include <fstream>
#include <filesystem>

#include <string>
#include <map>
#include <vector>

#include <Eigen/Sparse>

/// @brief Save matrix to a file. Matrix should be created as a vector of vectors.
/// @tparam T double or float
/// @param data Matrix data
/// @param name Name of the file, will be appended an "_out.csv"
template <typename T>
void saveMatrix(std::vector<std::vector<T>> data, std::string name = "matrix");

/// @brief Print matrix. Matrix should be created as a vector of vectors.
/// @tparam T double or float
/// @param data Matrix data
template <typename T>
void printMatrix(std::vector<std::vector<T>> data);

/// @brief Save vector to a file.
/// @tparam T double or float
/// @param data Vector to save
/// @param name Name of the file, will be appended an "_out.csv"
template <typename T>
void saveVector(std::vector<T> data, std::string name = "vector");

/// @brief Save matrix to a file. Matrix should be created as a vector of arrays.
/// @tparam T double or float
/// @param data Matrix data
/// @param name Name of the file, will be appended an "_out.csv"
template <typename T, size_t l>
void saveMatrix(std::vector<std::array<T,l>> data, std::string name = "matrix");

/// @brief Save matrix to a file. Matrix should be created as a map of (row*cols + col) to coefficient, and it should be a square matrix.
/// @tparam T double or float
/// @param data Matrix data
/// @param rows Number of rows or columns
/// @param name Name of the file, will be appended an "_out.csv"
template <typename T>
void saveMatrix(std::map<uint64_t, T> data, int rows, std::string name = "matrix");

/// @brief Save matrix to a file. Matrix should be created as a vector of Triplets from the Eigen package.
/// @tparam T double or float
/// @param triplets Matrix data
/// @param name Name of the file, will be appended an "_out.csv"
template <typename T>
void saveMatrix(std::vector<Eigen::Triplet<T>> triplets, std::string name);

/// @brief Save array to a file. Matrix should be created as an Array from the Eigen package.
/// @tparam T double or float
/// @tparam i Only -1 for undefined number of rows
/// @tparam j -1 for undefined number of rows, 8 or 27
/// @param array Array to save
/// @param name Name of the file, will be appended an "_out.csv"
template <typename T, int i, int j>
void saveArray(Eigen::Array<T, i, j> array, std::string name);

#endif // __OUT_HPP__