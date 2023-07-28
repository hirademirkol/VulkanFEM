#ifndef __OUT_HPP__
#define __OUT_HPP__

#include <iostream>
#include <fstream>
#include <filesystem>

#include <string>
#include <map>

template <typename T>
void saveMatrix(std::vector<std::vector<T>> data, std::string name = "matrix");

template <typename T>
void saveVector(std::vector<T> data, std::string name = "vector");

template <typename T, size_t l>
void saveMatrix(std::vector<std::array<T,l>> data, std::string name = "matrix");

template <typename T>
void saveMatrix(std::map<uint64_t, T> data, int rows, std::string name = "matrix");

template <typename T>
void printMatrix(std::vector<std::vector<T>> data);


template void saveMatrix<double>(std::vector<std::vector<double>> data, std::string name);
template void saveMatrix<float>(std::vector<std::vector<float>> data, std::string name);

template void saveVector<double>(std::vector<double> data, std::string name);
template void saveVector<float>(std::vector<float> data, std::string name);

template void saveMatrix<double>(std::map<uint64_t, double> data, int rows, std::string name);
template void saveMatrix<float>(std::map<uint64_t, float> data, int rows, std::string name);

template void printMatrix<double>(std::vector<std::vector<double>> data);
template void printMatrix<float>(std::vector<std::vector<float>> data);

#endif // __OUT_HPP__