#ifndef __OUT_HPP__
#define __OUT_HPP__

#include <iostream>
#include <fstream>
#include <filesystem>

#include <string>
#include <map>

template <typename T>
void outMatrix(std::vector<std::vector<T>> data, std::ostream &outStream)
{
	int n = data.size();
	int m = data[0].size();
	
	if(n == 1)
	{
		outStream << "[ " << data[0][0];
		for (int j = 1; j < m-1; j++)
		{
			outStream << "\t" << data[0][j];
		}
		outStream << "\t" << data[0][m-1] << " ]" << std::endl;
		return;
	}

	outStream << "⎡ " << data[0][0];
	for (int j = 1; j < m-1; j++)
	{
		outStream << "\t" << data[0][j];
	}
	outStream << "\t" << data[0][m-1] << " ⎤" << std::endl;

	for (int i = 1; i < n-1; i++)
	{
		outStream << "⎢ " << data[i][0];
		for (int j= 1; j < m-1; j++)
		{
			outStream << "\t" << data[i][j];
		}
		outStream << "\t" << data[i][m-1] << " ⎥" << std::endl;
	}

	outStream << "⎣ " << data[n-1][0];
	for (int j= 1; j < m-1; j++)
	{
		outStream << "\t" << data[n-1][j];
	}
	outStream << "\t" << data[n-1][m-1] << " ⎦" << std::endl;
}

template <typename T, size_t lineLength>
void outMatrix(std::vector<std::array<T, lineLength>> data, std::ostream &outStream)
{
	int n = data.size();
	int m = data[0].size();
	
	if(n == 1)
	{
		outStream << "[ " << data[0][0];
		for (int j = 1; j < m-1; j++)
		{
			outStream << "\t" << data[0][j];
		}
		outStream << "\t" << data[0][m-1] << " ]" << std::endl;
		return;
	}

	outStream << "⎡ " << data[0][0];
	for (int j = 1; j < m-1; j++)
	{
		outStream << "\t" << data[0][j];
	}
	outStream << "\t" << data[0][m-1] << " ⎤" << std::endl;

	for (int i = 1; i < n-1; i++)
	{
		outStream << "⎢ " << data[i][0];
		for (int j= 1; j < m-1; j++)
		{
			outStream << "\t" << data[i][j];
		}
		outStream << "\t" << data[i][m-1] << " ⎥" << std::endl;
	}

	outStream << "⎣ " << data[n-1][0];
	for (int j= 1; j < m-1; j++)
	{
		outStream << "\t" << data[n-1][j];
	}
	outStream << "\t" << data[n-1][m-1] << " ⎦" << std::endl;
}

template <typename T>
void saveMatrix(std::vector<std::vector<T>> data, std::string name = "matrix")
{
	std::filesystem::path path = std::filesystem::current_path();
    std::filesystem::path file_path = "\\..\\data\\" + name + ".txt";
	path += file_path;
    
    std::ofstream myfile (path.string());

	outMatrix(data, myfile);

	myfile.close();
}

template <typename T, size_t l>
void saveMatrix(std::vector<std::array<T,l>> data, std::string name = "matrix")
{
	std::filesystem::path path = std::filesystem::current_path();
    std::filesystem::path file_path = "\\..\\data\\" + name + ".txt";
	path += file_path;
    
    std::ofstream myfile (path.string());

	for(size_t i = 0; i < data.size(); i++)
	{
		for(size_t j = 0; j < l; j++)
			myfile << data[i][j] << ",";
		
		myfile << std::endl;
	}

	myfile.close();
}

template <typename T>
void saveMatrix(std::map<uint64_t, T> data, int rows, std::string name = "matrix")
{
	std::filesystem::path path = std::filesystem::current_path();
    std::filesystem::path file_path = "\\..\\data\\" + name + ".csv";
	path += file_path;
    
    std::ofstream myfile (path.string());

	for(auto line = data.begin(); line != data.end(); line++)
	{
		myfile << line->first/rows <<"," << line->first%rows << ","<< line->second << std::endl;
	}

	myfile.close();
}

template <typename T>
void printMatrix(std::vector<std::vector<T>> data)
{
	outMatrix(data, std::cout);
}

template void printMatrix<double>(std::vector<std::vector<double>> data);
template void printMatrix<float>(std::vector<std::vector<float>> data);

template<typename T>
void printModel(T* model, int numElements)
{
    std::filesystem::path path = std::filesystem::current_path();
    std::filesystem::path file_path = "\\..\\data\\result.bin";
	path += file_path;
    
    std::ofstream myfile (path.string());
	for(int i = 0; i < numElements; i++)
		myfile << model[i] << std::endl;

	myfile.close();
}

#endif // __OUT_HPP__