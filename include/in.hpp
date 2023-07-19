#ifndef __IN_HPP__
#define __IN_HPP__

#include <fstream>
#include <string>
#include <filesystem>

int* loadVoxel(int &sizeX, int &sizeY, int &sizeZ)
{
    std::string line;
    std::filesystem::path path = std::filesystem::current_path();
    std::filesystem::path file_path = "\\..\\data\\voxelizedModel.bin";
    path += file_path;
    
    std::ifstream myfile (path.string());

    std::getline(myfile, line, ' ');
    sizeX = std::stoi(line);
    std::getline(myfile, line, ' ');
    sizeY = std::stoi(line);
    std::getline(myfile, line);
    sizeZ = std::stoi(line);

    int* model = new int[sizeX*sizeY*sizeZ];

    // for (int i = 0; i++; i < sizeX)
    // {
    //     for (int j = 0; j++; j < sizeY)
    //     {
    //         for (int k = 0; k++; k < sizeZ)
    //         {
    //             std::getline(myfile, line);

    //         }
    //     }
    // }

    int i = 0;
    while(std::getline(myfile, line))
    {
        model[i++] = std::stoi(line);
    }

    myfile.close();

    return model;
}

template <typename scalar>
void getKe(scalar Ke[24][24])
{
    std::string line;
    std::filesystem::path path = std::filesystem::current_path();
    std::filesystem::path file_path = "\\..\\data\\Ke.txt";
    path += file_path;
    
    std::ifstream myfile (path.string());

    for(int i = 0; i < 24; i++)
    {
        for(int j = 0; j < 23; j++)
        {
            std::getline(myfile, line, ' ');            
            Ke[i][j] = (scalar)std::stod(line);
        }   
        std::getline(myfile, line);
        Ke[i][23] = (scalar)std::stod(line);
    }
    
    myfile.close();
}

template void getKe<double>(double Ke[24][24]);
template void getKe<float>(float Ke[24][24]);

template <typename scalar>
void getBoundaryConditions(std::vector<scalar>& f, std::vector<uint64_t>& fixedNodes)
{
    std::string line;
    
    //fixing conditions
    std::filesystem::path path = std::filesystem::current_path();
    std::filesystem::path file_path = "\\..\\data\\fixed.dat";
    path += file_path;

    std::ifstream myfile (path.string());

    std::getline(myfile, line);
    fixedNodes.resize(std::stoi(line));

    for(int i = 0; i < fixedNodes.size(); i++)
    {
        std::getline(myfile, line);
        fixedNodes[i] = std::stoull(line) - 1;
        // f[fixedNodes[i]] = 0.0;
    }

    myfile.close();

    //loading conditions
    path = std::filesystem::current_path();
    file_path = "\\..\\data\\loading.dat";
    path += file_path;

    myfile.open(path.string());

    std::getline(myfile, line);
    int numLoad = std::stoi(line);

    for (int i = 0; i < numLoad; i++)
    {
        std::getline(myfile, line, '\t');
        uint64_t ind = std::stoull(line) - 1;

        std::getline(myfile, line, '\t');
        f[ind*3] = std::stod(line);
        
        std::getline(myfile, line, '\t');
        f[ind*3 + 1] = std::stod(line);

        std::getline(myfile, line);
        f[ind*3 + 2] = std::stod(line);
    }

    myfile.close();
}

template void getBoundaryConditions<double>(std::vector<double>& f, std::vector<uint64_t>& fixedNodes);
template void getBoundaryConditions<float>(std::vector<float>& f, std::vector<uint64_t>& fixedNodes);

#endif // __IN_HPP__
