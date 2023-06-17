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
void GetKe(scalar Ke[24][24])
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
            std::getline(myfile, line, ' ');            Ke[i][j] = (scalar)std::stod(line);
        }   
        std::getline(myfile, line);
        Ke[i][23] = (scalar)std::stod(line);
    }
    
    myfile.close();
}

template void GetKe<double>(double Ke[24][24]);
template void GetKe<float>(float Ke[24][24]);

#endif // __IN_HPP__
