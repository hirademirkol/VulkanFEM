#ifndef __IN_HPP__
#define __IN_HPP__

#include <fstream>
#include <string>
#include <filesystem>

#include <set>
#include <map>

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
void getBoundaryConditions(std::set<uint64_t>& fixedNodes, std::map<uint64_t, Vec3<scalar>>& loadedNodes)
{
    std::string line;
    
    //fixing conditions
    std::filesystem::path path = std::filesystem::current_path();
    std::filesystem::path file_path = "\\..\\data\\fixed.dat";
    path += file_path;

    std::ifstream myfile (path.string());

    std::getline(myfile, line);
    int numFixed = std::stoi(line);

    for(int i = 0; i < numFixed; i++)
    {
        std::getline(myfile, line);
        fixedNodes.insert(std::stoull(line) - 1);
    }

    myfile.close();

    //loading conditions
    path = std::filesystem::current_path();
    file_path = "\\..\\data\\loading.dat";
    path += file_path;

    myfile.open(path.string());

    std::getline(myfile, line);
    int numLoads = std::stoi(line);

    for (int i = 0; i < numLoads; i++)
    {
        std::getline(myfile, line, '\t');
        uint64_t ind = std::stoull(line) - 1 - numFixed;

        std::getline(myfile, line, '\t');
        scalar x = std::stod(line);
        
        std::getline(myfile, line, '\t');
        scalar y = std::stod(line);

        std::getline(myfile, line);
        scalar z = std::stod(line);

        loadedNodes[ind] = Vec3<scalar>(x,y,z);
    }

    myfile.close();
}

// template void getBoundaryConditions<double>(std::vector<double>& f, std::vector<uint64_t>& fixedNodes);
// template void getBoundaryConditions<float>(std::vector<float>& f, std::vector<uint64_t>& fixedNodes);

#endif // __IN_HPP__
