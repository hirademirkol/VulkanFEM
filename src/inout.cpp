#include "in.hpp"
#include "out.hpp"

#pragma region input


int* loadVoxel(int &sizeX, int &sizeY, int &sizeZ)
{
    std::string line;
    std::filesystem::path path = std::filesystem::current_path();
    std::filesystem::path file_path = "/../data/voxelizedModel.bin";
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
    std::filesystem::path file_path = "/../data/Ke.txt";
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

template <typename scalar>
void getBoundaryConditions(std::set<uint64_t>& fixedNodes, std::map<uint64_t, Vec3<scalar>>& loadedNodes)
{
    std::string line;
    
    //fixing conditions
    std::filesystem::path path = std::filesystem::current_path();
    std::filesystem::path file_path = "/../data/fixed.dat";
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
    file_path = "/../data/loading.dat";
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

#pragma endregion input

#pragma region output

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
void saveMatrix(std::vector<std::vector<T>> data, std::string name)
{
	std::filesystem::path path = std::filesystem::current_path();
    std::filesystem::path file_path = "/../data/out/" + name + ".txt";
	path += file_path;
    
    std::ofstream myfile (path.string());

	outMatrix(data, myfile);

	myfile.close();
}

template <typename T>
void saveVector(std::vector<T> data, std::string name)
{
	std::filesystem::path path = std::filesystem::current_path();
    std::filesystem::path file_path = "/../data/out/" + name + ".csv";
	path += file_path;
    
    std::ofstream myfile (path.string());

	int n = data.size();
	for (int j = 0; j < n; j++)
	{
		myfile << "," << data[j];
	}
	myfile << std::endl;
	
	myfile.close();
}

template <typename T, size_t l>
void saveMatrix(std::vector<std::array<T,l>> data, std::string name)
{
	std::filesystem::path path = std::filesystem::current_path();
    std::filesystem::path file_path = "/../data/out/out/" + name + ".txt";
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
void saveMatrix(std::map<uint64_t, T> data, int rows, std::string name)
{
	std::filesystem::path path = std::filesystem::current_path();
    std::filesystem::path file_path = "/../data/out/" + name + ".csv";
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

template<typename T>
void printModel(T* model, int numElements)
{
    std::filesystem::path path = std::filesystem::current_path();
    std::filesystem::path file_path = "/../data/out/result.bin";
	path += file_path;
    
    std::ofstream myfile (path.string());
	for(int i = 0; i < numElements; i++)
		myfile << model[i] << std::endl;

	myfile.close();
}

#pragma endregion output