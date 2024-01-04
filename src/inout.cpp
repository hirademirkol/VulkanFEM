#include "in.hpp"
#include "out.hpp"

#pragma region input

int* loadVoxel(std::string& model_name, int &sizeX, int &sizeY, int &sizeZ)
{
    std::string line;
    std::filesystem::path path = std::filesystem::current_path();
    std::filesystem::path file_path = "/" + model_name + ".dat";
    path += file_path;
    
    std::ifstream myfile (path.string());

    // Read the dimensions first
    std::getline(myfile, line, ' ');
    sizeX = std::stoi(line);
    std::getline(myfile, line, ' ');
    sizeY = std::stoi(line);
    std::getline(myfile, line);
    sizeZ = std::stoi(line);

    int* model = new int[sizeX*sizeY*sizeZ];

    // Read each element into the model
    int element = 0;
    for(int i = 0; i < sizeZ; i++)
    {
        for(int j = 0; j < sizeX*sizeY - 1; j++)
        {
            std::getline(myfile, line, ',');
            model[element++] = std::stoi(line);
        }
        std::getline(myfile, line);
        model[element++] = std::stoi(line);
    }

    myfile.close();

    return model;
}

template <typename scalar>
void getKe(std::string& model_name, scalar Ke[24][24])
{
    std::string line;
    std::filesystem::path path = std::filesystem::current_path();
    std::filesystem::path file_path = "/" + model_name + "_Ke0.dat";
    path += file_path;

    if(!std::filesystem::exists(path))
    {
        path = std::filesystem::current_path();
        file_path = "/data/Ke0.dat";
        path += file_path;
    }
    
    std::ifstream myfile (path.string());

    for(int i = 0; i < 24; i++)
    {
        for(int j = 0; j < 23; j++)
        {
            std::getline(myfile, line, ',');            
            Ke[i][j] = (scalar)std::stod(line);
        }   
        std::getline(myfile, line);
        Ke[i][23] = (scalar)std::stod(line);
    }
    
    myfile.close();
}

template void getKe<double>(std::string& model_name, double Ke[24][24]);
template void getKe<float>(std::string& model_name, float Ke[24][24]);

template <typename scalar>
void getBoundaryConditions(std::string& model_name, std::set<uint64_t>& fixedNodes, std::map<uint64_t, Vec3<scalar>>& loadedNodes)
{
    std::string line;
    
    // Fixing conditions
    std::filesystem::path path = std::filesystem::current_path();
    std::filesystem::path file_path = "/" + model_name + "_fixed.dat";
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

    // Loading conditions
    path = std::filesystem::current_path();
    file_path = "/" + model_name + "_loading.dat";
    path += file_path;

    myfile.open(path.string());

    std::getline(myfile, line);
    int numLoads = std::stoi(line);

    for (int i = 0; i < numLoads; i++)
    {
        std::getline(myfile, line, '\t');
        uint64_t ind = std::stoull(line) - 1;

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

template void getBoundaryConditions<double>(std::string& model_name, std::set<uint64_t>& fixedNodes, std::map<uint64_t, Vec3<double>>& loadedNodes);
template void getBoundaryConditions<float>(std::string& model_name, std::set<uint64_t>& fixedNodes, std::map<uint64_t, Vec3<float>>& loadedNodes);

#pragma endregion input

#pragma region output

// Outputs the matrix in a formatted form to an output stream, not suitable for bigger matrices
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

template <typename T>
void saveMatrix(std::vector<std::vector<T>> data, std::string name)
{
	std::filesystem::path path = std::filesystem::current_path();
    std::filesystem::path file_path = "/" + name + "_out.csv";
	path += file_path;
    
    std::ofstream myfile (path.string());

	outMatrix(data, myfile);

	myfile.close();
}

template void saveMatrix<double>(std::vector<std::vector<double>> data, std::string name);
template void saveMatrix<float>(std::vector<std::vector<float>> data, std::string name);

template <typename T>
void printMatrix(std::vector<std::vector<T>> data)
{
	outMatrix(data, std::cout);
}

template void printMatrix<double>(std::vector<std::vector<double>> data);
template void printMatrix<float>(std::vector<std::vector<float>> data);

template <typename T>
void saveVector(std::vector<T> data, std::string name)
{
	std::filesystem::path path = std::filesystem::current_path();
    std::filesystem::path file_path = "/" + name + "_out.csv";
	path += file_path;
    
    std::ofstream myfile (path.string());

    myfile << std::to_string(data[0]);

	int n = data.size();
	for (int j = 1; j < n; j++)
	{
		myfile << "," << std::to_string(data[j]);
	}
	myfile << std::endl;
	
	myfile.close();
}

template void saveVector<double>(std::vector<double> data, std::string name);
template void saveVector<float>(std::vector<float> data, std::string name);

template <typename T>
void saveMatrix(std::map<uint64_t, T> data, int rows, std::string name)
{
	std::filesystem::path path = std::filesystem::current_path();
    std::filesystem::path file_path = "/" + name + "_out.csv";
	path += file_path;
    
    std::ofstream myfile (path.string());

	for(auto line = data.begin(); line != data.end(); line++)
	{
		myfile << line->first/rows <<"," << line->first%rows << ","<< line->second << std::endl;
	}

	myfile.close();
}

template void saveMatrix<double>(std::map<uint64_t, double> data, int rows, std::string name);
template void saveMatrix<float>(std::map<uint64_t, float> data, int rows, std::string name);

template <typename T>
void saveMatrix(std::vector<Eigen::Triplet<T>> triplets, std::string name)
{
	std::filesystem::path path = std::filesystem::current_path();
    std::filesystem::path file_path = "/" + name + "_out.csv";
	path += file_path;

    std::ofstream myfile (path.string());

    for(auto triplet : triplets)
    {
        myfile << triplet.row()   << ",";
        myfile << triplet.col()   << ",";
        myfile << triplet.value() << std::endl;
    }

    myfile.close();    
}

template void saveMatrix<double>(std::vector<Eigen::Triplet<double>> data, std::string name);
template void saveMatrix<float>(std::vector<Eigen::Triplet<float>> data, std::string name);

template <typename T, int i, int j>
void saveArray(Eigen::Array<T, i, j> array, std::string name)
{
    std::filesystem::path path = std::filesystem::current_path();
    std::filesystem::path file_path = "/" + name + "_out.csv";
	path += file_path;

    std::ofstream myfile (path.string());

    for(int x = 0; x < array.rows(); x++)
    {
        for(int y = 0; y < array.cols(); y++)
        {
            myfile << array(x,y) << ",";
        }
        myfile << std::endl;
    }

    myfile.close(); 
}

template void saveArray<int, -1, -1>(Eigen::Array<int, -1, -1> array, std::string name);
template void saveArray<int, -1, 8>(Eigen::Array<int, -1, 8> array, std::string name);
template void saveArray<int, -1, 27>(Eigen::Array<int, -1, 27> array, std::string name);

#pragma endregion output