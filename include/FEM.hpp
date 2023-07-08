#include "utils.hpp"

#include <utility>
#include <map>
#include <unordered_map>
#include <array>

#include <math.h>

#define MAX_ITER 10
#define index(i,j,n) (i)*(n)+(j)
#define get_symmetric(A,i,j) (i) <= (j) ? A[(i)][(j)] : A[(j)][(i)]

template <typename scalar> 
struct Mat{
	std::map<uint64_t, scalar> data;
	uint64_t rows;
	const scalar zero = (scalar)0.0;

	Mat<scalar>(uint64_t size)
	{
		rows = size;
	}

	const scalar& operator[](uint64_t index) const
	{
		if(auto a = data.find(index); a != data.end())
			return a->second;
		else
			return zero;
	}

	scalar& operator[](uint64_t index)
	{
		if(auto a = data.find(index); a != data.end())
			return a->second;
		else
			data[index] = 0; return data[index];
	}

	std::vector<uint64_t> indices()
	{
		std::vector<uint64_t> keys;

		for(auto it = data.begin(); it != data.end(); it++)
			keys.push_back(it->first);

		return keys;
	}

	std::vector<scalar> values()
	{
		std::vector<scalar> values;

		for(auto it = data.begin(); it != data.end(); it++)
			values.push_back(it->second);

		return values;
	}

	void clearZeros()
	{
		std::erase_if(data, [](const auto& item){
			auto const& [key, value] = item;
			return value == 0;
		});
	}
};

typedef Mat<double> Matd;
typedef Mat<float> Matf;

inline bool stationary(Vec3i element)
{
	//Return true if at staionary boundary
	return element.z == 0;
}

template<typename scalar>
#ifdef CSPARSE
cs* assembleSystemMatrix(int* voxelModel, Vec3i voxelGridDimensions, scalar elementStiffness[24][24])
#else
Mat<scalar> assembleSystemMatrix(int* voxelModel, Vec3i voxelGridDimensions, scalar elementStiffness[24][24])
#endif
{
	Vec3i vertexGridDimensions = voxelGridDimensions + Vec3i(1);

	std::vector<Vec3i> usedElements;
	Vec3i element;
	FOR3(element,  Vec3i(0), voxelGridDimensions)
	{
		if(voxelModel[Linearize(element, voxelGridDimensions)] == 1 /*&& !stationary(element)*/)
		{
			usedElements.push_back(element);
		}
	}

	Vec3i::MAX = std::max(vertexGridDimensions.x,vertexGridDimensions.y);
	Vec3i::MAX = std::max(vertexGridDimensions.z, Vec3i::MAX);

	std::map<uint64_t, uint64_t> usedVertices;
	// std::vector<int> usedVertices;
	Vec3i vert;

	for(auto element : usedElements)
	{
		FOR3(vert, Vec3i(0), Vec3i(2))
		{
			auto val = Linearize(element + vert, vertexGridDimensions);
			if(!usedVertices.contains(val))
			{
				usedVertices[val] = val;
			}
		}
	}

	uint64_t index = 0;
	for (auto line : usedVertices)
	{
		usedVertices[line.first] = index++;
	}

	std::vector<std::array<uint64_t, 8>> elementToGlobal;
	for (auto element : usedElements)
	{
		std::array<uint64_t, 8> verts;
		FOR3(vert, Vec3i(0), Vec3i(2))
		{
			verts[Linearize(vert, Vec3i(2))] = usedVertices[Linearize(element + vert, vertexGridDimensions)];
		}
		elementToGlobal.push_back(verts);
	}

	// saveMatrix<uint64_t, 8>(elementToGlobal);

	int size = usedVertices.size() * 3;

#ifdef CSPARSE
	cs* systemMatrix = new cs();
	systemMatrix->m = size;
	systemMatrix->n = size;
#else
	Mat<scalar> systemMatrix(size);
#endif
	for(auto line : elementToGlobal)
	{
		FOR3(vert, Vec3i(0), Vec3i(2))
		{
			Vec3i vert2; 
			FOR3(vert2, Vec3i(0), Vec3i(2))
			{
				auto i = Linearize(vert,  Vec3i(2));
				auto j = Linearize(vert2, Vec3i(2));

				auto iGlobal = line[i];
				auto jGlobal = line[j];

				if(iGlobal == 10412 && jGlobal == 10418)
				{
					std::cout << elementStiffness[i*3 + 2][j*3] << std::endl;
				}

				for(int c1 = 0; c1 < 3; c1++)
				{
					for(int c2 = 0; c2 < 3; c2++)
					{
#ifdef CSPARSE
						cs_entry(systemMatrix, iGlobal, jGlobal, get_symmetric(elementStiffness,i*3 + c1,j*3 + c2));
#else
						systemMatrix[index(iGlobal*3 + c1, jGlobal*3 + c2, size)] += get_symmetric(elementStiffness,i*3 + c1,j*3 + c2);
#endif
					}
				}
			}
		}
	}

#ifdef CSPARSE
	cs_dupl(systemMatrix);
	std::cout << "System Matrix : Size:" << size << "x" << size << ", Non-Zero:" << systemMatrix->nz << ", " << ((float)systemMatrix->nz/size) << " full per row" << std::endl;
#else
	std::cout << "System Matrix : Size:" << size << "x" << size << ", Non-Zero:" << systemMatrix.data.size() << ", " << ((float)systemMatrix.data.size()/size) << " full per row" << std::endl;
#endif
	return systemMatrix;
}

template <typename scalar>
std::vector<scalar> multiply(const Mat<scalar>& A, const std::vector<scalar>& b)
{
	uint64_t size = b.size();
	std::vector<scalar> f;
	f.resize(size, 0);

	//#pragma omp parallel
	for(auto it = A.data.begin(); it != A.data.end(); it++)
	{
		uint64_t i = it->first / size;
		uint64_t j = it->first % size;
		
		f[i] += it->second * b[j];
	}

	return f;
}

template <typename scalar>
std::vector<scalar> multiply(const std::vector<scalar>& a, const std::vector<scalar>& b)
{
	uint64_t size = a.size();
	std::vector<scalar> result(a);

	//#pragma omp parallel for
	for(uint64_t i = 0; i < size; i++)
	{
		result[i] *= b[i];
	}

	return result;
}

template <typename scalar>
void add(std::vector<scalar>& addee, const std::vector<scalar>& addend)
{
	uint64_t size = addee.size();

	//#pragma omp parallel for
	for(uint64_t i = 0; i < size; i++)
	{
		addee[i] += addend[i];
	}
}

template <typename scalar>
void subtract(std::vector<scalar>& subtractee, const std::vector<scalar>& subtracted)
{
	uint64_t size = subtractee.size();

	//#pragma omp parallel for
	for(uint64_t i = 0; i < size; i++)
	{
		subtractee[i] -= subtracted[i];
	}
}

template <typename scalar>
scalar multiplyTranspose(const std::vector<scalar>& a, const std::vector<scalar>& b)
{
	uint64_t size = a.size();
	scalar result = 0;
	
	//#pragma omp parallel for
	for(int64_t i = 0; i < size; i++)
	{
		result += a[i]*b[i];
	}

	return result;
}

template <typename scalar>
std::vector<scalar> multiply(const scalar& a, const std::vector<scalar>& b)
{
	uint64_t size = b.size();
	std::vector<scalar> result(b);
	
	//#pragma omp parallel for
	for(int64_t i = 0; i < size; i++)
	{
		result[i] *= a;
	}

	return result;
}

template <typename scalar>
scalar two_norm(const std::vector<scalar>& b)
{
	scalar norm = (scalar)0;

	for(auto val : b)
	{
		norm += val*val;
	}

	return sqrt(norm);
}

template <typename scalar>
void solveWithCG(const Mat<scalar>& A, const std::vector<scalar>& b, std::vector<scalar>& x)
{
	uint64_t size = x.size();
	std::vector<scalar> r(b);
	std::vector<scalar> r2(r);

	scalar lambda = (scalar)0;
	subtract(r, multiply(A, x));
	std::vector<scalar> p(r);

	for(int i = 0; i < MAX_ITER; i++)
	{
		scalar n1 = multiplyTranspose(r,r);
		lambda = n1 / multiplyTranspose(p, multiply(A, p));

		// lambda *= (scalar)0.5;

		add(x, multiply(lambda, p));

		subtract(r, multiply(lambda, multiply(A, p)));
		r2 = r;
		
		scalar n2 = multiplyTranspose(r2, r2);

		add(r2, multiply(n2/n1, p));
		p = r2;

		r.swap(r2);

		std::cout << "iteration: " << i << ", " << "Norm: " << two_norm(r2) << std::endl;
	}
}

template void solveWithCG<double>(const Matd& A, const std::vector<double>& b, std::vector<double>& x);
template void solveWithCG<float>(const Matf& A, const std::vector<float>& b, std::vector<float>& x);
