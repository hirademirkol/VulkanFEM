#include <iostream>


template <typename T>
void printMatrix(std::vector<T> data, int32_t n, int32_t m)
{
if(n == 1)
{
	std::cout << "[ ";
	for (int j= 0; j < m; j++)
	{
		std::cout << data[j] << "  ";
	}
	std::cout << "]" << std::endl;
	return;
}

	std::cout << "⎡ ";
	for (int j= 0; j < m; j++)
	{
		std::cout << data[j] << " ";
	}
	std::cout << "⎤" << std::endl;

	for (int i = 1; i < n-1; i++)
	{
		std::cout << "⎢ ";
		for (int j= 0; j < m; j++)
		{
			std::cout << data[m*i + j] << " ";
		}
		std::cout << "⎥" << std::endl;
	}

	std::cout << "⎣ ";
	for (int j= 0; j < m; j++)
	{
		std::cout << data[m*(n-1) + j] << " ";
	}
	std::cout << "⎦" << std::endl;
}
