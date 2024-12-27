#include <iostream>
#include <fstream>
#include <random>
#include <chrono>

int main()
{
	int n;
	std::cout << "Enter the dimension n for the n*n matrix: ";
	std::cin >> n;

	// 检查输入的维度是否有效
	if (n <= 0)
	{
		std::cerr << "Invalid matrix dimension. n must be positive." << std::endl;
		return 1;
	}

	// 动态分配二维数组
	double **matrix = new double *[n];
	for (int i = 0; i < n; ++i)
	{
		matrix[i] = new double[n + 1]; // 按行储存矩阵
	}

	// 初始化随机数生成器
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::mt19937 gen(seed);
	std::uniform_real_distribution<double> dist(0.0, 1.0);

	// 填充矩阵
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j <= n; ++j)
		{
			matrix[i][j] = dist(gen);
		}
	}

	// 写入矩阵到文件
	std::ofstream outfile("matrix.txt");
	if (!outfile)
	{
		std::cerr << "Error: Could not open file for writing." << std::endl;
		// 释放内存
		for (int i = 0; i < n; ++i)
		{
			delete[] matrix[i];
		}
		delete[] matrix;
		return 1;
	}

	// 写入矩阵维度
	outfile << n << std::endl;

	// 写入矩阵数据
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j <= n; ++j)
		{
			outfile << matrix[i][j] << " ";
		}
		outfile << std::endl;
	}

	outfile.close();
	std::cout << "Matrix has been written to matrix.txt" << std::endl;

	// 释放动态分配的内存
	for (int i = 0; i < n; ++i)
	{
		delete[] matrix[i];
	}
	delete[] matrix;

	return 0;
}