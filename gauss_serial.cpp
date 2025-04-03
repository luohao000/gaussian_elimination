﻿#include <iostream>
#include <fstream>
#include <chrono>
using namespace std;
using namespace chrono;

// 从文件 matrix.txt 中读取增广矩阵
// matrix.txt 中的第一行是矩阵 A 的维数 n
// 后面 n 行是增广矩阵，矩阵同行元素用空格隔开
void read_A(const char *fileName, double *&A, int &n)
{
	FILE *file = fopen(fileName, "r");
	if (!file)
	{
		cerr << "Error: Unable to open file " << fileName << " for reading." << endl;
		exit(EXIT_FAILURE);
	}

	fscanf(file, "%d", &n);
	A = new double[n * (n + 1)];
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j <= n; j++)
		{
			fscanf(file, "%lf", &A[i * (n + 1) + j]); // 矩阵 A 是按行储存的, A_ij = A[i * (n + 1) + j]
		}
	}
	fclose(file);
}

// 将方程组的解输出到文件 solution.txt 中
void export_x(const char *fileName, double *x, int n)
{
	ofstream os(fileName);
	if (!os)
	{
		cerr << "Error: Unable to open file " << fileName << " for writing." << endl;
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < n; i++)
	{
		os << x[i] << endl;
	}
	os.close();
	cout << "Solutions have been written to " << fileName << endl;
}

int main()
{
	auto start = system_clock::now();

	int n;
	int *loc;
	double *x, *A;
	double magnitude;
	double t;
	int picked;
	int tmp;

	read_A("matrix.txt", A, n);

	x = new double[n];
	loc = new int[n];

	for (int i = 0; i < n; i++)
	{
		loc[i] = i;
	}

	for (int i = 0; i < n; i++)
	{
		// 选主元
		magnitude = 0;
		for (int j = i; j < n; j++)
		{
			if (abs(A[loc[j] * (n + 1) + i]) > magnitude)
			{
				magnitude = abs(A[loc[j] * (n + 1) + i]);
				picked = j;
			}
		}
		if (magnitude == 0)
		{
			cerr << "Error: Singular matrix detected during pivoting." << endl;
			delete[] loc;
			delete[] x;
			delete[] A;
			exit(EXIT_FAILURE);
		}

		tmp = loc[i];
		loc[i] = loc[picked];
		loc[picked] = tmp;

		// 消元
		for (int j = i + 1; j < n; j++)
		{
			t = A[loc[j] * (n + 1) + i] / A[loc[i] * (n + 1) + i];
			for (int k = i + 1; k <= n; k++)
			{
				A[loc[j] * (n + 1) + k] = A[loc[j] * (n + 1) + k] - A[loc[i] * (n + 1) + k] * t;
			}
		}
	}

	// 回代求解
	for (int i = n - 1; i >= 0; i--)
	{
		x[i] = A[loc[i] * (n + 1) + n] / A[loc[i] * (n + 1) + i];
		for (int j = 0; j < i; j++)
		{
			A[loc[j] * (n + 1) + n] = A[loc[j] * (n + 1) + n] - x[i] * A[loc[j] * (n + 1) + i];
		}
	}

	export_x("solution.txt", x, n);

	auto end = system_clock::now();
	auto duration = duration_cast<microseconds>(end - start);
	cout << "Execution time: " << double(duration.count()) * microseconds::period::num / microseconds::period::den << " seconds" << endl;

	delete[] loc;
	delete[] x;
	delete[] A;

	return 0;
}