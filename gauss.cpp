#include <iostream>
#include <fstream>
using namespace std;

// 从文件 matrix.txt 中读取增广矩阵
// matrix.txt 中的第一行是矩阵 A 的维数 n
// 后面 n 行是增广矩阵，矩阵同行元素用空格隔开
// 可以从 matrix.txt 中输入矩阵，运行程序后会将解输出到文件 solution.txt 中
void read_A(const char* fileName, double*& A, int& n)
{
	FILE* file = NULL;
	file = fopen(fileName, "r");
	fscanf(file, "%d", &n);
	A = new double[n * (n + 1)];
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j <= n; j++)
		{
			fscanf(file, "%lf", &A[i + j * n]);
		}
	}
	fclose(file);
}

// 将方程组的解输出到文件 solution.txt 中
void export_x(const char* fileName, double* x, int n, int* loc)
{
	ofstream os(fileName);
	for (int i = 0; i < n - 1; i++)
	{
		os << x[i] << endl;
	}
	os << x[n - 1];
	os.close();
}

int main()
{
	int n;
	int* loc;
	double* x, * A;
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

	for (int i = 0; i < n - 1; i++)
	{
		// 选主元
		magnitude = 0;
		for (int j = i; j < n; j++)
		{
			if (abs(A[loc[i] + j * n]) > magnitude)
			{
				magnitude = A[loc[i] + j * n];
				picked = j;
			}
		}
		tmp = loc[i];
		loc[i] = loc[picked];
		loc[picked] = tmp;

		// 消元
		for (int j = i + 1; j < n; j++)
		{
			t = A[loc[j] + i * n] / A[loc[i] + i * n];
			for (int k = i + 1; k <= n; k++)
			{
				A[loc[j] + k * n] = A[loc[j] + k * n] - A[loc[i] + k * n] * t;
			}
		}
	}

	// 回代求解
	for (int i = n - 1; i >= 0; i--)
	{
		x[i] = A[loc[i] + n * n] / A[loc[i] + i * n];
		for (int j = 0; j < i; j++)
		{
			A[loc[j] + n * n] = A[loc[j] + n * n] - x[i] * A[loc[j] + i * n];
		}
	}

	export_x("solution.txt", x, n, loc);

	delete[] loc;
	delete[] x;
	delete[] A;

	return 0;
}