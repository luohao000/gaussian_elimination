#include <iostream>
#include <fstream>
#include <mpi.h>
using namespace std;

int get_local_rows(int n, int size, int rank) // 计算每个进程负责的行的数量
{
	int rows = n / size;
	if (n % size > rank)
	{
		rows += 1;
	}
	return rows;
}

void read_matrix(const char *fileName, double *&local_matrix, int &n, int &local_rows, int rank, int size)
{
	double *full_matrix = nullptr;

	if (rank == 0)
	{
		// 主进程读取整个矩阵
		ifstream file(fileName);
		if (!file.is_open())
		{
			cerr << "Error: Unable to open file " << fileName << " for reading." << endl;
			MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
		}

		file >> n; // 读取矩阵维数
		full_matrix = new double[n * (n + 1)];
		for (int i = 0; i < n * (n + 1); i++)
		{
			file >> full_matrix[i]; // full_matrix 按行读取
		}
		file.close();
	}

	// 广播矩阵维数
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// 计算本地行数
	local_rows = get_local_rows(n, size, rank); // 每个进程分配行数，向上取整
	int local_matrix_size = local_rows * (n + 1);

	// 分配本地矩阵存储空间
	local_matrix = new double[local_matrix_size];

	if (rank == 0)
	{
		// 分发矩阵行
		for (int i = 0; i < n; i++)
		{
			int dest_rank = i % size;
			if (dest_rank == 0)
			{
				// 主进程保留自己的行
				for (int j = 0; j <= n; j++)
				{
					local_matrix[(i / size) * (n + 1) + j] = full_matrix[i * (n + 1) + j]; // local_matrix 按行储存
				}
			}
			else
			{
				// 发送给目标进程
				MPI_Send(&full_matrix[i * (n + 1)], n + 1, MPI_DOUBLE, dest_rank, 0, MPI_COMM_WORLD);
			}
		}

		delete[] full_matrix; // 释放主进程的完整矩阵内存
	}
	else
	{
		// 接收本地的行
		MPI_Status status;
		for (int i = rank; i < n; i += size)
		{
			MPI_Recv(&local_matrix[(i / size) * (n + 1)], n + 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
		}
	}
}

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

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); // 获取当前进程号
	MPI_Comm_size(MPI_COMM_WORLD, &size); // 获取进程总数

	double time;
	if (rank == 0)
	{
		time = -MPI_Wtime();
	}

	double *local_matrix = nullptr; // 当前进程的本地矩阵
	int n = 0;						// 矩阵维数
	int local_rows = 0;				// 本地行数

	// 从文件读取矩阵并按行分配
	read_matrix("matrix.txt", local_matrix, n, local_rows, rank, size);

	double magnitude;
	double t;
	int picked;
	int tmp;
	double *pivot = new double[n + 1];

	int *label = new int[local_rows];
	for (int i = 0; i < local_rows; i++)
	{
		label[i] = -1;
	}

	// 高斯消元
	for (int i = 0; i < n; i++)
	{
		// 选主元
		magnitude = 0;
		for (int j = 0; j < local_rows; j++)
		{
			if (label[j] < 0)
			{
				if (abs(local_matrix[j * (n + 1) + i]) > magnitude)
				{
					magnitude = abs(local_matrix[j * (n + 1) + i]);
					picked = j;
				}
			}
		}
		struct
		{
			double m;
			int r;
		} local_max = {magnitude, rank}, global_max;
		MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
		if (rank == global_max.r)
		{
			if (global_max.m == 0)
			{
				cerr << "Error: Singular matrix detected during pivoting." << endl;
				MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
			}
			for (int k = i; k <= n; k++)
			{
				pivot[k] = local_matrix[picked * (n + 1) + k];
			}
			label[picked] = i;
		}
		MPI_Bcast(&pivot[i], n + 1 - i, MPI_DOUBLE, global_max.r, MPI_COMM_WORLD);

		// 消元
		for (int j = 0; j < local_rows; j++)
		{
			if (label[j] < 0)
			{
				t = local_matrix[j * (n + 1) + i] / pivot[i];
				for (int k = i + 1; k <= n; k++)
				{
					local_matrix[j * (n + 1) + k] -= pivot[k] * t;
				}
			}
		}
	}

	double *x = new double[n];

	// 回代求解
	for(int i = n - 1; i >= 0; i--)
	{
		int root;
		int b = 0;
		for (int j = 0; j < local_rows; j++)
		{
			if (label[j] == i)
			{
				x[i] = local_matrix[j * (n + 1) + n] / local_matrix[j * (n + 1) + i];
				b = rank;
				break;
			}
		}
		MPI_Allreduce(&b, &root, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
		MPI_Bcast(&x[i], 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
		for (int j = 0; j < local_rows; j++)
		{
			if (label[j] < i)
			{
				local_matrix[j * (n + 1) + n] -= x[i] * local_matrix[j * (n + 1) + i];
			}
		}
	}

	if (rank == 0)
	{
		export_x("p_solution.txt", x, n);
		time += MPI_Wtime();
		cout << "Execution time: " << time << " seconds" << endl;
	}

	delete[] local_matrix; // 释放本地矩阵内存
	delete[] label;
	delete[] pivot;
	delete[] x;

	MPI_Finalize();
	return 0;
}
