#include <iostream>
#include <fstream>
#include <mpi.h>
#include <cmath>
#include <cstdlib>

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
    local_rows = get_local_rows(n, size, rank);
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
                    local_matrix[(i / size) * (n + 1) + j] = full_matrix[i * (n + 1) + j];
                }
            }
            else
            {
                // 发送给目标进程
                MPI_Send(&full_matrix[i * (n + 1)], n + 1, MPI_DOUBLE, dest_rank, 0, MPI_COMM_WORLD);
            }
        }

        delete[] full_matrix;
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
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double time;
    if (rank == 0)
    {
        time = -MPI_Wtime();
    }

    double *local_matrix = nullptr; 
    int n = 0;
    int local_rows = 0;

    // 从文件读取矩阵并按行分配
    read_matrix("matrix.txt", local_matrix, n, local_rows, rank, size);

    double *pivot = new double[n + 1];

    // 使用 isPivoted 代替 label，逻辑更清晰：false表示该行未被用作主元
    bool *isPivoted = new bool[local_rows];
    for (int i = 0; i < local_rows; i++)
    {
        isPivoted[i] = false;
    }

    // 高斯消元
    for (int i = 0; i < n; i++)
    {
        double magnitude = 0.0;
        int picked = -1;

        // 在本地行中选主元
        for (int j = 0; j < local_rows; j++)
        {
            if (!isPivoted[j])
            {
                double val = fabs(local_matrix[j * (n + 1) + i]);
                if (val > magnitude)
                {
                    magnitude = val;
                    picked = j;
                }
            }
        }

        // MPI MAXLOC
        struct
        {
            double m;
            int r;
        } local_max = {magnitude, rank}, global_max;

        MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

        if (rank == global_max.r)
        {
            if (global_max.m == 0.0)
            {
                cerr << "Error: Singular matrix detected during pivoting." << endl;
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }
            for (int k = i; k <= n; k++)
            {
                pivot[k] = local_matrix[picked * (n + 1) + k];
            }
            isPivoted[picked] = true;
        }

        MPI_Bcast(&pivot[i], n + 1 - i, MPI_DOUBLE, global_max.r, MPI_COMM_WORLD);

        // 消元
        for (int j = 0; j < local_rows; j++)
        {
            if (!isPivoted[j])
            {
                double factor = local_matrix[j * (n + 1) + i] / pivot[i];
                for (int k = i + 1; k <= n; k++)
                {
                    local_matrix[j * (n + 1) + k] -= pivot[k] * factor;
                }
            }
        }
    }

    double *x = new double[n];

    // 回代求解
    for (int i = n - 1; i >= 0; i--)
    {
        double xi = 0.0;
        int root = -1;
        for (int j = 0; j < local_rows; j++)
        {
            // 在 pivot[i] 所在行找到对应的 x[i]
            // 此处通过判断 local_matrix[j*(n+1)+i] 的系数是否等于 pivot[i]
            // 也可以在前面记录具体被选中的主元行以简化这里的查找
            double val = local_matrix[j * (n + 1) + i];
            // 如果该行确实被用作第 i 个主元行，则 val 应该对应 pivot[i]
            // (允许一些浮动误差)
            if (fabs(val - pivot[i]) < 1e-14 && fabs(pivot[i]) > 1e-14)
            {
                xi = local_matrix[j * (n + 1) + n] / val;
                root = rank;
                break;
            }
        }

        // 广播 x[i] 值，找到的进程号进行广播
        int broadcastRank;
        MPI_Allreduce(&root, &broadcastRank, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        MPI_Bcast(&xi, 1, MPI_DOUBLE, broadcastRank, MPI_COMM_WORLD);
        x[i] = xi;

        // 更新本地行
        for (int j = 0; j < local_rows; j++)
        {
            local_matrix[j * (n + 1) + n] -= x[i] * local_matrix[j * (n + 1) + i];
        }
    }

    if (rank == 0)
    {
        export_x("p_solution.txt", x, n);
        time += MPI_Wtime();
        cout << "Execution time: " << time << " seconds" << endl;
    }

    delete[] local_matrix;
    delete[] pivot;
    delete[] isPivoted;
    delete[] x;

    MPI_Finalize();
    return 0;
}