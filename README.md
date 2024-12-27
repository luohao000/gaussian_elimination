## 程序说明

### gauss_parallel.exe

部分选主元高斯消元法的 mpi 并行程序，从 "matrix.txt" 中读取增广矩阵，解输出到 "p_solution.txt" 中。

### gauss_serial.exe

部分选主元高斯消元法的串行程序，从 "matrix.txt" 中读取增广矩阵，解输出到 "solution.txt" 中。

### matrix_generator.exe

用于生成随机 n*(n+1) 维矩阵, 矩阵的元素服从 (0,1) 上的均匀分布, 生成结果会输出到文件 "matrix.txt" 中。

### matrix.txt

n*(n+1) 维矩阵，文件的第一行是 n，后面 n 行是增广矩阵，矩阵同行元素用空格隔开。

## 如何运行程序?

在 cmd 或者 powershell 中进行如下操作

1. 运行 matrix_generator
  - `make matrix_generator`
  - `./matrix_generator`
  - `2000`(这里输入的是要生成的增广矩阵的行数)
2. 运行 gauss_serial
  - `make gauss_serial`
  - `./gauss_serial`
3. 运行 gauss_parallel
  - `g++ gauss_parallel.cpp -o gauss_parallel.exe -I "C:\Program Files (x86)\Microsoft SDKs\MPI\Include" -L "C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64" -lmsmpi` (这里输入自己的Include和Lib路径)
  - `mpiexec -n 8 gauss_parallel` (进程数可以更改)

### 注意事项

只有当矩阵维数较大时（大约上千维时），并行程序才会明显快于串行程序。

### 运行结果示例

![image](运行结果.png)
