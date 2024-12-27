## gauss_parallel.exe

部分选主元高斯消元法的 mpi 并行程序，从 "matrix.txt" 中读取增广矩阵，解输出到 "p_solution.txt" 中。

## gauss_serial.exe

部分选主元高斯消元法的串行程序，从 "matrix.txt" 中读取增广矩阵，解输出到 "solution.txt" 中。

## matrix_generator.exe

用于生成随机 $ n \time (n+1) $ 维矩阵, 矩阵的元素服从 \( (0,1) \) 上的均匀分布, 生成结果会输出到文件 "matrix.txt" 中。

## matrix.txt

\( n \time (n+1) \) 维矩阵，文件的第一行是 \(n\)，后面 \(n\) 行是增广矩阵，矩阵同行元素用空格隔开。
