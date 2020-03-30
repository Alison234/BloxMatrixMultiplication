#include <iostream>
#include <omp.h>
#include <cstdlib>

long** A_matrix, ** B_matrix, ** result_usual;
long* A, * B, ** result_block_parallel;
long* indexes_matrix, N;

void init_matrix()
{
	A_matrix = new long* [N];
	for (int i = 0; i < N; ++i)
		A_matrix[i] = new long[N];

	B_matrix = new long* [N];
	for (int i = 0; i < N; ++i)
		B_matrix[i] = new long[N];

	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
		{
			if (j > i)
			{
				A_matrix[i][j] = 0;
				B_matrix[i][j] = 0;
			}
			else
			{
				A_matrix[i][j] = rand() % 10+1;
				B_matrix[i][j] = rand() % 10+1;
			}
		}
}

void init_matrix_array(int block_size)
{
	int number_of_blocks = N / block_size;
	long number_of_element = ((number_of_blocks + 1) * number_of_blocks / 2) * block_size * block_size;

	A = new long[number_of_element];
	B = new long[number_of_element];
	indexes_matrix = new long[N * N];

	int index = 0;
	for (int i = 0; i < number_of_blocks; ++i)
		for (int j = 0; j <= i; ++j)
			for (int iInto = 0; iInto < block_size; ++iInto)
				for (int jInto = 0; jInto < block_size; ++jInto)
				{
					B[index] = B_matrix[i * block_size + iInto][j * block_size + jInto];
					index++;
				}
	 index = 0;
	for (int i = 0; i < number_of_blocks; ++i)
		for (int j = i; j < number_of_blocks;++j)
			for (int iInto = 0; iInto < block_size; ++iInto)
				for (int jInto = 0; jInto < block_size; ++jInto)
				{
					A[index] = A_matrix[j * block_size + iInto][i * block_size + jInto];
					index++;
				}

}
void print_matrix(long** matrix)
{
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
			std::cout << matrix[i][j] << " ";
		std::cout << std::endl;
	}
}

void print_array(long* matrix,int block_size)
{
	int number_of_blocks = N / block_size;
	int number_of_element = ((number_of_blocks + 1) * number_of_blocks / 2) * block_size * block_size;

	for (int i = 0; i < number_of_element; ++i)
	{
		std::cout << matrix[i] << " ";
	}
}

void multiply_matrixUsual()
{
	result_usual = new long* [N];
	for (int i = 0; i < N; ++i)
		result_usual[i] = new long[N];

	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
			result_usual[i][j] = 0;

	double start_time = omp_get_wtime();

	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
			for (int k = 0; k < N; ++k)
				result_usual[i][j] += A_matrix[i][k] * B_matrix[k][j];

	double end_time = omp_get_wtime();

	std::cout << "\n Non parallel multiplication \n" << std::endl;
	std::cout << "Run time: " << end_time - start_time << std::endl;

}
void printBloxA(int index, int blocksize) {
	for (int i = 0; i < blocksize * blocksize; i++) {
		std::cout << A[index* blocksize * blocksize + i] << " ";
	}
}

void printBloxB(int index, int blocksize) {
	for (int i = 0; i < blocksize * blocksize; i++) {
		std::cout << B[index * blocksize * blocksize + i] << " ";
	}
}
void parallelMultiplyOne(int sizeOfBlock, int threadsAmount,int n)
{
	int Num = n / sizeOfBlock;
	result_block_parallel = new long* [n];
	for (int i = 0; i < n; ++i)
		result_block_parallel[i] = new long[n];

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			result_block_parallel[i][j] = 0;

	double start_time = omp_get_wtime();
#pragma omp parallel num_threads(threadsAmount)
	for (int i = 0; i < Num; ++i)
	{
#pragma omp	for
		for (int j = 0; j < Num; ++j)
		{
			for (int k = j; k <= i; ++k)
			{
				int beginA = ((2 * Num - k + 1) * k / 2 + i - k)* sizeOfBlock* sizeOfBlock;
				int beginB = ((1 + k) * k / 2 + j) * sizeOfBlock * sizeOfBlock;
				for (int ii = 0; ii < sizeOfBlock; ++ii)
					for (int jj = 0; jj < sizeOfBlock; ++jj)
						for (int kk = 0; kk < sizeOfBlock; ++kk)
						{
							result_block_parallel[i * sizeOfBlock + ii][j * sizeOfBlock + jj] +=
								A[beginA + ii * sizeOfBlock + kk]
								* B[beginB + kk * sizeOfBlock + jj];
						}
			}
		}
	}

	double end_time = omp_get_wtime();
	std::cout << "\n\t Parallel multiplication 1 \t\n" << std::endl;
	std::cout << "Block_size " << sizeOfBlock  <<std::endl;
	std::cout << "Run time: " << end_time - start_time << std::endl;
}

void parallelMultiplyTwo(int sizeOfBlock, int threadsAmount, int n)
{
	int Num = n / sizeOfBlock;
	result_block_parallel = new long* [n];
	for (int i = 0; i < n; ++i)
		result_block_parallel[i] = new long[n];

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			result_block_parallel[i][j] = 0;


	double start_time = omp_get_wtime();
#pragma omp parallel num_threads(threadsAmount)
	for (int i = 0; i < Num; ++i)
	{
		for (int j = 0; j < Num; ++j)
		{
#pragma omp	for
			for (int k = j; k <= i; ++k)
			{
				int beginA = ((2 * Num - k + 1) * k / 2 + i - k) * sizeOfBlock * sizeOfBlock;
				int beginB = ((1 + k) * k / 2 + j) * sizeOfBlock * sizeOfBlock;
				for (int ii = 0; ii < sizeOfBlock; ++ii)
					for (int jj = 0; jj < sizeOfBlock; ++jj)
						for (int kk = 0; kk < sizeOfBlock; ++kk)
						{
							result_block_parallel[i * sizeOfBlock + ii][j * sizeOfBlock + jj] +=
								A[beginA + ii * sizeOfBlock + kk]
								* B[beginB + kk * sizeOfBlock + jj];
						}
			}
		}
	}

	double end_time = omp_get_wtime();
	std::cout << "\n\t Parallel multiplication 2 \t\n" << std::endl;
	std::cout << "Block_size " << sizeOfBlock << std::endl;
	std::cout << "Run time: " << end_time - start_time << std::endl;
}



void checkAlgorithm(int n)
{
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			if (result_usual[i][j] != result_block_parallel[i][j])
			{
				std::cout << "Result is not correct!" << std::endl;
				return;
			}
	std::cout << "Result is correct!" << std::endl;
}



int main() {
	N = 16;
	init_matrix();
	multiply_matrixUsual();
	
	//std::cout << std::endl;
	//std::cout << std::endl;
	//print_array(A, 2);
	int BlockSize[] = { 10, 15, 20, 24, 30, 36, 40, 60, 72, 80,
		96, 120, 144, 160, 180, 240, 360, 480, 720 };

	//print_matrix(A_matrix);

	//init_matrix_array(144);
	//parallelMultiplyOne(144, 3, N);
	//for (int i = 0; i < 18; i++) {
		init_matrix_array(4);
		parallelMultiplyOne(4, 8, N);
		//parallelMultiplyTwo(BlockSize[i], 8, N);
		checkAlgorithm(N);
	//}
	//std::cout << "Parallel matrix";
	//std::cout << std::endl;
	//print_matrix(result_block_parallel);
	//std::cout << "usual matrix";
	//std::cout << std::endl;
	//print_matrix(result_usual);
	

	return 1;
}