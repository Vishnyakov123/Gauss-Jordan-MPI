#include "header.h"

int fill_from_file(double *matrix_buffer, int n, char* file, int rank, int size, double *b){

	int rows;
	FILE* f = fopen(file, "r");

	if (rank != size - 1)
		rows = (rank + 1) * (n / size) - rank * (n / size);
	else
		rows = (rank + 1) * (n / size) - rank * (n / size) + n % size;

	if (f == NULL) {
		fprintf(stderr, "Can not open file\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
		return -1;
	}

	if (rank == 0) {
		
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < n; j++)
			{
				if (fscanf(f, "%lf", &matrix_buffer[i * n + j]) <= 0) {
					fclose(f);
					fprintf(stderr, "Error in file scanning\n");
					MPI_Abort(MPI_COMM_WORLD, 1);
				}
			}
		}

		for (int i = 0; i < rows; i++)
		{	
			b[i] = 0.; 

			for (int l = 0; l < (n + 1) / 2; l++)
			{
				b[i] = b[i] + matrix_buffer[i * n + 2 * l];
			}
		}
	}

	if (rank != 0) {

		int counter = 0;
		char c;

		while (1)
		{
			fscanf(f, "%c", &c);
			if (c == '\n')
				counter += 1;
			if (counter == n / size * rank || c == EOF)
				break;
		}

		for (int i = 0; i < rows; i++) 
		{
			for (int j = 0; j < n; j++)
			{
				if (fscanf(f, "%lf", &matrix_buffer[i * n + j]) <= 0) {
					fclose(f);
					fprintf(stderr, "Error in file scanning\n");
					MPI_Abort(MPI_COMM_WORLD, 1);
					return -1;
				}
			}
		}

		for (int i = 0; i < rows; i++)
		{
			b[i] = 0.;

			for (int l = 0; l < (n + 1) / 2; l++)
			{
				b[i] = b[i] + matrix_buffer[i * n + 2 * l];
			}
		}
	}

	fclose(f);
	return 0;
}

void print_matrix(double *matrix_buffer, int n, int size, int rank, double *print_buf,int m, int stop) {

	int rows;
	MPI_Status status;

	if (m != 0 && m != -1 && stop!=0)
		printf("Matrix:\n");
	else
		printf("Matrix after method:\n");

	rows = (rank + 1) * (n / size) - rank * (n / size);

	for (int i = 0; i < rows; i++)
	{
		if (i == m)
			return;/*чтобы  выводило не более m строк и столбцов*/
		for (int j = 0; j < n; j++)
		{
			if (j == m)
				break;/*чтобы  выводило не более m строк и столбцов*/
			printf(" %10.3e ", matrix_buffer[i * n + j]);
		}
		printf("\n");
	}

	for (int i = 1; i < size; i++)
	{
		if (i != size - 1)
			rows = (rank + 1) * (n / size) - rank * (n / size);
		else
			rows = (rank + 1) * (n / size) - rank * (n / size) + n % size;

		MPI_Recv(print_buf, rows * n, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);

		for (int k = 0; k < rows; k++)
		{
			if (k + i * (n / size) == m)
				return;

			for (int j = 0; j < n; j++)
			{
				if (j == m)
					break;
				printf(" %10.3e ", print_buf[k * n + j]);
			}
			printf("\n");
		}			
	}		
}

void fill_matrix(int n, double *matrix_buffer, int size, int rank, double *b,int k){

	int rows = (rank + 1) * (n / size) - rank * (n / size);
	
	if (rank == size - 1)
		rows = (rank + 1) * (n / size) - rank * (n / size) + n % size;

	for (int i = 0; i < rows; i++)
	{
		b[i] = 0.;
		for (int j = 0; j < n; j++)
		{	
			switch (k)
			{
			case 1:
				matrix_buffer[i * n + j] = n - max(i + rank * (n / size) + 1, j + 1) + 1;
				break;
			case 2:
				matrix_buffer[i * n + j] = max(i + rank * (n / size) + 1, j + 1);
				break;
			case 3:
				matrix_buffer[i * n + j] = fabs((double)(i + rank * (n / size) + 1 - j - 1));
				break;
			case 4:
				matrix_buffer[i * n + j] = 1.0 / ((double)(i + rank * (n / size) + 1.0 + j));
			default:
				break;
			}
		}

		for (int l = 0; l < (n + 1) / 2; l++)
		{
			b[i] = b[i] + matrix_buffer[i * n + 2 * l];
		}
	}
}

int max(int i, int j)
{
	if (i >= j)
		return i;
	return j;
}

