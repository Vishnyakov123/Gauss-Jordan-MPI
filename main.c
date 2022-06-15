#include "header.h"

int main(int argc, char **argv) {

	int rank;
	int size;
	double alltime;
	double *b;
	double *x;
	int* perest;
	FILE *f;
	double temp1;
	int temp4;
	char* file;

	if (argc < 3 || argc > 5) {
		fprintf(stderr, "Incorrect parametrs 1\n");
		return 0;
	}
	
	double * print_buf;
	double *matrix_buffer;
	double* output;
	double *row_buffer;
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	 
	int matrix_size = atoi(argv[1]);
	int output_number = atoi(argv[2]);
	int formula_number = atoi(argv[3]);

	if (argc < 3 || argc > 5 || output_number < -1 || formula_number < 0 || formula_number > 4) {
		fprintf(stderr, "incorrect parametrs 1\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	file = NULL;

	if (formula_number == 0)
	{
		file = argv[4];
		f = fopen(file, "r");

		if (f == NULL)
		{
			fprintf(stderr, "Can not open file\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
	}
	
	if (rank == 0)
	{
		int row_size = matrix_size / size;

		if (!(print_buf = (double *)malloc(matrix_size*(row_size + matrix_size % size) * sizeof(double))))
		{
			fprintf(stderr, "Memory allocation error\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		if (!(x = (double*)malloc(matrix_size * sizeof(double))))
		{
			fprintf(stderr, "Memory allocation error\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		if (!(perest = (int*)malloc(matrix_size * sizeof(int))))
		{
			fprintf(stderr, "Memory allocation error\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		for (int i = 0; i < matrix_size; i++)
			perest[i] = i;
	}

	if (rank != size - 1)
	{
		int row_size = matrix_size / size;

		if (!(row_buffer = (double *)malloc(matrix_size * sizeof(double))))
		{
			fprintf(stderr, "Memory allocation error\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		if (!(matrix_buffer = (double *)malloc(matrix_size*row_size * sizeof(double))))
		{
			fprintf(stderr, "Memory allocation error\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		if (!(b = (double *)malloc(row_size * sizeof(double))))
		{
			fprintf(stderr, "Memory allocation error\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		
		if (!(output = (double *)malloc(row_size * sizeof(double))))
		{
			fprintf(stderr, "Memory allocation error\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
	}

	if (rank == size - 1)
	{
		int row_size = matrix_size / size;

		if (!(row_buffer = (double *)malloc(matrix_size * sizeof(double))))
		{
			fprintf(stderr, "Memory allocation error\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		if (!(matrix_buffer = (double *)malloc(matrix_size*(row_size + matrix_size % size) * sizeof(double))))
		{
			fprintf(stderr, "Memory allocation error\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		if (!(b = (double *)malloc((row_size + matrix_size%size) * sizeof(double))))
		{
			fprintf(stderr, "Memory allocation error\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		if (!(output = (double *)malloc((row_size + matrix_size%size)* sizeof(double))))
		{
			fprintf(stderr, "Memory allocation error\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
	}

	if (output_number > matrix_size)
		output_number = matrix_size;

	if (formula_number != 0)
		fill_matrix(matrix_size, matrix_buffer, size, rank, b, formula_number);
	if (formula_number == 0)
		fill_from_file(matrix_buffer, matrix_size, file, rank, size, b);

	int currentproc;
	int countrows;

	if (rank != size - 1)
		countrows = (rank + 1) * (matrix_size / size) - rank * (matrix_size / size);
	else
		countrows = (rank + 1) * (matrix_size / size) - rank * (matrix_size / size) + matrix_size % size;

	currentproc = output_number/(matrix_size/size);

	if (output_number != 0 && output_number != -1)
	{
		if (output_number == matrix_size) {
			if (rank != 0)
				MPI_Send(matrix_buffer, countrows * matrix_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			if (rank == 0)
				print_matrix(matrix_buffer, matrix_size, size, rank, print_buf, output_number, 1);
		}
		else {
			if (rank != 0 && rank <= currentproc)
				MPI_Send(matrix_buffer, countrows * matrix_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			if (rank == 0)
				print_matrix(matrix_buffer, matrix_size, size, rank, print_buf, output_number, 1);
		}
	}

	alltime = MPI_Wtime();
	solving_func(matrix_buffer, matrix_size, rank, size, b, row_buffer, perest);
	alltime = MPI_Wtime() - alltime;

	if (output_number != 0 && output_number != -1)
	{
		if (output_number == matrix_size) {
			if (rank != 0)
				MPI_Send(matrix_buffer, countrows * matrix_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			if (rank == 0)
				print_matrix(matrix_buffer, matrix_size, size, rank, print_buf, output_number, 0);
		}
		else {
			if (rank != 0 && rank <= currentproc)
				MPI_Send(matrix_buffer, countrows * matrix_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			if (rank == 0)
				print_matrix(matrix_buffer, matrix_size, size, rank, print_buf, output_number, 0);
		}
	}

	if (rank == 0)
	{
		int rows;

		if (rank != size - 1)
			rows = (rank + 1) * (matrix_size / size) - rank * (matrix_size / size);
		else
			rows = (rank + 1) * (matrix_size / size) - rank * (matrix_size / size) + matrix_size % size;

		for (int i = 0; i < rows; i++)
		{
			x[i] = b[i] / matrix_buffer[(i)*matrix_size + i];
		}
	}
	
	if (rank == 0 && output_number != 0 && output_number != -1)
		printf("Result:\n");

	for (int w = 0; w < size; w++)
	{
		int rows;

		if (rank != size - 1)
			rows = (rank + 1) * (matrix_size / size) - rank * (matrix_size / size);
		else
			rows = (rank + 1) * (matrix_size / size) - rank * (matrix_size / size) + matrix_size % size;

		if (rank == w)
		{
			for (int j = 0; j < rows; j++)
			{
				double send = b[j] / matrix_buffer[(j)*matrix_size + j + (matrix_size / size) * rank];
				if (w != 0)
					MPI_Send(&send, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			}
		}
	}

	if (rank == 0)
	{	
		double tmp;
		for (int j = 1; j < size; j++)
		{
			int rows;

			if (j != size - 1)
				rows = (j + 1) * (matrix_size / size) - j * (matrix_size / size);
			else
				rows = (j + 1) * (matrix_size / size) - j * (matrix_size / size) + matrix_size % size;

			for (int k = 0; k < rows; k++)
			{
				MPI_Recv(&tmp, 1, MPI_DOUBLE, j, 0, MPI_COMM_WORLD, &status);
				x[k + (matrix_size / size) * j] = tmp;
			}
		}
			
		for (int i = 0; i < matrix_size; i++)
		{
			for (int j = 0; j < matrix_size - 1; j++)
			{
				if (perest[j] > perest[j + 1])
				{
					temp1 = x[j];
					x[j] = x[j + 1];
					x[j + 1] = temp1;

					temp4 = perest[j];
					perest[j] = perest[j + 1];
					perest[j + 1] = temp4;
				}
			}
		}
			
		if (output_number == -1)
		{
			for (int w = 0; w < matrix_size; w++)
				printf(" %10.3e ", x[w]);
			printf("\n");
		}
		else
		{
			for (int w = 0; w <output_number; w++)
				printf(" %10.3e ", x[w]);
			printf("\n");
		}
	}

	if (formula_number != 0)
		fill_matrix(matrix_size, matrix_buffer, size, rank, b, formula_number);
	if (formula_number == 0)
		fill_from_file(matrix_buffer, matrix_size, file, rank, size, b);
	
	residual(matrix_buffer, matrix_size, x, b, rank, size, output_number, output, row_buffer);

	if (rank == 0 && output_number != -1)
		printf("Time: %.2f sec\n", alltime);

	if (rank == 0)
	{
		free(print_buf);
		free(x);
		free(perest);
	}

	free(matrix_buffer);
	free(row_buffer);
	free(b);
	free(output);

	MPI_Finalize();
	return 0;
}
