#include "header.h"

void residual(double *matrix_buffer, int n, double *x, double *b, int rank, int size, int m, double *output, double *row_buffer) {

	int rows;

	if (rank != size - 1)
		rows = (rank + 1) * (n / size) - rank * (n / size);
	else
		rows = (rank + 1) * (n / size) - rank * (n / size) + n % size;

	for (int i = 0; i < size; i++)
	{
		MPI_Status status;
		int rows;

		if (rank != size - 1)
			rows = (rank + 1) * (n / size) - rank * (n / size);
		else
			rows = (rank + 1) * (n / size) - rank * (n / size) + n % size;

		if (rank == i)
		{				
			if (rank != 0)
			{
				MPI_Recv(row_buffer, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
				
				for (int w = 0; w < rows; w++)
				{
					double sum = 0.0;

					for (int r = 0; r < n; r++)
						sum = sum + matrix_buffer[w * n + r] * row_buffer[r];
					output[w] = sum;
				}
			}

			for (int j = 0; j < rows; j++)
			{
				output[j] = output[j] - b[j];
			}

			if (rank == 0)
			{
				for (int i = 1; i < size; i++)
				{
					MPI_Send(x, n, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
				}
			}
		}
	}	

	for (int i = 1; i < size; i++)
	{
		if (rank == i)
		{
			double sumb = 0.0;
			double summ = 0.0;

			for (int j = 0; j < rows; j++)
			{
				summ = summ + output[i] * output[i];
				sumb = sumb + b[i] * b[i];
			}

			MPI_Send(&summ, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			MPI_Send(&sumb, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
		}
	}	

	if (rank == 0)
	{
		double normb;
		MPI_Status status;
		double temp;
		double sum;

		normb = 0.0;
		sum = 0.0;

		for (int w = 0; w < rows; w++)
		{
			double sum1 = 0.0;
			for (int r = 0; r < n; r++)
				sum1 = sum1 + matrix_buffer[w * n + r] * x[r];
			output[w] = sum1;
		}

		for (int i = 0; i < rows; i++)
		{
			normb = normb + b[i] * b[i];
			output[i] = output[i] - b[i];
		}

		for (int i = 0; i < rows; i++)
		{
			sum = sum + output[i] * output[i];
		}

		for (int i = 1; i < size; i++)
		{			
			MPI_Recv(&temp, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
			sum = sum + temp;		
		}

		sum = sqrt(sum);/* ||Ax-b|| */

		for (int i = 1; i < size; i++)
		{			
			MPI_Recv(&temp, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &status);
			normb = normb + temp;
		}

		double error = 0.0;

		for (int w = 0; w < n; w++)
		{
			row_buffer[w] = 0.0;
			if (w % 2 == 0)
				row_buffer[w] = x[w] - 1.0;
			else
				row_buffer[w] = x[w] - 0.0;
			error = error + row_buffer[w] * row_buffer[w];
		}

		normb = sqrt(normb);/* ||b|| */

		if (m != -1) {
			printf("Residual: %10.3e \n", sum / normb);
			printf("Error: %10.3e \n", sqrt(error));
		}
	}		
}

void multiply(double *matrix_buffer, int n, int size, int rank, double *x,double *row_buffer) {

	int rows;
	double sum;

	sum = 0.;

	if (rank != size - 1)
		rows = (rank + 1) * (n / size) - rank * (n / size);
	else
		rows = (rank + 1) * (n / size) - rank * (n / size) + n % size;

	for (int i = 0; i < rows;i++)
	{
		for (int j = 0; j < n; j++)
		{
			sum = sum + matrix_buffer[i * n + j] * x[j];
		}
		row_buffer[i] = sum;
	}
}

void sum1(double *a, double* b, int n, int i, int j, double lambda) {


	for (int k = j; k < n; k++)
	{
		a[i * n + k] = a[i * n + k] + a[j * n + k] * lambda;
	}

	b[i] = b[i] + b[j] * lambda;
}

void solving_func(double *matrix_buffer, int n, int rank, int size, double *b, double *row_buffer, int *perest) {

	int temp5;
	int column_number;
	int current;
	int rows;

	if (rank != size - 1)
		rows = (rank + 1) * (n / size) - rank * (n / size);
	else
		rows = (rank + 1) * (n / size) - rank * (n / size) + n % size;

	for (int i = 0; i < n-1; i++)
	{
		for (int j = 0; j < n; j++)
		{		
			if ((i >= (n / size) * j && i < (n / size) * (j + 1) && (j < size - 1)) || ((j == size - 1) && (i >= (n / size) * j && i < (n / size) * (j + 1) + n % size)))
			{
				current = j;
				break;
			}
		}
		
		if (rank == current)
		{	
			int number = i;
			int local_row;
		
			for (int t = 0; t < rows; t++)
			{
				if (i == current * (n / size) + t)
				{
					local_row = t;
					break;
				}
			}
					
			double max = fabs(matrix_buffer[local_row * n + i]);
			for (int t = i + 1; t < n; t++)
			{
				if (max < fabs(matrix_buffer[n * local_row + t]))
				{
					max = fabs(matrix_buffer[n * local_row + t]);
					number = t;
				}
			}
	
			for (int k = 0; k < rows; k++)
			{
				double temp = matrix_buffer[k * n + number];
				matrix_buffer[k * n + number] = matrix_buffer[k * n + i];
				matrix_buffer[k * n + i] = temp;
			}

			if (rank == 0)
			{
				temp5 = perest[i];
				perest[i] = perest[number];
				perest[number] = temp5;
			}
			
			for (int t = 0; t < rows; t++)
			{
				if (t == local_row)
					continue;
				if (fabs(matrix_buffer[local_row * n + i]) < 1e-20)
				{
					fprintf(stderr, "Error a_ii < 1e-20 \n");
					MPI_Abort(MPI_COMM_WORLD, 1);
				}
				sum1(matrix_buffer, b, n, t, local_row, -matrix_buffer[t * n + i] / matrix_buffer[local_row * n + i]);
			}

			for (int w = 0; w < size; w++)
			{
				if (w == current)
					continue;
				MPI_Send(&number, 1, MPI_INT, w, 0, MPI_COMM_WORLD);
				MPI_Send(matrix_buffer + n * local_row, n, MPI_DOUBLE, w, 1, MPI_COMM_WORLD);
				MPI_Send(&b[local_row], 1, MPI_DOUBLE, w, 2, MPI_COMM_WORLD);
			}
		}

		if (rank != current)
		{
			double btemp;
			MPI_Status status;
			MPI_Recv(&column_number, 1, MPI_INT, current, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(row_buffer, n, MPI_DOUBLE, current, 1, MPI_COMM_WORLD, &status);
			MPI_Recv(&btemp, 1, MPI_DOUBLE, current, 2, MPI_COMM_WORLD, &status);

			for (int k = 0; k < rows; k++)
			{
				double temp = matrix_buffer[k * n + column_number];
				matrix_buffer[k * n + column_number] = matrix_buffer[k * n + i];
				matrix_buffer[k * n + i] = temp;
			}

			if (rank == 0)
			{
				temp5 = perest[i];
				perest[i] = perest[column_number];
				perest[column_number] = temp5;
			}

			for (int k = 0; k < rows; k++)
			{
				sum2(matrix_buffer, b, n, k, i, -matrix_buffer[k * n + i] / row_buffer[0 * n + i], row_buffer, btemp);
			}
		}
	}
	
	if (rank == size - 1)
	{
		
		for (int i = 0; i < rows-1; i++)
		{
			b[i] = b[i] + (-1) * b[rows - 1] * matrix_buffer[i * n + n - 1] / matrix_buffer[(rows - 1) * n + n - 1];
			matrix_buffer[i * n + n - 1] = matrix_buffer[i * n + n - 1] + matrix_buffer[(rows - 1) * n + n - 1] * (-1) * matrix_buffer[i * n + n - 1] / matrix_buffer[(rows - 1) * n + n - 1];
		}

		for (int w = 0; w < size - 1; w++)
		{
			MPI_Send(&matrix_buffer[(rows - 1) * n + n - 1], 1, MPI_DOUBLE, w, 0, MPI_COMM_WORLD);
			MPI_Send(&b[rows - 1], 1, MPI_DOUBLE, w, 1, MPI_COMM_WORLD);
		}
	}
	
	if (rank != size - 1)
	{
		double btemp;
		double temp;
		MPI_Status status;

		MPI_Recv(&temp, 1, MPI_DOUBLE, size -1 , 0,MPI_COMM_WORLD, &status);
		MPI_Recv(&btemp, 1, MPI_DOUBLE, size - 1,1, MPI_COMM_WORLD, &status);

		for (int i = 0; i < rows; i++)
		{
			b[i] = b[i] + (-1.0) * btemp * matrix_buffer[i * n + n - 1] / temp;
			matrix_buffer[i * n + n - 1] = matrix_buffer[i * n + n - 1] + (-1.0) * temp * matrix_buffer[i * n + n - 1] / temp;
		}
	}
}

void sum2(double *a, double* b, int n, int i, int j, double lambda, double *row_buffer, double btemp) {

	for (int k = j; k < n; k++)
	{
		a[i * n + k] = a[i * n + k] + row_buffer[k] * lambda;
	}
	b[i] = b[i] + btemp * lambda;
}
