#ifndef HEADER_H
#define HEADER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int max(int i, int j);
void print_row(int n, int row, double *row_buffer);
void read_row(int n, int row, double *row_buffer, int k);
void fill_matrix(int n, double *matrix_buffer, int size, int rank, double *b, int k);
void print_matrix(double *matrix_buffer, int n, int size, int rank, double *print_buf, int m,int l);
int fill_from_file(double *matrix_buffer, int n, char* file, int rank,int size,double *b);
void sum1(double* a, double* b, int n, int i, int j, double lambda);
void sum2(double* a, double* b, int n, int i, int j, double lambda, double* row_buffer, double btemp);
void solving_func(double* matrix_buffer, int n, int rank, int size, double* b, double* row_buffer, int* perest);
void residual(double* matrix_buffer, int n, double* x, double* b, int rank, int size, int m, double* output, double* row_buffer);
void multiply(double* matrix_buffer, int n, int size, int rank, double* x, double* row_buffer);

#endif 