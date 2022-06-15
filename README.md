# Gauss-Jordan-MPI
Open MPI 4.0.3

My parallel implementation of the Gauss-Jordan algorithm for solving a linear system Ax=b, where the matrix A is square non-degenerate.
The vector b is determined by the matrix A as follows bi is the sum of k from 0 to (n-2)/2 a_{i, 2k+1}, where n is the size of the matrix. It follows from this that the theoretical solution is a vector (1, 0, 1, 0 ...)

To compile the project

~$ make

To run the project 

~$ mpirun -np number_of_processes ./prog matrix_size matrix_output formula_number

for examle:

~$ mpirun -np 3 ./prog 5 5 0 matrix.txt

~$ mpirun -np 2 ./prog 2000 5 1


formula_number takes values from 0 to 4, if it is equal to 0, then the following argument specifies the name of the file (.txt) containing the elements of the matrix.

If matrix_output is equal to 0, then only the residual ( ‖Ax − b‖ / ‖b‖ ), error ( ‖x − (1,0,1,0,...)‖ ) and time will be output
If matrix_output is equal to -1, then only the vector of solutions will be output in the terminal. if >0, then a submatrix of this size will be output, before and after the algorithm.
