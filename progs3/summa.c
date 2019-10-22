/******************************************************************************************
*
*	Filename:	summa.c
*	Purpose:	A paritally implemented program for MSCS6060 HW. Students will complete 
*			the program by adding SUMMA implementation for matrix multiplication C = A * B.  
*	Assumptions:    A, B, and C are square matrices n by n; 
*			the total number of processors (np) is a square number (q^2).
*	To compile, use 
*	    mpicc -o summa summa.c
*       To run, use
*	    mpiexec -n $(NPROCS) ./summa
*********************************************************************************************/


#include <stdio.h>
#include <time.h>	
#include <stdlib.h>	
#include <math.h>	
#include "mpi.h"
#include <string.h>

#define min(a, b) ((a < b) ? a : b)
// #define SZ 4000		//Each matrix of entire A, B, and C is SZ by SZ. Set a small value for testing, and set a large value for collecting experimental data.


/**
*   Allocate space for a two-dimensional array
*/
double **alloc_2d_double(int n_rows, int n_cols) {
	int i;
	double **array;
	array = (double **)malloc(n_rows * sizeof (double *));
        array[0] = (double *) malloc(n_rows * n_cols * sizeof(double));
        for (i=1; i<n_rows; i++){
                array[i] = array[0] + i * n_cols;
        }
        return array;
}

/**
*	Initialize arrays A and B with random numbers, and array C with zeros. 
*	Each array is setup as a square block of blck_sz.
*	ii,jj is the index of each block
*	And i,j is the index of the whole matrix
**/

void init_test(double **lA, double **lB, double **lC, int blck_sz, int coordinates[2]){
	int i, j, ii, jj;

	for (ii=0; ii<blck_sz; ii++){
		for (jj=0; jj<blck_sz; jj++){		
			i = ii + blck_sz * coordinates[0];
			j = jj + blck_sz * coordinates[1];
			if (i == j) {
				lA[ii][jj] = 1.0;
				lB[ii][jj] = 1.0;
			}else if((j-1) == i){
				lB[ii][jj] = 1.0;
				lA[ii][jj] = 0.0;
			}else if((i-1) == j){
				lA[ii][jj] = 1.0;
				lB[ii][jj] = 0.0;
			}else{
				lA[ii][jj] =0.0;
				lB[ii][jj] =0.0;
			}
			lC[ii][jj] = 0.0;
		}
	}
}

/**
*	Initialize arrays A and B with random numbers, and array C with zeros. 
*	Each array is setup as a square block of blck_sz.
**/
void initialize(double **lA, double **lB, double **lC, int blck_sz){
	int i, j;
	double value;
	// Set random values...technically it is already random and this is redundant
	for (i=0; i<blck_sz; i++){
		for (j=0; j<blck_sz; j++){
			lA[i][j] = (double)rand() / (double)RAND_MAX;
			lB[i][j] = (double)rand() / (double)RAND_MAX;
			lC[i][j] = 0.0;
		}
	}
}

//basic SUMMA algorithm 
void multi_SUMMA(double **my_C, double **my_A, double **my_B, int block_sz) {
	for (int k = 0; k < block_sz; ++k) {
		for (int i = 0; i < block_sz; ++i) {
			for (int j = 0; j < block_sz; ++j) {
				my_C[i][j] += my_A[i][k] * my_B[k][j];
			}
		}
	}
}

/**
*	Perform the SUMMA matrix multiplication. 
*       Follow the pseudo code in lecture slides.
*/
void matmul(int my_rank, int proc_grid_sz, int block_sz, double **my_A,
						double **my_B, double **my_C, int coordinates[2], MPI_Comm grid_comm){

	//Add your implementation of SUMMA algorithm
	int remainDims[2];
	double **buff_A;
	buff_A =  alloc_2d_double(block_sz, block_sz);
	double **buff_B; 
	buff_B = alloc_2d_double(block_sz, block_sz);

	MPI_Comm row_comm;
	MPI_Comm col_comm;

	// Partitions a communicator into subgroups 
	// which form lower-dimensional cartesian subgrids
	remainDims[0] = 0;
	remainDims[1] = 1;
	MPI_Cart_sub(grid_comm, remainDims, &row_comm);

	remainDims[0] = 1;
	remainDims[1] = 0;
	MPI_Cart_sub(grid_comm, remainDims, &col_comm);

	for (int k = 0; k < proc_grid_sz; ++k) {
		if (coordinates[1] == k) {
			memcpy(buff_A[0], my_A[0], block_sz*block_sz*sizeof(double));
		}
		MPI_Bcast(*buff_A, block_sz*block_sz, MPI_DOUBLE, k, row_comm);

		if (coordinates[0] == k) {
			memcpy(buff_B[0], my_B[0], block_sz*block_sz*sizeof(double));
		}
		MPI_Bcast(*buff_B, block_sz*block_sz, MPI_DOUBLE, k, col_comm);
		
		if (coordinates[0] == k && coordinates[1] == k) {
			multi_SUMMA(my_C, my_A, my_B, block_sz);
		}else if(coordinates[0] == k){
			multi_SUMMA(my_C, buff_A, my_B, block_sz);
		}else if(coordinates[1] == k){
			multi_SUMMA(my_C, my_A, buff_B, block_sz);
		}else{
			multi_SUMMA(my_C, buff_A, buff_B, block_sz);
		}
	}
}


int main(int argc, char *argv[]) {
	int p_rank, p_num;							//process rank and total number of processes
	double start_time, end_time, total_time;	// for timing
	int block_sz;								// Block size length for each processor to handle
	int proc_grid_sz;							// 'q' from the slides
		
	srand(time(NULL));							// Seed random numbers
	
	int SZ = 0;
	
	if (argc == 2) {
		SZ = atoi(argv[1]);
	}else{
		SZ = 4000;
	}

	int wrap[2];
	int coordinates[2];
	int dim_sz[2];
	int reorder = 1;

/* insert MPI functions to 1) start process, 2) get total number of processors and 3) process rank*/

	MPI_Init(&argc, &argv);
	MPI_Comm gridComm;
	MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p_num);

/* assign values to 1) proc_grid_sz and 2) block_sz*/
	proc_grid_sz = (int)sqrt((double)p_num);
	block_sz = SZ/proc_grid_sz;
	dim_sz[0] = dim_sz[1] = proc_grid_sz;


	wrap[0] = wrap[1] = 1;

	int my_rank;

	MPI_Cart_create(MPI_COMM_WORLD, 2, dim_sz, wrap, reorder, &gridComm);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Cart_coords(gridComm, my_rank, 2, coordinates);
	if (SZ % proc_grid_sz != 0){
		printf("The matrix size cannot be evenly split!\n");
		exit(-1);
	}

	// Create the local matrices on each process

	double **A, **B, **C;
	A = alloc_2d_double(block_sz, block_sz);
	B = alloc_2d_double(block_sz, block_sz);
	C = alloc_2d_double(block_sz, block_sz);

	//Call this function to create Bi/tridigonal for testing
	init_test(A, B, C, block_sz, coordinates);
	//initialize(A, B, C, block_sz);

	// Use MPI_Wtime to get the starting time
	start_time = MPI_Wtime();


	// Use SUMMA algorithm to calculate product C
	matmul(p_rank, proc_grid_sz, block_sz, A, B, C, coordinates, gridComm);


	// Use MPI_Wtime to get the finishing time
	end_time = MPI_Wtime();


	// Obtain the elapsed time and assign it to total_time
	total_time = end_time - start_time; 
	printf("Computation Time: %f ms. \n", total_time);

	//Insert statements for testing
	int ii,jj,i,j;
	for(ii = 0; ii < block_sz; ii++) {
		for(jj = 0; jj < block_sz; jj++) {
			i = ii + block_sz * coordinates[0];
			j = jj + block_sz * coordinates[1];
			if (i == 0 && j==0 ) {	
				if (C[ii][jj]!=1) {
					printf("C[%d][%d] is incorrect\n", ii,jj);
				}
			}else if(i == j){
				if (C[ii][jj]!=2) {
					printf("C[%d][%d] is incorrect\n", ii,jj);
				}
			}else if( (i-1) == j){
				if (C[ii][jj]!=1) {
					printf("C[%d][%d] is incorrect\n", ii,jj);
				}
			}else if(i == (j-1) ){
				if (C[ii][jj]!=1) {
					printf("C[%d][%d] is incorrect\n", ii,jj);
				}
			}	
		}
	}
	printf("Test pass\n");


	if (p_rank == 0){
		// Print in pseudo csv format for easier results compilation
		printf("squareMatrixSideLength,%d,numMPICopies,%d,walltime,%lf\n",
			SZ, p_num, total_time);
	}

	// Destroy MPI processes

	MPI_Finalize();
	printf("The Program Compeleted....\n");
	return 0;
}
