#include "math.h"
#include "blas_sparse.h"
#include "stdio.h"

int BASE = 1;
void Function_PCG_WoodScratch(
		const int numRowA, const int numColA,const int numNZA,
		double A_val[29], int A_indx[29], int A_jndx[29],
		double b[17], double c[10], double gamma,
		double x[10], int *NumVar, int *IPMit, double *CondNum, int *PCGit){

/*		const int N = 4;
		const int nx = 6;
		double val[] = {1.1,2.2,2.4,3.3,4.1,4.4};
		int indx[] = {0,1,1,2,3,3};
		int jndx[] = {0,1,3,2,0,3};
		double x[] = {1.0,1.0,1.0,1.0};
		double y[] = {0.0,0.0,0.0,0.0}; // this is the solution

		blas_sparse_matrix A;
		int i;
		double alpha = 1.0;

		// Step 1: Create Sparse BLAS Handle
		A = BLAS_duscr_begin(N,N);

		//Step 2: Insert Entries
		for (i=0;i<nx;i++)
			BLAS_duscr_insert_entry(A,val[i],indx[i],jndx[i]);

		//Step 3: Complete construction of sparse matrix
		BLAS_duscr_end(A);

		//Step 4: Compute matrix vector product y=A*x
		BLAS_dusmv(blas_no_trans,alpha,A,x,1,y,1);

		//Step 5: release matrix handle
		BLAS_usds(A);

		printf("The solution is:\n");
		for (i=0;i<N;i++)
		printf("%f, ",y[i]);*/

	/*	INITIALIZE	*/
	blas_sparse_matrix Minv;
	int test_true, bool_row_v;
	double row[numColA], dx_1[numColA];
	double sum=0, epsilon;


	//CONSTRUCT A
	blas_sparse_matrix A;
	int i;

	/*  Step 1: Create Spares BLAS Handle  */
	A=BLAS_duscr_begin(numRowA,numColA);

	/*  Step 2: Insert entries  */
	for (i=0;i<numNZA;i++){
		BLAS_duscr_insert_entry(A,A_val[i],A_indx[i],A_jndx[i]);
		//printf("At i = %d, and j = %d ==> A[%d] = %f \n",A_indx[i], A_jndx[i], i,A_val[i]);
	}

	/*  Step 3: Complete construction of sparse matrix  */
	BLAS_duscr_end(A);


	/* generate "row" and sum row for creation of "epsilon"*/
	for(i=0;i<numColA;i++)
		row[i]=0.0;


	test_true = 0+BASE;
	for(i=0;i<numNZA;i++){
		if(A_indx[i]==test_true){
			row[A_jndx[i]-BASE] = A_val[i];
			sum+=A_val[i];
		}
	}
	epsilon = 1/(2*sum);
	/*-----------------*/

	for(i=0;i<numColA;i++){
		x[i] = epsilon;
		dx_1[i] = epsilon;
	}

	*NumVar = numColA;

	//Find where the components of A begin and end [A_tree;A_edge;A_nonneg]
	//1) The Tree constraints end
	test_true = numColA+BASE;
	bool_row_v = numRowA;
	for(i=0;i<numNZA;i++){
		// If the element is in the last column and has the smallest number row
		if(A_jndx[i]==test_true && A_indx[i]<bool_row_v){
			bool_row_v = A_indx[i];
		}
	}
	printf()
	//array123_row_v=1:1:numel(bool_row_v);
	//row_v=array123_row_v(bool_row_v);//Find the last row of net constratints (A must be constructed such that there are blocks of cascading ones)

//	/*  Step 4: Compute matrix vector product y=A*x  */
//	BLAS_dusmv(blas_no_trans,alpha,A,x,1,y,1);
//
//	/*  Step 5: release matrix handle  */
//	BLAS_usds(A);

}
