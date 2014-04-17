#include "math.h"
#include "blas_sparse.h"
#include "stdio.h"

int BASE = 1;
void Function_PCG_WoodScratch(
		const int numRowA_tree,const int numRowA_track,const int numRowA_nonneg,
		const int numColA,const int numNZA_tree,const int numNZA_track,const int numNZA_nonneg,
		double A_tree_val[29], int A_tree_indx[29], int A_tree_jndx[29],
		double A_track_val[29], int A_track_indx[29], int A_track_jndx[29],
		double A_nonneg_val[29], int A_nonneg_indx[29], int A_nonneg_jndx[29],
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
	blas_sparse_matrix A_tree;
	blas_sparse_matrix A_track;
	blas_sparse_matrix A_nonneg;
	blas_sparse_matrix Minv;
	int test_true;
	double row[numColA], dx_1[numColA];
	double sum=0, epsilon;
	int bool_row_v[numRowA], array123_row_v[numRowA];
	int row_v[numRowA], row_v_ind;

	//CONSTRUCT A's
	blas_sparse_matrix A_tree;
	blas_sparse_matrix A_track;
	blas_sparse_matrix A_nonneg;
	int i;

	/*  Step 1: Create Spares BLAS Handle  */
	A_tree=BLAS_duscr_begin(numRowA_tree,numColA);
	A_track=BLAS_duscr_begin(numRowA_track,numColA);
	A_nonneg=BLAS_duscr_begin(numRowA_nonneg,numColA);

	/*  Step 2: Insert entries  */
	for (i=0;i<numNZA_tree;i++){
		BLAS_duscr_insert_entry(A_tree,A_tree_val[i],A_tree_indx[i],A_tree_jndx[i]);
		//printf("At i = %d, and j = %d ==> A[%d] = %f \n",A_indx[i], A_jndx[i], i,A_val[i]);
	}
	for (i=0;i<numNZA_track;i++){
		BLAS_duscr_insert_entry(A_track,A_track_val[i],A_track_indx[i],A_track_jndx[i]);
		//printf("At i = %d, and j = %d ==> A[%d] = %f \n",A_indx[i], A_jndx[i], i,A_val[i]);
	}
	for (i=0;i<numNZA_nonneg;i++){
		BLAS_duscr_insert_entry(A_nonneg,A_nonneg_val[i],A_nonneg_indx[i],A_nonneg_jndx[i]);
		//printf("At i = %d, and j = %d ==> A[%d] = %f \n",A_indx[i], A_jndx[i], i,A_val[i]);
	}

	/*  Step 3: Complete construction of sparse matrix  */
	BLAS_duscr_end(A_tree);
	BLAS_duscr_end(A_track);
	BLAS_duscr_end(A_nonneg);


	/* generate "row" and sum row for creation of "epsilon"*/
	for(i=0;i<numColA;i++)
		row[i]=0.0;


	test_true = 0+BASE;
	for(i=0;i<numNZA_tree;i++){
		if(A_tree_indx[i]==test_true){
			row[A_tree_jndx[i]-BASE] = A_tree_val[i];
			sum+=A_tree_val[i];
		}
	}
	epsilon = 1/(2*sum);
	/*-----------------*/

	for(i=0;i<numColA;i++){
		x[i] = epsilon;
		dx_1[i] = epsilon;
	}

	*NumVar = numColA;

	//------GOOD TILL HERE!!

	//Find where the components of A begin and end [A_tree;A_edge;A_nonneg]
	//1) The Tree constraints end
	test_true = numColA-1+BASE;
//	for(i=0;i<numRowA;i++){
//		bool_row_v[i]=0;
//		array123_row_v[i]=i;
//	}
//
//	for(i=0;i<numNZA;i++){
//		// If the element is in the last column and has the smallest number row
//		if(A_jndx[i]==test_true && A_val[i]==1){
//			printf("A_jndx[%d] = %d, A_val[%d] = %f, \n",i,A_jndx[i],i,A_val[i]);
//			bool_row_v[A_indx[i]-BASE] = 1;
//			printf("bool_row_v[%d] = %d, \n",A_indx[i],bool_row_v[A_indx[i]-BASE]);
//		}
//	}
//	printf("bool_row_v:\n");
//	for(i=0;i<numRowA;i++)
//		printf("%d ",bool_row_v[i]);

	for(i=0;i<numRowA;i++){
		row_v[i]=-1; //THIS IS a FLAG SO THAT IT WILL NOT BE USED AN INDICES
	}
	row_v_ind=0;
	for(i=0;i<numNZA;i++){
		// If the element is in the last column and has the smallest number row
		if(A_jndx[i]==test_true && A_val[i]==1){
			row_v[row_v_ind] = A_indx[i]-1;
			printf("A_jndx[%d] = %d, A_val[%d] = %f, \n",i,A_jndx[i],i,A_val[i]);
			printf("bool_row_v[%d] = %d, \n",row_v_ind,row_v[row_v_ind]);
			row_v_ind++;
		}
	}

	printf("row_v:\n");
	for(i=0;i<numRowA;i++)
		printf("%d ",row_v[i]);

	A_tree=BLAS_duscr_begin(row_v[1]+1,numColA);
	for (i=0;i<numNZA_tree;i++){
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


//	/*  Step 4: Compute matrix vector product y=A*x  */
//	BLAS_dusmv(blas_no_trans,alpha,A,x,1,y,1);
//
//	/*  Step 5: release matrix handle  */
//	BLAS_usds(A);

}
