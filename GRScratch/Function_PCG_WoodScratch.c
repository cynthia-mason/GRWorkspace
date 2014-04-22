#include "math.h"
#include "blas_sparse.h"
#include "stdio.h"

int BASE = 0;
void Function_PCG_WoodScratch(
		const int numRowA_tree,const int numRowA_track,const int numRowA_nonneg,
		const int numColA,const int numNZA_tree,const int numNZA_track,const int numNZA_nonneg,
		double A_val[29], int A_indx[29], int A_jndx[29],
		double A_tree_val[10], int A_tree_indx[10], int A_tree_jndx[10],
		double A_track_val[9], int A_track_indx[9], int A_track_jndx[9],
		double A_nonneg_val[10], int A_nonneg_indx[10], int A_nonneg_jndx[10],
		double b_tree[3],double b_track[4],double b_nonneg[10],
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
	blas_sparse_matrix A;
	blas_sparse_matrix A_tree;
	blas_sparse_matrix A_track;
	blas_sparse_matrix A_nonneg;
/*	blas_sparse_matrix D_treeBLAS;
	blas_sparse_matrix D_trackBLAS;
	blas_sparse_matrix D_nonnegBLAS;*/
	double D_treeBLAS[numRowA_tree][numRowA_tree];
	double D_trackBLAS[numRowA_track][numRowA_track];
	double D_nonnegBLAS[numRowA_nonneg][numRowA_nonneg];
	blas_sparse_matrix Minv;
	int test_true;
	int numRowA=numRowA_tree+numRowA_track+numRowA_nonneg;
	double row[numColA], dx_1[numColA],w_temp[numRowA],w[numRowA];
	double w_temp_tree[numRowA_tree],w_temp_track[numRowA_track],w_temp_nonneg[numRowA_nonneg];
	double w_tree[numRowA_tree], D_tree[numRowA_tree], w_track[numRowA_track], D_track[numRowA_track], w_nonneg[numRowA_nonneg], D_nonneg[numRowA_nonneg];
	int D_tree_indices[numRowA_tree],D_track_indices[numRowA_track],D_nonneg_indices[numRowA_nonneg];
	double sum=0, epsilon;
	int check;
	int i;
	double D[numColA];

	//double x_test[9]={1,1,1,1,1,1,1,1,1}, w_test[9]={0,0,0,0,0,0,0,0,0};


	// ------------- CONSTRUCT A's -------------


	/*  Step 1: Create Spares BLAS Handle  */
	A=BLAS_duscr_begin(numRowA,numColA);
	A_tree=BLAS_duscr_begin(numRowA_tree,numColA);
	A_track=BLAS_duscr_begin(numRowA_track,numColA);
	A_nonneg=BLAS_duscr_begin(numRowA_nonneg,numColA);

	//Wanted to have sparse diagonal matrices, but BLAS won't allow multiplication between 2 sparse matrices.
/*	D_treeBLAS=BLAS_duscr_begin(numRowA_tree,numRowA_tree);
	D_trackBLAS=BLAS_duscr_begin(numRowA_track,numRowA_track);
	D_nonnegBLAS=BLAS_duscr_begin(numRowA_nonneg,numRowA_nonneg);*/

	/*  Step 2: Insert entries  */
	for (i=0;i<(numNZA_tree+numNZA_track+numNZA_nonneg);i++){
			BLAS_duscr_insert_entry(A,A_val[i],A_indx[i],A_jndx[i]);
			//printf("At i = %d, and j = %d ==> A[%d] = %f \n",A_indx[i], A_jndx[i], i,A_val[i]);
	}
	for (i=0;i<numNZA_tree;i++){
		BLAS_duscr_insert_entry(A_tree,A_tree_val[i],A_tree_indx[i],A_tree_jndx[i]);
		//printf("At i = %d, and j = %d ==> A[%d] = %f \n",A_indx[i], A_jndx[i], i,A_val[i]);
	}
	for (i=0;i<numNZA_track;i++){
		BLAS_duscr_insert_entry(A_track,A_track_val[i],A_track_indx[i],A_track_jndx[i]);
	}
	for (i=0;i<numNZA_nonneg;i++){
		BLAS_duscr_insert_entry(A_nonneg,A_nonneg_val[i],A_nonneg_indx[i],A_nonneg_jndx[i]);
	}

	/*  Step 3: Complete construction of sparse matrix  */
	BLAS_duscr_end(A);
	BLAS_duscr_end(A_tree);
	BLAS_duscr_end(A_track);
	BLAS_duscr_end(A_nonneg);
	// -------------               -------------

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

	printf("Epsilon: %f",epsilon);


	// ----------- CALCULATIONS ------------

//	for(check=0;check<20;check++){
		printf("\n--- Interior point iteration: %d ---\n",check);

		//Step 2
		for(i=0;i<numRowA;i++){
			w_temp[i] = 0.0;
		}

		BLAS_dusmv(blas_no_trans,1.0,A_tree,x,1,w_temp_tree,1);
		BLAS_dusmv(blas_no_trans,1.0,A_track,x,1,w_temp_track,1);
		BLAS_dusmv(blas_no_trans,1.0,A_nonneg,x,1,w_temp_nonneg,1);
		//int BLAS_dusmv( enum blas_trans_type transa, double alpha,blas_sparse_matrix A, const double *x, int incx, double *y, int incy );

		for(i=0;i<numRowA_tree;i++){
			w_tree[i]=b_tree[i]-w_temp_tree[i];
			D_tree[i]=1/w_tree[i];
			D_tree_indices[i]=i;
		}
		for(i=0;i<numRowA_track;i++){
			w_track[i]=b_track[i]-w_temp_track[i];
			D_track[i]=1/w_track[i];
			D_track_indices[i]=i;
		}
		for(i=0;i<numRowA_nonneg;i++){
			w_nonneg[i]=b_nonneg[i]-w_temp_nonneg[i];
			D_nonneg[i]=1/w_nonneg[i];
			D_nonneg_indices[i]=i;
		}

	/*	//Construct Diagonal matrices
		for (i=0;i<numRowA_tree;i++){
				BLAS_duscr_insert_entry(D_treeBLAS,D_tree[i],D_tree_indices[i],D_tree_indices[i]);
				//printf("At i = %d, and j = %d ==> A[%d] = %f \n",A_indx[i], A_jndx[i], i,A_val[i]);
		}
		for (i=0;i<numRowA_track;i++){
			BLAS_duscr_insert_entry(D_trackBLAS,D_track[i],D_track_indices[i],D_track_indices[i]);
		}
		for (i=0;i<numRowA_nonneg;i++){
			BLAS_duscr_insert_entry(D_nonnegBLAS,D_nonneg[i],D_nonneg_indices[i],D_nonneg_indices[i]);
		}

		BLAS_duscr_end(D_treeBLAS);
		BLAS_duscr_end(D_trackBLAS);
		BLAS_duscr_end(D_nonnegBLAS);*/

		for (i=0;i<numRowA_tree;i++){
			for (i=0;i<numRowA_tree;i++){

			}
		}
		for (i=0;i<numRowA_track;i++){
		}
		for (i=0;i<numRowA_nonneg;i++){
		}

		//ADDA=(A_tree')*D_tree*D_tree*A_tree
		//			+(A_track')*D_track*D_track*A_track
		//			+(A_nonneg')*D_nonneg*D_nonneg*A_nonneg;

		/*------CONSTRUCT ADDA------*/

		//(A_tree')*D_tree*D_tree*A_tree
/*
		printf("w:     w_temp:       D:\n");
		for(i=0;i<numRowA;i++){
			w[i]=b[i]-w_temp[i];
			D[i]=1/w[i];
			printf("%lf , %lf , %lf\n",w[i], w_temp[i],D[i]);
		}
*/



//	}

//	/*  Step 4: Compute matrix vector product y=A*x  */
//	BLAS_dusmv(blas_no_trans,alpha,A,x,1,y,1);
//
//	/*  Step 5: release matrix handle  */
//	BLAS_usds(A);

}
