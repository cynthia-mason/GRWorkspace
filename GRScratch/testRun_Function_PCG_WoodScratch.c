
#include <stdlib.h>
#include <string.h>
#include <stdio.h>


  int main(){
  
  const int numRowA_tree=3;	//Pass in
  const int numRowA_track=4;	//Pass in
  const int numRowA_nonneg=10;	//Pass in

  const int numColA=10;	//Pass in

  const int numNZA_tree=10;	//Pass in
  const int numNZA_track=9;	//Pass in
  const int numNZA_nonneg=10;	//Pass in


  double A_tree_val[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};	//Pass in
  int A_tree_indx[] ={1, 1, 1, 1, 2, 2, 2, 3, 3, 3};	//Pass in
  int A_tree_jndx[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};	//Pass in

  double A_track_val[] = {1, 1, 1, 1, 1, 1, 1, 1, 1};	//Pass in
  int A_track_indx[] ={1, 2, 4, 1, 2, 3, 1, 4, 3};	//Pass in
  int A_track_jndx[] = {1, 2, 3, 5, 6, 7, 8, 9, 10};	//Pass in

  double A_nonneg_val[] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1};	//Pass in
  int A_nonneg_indx[] ={1, 2, 3, 4, 5, 6, 7, 8, 9, 10};	//Pass in
  int A_nonneg_jndx[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};	//Pass in

  double A_val[] = {1,1,-1,1,1,-1,1,1,-1,1,-1,1,1,-1,1,1,-1,1,1,-1,1,1,-1,1,1,-1,1,1,-1};	//Pass in
  int A_indx[] ={1,4,8,1,5,9,1,7,10,1,11,2,4,12,2,5,13,2,6,14,3,4,15,3,7,16,3,6,17};	//Pass in
  int A_jndx[] = {1,1,1,2,2,2,3,3,3,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10};	//Pass in

  int i=0;

  for(i=0;i<10;i++){
	  A_tree_indx[i]-=1;
	  A_tree_jndx[i]-=1;
	  A_nonneg_indx[i]-=1;
	  A_nonneg_jndx[i]-=1;
  }
  for(i=0;i<9;i++){
	  A_track_indx[i]-=1;
	  A_track_jndx[i]-=1;
  }
  for(i=0;i<29;i++){
	  A_indx[i]-=1;
	  A_jndx[i]-=1;
  }

/*
    double A_val[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};	//Pass in
    int A_indx[9] ={1, 2, 3, 4, 5, 6, 7, 8, 9};	//Pass in
    int A_jndx[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};	//Pass in
*/



  double b[] = {1,1,1,2,2,1,1,0,0,0,0,0,0,0,0,0,0};	//Pass in
  double b_tree[] = {1,1,1};	//Pass in
  double b_track[] = {2,2,1,1};	//Pass in
  double b_nonneg[] = {0,0,0,0,0,0,0,0,0,0};	//Pass in
  double c[] ={10,8,6,4,5,4,13,16,7,8};
  double gamma = 0.9;	//Pass in
  double x[10];
  int* NumVar;
  int* IPMit;
  int* CondNum;
  int* PCGit;


  
  //Allocate an int pointee
  NumVar = malloc(sizeof(int));
  IPMit = malloc(sizeof(int));
  CondNum = malloc(sizeof(int));
  PCGit = malloc(sizeof(int));

  Function_PCG_WoodScratch(numRowA_tree,numRowA_track,numRowA_nonneg, numColA,
		  numNZA_tree,numNZA_track,numNZA_nonneg,
		  A_val, A_indx, A_jndx,
		  A_tree_val, A_tree_indx, A_tree_jndx,
		  A_track_val, A_track_indx, A_track_jndx,
		  A_nonneg_val, A_nonneg_indx, A_nonneg_jndx,
		  b_tree,b_track,b_nonneg,
		  b, c, gamma, x, NumVar, IPMit, CondNum, PCGit);
  //Function_PCG_WoodScratch(const int numRowA, const int numColA, const int numNZA, double[] Aval, int[] A_indx, int[] A_jndx, double[] b, double[] c, double gamma, double[] x, int *NumVar, int *IPMit, double *CondNum, int *PCGit)

  printf("\n\n x = \n");
  for(i=0;i<10;i++){
      printf("element %d = %f\n", i+1,x[i]);
  }

  printf("This worked!\n\n");
	return 0;
	}
