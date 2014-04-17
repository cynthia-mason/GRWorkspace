
#include <stdlib.h>
#include <string.h>
#include <stdio.h>


  int main(){
  
  const int numRowA=17;	//Pass in
  const int numColA=10;	//Pass in
  const int numNZA=29;	//Pass in
  double A_val[29] = {1, 1, -1, 1, 1, -1, 1, 1, -1, 1, -1, 1, 1, -1, 1, 1, -1, 1, 1, -1, 1, 1, -1, 1, 1, -1, 1, 1, -1};	//Pass in
  int A_indx[29] ={1, 4, 8, 1, 5, 9, 1, 7, 10, 1, 11, 2, 4, 12, 2, 5, 13, 2, 6, 14, 3, 4, 15, 3, 7, 16, 3, 6, 17};	//Pass in
  int A_jndx[29] = {1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10};	//Pass in
  double b[17] = {1,1,1,2,2,1,1,0,0,0,0,0,0,0,0,0,0};	//Pass in
  double c[10] ={10,8,6,4,5,4,13,16,7,8};
  double gamma = 0.9;	//Pass in
  double x[10];
  int* NumVar;
  int* IPMit;
  int* CondNum;
  int* PCGit;

  int i;
  
  //Allocate an int pointee
  NumVar = malloc(sizeof(int));
  IPMit = malloc(sizeof(int));
  CondNum = malloc(sizeof(int));
  PCGit = malloc(sizeof(int));

  Function_PCG_WoodScratch(numRowA, numColA, numNZA, A_val, A_indx, A_jndx, b, c, gamma, x, NumVar, IPMit, CondNum, PCGit);
  //Function_PCG_WoodScratch(const int numRowA, const int numColA, const int numNZA, double[] Aval, int[] A_indx, int[] A_jndx, double[] b, double[] c, double gamma, double[] x, int *NumVar, int *IPMit, double *CondNum, int *PCGit)

  printf("\n\n x = \n");
  for(i=0;i<10;i++){
      printf("element %d = %f\n", i+1,x[i]);
  }

  printf("This worked!\n\n");
	return 0;
	}
