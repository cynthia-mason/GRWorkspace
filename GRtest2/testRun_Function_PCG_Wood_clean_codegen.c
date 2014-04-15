
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "rt_nonfinite.h"
#include "Function_PCG_Wood_clean_codegen2.h"
#include "Function_PCG_Wood_clean_codegen2_emxutil.h"
#include "Function_PCG_Wood_clean_codegen2_types.h"
#include "Function_PCG_Wood_clean_codegen2_terminate.h"
#include "Function_PCG_Wood_clean_codegen2_initialize.h"

  //testRun_Function_PCG_Wood_clean_codegen
  int main(){
  
  const real_T A[170] = {1,0,0,1,0,0,0,-1,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,-1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,-1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,-1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,-1,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,-1};
  const real_T b[17] = {1,1,1,2,2,1,1,0,0,0,0,0,0,0,0,0,0};
  const real_T c[10] ={10,8,6,4,5,4,13,16,7,8};
  real_T b_gamma = 0.9;
  real_T x[10];
  real_T *NumVar;
  real_T *IPMit;
  real_T *CondNum;
  real_T *PCGit;
  int i;
  
  Function_PCG_Wood_clean_codegen2(A,b,c,b_gamma,x,&NumVar, &IPMit, &CondNum, &PCGit);
 // void Function_PCG_Wood_clean_codegen2(const real_T A[170], const real_T b[17],const real_T c[10], real_T b_gamma, real_T x[10], real_T *NumVar, real_T*IPMit, real_T *CondNum, real_T *PCGit)
  printf("x = ");
  for(i=0;i<10;i++){
      printf("element %d = %f\n", i+1,x[i]);
  }

  printf("This worked!\n\n");
	return 0;
	}
