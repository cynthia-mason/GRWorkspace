/* hsl_ma97cs1.c */
/* Simple code to illustrate coordinate entry for hsl_ma97 */

#include <stdio.h>
#include <stdlib.h>
#include "hsl_ma97c.h"

int main(void) {
   typedef float complex pkgtype;
   typedef float pkgrealtype;
   
   void *akeep, *fkeep;
   struct ma97_control control;
   struct ma97_info info;

   int *row, *col;
   pkgtype *val, *x;
   pkgrealtype re, im;

   int i,matrix_type,n,ne;

   /* Read in the order n of the matrix and number of entries in lwr triangle */
   scanf("%d %d", &n, &ne);

   /* Allocate arrays for matrix data and arrays for hsl_ma97 */
   row = (int *) malloc(ne*sizeof(int));
   col = (int *) malloc(ne*sizeof(int));
   val = (pkgtype *) malloc(ne*sizeof(pkgtype));
   x = (pkgtype *) malloc(2*n*sizeof(pkgtype));

   for(i=0; i<ne; i++) scanf("%d", &(row[i]));
   for(i=0; i<ne; i++) scanf("%d", &(col[i]));
   for(i=0; i<ne; i++) {
      scanf(" (%f,%f)", &re, &im);
      val[i] = re + im*I;
   }

   ma97_default_control(&control);

   /* Perform analyse and factorise with data checking */
   ma97_analyse_coord(n,ne,row,col,NULL,&akeep,&control,&info,NULL);
   if (info.flag < 0) {
      ma97_free_akeep(&akeep);
      free(row); free(col); free(val); free(x);
      return 1;
   }
   matrix_type = -4; /* Complex, Hermitian indefinite */
   fkeep = NULL; /* important that this is initialised to NULL on first call */
   ma97_factor(matrix_type,NULL,NULL,val,&akeep,&fkeep,&control,&info,NULL);
   if (info.flag < 0) {
      ma97_finalise(&akeep,&fkeep);
      free(row); free(col); free(val); free(x);
      return 1;
   }

   /* Read in the right-hand sides. */
   for(i=0; i<2*n; i++) {
      scanf(" (%f,%f)", &re, &im);
      x[i] = re + im*I;
   }

   /* Solve */
   ma97_solve(0,2,x,n,&akeep,&fkeep,&control,&info);
   if (info.flag < 0) {
      ma97_finalise(&akeep,&fkeep);
      free(row); free(col); free(val); free(x);
      return 1;
   }
   printf("\nThe computed solution is:\n");
   for(i=0; i<n; i++) printf("(%5.3f,%5.3f)", creal(x[i]), cimag(x[i]));
   printf("\n");
   for(i=0; i<n; i++) printf("(%5.3f,%5.3f)", creal(x[n+i]), cimag(x[n+i]));

   /* Read second matrix with same pattern */
   for(i=0; i<ne; i++) {
      scanf(" (%f,%f)", &re, &im);
      val[i] = re + im*I;
   }

   /* Read another right hand side */
   for(i=0; i<n; i++) {
      scanf(" (%f,%f)", &re, &im);
      x[i] = re + im*I;
   }

   /* Perform combined factor and solve */
   /* Note: no need to set fkeep to NULL, we will overwrite existing factors */
   ma97_factor_solve(matrix_type,NULL,NULL,val,1,x,n,&akeep,&fkeep,&control,
         &info,NULL);
   printf("\nNext solution is:\n");
   for(i=0; i<n; i++) printf("(%5.3f,%5.3f)", creal(x[i]), cimag(x[i]));
   printf("\n");

   ma97_finalise(&akeep,&fkeep);

   /* Deallocate all arrays */
   free(row); free(col); free(val); free(x);

   return 0;
}
