/* hsl_ma97cs.c */
/* Simple code to illustrate entry by columns to hsl_ma97 */

#include <stdio.h>
#include <stdlib.h>
#include "hsl_ma97c.h"

int main(void) {
   typedef float complex pkgtype;
   typedef float pkgrealtype;
   
   void *akeep, *fkeep;
   struct ma97_control control;
   struct ma97_info info;

   int *ptr, *row, *piv_order;
   pkgtype *val, *x;
   pkgrealtype re, im;

   int i,matrix_type,n,ne,check;

   /* Read in the order n of the matrix and number of entries in lwr triangle */
   scanf("%d %d", &n, &ne);

   /* Allocate arrays for matrix data and arrays for hsl_ma97 */
   ptr = (int *) malloc((n+1)*sizeof(int));
   row = (int *) malloc(ne*sizeof(int));
   val = (pkgtype *) malloc(ne*sizeof(pkgtype));
   x = (pkgtype *) malloc(n*sizeof(pkgtype));
   piv_order = (int *) malloc(n*sizeof(int));

   for(i=0; i<n+1; i++) scanf("%d", &(ptr[i]));
   for(i=0; i<ne; i++) scanf("%d", &(row[i]));
   for(i=0; i<ne; i++) {
      scanf(" (%f,%f)", &re, &im);
      val[i] = re + im*I;
   }

   ma97_default_control(&control);

   /* Perform analyse and factorise with data checking */
   check = 1; /* true */
   ma97_analyse(check,n,ptr,row,NULL,&akeep,&control,&info, NULL);
   if (info.flag < 0) {
      ma97_free_akeep(&akeep);
      free(ptr); free(row); free(val); free(x); free(piv_order);
      return 1;
   }
   matrix_type = -4;
   fkeep = NULL; /* important that this is initialised to NULL on first call */
   ma97_factor(matrix_type,ptr,row,val,&akeep,&fkeep,&control,&info,NULL);
   if (info.flag < 0) {
      ma97_finalise(&akeep,&fkeep);
      free(ptr); free(row); free(val); free(x); free(piv_order);
      return 1;
   }

   /* Read in the right-hand side. */
   for(i=0; i<n; i++) {
      scanf(" (%f,%f)", &re, &im);
      x[i] = re + im*I;
   }

   /* Solve */
   ma97_solve(0,1,x,n,&akeep,&fkeep,&control,&info);
   if (info.flag < 0) {
      ma97_finalise(&akeep,&fkeep);
      free(ptr); free(row); free(val); free(x); free(piv_order);
      return 1;
   }
   printf("\nThe computed solution is:\n");
   for(i=0; i<n; i++) printf("(%5.3f,%5.3f)", creal(x[i]), cimag(x[i]));

   /* Determine the pivot order used */
   printf("\nPivot order:");
   ma97_enquire_indef(&akeep,&fkeep,&control,&info,piv_order,NULL);
   for(i=0; i<n; i++) printf(" %d", piv_order[i]);
   printf("\n");

   ma97_finalise(&akeep,&fkeep);

   /* Deallocate all arrays */
   free(ptr); free(row); free(val); free(x); free(piv_order);

   return 0;
}
