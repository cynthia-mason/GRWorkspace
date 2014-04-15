/* hsl_ma97ds2.c */
/* To illustrate use of advanced features of hsl_ma97 */

#include <stdio.h>
#include <stdlib.h>
#include "hsl_ma97s.h"

int main(void) {
   typedef float pkgtype;
   
   void *akeep, *fkeep;
   struct ma97_control control;
   struct ma97_info info;

   int *ptr, *row, *order, *flag_out, *bindex, *xindex;
   pkgtype *val, *x, *y, *b;

   int i,j,matrix_type,n,ne,check,nrhs,nbi,nxi;

   /* Read in the order n of the matrix and number of entries in lwr triangle */
   scanf("%d %d %d", &n, &ne, &nrhs);

   /* Allocate arrays for matrix data and arrays for hsl_ma97 */
   ptr = (int *) malloc((n+1)*sizeof(int));
   row = (int *) malloc(ne*sizeof(int));
   val = (pkgtype *) malloc(ne*sizeof(pkgtype));
   x = (pkgtype *) malloc(2*nrhs*n*sizeof(pkgtype));
   y = (pkgtype *) malloc(nrhs*n*sizeof(pkgtype));
   order = (int *) malloc(n*sizeof(int));
   flag_out = (int *) malloc(nrhs*sizeof(int));

   for(i=0; i<n+1; i++) scanf("%d", &(ptr[i]));
   for(i=0; i<ne; i++) scanf("%d", &(row[i]));
   for(i=0; i<ne; i++) scanf("%f", &(val[i]));

   ma97_default_control(&control);

   /* Perform analyse and factorise with data checking */
   check = 1; /* true */
   ma97_analyse(check,n,ptr,row,NULL,&akeep,&control,&info,order);
   if (info.flag < 0) {
      ma97_free_akeep(&akeep);
      free(ptr); free(row); free(val); free(x); free(y); free(flag_out);
      free(order);
      return 1;
   }
   matrix_type = 4; /* Real, symmetric indefinite */
   fkeep = NULL; /* important that this is initialised to NULL on first call */
   ma97_factor(matrix_type,ptr,row,val,&akeep,&fkeep,&control,&info,NULL);
   if (info.flag < 0) {
      ma97_finalise(&akeep,&fkeep);
      free(ptr); free(row); free(val); free(x); free(y); free(flag_out);
      free(order);
      return 1;
   }

   /* Read in the right-hand sides. */
   for(i=0; i<n*nrhs; i++) scanf("%f", &(x[i]));

   /* Solve with Fredholm alternative if rhs is inconsistent */
   ma97_solve_fredholm(nrhs,flag_out,x,n,&akeep,&fkeep,&control,&info);
   if (info.flag < 0) {
      ma97_finalise(&akeep,&fkeep);
      free(ptr); free(row); free(val); free(x); free(y); free(flag_out);
      free(order);
      return 1;
   }
   for(j=0; j<nrhs; j++) {
      if(flag_out[j]) {
         printf("Right-hand side %d is consistent with solution:\n", j);
         for(i=0; i<n; i++) printf(" %f", x[i+j*n]);
      } else {
         printf("Right-hand side %d inconsistent. Ax=0, x^Tb/=0 given by:\n",
            j);
         for(i=0; i<n; i++) printf(" %f", x[nrhs*n+i+j*n]);
      }
      printf("\n");
   }

   /* Form (S^{-1}PL)X */
   ma97_lmultiply(0,2,x,n,y,n,&akeep,&fkeep,&control,&info);
   if (info.flag < 0) {
      ma97_finalise(&akeep,&fkeep);
      free(ptr); free(row); free(val); free(x); free(y); free(flag_out);
      free(order);
      return 1;
   }
   for(j=0; j<nrhs; j++) {
      printf("S^{-1}PLX_%d = ", j);
      for(i=0; i<n; i++) printf(" %f", y[i+j*n]);
      printf("\n");
   }

   /* Read sparse right-hand side */
   scanf("%d", &nbi);
   bindex = (int *) malloc(nbi*sizeof(int));
   b = (pkgtype *) malloc(n*sizeof(pkgtype));
   xindex = (int *) malloc(n*sizeof(int));
   for(i=0; i<nbi; i++) scanf("%d", &(bindex[i]));
   for(i=0; i<nbi; i++) scanf("%f", &(b[bindex[i]]));

   /* Perform sparse fwd solve */
   ma97_sparse_fwd_solve(nbi, bindex, b, order, &nxi, xindex, x, &akeep,
      &fkeep, &control, &info);
   if (info.flag < 0) {
      ma97_finalise(&akeep,&fkeep);
      free(ptr); free(row); free(val); free(x); free(y); free(flag_out);
      free(order); free(bindex); free(b); free(xindex);
      return 1;
   }
   printf("Sparse solution has entries:\n");
   for(i=0; i<nxi; i++) printf("%d %f\n", xindex[i], x[xindex[i]]);

   ma97_finalise(&akeep,&fkeep);

   /* Deallocate all arrays */
   free(ptr); free(row); free(val); free(x); free(y); free(flag_out);
   free(order); free(bindex); free(b); free(xindex);

   return 0;
}
