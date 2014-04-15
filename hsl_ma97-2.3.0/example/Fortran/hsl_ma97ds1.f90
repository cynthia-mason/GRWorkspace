! Simple code to illustrate use of hsl_ma97
program hsl_ma97ds1
   use hsl_ma97_double
   implicit none

   ! Derived types
   type (ma97_akeep)   :: akeep
   type (ma97_fkeep)   :: fkeep
   type (ma97_control) :: control
   type (ma97_info)    :: info

   ! Parameters
   integer, parameter :: wp = kind(0.0d0)

   integer, dimension (:),  allocatable :: row
   integer, dimension (:),  allocatable :: col
   real(wp), dimension (:), allocatable :: val
   real(wp), dimension (:,:), allocatable :: x

   integer :: matrix_type,n,ne

   ! Read in the order n of the matrix and number of entries in lower triangle
   read (*,*) n,ne

   ! Allocate arrays for matrix data and arrays for hsl_ma97
   allocate (row(ne),col(ne),val(ne))
   allocate (x(n,2))

   read (*,*) row(1:ne)
   read (*,*) col(1:ne)
   read (*,*) val(1:ne)

   ! Perform analyse and factorise with coordinate input
   call ma97_analyse_coord(n,ne,row,col,akeep,control,info)
   if (info%flag < 0) go to 100
   matrix_type = 4 ! Real, symmetric indefinite
   call ma97_factor(matrix_type,val,akeep,fkeep,control,info)
   if (info%flag < 0) go to 100

   ! Read in the right-hand side
   read (*,*) x(1:n,1:2)

   ! Solve
   call ma97_solve(2,x,n,akeep,fkeep,control,info)
   if (info%flag < 0) go to 100
   write (*,'(/a,/,(3es18.10))') ' The computed solution is:', x(1:n,1)
   write (*,'(3es18.10)') x(1:n,2)

   ! Read values of second matrix with same pattern
   read (*,*) val(1:ne)

   ! Read another right hand side
   read (*,*) x(1:n,1)

   ! Perform combined factor and solve
   call ma97_factor_solve(matrix_type,val,x(1:n,1),akeep,fkeep,control,info)
   write (*,'(/a,/,(3es18.10))') ' Next solution is:', x(1:n,1)

   100 continue
   call ma97_finalise(akeep, fkeep)

end program hsl_ma97ds1
