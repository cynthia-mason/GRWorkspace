! Simple code to illustrate use of hsl_ma97
program hsl_ma97zs
   use hsl_ma97_double_complex
   implicit none

   ! Derived types
   type (ma97_akeep)   :: akeep
   type (ma97_fkeep)   :: fkeep
   type (ma97_control) :: control
   type (ma97_info)    :: info

   ! Parameters
   integer, parameter :: wp = kind(0.0d0)

   integer, dimension (:),  allocatable :: ptr
   integer, dimension (:),  allocatable :: piv_order
   integer, dimension (:),  allocatable :: row
   complex(wp), dimension (:), allocatable :: val
   complex(wp), dimension (:), allocatable :: x

   integer :: matrix_type,n,ne
   logical :: check

   ! Read in the order n of the matrix and number of entries in lower triangle
   read (*,*) n,ne

   ! Allocate arrays for matrix data and arrays for hsl_ma97
   allocate (ptr(n+1),row(ne),val(ne))
   allocate (x(n))

   read (*,*) ptr(1:n+1)
   read (*,*) row(1:ne)
   read (*,*) val(1:ne)

   ! Perform analyse and factorise with data checking
   check = .true.
   call ma97_analyse(check,n,ptr,row,akeep,control,info)
   if (info%flag < 0) go to 100
   matrix_type = -4 ! Complex, Hermitian indefinite
   call ma97_factor(matrix_type,val,akeep,fkeep,control,info)
   if (info%flag < 0) go to 100

   ! Read in the right-hand side
   read (*,*) x(1:n)

   ! Solve
   call ma97_solve(x,akeep,fkeep,control,info)
   if (info%flag < 0) go to 100
   write (*,"(/a,/,(2('(',2es18.10,')')))") ' The computed solution is:', x(1:n)

   ! Determine the pivot order used
   allocate (piv_order(1:n))
   call ma97_enquire_indef(akeep, fkeep, control, info, piv_order=piv_order)
   write (6,"(a,10i5)") 'piv_order', piv_order(1:n)

   100 continue
   call ma97_finalise(akeep, fkeep)

end program hsl_ma97zs
