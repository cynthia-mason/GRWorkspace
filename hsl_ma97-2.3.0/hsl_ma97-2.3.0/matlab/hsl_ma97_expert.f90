!
! This routine offers direct access to hsl_ma97 functionality
!     [handle, info] = ma97_expert('factor', A[, control, P])
!     [x, info] = ma97_expert('solve', handle, b[, control])
!     [x, info, consist] = ma97_expert('solve', handle, b[, control])
!     [x, info, handle] = ma97_expert('backslash', A, b[, control, P])
!     [y] = ma97_expert('multiply', handle, transpose, inverse, x)
!     ma97_expert('destroy', handle)
! P is an optional ordering.
!

! Structure of code:
! The module ma97_matlab_main provides a series of functions, one for each
!    possible 'action' in real or complex data.
! The module ma97_handles provides a series of functions for dealing with
!    integer handles associated on the Fortran side with keep and order
! The mexFunction unwraps handles and calls the corrected subroutine
!    for each 'action'.

module ma97_matlab_main
!$ use omp_lib
   use hsl_matlab
   use hsl_mc68_double
   use hsl_mc69_double, only: mc69_cscl_clean_real => mc69_cscl_clean, &
                              HSL_MATRIX_REAL_SYM_PSDEF, &
                              HSL_MATRIX_REAL_SYM_INDEF, &
                              HSL_MATRIX_CPLX_SYM, &
                              HSL_MATRIX_CPLX_HERM_PSDEF, &
                              HSL_MATRIX_CPLX_HERM_INDEF
   use hsl_mc69_double_complex, only: mc69_cscl_clean_complex => mc69_cscl_clean
   use hsl_ma97_double, only: ma97_rakeep => ma97_akeep,          &
                              ma97_rfkeep => ma97_fkeep,          &
                              ma97_rcontrol => ma97_control,      &
                              ma97_rinfo => ma97_info,            &
                              ma97_analyse,                       &
                              ma97_factor,                        &
                              ma97_factor_solve,                  &
                              ma97_solve,                         &
                              ma97_solve_fredholm,                &
                              ma97_sparse_fwd_solve,              &
                              ma97_lmultiply,                     &
                              ma97_finalise,                      &
                              ma97_get_n__
   use hsl_ma97_double_complex, only: ma97_cakeep => ma97_akeep,  &
                              ma97_cfkeep => ma97_fkeep,          &
                              ma97_ccontrol => ma97_control,      &
                              ma97_cinfo => ma97_info,            &
                              ma97_analyse,                       &
                              ma97_factor,                        &
                              ma97_factor_solve,                  &
                              ma97_solve,                         &
                              ma97_solve_fredholm,                &
                              ma97_sparse_fwd_solve,              &
                              ma97_lmultiply,                     &
                              ma97_finalise,                      &
                              ma97_get_n__
   implicit none

   integer, parameter :: wp = kind(0d0)
   integer, parameter :: long = selected_int_kind(18)

   interface ma97_matlab_finalise
      module procedure ma97_matlab_finalise_real
      module procedure ma97_matlab_finalise_complex
   end interface ma97_matlab_finalise
   interface copy_control_in
      module procedure copy_control_in_real
      module procedure copy_control_in_complex
   end interface copy_control_in

contains
   subroutine copy_control_in_real(pm, control, num_threads, matrix_type)
      integer(mwp_) :: pm
      type(ma97_rcontrol), intent(inout) :: control
      integer, intent(out) :: num_threads
      integer, intent(out) :: matrix_type

      logical :: pos_def
      integer(mwp_) :: pc
      integer(int4_) :: fnum
      character(80) :: fname
      character(200) :: warnmsg

      num_threads = 0

      matrix_type = HSL_MATRIX_REAL_SYM_INDEF
      do fnum = 1, MATLAB_get_no_fields(pm)
         fname = MATLAB_get_field_name_by_no(pm, fnum)
         select case(trim(fname))
         case("nemin")
            call MATLAB_get_value(pm, fname, pc, control%nemin)
         case("num_threads")
            call MATLAB_get_value(pm, fname, pc, num_threads)
         case("ordering")
            call MATLAB_get_value(pm, fname, pc, control%ordering)
         case("pos_def")
            call MATLAB_get_value(pm, fname, pc, pos_def)
            if(pos_def) matrix_type = HSL_MATRIX_REAL_SYM_PSDEF
         case("scaling")
            call MATLAB_get_value(pm, fname, pc, control%scaling)
         case("small")
            call MATLAB_get_value(pm, fname, pc, control%small)
         case("u")
            call MATLAB_get_value(pm, fname, pc, control%u)
         case default
            write(warnmsg, "(3a)") "Ignored unrecognised control parameter '", &
               trim(fname), "'"
            call MATLAB_warning(warnmsg)
         end select
      end do
   end subroutine copy_control_in_real

   subroutine copy_control_in_complex(pm, control, num_threads, matrix_type)
      integer(mwp_) :: pm
      type(ma97_ccontrol), intent(inout) :: control
      integer, intent(out) :: num_threads
      integer, intent(out) :: matrix_type

      integer(mwp_) :: pc
      integer(int4_) :: fnum
      character(80) :: fname
      character(200) :: warnmsg
      logical :: herm, pos_def

      num_threads = 0
      pos_def = .false.
      herm = .false.

      do fnum = 1, MATLAB_get_no_fields(pm)
         fname = MATLAB_get_field_name_by_no(pm, fnum)
         select case(trim(fname))
         case("hermitian")
            call MATLAB_get_value(pm, fname, pc, herm)
         case("nemin")
            call MATLAB_get_value(pm, fname, pc, control%nemin)
         case("num_threads")
            call MATLAB_get_value(pm, fname, pc, num_threads)
         case("ordering")
            call MATLAB_get_value(pm, fname, pc, control%ordering)
         case("pos_def")
            call MATLAB_get_value(pm, fname, pc, pos_def)
         case("scaling")
            call MATLAB_get_value(pm, fname, pc, control%scaling)
         case("small")
            call MATLAB_get_value(pm, fname, pc, control%small)
         case("u")
            call MATLAB_get_value(pm, fname, pc, control%u)
         case default
            write(warnmsg, "(3a)") "Ignored unrecognised control parameter '", &
               trim(fname), "'"
            call MATLAB_warning(warnmsg)
         end select
      end do

      if(pos_def .and. .not.herm) then
         write(warnmsg, "(a)") "Specified non-Hemitian and positive-definite. &
         &These options conflict."
         call MATLAB_error(warnmsg)
      endif
 
      matrix_type = HSL_MATRIX_CPLX_SYM
      if(herm) then
         if(pos_def) then
            matrix_type = HSL_MATRIX_CPLX_HERM_PSDEF
         else
            matrix_type = HSL_MATRIX_CPLX_HERM_INDEF
         endif
      endif

   end subroutine copy_control_in_complex

   subroutine unknownError(routine, flag)
      character(len=*) :: routine
      integer :: flag
      character(len=200) :: errmsg

      write(errmsg, "(3a,i3)") "Error return from ", routine, ". flag = ", flag
      call MATLAB_error(errmsg)
   end subroutine unknownError

   ! [handle, info] = ma97_expert('factor', A[, P])
   ! handle and 'factor' already dealt with
   subroutine ma97_matlab_analyse_factor_real(nlhs_in, plhs, nrhs_in, prhs, &
         akeep, fkeep, order)
      integer*4 :: nlhs_in, nrhs_in
      integer(mwp_) :: plhs(*), prhs(*)
      type(ma97_rakeep), intent(inout) :: akeep
      type(ma97_rfkeep), intent(inout) :: fkeep
      integer, dimension(:), pointer, intent(out) :: order

      integer, parameter :: A_in = 1, &
                            control_in = 2, &
                            P_in = 3
      integer, parameter :: info_out = 1

      integer(mws_) :: mwm, mwn, mwnrhs, nz, temp
      integer(mwp_) :: ml_ptr, ml_row, ml_val
      integer :: m, n, nrhs, i, matrix_type
      integer, dimension(:), allocatable :: ptr, row, invp
      real(wp), dimension(:), allocatable :: val
      real(wp), dimension(:,:), allocatable :: rhs
      integer :: flag
      integer :: st
      integer :: num_threads, orig_num_threads

      type(mc68_control) :: control68
      type(mc68_info) :: info68

      type(ma97_rcontrol) :: control
      type(ma97_rinfo) :: info

      character(len=200) :: errmsg

      real(wp) :: atime, ftime
      integer :: t_start, t_stop, t_rate
      character(len=10) :: orderused

      ! Check number of arguments
      if(nrhs_in.lt.1 .or. nrhs_in .gt.3) &
         call MATLAB_error("Wrong number of input arguments")
      if(nlhs_in.ne.0 .and. nlhs_in.ne.1) &
         call MATLAB_error("Wrong number of output arguments")

      ! Copy matrix in
      call matlab_to_fortran(prhs(A_in), mwm, mwn, ptr, row, val, 'A')
      if(mwm.ne.mwn) call MATLAB_error("The matrix must be square")
      n = mwn

      ! Read control in
      num_threads = 0
      matrix_type = HSL_MATRIX_REAL_SYM_INDEF ! default if control not present
      if(nrhs_in.ge.control_in) then
         ! User supplied control structure
         if(MATLAB_is_structure(prhs(control_in))) then
            call copy_control_in(prhs(control_in), control, num_threads, &
               matrix_type)
         endif
      endif

      ! Convert to HSL standard form
      call mc69_cscl_clean_real(matrix_type, n, n, ptr, row, &
         flag, val=val)
      select case(flag)
      case(0,1,4,5)
         ! Sucess (expect oor entries as we are throwing away upr triangle!)
      case(-12)
         ! Hermitian but non-zero imaginary component on diagonal
         call MATLAB_context_error("hsl_ma97:ExpectsHermitian", &
            "HSL_MA97 Error: Hermitian factorization requested, but matrix has &
            &non-zero imaginary component on diagonal.")
      case default
         call unknownError("mc69_cscl_clean", flag)
      end select

      ! Set number of threads
      orig_num_threads = 0
      if(num_threads.ne.0) then
         ! User specified a specific number of threads
         if(num_threads.le.0) &
            call MATLAB_error("control.num_threads must be positive")
!$       orig_num_threads = omp_get_max_threads()
         if(orig_num_threads.eq.0 .and. num_threads.gt.1) then
            call MATLAB_warning( &
               "Parallel run requested, but compiled without OpenMP. Running in&
               & serial")
         endif
!$       if(num_threads .gt. orig_num_threads) &
!$          call MATLAB_warning( &
!$             "control.num_threads exceeds omp_get_max_threads(). Proceeding &
!$             &with control.num_threads anyway.")
!$       call omp_set_num_threads(num_threads)
      endif

      ! Obtain ordering
      allocate(order(n), stat=st)
      if(st.ne.0) goto 100
      if(nrhs_in.ge.P_in) then
         ! User supplied (note: input as double, inverse of order)
         call matlab_to_fortran(prhs(P_in), invp, temp, 'P')
         if(mwn .ne. temp) then
            write(errmsg, "(4(a,i5))") &
               "Dimensions of A and P inconsistent: A=", &
               n, "x", n, ", size(P)=", temp
            call MATLAB_error(errmsg)
         endif
         order(:) = (n+1) ! initialise to invalid value
         do i = 1, n
            order(invp(i)) = i
         end do
         control%ordering = 0 ! user defined
      endif

      ! Call analyse
      call system_clock(t_start)
      call ma97_analyse(.false., n, ptr, row, akeep, control, info, &
         order=order, val=val)
      select case (info%flag)
      case(0)
         ! Success. Do nothing
      case(6)
         ! Structually singular. Ignore as will be picked up in factorize.
      case(-11)
         ! Problem with order
         if(allocated(invp)) then
            write(errmsg, "(a, 10i4)") "Problem with P. P = ", &
               invp(1:min(n,10))
         else
            write(errmsg, "(a, 10i4)") "Bad value of control.ordering = ", &
               control%ordering
         endif
         call MATLAB_error(errmsg)
      case default
!$       if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)
         call unknownError("ma97_analyse", info%flag)
      end select
      select case(info%ordering)
      case(0)
         orderused = 'user'
      case(1)
         orderused = 'AMD'
      case(2)
         orderused = 'MD'
      case(3)
         orderused = 'MeTiS'
      case(4)
         orderused = 'MC47'
      case(7)
         orderused = 'MC80:AMD'
      case(8)
         orderused = 'MC80:METIS'
      end select
      call system_clock(t_stop, t_rate)
      atime = real(t_stop-t_start)/t_rate

      ! Use free MC64 scaling if available and desired
      if(control%scaling.eq.1 .and. &
         (info%ordering.eq.7 .or. info%ordering.eq.8)) control%scaling = 3

      ! Call factor
      call system_clock(t_start)
      call ma97_factor(matrix_type, val, akeep, fkeep, control, info, &
         ptr=ptr, row=row)
      call system_clock(t_stop, t_rate)
      select case (info%flag)
      case(0:1)
         ! Success. Do nothing.
      case(-8)
         ! Not positive definite
         call MATLAB_context_error("hsl_ma97:ExpectsPosdef", &
            "HSL_MA97 Error: control.pos_def=true, but matrix is not positive-&
            &definite.")
      case(-9)
         ! IEEE infinities in factorization
         call MATLAB_context_error("hsl_ma97:Infs", &
            "HSL_MA97 Error: IEEE infinities encountered during factorization. &
            &Try larger control.u or control.small.")
      case(-11)
         ! Bad value of ordering.
         call MATLAB_context_error("hsl_ma97:BadOrdering", &
            "HSL_MA97 Error: control.ordering has invalid value")
      case(6,7)
         ! Singular matrix
         call MATLAB_warning("HSL_MA97: Matrix was found to be singular.")
      case(8)
         ! Used Matching-based ordering but not scaling
         call MATLAB_warning("HSL_MA97: Used matching-based ordering but &
            &ignored free scaling.")
      case default
!$       if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)
         call unknownError("ma97_factor", info%flag)
      end select
      ftime = real(t_stop-t_start)/t_rate

      ! Setup info structure if required
      if(nlhs_in.ge.info_out) then
         call copy_info_out(plhs(info_out), info, atime, ftime, &
            orderused)
      endif

!$    if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)

      return
      100 continue
      call MATLAB_error("Insufficient memory")

   contains
      subroutine copy_info_out(ml_info, info, atime, ftime, order)
         integer(mwp_), intent(inout) :: ml_info
         type(ma97_rinfo), intent(in) :: info
         real(wp), intent(in) :: atime, ftime
         character(len=*) :: order

         integer(mwp_) :: comp
         integer*4 :: szf

         character(len=50), dimension(9) :: fields = (/ &
            "matrix_rank   ", &
            "num_delay     ", &
            "num_factor    ", &
            "num_flops     ", &
            "num_neg       ", &
            "num_two       ", &
            "order         ", &
            "analyse_time  ", &
            "factor_time   " /)

         ml_info = MATLAB_create_structure(fields)

         call MATLAB_set_field(ml_info, 'matrix_rank', info%matrix_rank)
         call MATLAB_set_field(ml_info, 'num_delay', info%num_delay)
         call MATLAB_set_field(ml_info, 'num_factor', info%num_factor)
         call MATLAB_set_field(ml_info, 'num_flops', info%num_flops)
         call MATLAB_set_field(ml_info, 'num_neg', info%num_neg)
         call MATLAB_set_field(ml_info, 'num_two', info%num_two)
         call MATLAB_set_field(ml_info, 'order', order)
         call MATLAB_set_field(ml_info, 'analyse_time', atime)
         call MATLAB_set_field(ml_info, 'factor_time', ftime)
      end subroutine copy_info_out

   end subroutine ma97_matlab_analyse_factor_real

   ! [handle, info] = ma97_expert('factor', A[, P])
   ! handle and 'factor' already dealt with
   subroutine ma97_matlab_analyse_factor_complex(nlhs_in, plhs, nrhs_in, prhs, &
         akeep, fkeep, order)
      integer*4 :: nlhs_in, nrhs_in
      integer(mwp_) :: plhs(*), prhs(*)
      type(ma97_cakeep), intent(inout) :: akeep
      type(ma97_cfkeep), intent(inout) :: fkeep
      integer, dimension(:), pointer, intent(out) :: order

      integer, parameter :: A_in = 1, &
                            control_in = 2, &
                            P_in = 3
      integer, parameter :: info_out = 1

      integer(mws_) :: mwm, mwn, mwnrhs, nz
      integer(mwp_) :: ml_ptr, ml_row, ml_val
      integer :: m, n, nrhs, i
      integer, dimension(:), allocatable :: ptr, row, invp
      complex(wp), dimension(:), allocatable :: val
      integer :: flag
      integer :: st
      integer(mws_) :: temp
      integer :: num_threads, orig_num_threads
      character(len=10) :: orderused

      type(mc68_control) :: control68
      type(mc68_info) :: info68

      type(ma97_ccontrol) :: control
      type(ma97_cinfo) :: info

      character(len=200) :: errmsg

      real(wp) :: atime, ftime
      integer :: t_start, t_stop, t_rate

      integer :: matrix_type

      ! Check number of arguments
      if(nrhs_in.lt.1 .or. nrhs_in .gt.3) &
         call MATLAB_error("Wrong number of input arguments")
      if(nlhs_in.ne.0 .and. nlhs_in.ne.1) &
         call MATLAB_error("Wrong number of output arguments")
         
      ! Copy matrix in
      call matlab_to_fortran(prhs(A_in), mwm, mwn, ptr, row, val, 'A')
      if(mwm.ne.mwn) call MATLAB_error("The matrix must be square")
      n = mwn

      ! Copy control in
      num_threads = 0
      matrix_type = HSL_MATRIX_CPLX_SYM
      if(nrhs_in.ge.control_in) then
         ! User supplied control structure
         if(MATLAB_is_structure(prhs(control_in))) then
            call copy_control_in(prhs(control_in), control, num_threads, &
               matrix_type)
         endif
      endif

      ! Convert to HSL standard form
      call mc69_cscl_clean_complex(matrix_type, n, n, ptr, row, flag, val=val)
      select case(flag)
      case(0,1,4,5)
         ! Sucess (expect oor entries as we are throwing away upr triangle!)
      case(-12)
         ! Hermitian but non-zero imaginary component on diagonal
         call MATLAB_context_error("hsl_ma97:ExpectsHermitian", &
            "HSL_MA97 Error: Hermitian factorization requested, but matrix has &
            &non-zero imaginary component on diagonal.")
      case default
         call unknownError("mc69_cscl_clean", flag)
      end select

      ! Set number of threads
      orig_num_threads = 0
      if(num_threads.ne.0) then
         ! User specified a specific number of threads
         if(num_threads.le.0) &
            call MATLAB_error("control.num_threads must be positive")
!$       orig_num_threads = omp_get_max_threads()
         if(orig_num_threads.eq.0 .and. num_threads.gt.1) then
            call MATLAB_warning( &
               "Parallel run requested, but compiled without OpenMP. Running in&
               & serial")
         endif
!$       if(num_threads .gt. orig_num_threads) &
!$          call MATLAB_warning( &
!$             "control.num_threads exceeds omp_get_max_threads(). Proceeding &
!$             &with control.num_threads anyway.")
!$       call omp_set_num_threads(num_threads)
      endif

      ! Obtain ordering
      allocate(order(n), stat=st)
      if(st.ne.0) goto 100
      if(nrhs_in.ge.P_in) then
         ! User supplied (note: input as double, inverse of order)
         call matlab_to_fortran(prhs(P_in), invp, temp, 'P')
         if(mwn .ne. temp) then
            write(errmsg, "(4(a,i5))") &
               "Dimensions of A and P inconsistent: A=", &
               n, "x", n, ", size(P)=", temp
            call MATLAB_error(errmsg)
         endif
         order(:) = (n+1) ! initialise to invalid value
         do i = 1, n
            order(invp(i)) = i
         end do
         control%ordering = 0 ! user defined
      endif

      ! Call analyse
      call system_clock(t_start)
      call ma97_analyse(.false., n, ptr, row, akeep, control, info, &
         order=order)
      call system_clock(t_stop, t_rate)
      select case (info%flag)
      case(0)
         ! Success. Do nothing
      case(6)
         ! Structually singular. Ignore as will be picked up in factorize.
      case(-11)
         ! Problem with order
         if(allocated(invp)) then
            write(errmsg, "(a, 10i4)") "Problem with P. P = ", &
               invp(1:min(n,10))
         else
            write(errmsg, "(a, 10i4)") "Bad value of control.ordering = ", &
               control%ordering
         endif
         call MATLAB_error(errmsg)
      case default
!$       if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)
         call unknownError("ma97_analyse", info%flag)
      end select
      select case(info%ordering)
      case(0)
         orderused = 'user'
      case(1)
         orderused = 'AMD'
      case(2)
         orderused = 'MD'
      case(3)
         orderused = 'MeTiS'
      case(4)
         orderused = 'MC47'
      case(7)
         orderused = 'MC80:AMD'
      case(8)
         orderused = 'MC80:METIS'
      end select
      atime = real(t_stop-t_start)/t_rate

      ! Use free MC64 scaling if available and desired
      if(control%scaling.eq.1 .and. &
         (info%ordering.eq.7 .or. info%ordering.eq.8)) control%scaling = 3

      ! Call factor
      call system_clock(t_start)
      call ma97_factor(matrix_type, val, akeep, fkeep, control, info, &
         ptr=ptr, row=row)
      call system_clock(t_stop, t_rate)
      select case (info%flag)
      case(0:1)
         ! Success. Do nothing.
      case(-8)
         ! Not positive definite
         call MATLAB_context_error("hsl_ma97:ExpectsPosdef", &
            "HSL_MA97 Error: control.pos_def=true, but matrix is not positive-&
            &definite.")
      case(-9)
         ! IEEE infinities in factorization
         call MATLAB_context_error("hsl_ma97:Infs", &
            "HSL_MA97 Error: IEEE infinities encountered during factorization. &
            &Try larger control.u or control.small.")
      case(-11)
         ! Bad value of ordering.
         call MATLAB_context_error("hsl_ma97:BadOrdering", &
            "HSL_MA97 Error: control.ordering has invalid value")
      case(6,7)
         ! Singular matrix
         call MATLAB_warning("HSL_MA97: Matrix was found to be singular.")
      case(8)
         ! Used Matching-based ordering but not scaling
         call MATLAB_warning("HSL_MA97: Used matching-based ordering but &
            &ignored free scaling.")
      case default
!$       if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)
         call unknownError("ma97_factor", info%flag)
      end select
      ftime = real(t_stop-t_start)/t_rate

      ! Setup info structure if required
      if(nlhs_in.ge.info_out) then
         call copy_info_out(plhs(info_out), info, atime, ftime, &
            orderused)
      endif

!$    if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)

      return
      100 continue
      call MATLAB_error("Insufficient memory")

   contains
      subroutine copy_info_out(ml_info, info, atime, ftime, order)
         integer(mwp_), intent(inout) :: ml_info
         type(ma97_cinfo), intent(in) :: info
         real(wp), intent(in) :: atime, ftime
         character(len=*) :: order

         integer(mwp_) :: comp
         integer*4 :: szf

         character(len=50), dimension(9) :: fields = (/ &
            "matrix_rank   ", &
            "num_delay     ", &
            "num_factor    ", &
            "num_flops     ", &
            "num_neg       ", &
            "num_two       ", &
            "order         ", &
            "analyse_time  ", &
            "factor_time   " /)

         ml_info = MATLAB_create_structure(fields)

         call MATLAB_set_field(ml_info, 'matrix_rank', info%matrix_rank)
         call MATLAB_set_field(ml_info, 'num_delay', info%num_delay)
         call MATLAB_set_field(ml_info, 'num_factor', info%num_factor)
         call MATLAB_set_field(ml_info, 'num_flops', info%num_flops)
         call MATLAB_set_field(ml_info, 'num_neg', info%num_neg)
         call MATLAB_set_field(ml_info, 'num_two', info%num_two)
         call MATLAB_set_field(ml_info, 'order', order)
         call MATLAB_set_field(ml_info, 'analyse_time', atime)
         call MATLAB_set_field(ml_info, 'factor_time', ftime)
      end subroutine copy_info_out

   end subroutine ma97_matlab_analyse_factor_complex

   ! [x, info] = ma97_expert('solve', handle, b)
   ! 'solve' and handle already dealt with
   subroutine ma97_matlab_solve_real(nlhs_in, plhs, nrhs_in, prhs, akeep, &
         fkeep)
      integer*4 :: nlhs_in, nrhs_in
      integer(mwp_) :: plhs(*), prhs(*)
      type(ma97_rakeep), intent(inout) :: akeep
      type(ma97_rfkeep), intent(inout) :: fkeep

      integer, parameter :: b_in = 1, &
                            control_in = 2
      integer, parameter :: x_out = 1, &
                            info_out = 2, &
                            consist_out = 3

      integer(mws_) :: mwn, mwnrhs
      integer(mwp_) :: ml_ptr, ml_row, ml_val
      real(wp), dimension(:,:), allocatable :: rhs, rhs2
      integer :: flag
      integer :: st
      integer :: temp
      integer :: matrix_type

      type(ma97_rcontrol) :: control
      type(ma97_rinfo) :: info

      character(len=200) :: errmsg
      integer :: n, nrhs
      integer :: num_threads, orig_num_threads
      logical, dimension(:), allocatable :: consist

      real(wp) :: stime
      integer :: t_start, t_stop, t_rate

      ! Check number of arguments
      if(nrhs_in.lt.1) call MATLAB_error("Insufficient input arguments")
      if(nrhs_in.gt.2) call MATLAB_error("Too many input arguments")
      if(nlhs_in.lt.1) call MATLAB_error("Insufficient output arguments")
      if(nlhs_in.gt.3) call MATLAB_error("Too many output arguments")

      n = ma97_get_n__(akeep)

      ! Get rhs
      call matlab_to_fortran(prhs(b_in), rhs, mwn, mwnrhs, 'b')
      if(mwn .ne. n) then
         write(errmsg, "(4(a,i5))") "Dimensions of A and b inconsistent: A=", &
            n, "x", n, ", b=", mwn, "x", mwnrhs
         call MATLAB_error(errmsg)
      endif
      nrhs = mwnrhs

      ! Copy control in
      num_threads = 0
      matrix_type = HSL_MATRIX_REAL_SYM_INDEF ! default if control not present
      if(nrhs_in.ge.control_in) then
         ! User supplied control structure
         if(MATLAB_is_structure(prhs(control_in))) then
            call copy_control_in(prhs(control_in), control, num_threads, &
               matrix_type)
         endif
      endif

      ! Set number of threads
      orig_num_threads = 0
      if(num_threads.ne.0) then
         ! User specified a specific number of threads
         if(num_threads.le.0) &
            call MATLAB_error("control.num_threads must be positive")
!$          orig_num_threads = omp_get_max_threads()
         if(orig_num_threads.eq.0 .and. num_threads.gt.1) then
            call MATLAB_warning( &
               "Parallel run requested, but compiled without OpenMP. Running in&
               & serial")
         endif
!$       if(num_threads .gt. orig_num_threads) &
!$          call MATLAB_warning( &
!$             "control.num_threads exceeds omp_get_max_threads(). Proceeding &
!$             &with control.num_threads anyway.")
!$       call omp_set_num_threads(num_threads)
      endif

      ! Setup flag array if required
      if(nlhs_in.eq.consist_out) allocate(consist(nrhs))
      
      ! Call solve
      call system_clock(t_start)
      if(nlhs_in.ge.consist_out) then
         allocate(rhs2(n,2*nrhs))
         rhs2(:,1:nrhs) = rhs(:,:)
         deallocate(rhs)
         call ma97_solve_fredholm(nrhs, consist, rhs2, n, &
            akeep, fkeep, control, info)
      else
         call ma97_solve(nrhs, rhs, n, akeep, fkeep, control, info)
      endif
      call system_clock(t_stop, t_rate)
      select case (info%flag)
      case(0)
         ! Success. Do nothing.
      case(7)
         ! Original matrix was singular
         if(nlhs_in.lt.consist_out) &
            call MATLAB_warning("HSL_MA97: Matrix was singular and consistency &
               &of rhs has not been checked.")
      case default
!$       if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)
         call unknownError("ma97_solve", info%flag)
      end select
      stime = real(t_stop-t_start)/t_rate

      ! Setup info structure if required
      if(nlhs_in.ge.info_out) then
         call copy_info_out(plhs(info_out), info, stime)
      endif

      ! Copy solution out
      if(nlhs_in.ge.consist_out) then
         plhs(x_out) = fortran_to_matlab(rhs2)
      else
         plhs(x_out) = fortran_to_matlab(rhs)
      endif

      ! Copy consist out
      if(nlhs_in.ge.consist_out) then
         plhs(consist_out) = fortran_to_matlab(consist)
      endif

!$    if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)

      return
      100 continue
      call MATLAB_error("Insufficient memory")

   contains
      subroutine copy_info_out(ml_info, info, stime)
         integer(mwp_), intent(inout) :: ml_info
         type(ma97_rinfo), intent(in) :: info
         real(wp), intent(in) :: stime

         integer(mwp_) :: comp
         integer*4 :: szf

         character(len=50), dimension(1) :: fields = (/ &
            "solve_time    " /)

         ml_info = MATLAB_create_structure(fields)

         call MATLAB_set_field(ml_info, 'solve_time', stime)
      end subroutine copy_info_out

   end subroutine ma97_matlab_solve_real

   ! [x, info] = ma97_expert('solve', handle, b)
   ! 'solve' and handle already dealt with
   subroutine ma97_matlab_solve_complex(nlhs_in, plhs, nrhs_in, prhs, akeep, &
         fkeep)
      integer*4 :: nlhs_in, nrhs_in
      integer(mwp_) :: plhs(*), prhs(*)
      type(ma97_cakeep), intent(inout) :: akeep
      type(ma97_cfkeep), intent(inout) :: fkeep

      integer, parameter :: b_in = 1, &
                            control_in = 2
      integer, parameter :: x_out = 1, &
                            info_out = 2, &
                            consist_out = 3

      integer(mws_) :: mwn, mwnrhs
      integer(mwp_) :: ml_ptr, ml_row, ml_val
      complex(wp), dimension(:,:), allocatable :: rhs, rhs2
      integer :: flag
      integer :: st
      integer :: temp

      type(ma97_ccontrol) :: control
      type(ma97_cinfo) :: info

      character(len=200) :: errmsg
      integer :: n, nrhs
      integer :: num_threads, orig_num_threads
      logical, dimension(:), allocatable :: consist

      real(wp) :: stime
      integer :: t_start, t_stop, t_rate
      integer :: matrix_type

      ! Check number of arguments
      if(nrhs_in.lt.1) call MATLAB_error("Insufficient input arguments")
      if(nrhs_in.gt.2) call MATLAB_error("Too many input arguments")
      if(nlhs_in.lt.1) call MATLAB_error("Insufficient output arguments")
      if(nlhs_in.gt.3) call MATLAB_error("Too many output arguments")

      n = ma97_get_n__(akeep)

      ! Obtain rhs vector
      call matlab_to_fortran(prhs(b_in), rhs, mwn, mwnrhs, 'b')
      if(mwn .ne. n) then
         write(errmsg, "(4(a,i5))") "Dimensions of A and b inconsistent: A=", &
           n, "x", n, ", b=", mwn, "x", mwnrhs
         call MATLAB_error(errmsg)
      endif
      nrhs = mwnrhs

      ! Copy control in
      num_threads = 0
      matrix_type = HSL_MATRIX_CPLX_SYM
      if(nrhs_in.ge.control_in) then
         ! User supplied control structure
         if(MATLAB_is_structure(prhs(control_in))) then
            call copy_control_in(prhs(control_in), control, num_threads, &
               matrix_type)
         endif
      endif

      ! Set number of threads
      orig_num_threads = 0
      if(num_threads.ne.0) then
         ! User specified a specific number of threads
         if(num_threads.le.0) &
            call MATLAB_error("control.num_threads must be positive")
!$          orig_num_threads = omp_get_max_threads()
         if(orig_num_threads.eq.0 .and. num_threads.gt.1) then
            call MATLAB_warning( &
               "Parallel run requested, but compiled without OpenMP. Running in&
               & serial")
         endif
!$       if(num_threads .gt. orig_num_threads) &
!$          call MATLAB_warning( &
!$             "control.num_threads exceeds omp_get_max_threads(). Proceeding &
!$             &with control.num_threads anyway.")
!$       call omp_set_num_threads(num_threads)
      endif

      ! Setup flag array if required
      if(nlhs_in.ge.consist_out) allocate(consist(nrhs))

      ! Call solve
      call system_clock(t_start)
      if(nlhs_in.ge.consist_out) then
         allocate(rhs2(n,2*nrhs))
         rhs2(:,1:nrhs) = rhs(:,:)
         deallocate(rhs)
         call ma97_solve_fredholm(nrhs, consist, rhs2, n, &
            akeep, fkeep, control, info)
      else
         call ma97_solve(nrhs, rhs, n, akeep, fkeep, control, info)
      endif
      call system_clock(t_stop, t_rate)
      select case (info%flag)
      case(0)
         ! Success. Do nothing.
      case(7)
         ! Original matrix was singular
         if(nlhs_in.lt.consist_out) &
            call MATLAB_warning("HSL_MA97: Matrix was singular and consistency &
               &of rhs has not been checked.")
      case default
!$       if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)
         call unknownError("ma97_solve", info%flag)
      end select
      stime = real(t_stop-t_start)/t_rate

      ! Setup info structure if required
      if(nlhs_in.ge.info_out) then
         call copy_info_out(plhs(info_out), info, stime)
      endif

      ! Copy solution out
      if(nlhs_in.ge.consist_out) then
         plhs(x_out) = fortran_to_matlab(rhs2)
      else
         plhs(x_out) = fortran_to_matlab(rhs)
      endif

      ! Copy consist out
      if(nlhs_in.ge.consist_out) then
         plhs(consist_out) = fortran_to_matlab(consist)
      endif

      ! Restore number of threads
!$    if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)

      return
      100 continue
      call MATLAB_error("Insufficient memory")
   contains
      subroutine copy_info_out(ml_info, info, stime)
         integer(mwp_), intent(inout) :: ml_info
         type(ma97_cinfo), intent(in) :: info
         real(wp), intent(in) :: stime

         integer(mwp_) :: comp
         integer*4 :: szf

         character(len=50), dimension(1) :: fields = (/ &
            "solve_time    " /)

         ml_info = MATLAB_create_structure(fields)

         call MATLAB_set_field(ml_info, 'solve_time', stime)
      end subroutine copy_info_out


   end subroutine ma97_matlab_solve_complex

   ! [x, info, handle] = ma97_expert('backslash', A, b[, control, P])
   ! handle and 'backslash' already dealt with
   subroutine ma97_matlab_backslash_real(nlhs_in, plhs, nrhs_in, prhs, &
         akeep, fkeep, order)
      integer*4 :: nlhs_in, nrhs_in
      integer(mwp_) :: plhs(*), prhs(*)
      type(ma97_rakeep), intent(inout) :: akeep
      type(ma97_rfkeep), intent(inout) :: fkeep
      integer, dimension(:), pointer, intent(out) :: order

      integer, parameter :: A_in = 1, &
                            b_in = 2, &
                            control_in = 3, &
                            P_in = 4
      integer, parameter :: x_out = 1, &
                            info_out = 2

      integer(mws_) :: mwm, mwn, mwnrhs, nz, temp
      integer(mwp_) :: ml_ptr, ml_row, ml_val
      integer :: m, n, nrhs, i, matrix_type
      integer, dimension(:), allocatable :: ptr, row, invp
      real(wp), dimension(:), allocatable :: val
      real(wp), dimension(:,:), allocatable :: rhs
      integer :: flag
      integer :: st
      integer :: num_threads, orig_num_threads

      type(mc68_control) :: control68
      type(mc68_info) :: info68

      type(ma97_rcontrol) :: control
      type(ma97_rinfo) :: info

      character(len=200) :: errmsg

      real(wp) :: atime, ftime
      integer :: t_start, t_stop, t_rate
      character(len=10) :: orderused

      ! Check number of arguments
      if(nrhs_in.lt.2) call MATLAB_error("Insufficient input arguments")
      if(nrhs_in.gt.4) call MATLAB_error("Too many input arguments")
      if(nlhs_in.lt.1) call MATLAB_error("Insufficient output arguments")
      if(nlhs_in.gt.3) call MATLAB_error("Too many output arguments")

      ! Get sparse matrix
      call matlab_to_fortran(prhs(A_in), mwm, mwn, ptr, row, val, 'A')
      if(mwm.ne.mwn) call MATLAB_error("The matrix must be square")
      n = mwn

      ! Read control in
      num_threads = 0
      matrix_type = HSL_MATRIX_REAL_SYM_INDEF ! default if control not present
      if(nrhs_in.ge.control_in) then
         ! User supplied control structure
         if(MATLAB_is_structure(prhs(control_in))) then
            call copy_control_in(prhs(control_in), control, num_threads, &
               matrix_type)
         endif
      endif

      ! Convert to HSL standard form
      call mc69_cscl_clean_real(matrix_type, n, n, ptr, row, flag, val=val)
      select case(flag)
      case(0,1,4,5)
         ! Sucess (expect oor entries as we are throwing away upr triangle!)
      case(-12)
         ! Hermitian but non-zero imaginary component on diagonal
         call MATLAB_context_error("hsl_ma97:ExpectsHermitian", &
            "HSL_MA97 Error: Hermitian factorization requested, but matrix has &
            &non-zero imaginary component on diagonal.")
      case default
         call unknownError("mc69_cscl_clean", flag)
      end select

      ! Obtain rhs vector
      call matlab_to_fortran(prhs(b_in), rhs, temp, mwnrhs, 'b')
      if(mwn .ne. temp) then
         write(errmsg, "(4(a,i5))") "Dimensions of A and b inconsistent: A=", &
            mwn, "x", mwn, ", b=", temp, "x", mwnrhs
         call MATLAB_error(errmsg)
      endif
      nrhs = mwnrhs

      ! Set number of threads
      orig_num_threads = 0
      if(num_threads.ne.0) then
         ! User specified a specific number of threads
         if(num_threads.le.0) &
            call MATLAB_error("control.num_threads must be positive")
!$       orig_num_threads = omp_get_max_threads()
         if(orig_num_threads.eq.0 .and. num_threads.gt.1) then
            call MATLAB_warning( &
               "Parallel run requested, but compiled without OpenMP. Running in&
               & serial")
         endif
!$       if(num_threads .gt. orig_num_threads) &
!$          call MATLAB_warning( &
!$             "control.num_threads exceeds omp_get_max_threads(). Proceeding &
!$             &with control.num_threads anyway.")
!$       call omp_set_num_threads(num_threads)
      endif

      ! Obtain ordering
      allocate(order(n), stat=st)
      if(st.ne.0) goto 100
      if(nrhs_in.ge.P_in) then
         ! User supplied (note: input as double, inverse of order)
         call matlab_to_fortran(prhs(P_in), invp, temp, 'P')
         if(mwn .ne. temp) then
            write(errmsg, "(4(a,i5))") &
               "Dimensions of A and P inconsistent: A=", &
               n, "x", n, ", size(P)=", temp
            call MATLAB_error(errmsg)
         endif
         order(:) = (n+1) ! initialise to invalid value
         do i = 1, n
            order(invp(i)) = i
         end do
         control%ordering = 0 ! user defined
      endif

      ! Call analyse
      call system_clock(t_start)
      call ma97_analyse(.false., n, ptr, row, akeep, control, info, &
         order=order)
      call system_clock(t_stop, t_rate)
      select case (info%flag)
      case(0)
         ! Success. Do nothing
      case(6)
         ! Structually singular, do nothing as will be picked up in factor_solve
      case(-11)
         ! Problem with order
         if(allocated(invp)) then
            write(errmsg, "(a, 10i4)") "Problem with P. P = ", &
               invp(1:min(n,10))
         else
            write(errmsg, "(a, 10i4)") "Bad value of control.ordering = ", &
               control%ordering
         endif
         call MATLAB_error(errmsg)
      case default
!$       if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)
         call unknownError("ma97_analyse", info%flag)
      end select
      select case(info%ordering)
      case(0)
         orderused = 'user'
      case(1)
         orderused = 'AMD'
      case(2)
         orderused = 'MD'
      case(3)
         orderused = 'MeTiS'
      case(4)
         orderused = 'MC47'
      case(7)
         orderused = 'MC80:AMD'
      case(8)
         orderused = 'MC80:METIS'
      end select
      atime = real(t_stop-t_start)/t_rate

      ! Use free MC64 scaling if available and desired
      if(control%scaling.eq.1 .and. &
         (info%ordering.eq.7 .or. info%ordering.eq.8)) control%scaling = 3

      ! Call factor_solve
      call system_clock(t_start)
      call ma97_factor_solve(matrix_type, val, nrhs, rhs, n, akeep, fkeep, &
         control, info, ptr=ptr, row=row)
      call system_clock(t_stop, t_rate)
      select case (info%flag)
      case(0:1)
         ! Success. Do nothing.
      case(-8)
         ! Not positive definite
         call MATLAB_context_error("hsl_ma97:ExpectsPosdef", &
            "HSL_MA97 Error: control.pos_def=true, but matrix is not positive-&
            &definite.")
      case(-9)
         ! IEEE infinities in factorization
         call MATLAB_context_error("hsl_ma97:Infs", &
            "HSL_MA97 Error: IEEE infinities encountered during factorization. &
            &Try larger control.u or control.small.")
      case(-11)
         ! Bad value of ordering.
         call MATLAB_context_error("hsl_ma97:BadOrdering", &
            "HSL_MA97 Error: control.ordering has invalid value")
      case(6,7)
         ! Singular matrix
         call MATLAB_warning("HSL_MA97: Matrix was found to be singular. &
            &RHS not checked for consistency.")
      case(8)
         ! Used Matching-based ordering but not scaling
         call MATLAB_warning("HSL_MA97: Used matching-based ordering but &
            &ignored free scaling.")
      case default
!$       if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)
         call unknownError("ma97_factor_solve", info%flag)
      end select
      ftime = real(t_stop-t_start)/t_rate

      ! Setup info structure if required
      if(nlhs_in.ge.info_out) then
         call copy_info_out(plhs(info_out), info, atime, ftime, &
            orderused)
      endif

      ! Copy solution out
      plhs(x_out) = fortran_to_matlab(rhs)

      ! Restore number of threads
!$    if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)

      return
      100 continue
      call MATLAB_error("Insufficient memory")

   contains
      subroutine copy_info_out(ml_info, info, atime, ftime, order)
         integer(mwp_), intent(inout) :: ml_info
         type(ma97_rinfo), intent(in) :: info
         real(wp), intent(in) :: atime, ftime
         character(len=*) :: order

         integer(mwp_) :: comp
         integer*4 :: szf

         character(len=50), dimension(9) :: fields = (/ &
            "matrix_rank      ", &
            "num_delay        ", &
            "num_factor       ", &
            "num_flops        ", &
            "num_neg          ", &
            "num_two          ", &
            "order            ", &
            "analyse_time     ", &
            "factor_solve_time" /)

         ml_info = MATLAB_create_structure(fields)

         call MATLAB_set_field(ml_info, 'matrix_rank', info%matrix_rank)
         call MATLAB_set_field(ml_info, 'num_delay', info%num_delay)
         call MATLAB_set_field(ml_info, 'num_factor', info%num_factor)
         call MATLAB_set_field(ml_info, 'num_flops', info%num_flops)
         call MATLAB_set_field(ml_info, 'num_neg', info%num_neg)
         call MATLAB_set_field(ml_info, 'num_two', info%num_two)
         call MATLAB_set_field(ml_info, 'order', order)
         call MATLAB_set_field(ml_info, 'analyse_time', atime)
         call MATLAB_set_field(ml_info, 'factor_solve_time', ftime)
      end subroutine copy_info_out

   end subroutine ma97_matlab_backslash_real

   ! [x, info, handle] = ma97_expert('backslash', A, b[, control, P])
   ! handle and 'backslash' already dealt with
   subroutine ma97_matlab_backslash_complex(nlhs_in, plhs, nrhs_in, prhs, &
         akeep, fkeep, order)
      integer*4 :: nlhs_in, nrhs_in
      integer(mwp_) :: plhs(*), prhs(*)
      type(ma97_cakeep), intent(inout) :: akeep
      type(ma97_cfkeep), intent(inout) :: fkeep
      integer, dimension(:), pointer, intent(out) :: order

      integer, parameter :: A_in = 1, &
                            b_in = 2, &
                            control_in = 3, &
                            P_in = 4
      integer, parameter :: x_out = 1, &
                            info_out = 2

      integer(mws_) :: mwm, mwn, mwnrhs, nz, temp
      integer(mwp_) :: ml_ptr, ml_row, ml_val
      integer :: m, n, nrhs, i
      integer, dimension(:), allocatable :: ptr, row, invp
      complex(wp), dimension(:), allocatable :: val
      complex(wp), dimension(:,:), allocatable :: rhs
      integer :: flag
      integer :: st
      integer :: num_threads, orig_num_threads

      type(mc68_control) :: control68
      type(mc68_info) :: info68

      type(ma97_ccontrol) :: control
      type(ma97_cinfo) :: info

      character(len=200) :: errmsg

      real(wp) :: atime, ftime
      integer :: t_start, t_stop, t_rate
      character(len=10) :: orderused
      integer :: matrix_type

      ! Check number of arguments
      if(nrhs_in.lt.2) call MATLAB_error("Insufficient input arguments")
      if(nrhs_in.gt.4) call MATLAB_error("Too many input arguments")
      if(nlhs_in.lt.1) call MATLAB_error("Insufficient output arguments")
      if(nlhs_in.gt.3) call MATLAB_error("Too many output arguments")

      ! Get sparse matrix
      call matlab_to_fortran(prhs(A_in), mwm, mwn, ptr, row, val, 'A')
      if(mwm.ne.mwn) call MATLAB_error("The matrix must be square")
      n = mwn

      ! Read control in
      num_threads = 0
      matrix_type = HSL_MATRIX_CPLX_SYM
      if(nrhs_in.ge.control_in) then
         ! User supplied control structure
         if(MATLAB_is_structure(prhs(control_in))) then
            call copy_control_in(prhs(control_in), control, num_threads, &
               matrix_type)
         endif
      endif

      ! Convert to HSL standard form
      call mc69_cscl_clean_complex(matrix_type, n, n, ptr, row, &
         flag, val=val)
      select case(flag)
      case(0,1,4,5)
         ! Sucess (expect oor entries as we are throwing away upr triangle!)
      case(-12)
         ! Hermitian but non-zero imaginary component on diagonal
         call MATLAB_context_error("hsl_ma97:ExpectsHermitian", &
            "HSL_MA97 Error: Hermitian factorization requested, but matrix has &
            &non-zero imaginary component on diagonal.")
      case default
         call unknownError("mc69_cscl_clean", flag)
      end select

      ! Obtain rhs vector
      call matlab_to_fortran(prhs(b_in), rhs, temp, mwnrhs, 'b')
      if(mwn .ne. temp) then
         write(errmsg, "(4(a,i5))") "Dimensions of A and b inconsistent: A=", &
            mwn, "x", mwn, ", b=", temp, "x", mwnrhs
         call MATLAB_error(errmsg)
      endif
      nrhs = mwnrhs

      ! Set number of threads
      orig_num_threads = 0
      if(num_threads.ne.0) then
         ! User specified a specific number of threads
         if(num_threads.le.0) &
            call MATLAB_error("control.num_threads must be positive")
!$       orig_num_threads = omp_get_max_threads()
         if(orig_num_threads.eq.0 .and. num_threads.gt.1) then
            call MATLAB_warning( &
               "Parallel run requested, but compiled without OpenMP. Running in&
               & serial")
         endif
!$       if(num_threads .gt. orig_num_threads) &
!$          call MATLAB_warning( &
!$             "control.num_threads exceeds omp_get_max_threads(). Proceeding &
!$             &with control.num_threads anyway.")
!$       call omp_set_num_threads(num_threads)
      endif

      ! Obtain ordering
      allocate(order(n), stat=st)
      if(st.ne.0) goto 100
      if(nrhs_in.ge.P_in) then
         ! User supplied (note: input as double, inverse of order)
         call matlab_to_fortran(prhs(P_in), invp, temp, 'P')
         if(mwn .ne. temp) then
            write(errmsg, "(4(a,i5))") &
               "Dimensions of A and P inconsistent: A=", &
               n, "x", n, ", size(P)=", temp
            call MATLAB_error(errmsg)
         endif
         order(:) = (n+1) ! initialise to invalid value
         do i = 1, n
            order(invp(i)) = i
         end do
         control%ordering = 0 ! user defined
      endif

      ! Call analyse
      call system_clock(t_start)
      call ma97_analyse(.false., n, ptr, row, akeep, control, info, &
         order=order)
      call system_clock(t_stop, t_rate)
      select case (info%flag)
      case(0)
         ! Success. Do nothing
      case(6)
         ! Structually singular, do nothing as will be picked up in factor_solve
      case(-11)
         ! Problem with order
         if(allocated(invp)) then
            write(errmsg, "(a, 10i4)") "Problem with P. P = ", &
               invp(1:min(n,10))
         else
            write(errmsg, "(a, 10i4)") "Bad value of control.ordering = ", &
               control%ordering
         endif
         call MATLAB_error(errmsg)
      case default
!$       if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)
         call unknownError("ma97_analyse", info%flag)
      end select
      select case(info%ordering)
      case(0)
         orderused = 'user'
      case(1)
         orderused = 'AMD'
      case(2)
         orderused = 'MD'
      case(3)
         orderused = 'MeTiS'
      case(4)
         orderused = 'MC47'
      case(7)
         orderused = 'MC80:AMD'
      case(8)
         orderused = 'MC80:METIS'
      end select
      atime = real(t_stop-t_start)/t_rate

      ! Use free MC64 scaling if available and desired
      if(control%scaling.eq.1 .and. &
         (info%ordering.eq.7 .or. info%ordering.eq.8)) control%scaling = 3

      ! Call factor_solve
      call system_clock(t_start)
      call ma97_factor_solve(matrix_type, val, nrhs, rhs, n, &
         akeep, fkeep, control, info, ptr=ptr, row=row)
      call system_clock(t_stop, t_rate)
      select case (info%flag)
      case(0:1)
         ! Success. Do nothing.
      case(-8)
         ! Not positive definite
         call MATLAB_context_error("hsl_ma97:ExpectsPosdef", &
            "HSL_MA97 Error: control.pos_def=true, but matrix is not positive-&
            &definite.")
      case(-9)
         ! IEEE infinities in factorization
         call MATLAB_context_error("hsl_ma97:Infs", &
            "HSL_MA97 Error: IEEE infinities encountered during factorization. &
            &Try larger control.u or control.small.")
      case(-11)
         ! Bad value of ordering.
         call MATLAB_context_error("hsl_ma97:BadOrdering", &
            "HSL_MA97 Error: control.ordering has invalid value")
      case(6,7)
         ! Singular matrix
         call MATLAB_warning("HSL_MA97: Matrix was found to be singular. &
            &RHS not checked for consistency.")
      case(8)
         ! Used Matching-based ordering but not scaling
         call MATLAB_warning("HSL_MA97: Used matching-based ordering but &
            &ignored free scaling.")
      case default
!$       if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)
         call unknownError("ma97_factor_solve", info%flag)
      end select
      ftime = real(t_stop-t_start)/t_rate

      ! Setup info structure if required
      if(nlhs_in.ge.info_out) then
         call copy_info_out(plhs(info_out), info, atime, ftime, &
            orderused)
      endif

      ! Copy solution out
      plhs(x_out) = fortran_to_matlab(rhs)

      ! Restore number of threads
!$    if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)

      return
      100 continue
      call MATLAB_error("Insufficient memory")

   contains
      subroutine copy_info_out(ml_info, info, atime, ftime, order)
         integer(mwp_), intent(inout) :: ml_info
         type(ma97_cinfo), intent(in) :: info
         real(wp), intent(in) :: atime, ftime
         character(len=*) :: order

         integer(mwp_) :: comp
         integer*4 :: szf

         character(len=50), dimension(9) :: fields = (/ &
            "matrix_rank      ", &
            "num_delay        ", &
            "num_factor       ", &
            "num_flops        ", &
            "num_neg          ", &
            "num_two          ", &
            "order            ", &
            "analyse_time     ", &
            "factor_solve_time" /)

         ml_info = MATLAB_create_structure(fields)

         call MATLAB_set_field(ml_info, 'matrix_rank', info%matrix_rank)
         call MATLAB_set_field(ml_info, 'num_delay', info%num_delay)
         call MATLAB_set_field(ml_info, 'num_factor', info%num_factor)
         call MATLAB_set_field(ml_info, 'num_flops', info%num_flops)
         call MATLAB_set_field(ml_info, 'num_neg', info%num_neg)
         call MATLAB_set_field(ml_info, 'num_two', info%num_two)
         call MATLAB_set_field(ml_info, 'order', order)
         call MATLAB_set_field(ml_info, 'analyse_time', atime)
         call MATLAB_set_field(ml_info, 'factor_solve_time', ftime)
      end subroutine copy_info_out

   end subroutine ma97_matlab_backslash_complex

   subroutine ma97_matlab_multiply_real(nlhs_in, plhs, nrhs_in, prhs, akeep, &
         fkeep, order, lflag, work)
      integer*4 :: nlhs_in, nrhs_in
      integer(mwp_) :: plhs(*), prhs(*)
      type(ma97_rakeep), intent(inout) :: akeep
      type(ma97_rfkeep), intent(inout) :: fkeep
      integer, dimension(:), intent(in) :: order
      logical, dimension(:), pointer, intent(inout) :: lflag
      real(wp), dimension(:), pointer, intent(inout) :: work

      integer, parameter :: transpose_in = 1, &
                            inverse_in = 2, &
                            x_in = 3, &
                            control_in = 4
      integer, parameter :: y_out = 1, &
                            info_out = 2

      integer :: i, j, nxi
      integer(mws_) :: mwn, mwnrhs
      integer(mwp_) :: ml_ptr, ml_row, ml_val
      integer, dimension(:), allocatable :: ptr, row
      real(wp), dimension(:,:), allocatable :: xvec, yvec
      integer, dimension(:), allocatable :: xindex
      real(wp), dimension(:), allocatable :: val
      integer :: flag
      integer :: st
      integer :: temp

      type(ma97_rcontrol) :: control
      type(ma97_rinfo) :: info

      character(len=200) :: errmsg
      integer :: n, nrhs, job
      logical :: trans, inv
      integer :: num_threads, orig_num_threads

      real(wp) :: stime
      integer :: t_start, t_stop, t_rate
      integer :: matrix_type

      ! Check number of arguments
      if(nrhs_in.lt.1) call MATLAB_error("Insufficient input arguments")
      if(nrhs_in.gt.4) call MATLAB_error("Too many input arguments")
      if(nlhs_in.lt.1) call MATLAB_error("Insufficient output arguments")
      if(nlhs_in.gt.2) call MATLAB_error("Too many output arguments")

      n = ma97_get_n__(akeep)

      ! Obtain transpose and inverse
      call matlab_to_fortran(prhs(transpose_in), trans, 'transpose')
      call matlab_to_fortran(prhs(inverse_in), inv, 'inverse')

      ! Obtain X vector
      if(matlab_is_sparse(prhs(x_in))) then
         call matlab_to_fortran(prhs(x_in), mwn, mwnrhs, ptr, row, val, 'X')
         if(mwn .ne. n) then
            write(errmsg,"(4(a,i5))") "Dimensions of L and X inconsistent: L=",&
              n, "x", n, ", X=", mwn, "x", mwnrhs
            call MATLAB_error(errmsg)
         endif
         nrhs = mwnrhs
         if(trans .or. .not.inv .or. nrhs.gt.1) then
            call matlab_warning("HSL_MA97: Converted sparse X vector to dense")
            ! Convert to dense
            allocate(xvec(n, nrhs))
            xvec(:,:) = 0
            do i = 1, nrhs
               do j = ptr(i), ptr(i+1)-1
                  xvec(row(j), i) = val(j)
               end do
            end do
            deallocate(ptr, row, val)
         endif
      else
         call matlab_to_fortran(prhs(x_in), xvec, mwn, mwnrhs, 'X')
         if(mwn .ne. n) then
            write(errmsg,"(4(a,i5))") "Dimensions of L and X inconsistent: L=",&
              n, "x", n, ", X=", mwn, "x", mwnrhs
            call MATLAB_error(errmsg)
         endif
         nrhs = mwnrhs
      endif

      ! Copy control in
      num_threads = 0
      matrix_type = HSL_MATRIX_CPLX_SYM
      if(nrhs_in.ge.control_in) then
         ! User supplied control structure
         if(MATLAB_is_structure(prhs(control_in))) then
            call copy_control_in(prhs(control_in), control, num_threads, &
               matrix_type)
         endif
      endif

      ! Set number of threads
      orig_num_threads = 0
      if(num_threads.ne.0) then
         ! User specified a specific number of threads
         if(num_threads.le.0) &
            call MATLAB_error("control.num_threads must be positive")
!$          orig_num_threads = omp_get_max_threads()
         if(orig_num_threads.eq.0 .and. num_threads.gt.1) then
            call MATLAB_warning( &
               "Parallel run requested, but compiled without OpenMP. Running in&
               & serial")
         endif
!$       if(num_threads .gt. orig_num_threads) &
!$          call MATLAB_warning( &
!$             "control.num_threads exceeds omp_get_max_threads(). Proceeding &
!$             &with control.num_threads anyway.")
!$       call omp_set_num_threads(num_threads)
      endif

      if(inv) then
         if(allocated(xvec)) then
            ! dense solve
            if(trans) then
               job = 3 ! bwd solve
            else
               job = 1 ! fwd solve
            endif
            ! Call solve
            call system_clock(t_start)
            call ma97_solve(nrhs, xvec, n, akeep, fkeep, control, info, job=job)
            call system_clock(t_stop, t_rate)
            select case (info%flag)
            case(0)
               ! Success. Do nothing.
            case(7)
               ! Original matrix was singular
               call MATLAB_warning("HSL_MA97: Matrix was singular and &
                  &consistency of rhs has not been checked.")
            case default
!$             if(orig_num_threads .ne. 0) &
!$                call omp_set_num_threads(orig_num_threads)
               call unknownError("ma97_solve", info%flag)
            end select
            stime = real(t_stop-t_start)/t_rate

            ! Copy solution out
            plhs(y_out) = fortran_to_matlab(xvec)
         else
            ! sparse forward solve

            ! First time only initialisation
            if(.not. associated(lflag)) then
               allocate(lflag(n), work(n))
               lflag(:) = .false.
               work(:) = 0
            endif

            ! allocate vectors
            allocate(xvec(n,1), xindex(n))

            ! Call solve
            do j = ptr(1), ptr(1+1)-1
               xvec(row(j), 1) = val(j)
            end do
            call system_clock(t_start)
            call ma97_sparse_fwd_solve(ptr(2)-ptr(1), &
               row(ptr(1):ptr(2)-1), xvec(:,1), order, lflag, nxi, xindex, &
               work, akeep, fkeep, control, info)
            call system_clock(t_stop, t_rate)
            select case (info%flag)
            case(0)
               ! Success. Do nothing.
            case(7)
               ! Original matrix was singular
               call MATLAB_warning("HSL_MA97: Matrix was singular and &
                  &consistency of rhs has not been checked.")
            case default
!$             if(orig_num_threads .ne. 0) &
!$                call omp_set_num_threads(orig_num_threads)
               call unknownError("ma97_solve", info%flag)
            end select
            ! extract soln, reset work arrays
            do j = 1, nxi
               lflag(xindex(j)) = .false.
               xvec(j,1) = work(xindex(j)) 
               work(xindex(j)) = 0
            end do
            stime = real(t_stop-t_start)/t_rate

            ! Copy solution out
            ptr(1) = 1; ptr(2) = nxi+1
            mwn = n
            plhs(y_out) = fortran_to_matlab(mwn, 1_mws_, ptr, xindex, xvec(:,1))
         endif
      else ! .not. inv
         allocate(yvec(n,nrhs))
         call system_clock(t_start)
         call ma97_lmultiply(trans, nrhs, xvec, n, yvec, n, akeep, fkeep, &
            control, info)
         select case (info%flag)
         case(0)
            ! Success. Do nothing.
         case(7)
            ! Original matrix was singular. Do nothing.
         case default
!$          if(orig_num_threads .ne. 0) &
!$             call omp_set_num_threads(orig_num_threads)
            call unknownError("ma97_lmultiply", info%flag)
         end select
         stime = real(t_stop-t_start)/t_rate
         call system_clock(t_stop, t_rate)

         ! Copy product out
         plhs(y_out) = fortran_to_matlab(yvec)
      endif

      ! Setup info structure if required
      if(nlhs_in.ge.info_out) then
         call copy_info_out(plhs(info_out), info, stime)
      endif

      ! Restore number of threads
!$    if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)

      return
      100 continue
      call MATLAB_error("Insufficient memory")
   contains
      subroutine copy_info_out(ml_info, info, stime)
         integer(mwp_), intent(inout) :: ml_info
         type(ma97_rinfo), intent(in) :: info
         real(wp), intent(in) :: stime

         integer(mwp_) :: comp
         integer*4 :: szf

         character(len=50), dimension(1) :: fields = (/ &
            "solve_time    " /)

         ml_info = MATLAB_create_structure(fields)

         call MATLAB_set_field(ml_info, 'solve_time', stime)
      end subroutine copy_info_out
   end subroutine ma97_matlab_multiply_real
   subroutine ma97_matlab_multiply_complex(nlhs_in, plhs, nrhs_in, prhs, akeep,&
         fkeep, order, lflag, work)
      integer*4 :: nlhs_in, nrhs_in
      integer(mwp_) :: plhs(*), prhs(*)
      type(ma97_cakeep), intent(inout) :: akeep
      type(ma97_cfkeep), intent(inout) :: fkeep
      integer, dimension(:), intent(in) :: order
      logical, dimension(:), pointer, intent(inout) :: lflag
      complex(wp), dimension(:), pointer, intent(inout) :: work

      integer, parameter :: transpose_in = 1, &
                            inverse_in = 2, &
                            x_in = 3, &
                            control_in = 4
      integer, parameter :: y_out = 1, &
                            info_out = 2

      integer :: i, j, nxi
      integer(mws_) :: mwn, mwnrhs
      integer(mwp_) :: ml_ptr, ml_row, ml_val
      integer, dimension(:), allocatable :: ptr, row
      complex(wp), dimension(:,:), allocatable :: xvec, yvec
      integer, dimension(:), allocatable :: xindex
      complex(wp), dimension(:), allocatable :: val
      integer :: flag
      integer :: st
      integer :: temp

      type(ma97_ccontrol) :: control
      type(ma97_cinfo) :: info

      character(len=200) :: errmsg
      integer :: n, nrhs, job
      logical :: trans, inv
      integer :: num_threads, orig_num_threads

      real(wp) :: stime
      integer :: t_start, t_stop, t_rate
      integer :: matrix_type

      ! Check number of arguments
      if(nrhs_in.lt.1) call MATLAB_error("Insufficient input arguments")
      if(nrhs_in.gt.4) call MATLAB_error("Too many input arguments")
      if(nlhs_in.lt.1) call MATLAB_error("Insufficient output arguments")
      if(nlhs_in.gt.2) call MATLAB_error("Too many output arguments")

      n = ma97_get_n__(akeep)

      ! Obtain transpose and inverse
      call matlab_to_fortran(prhs(transpose_in), trans, 'transpose')
      call matlab_to_fortran(prhs(inverse_in), inv, 'inverse')

      ! Obtain X vector
      if(matlab_is_sparse(prhs(x_in))) then
         call matlab_to_fortran(prhs(x_in), mwn, mwnrhs, ptr, row, val, 'X')
         if(mwn .ne. n) then
            write(errmsg,"(4(a,i5))") "Dimensions of L and X inconsistent: L=",&
              n, "x", n, ", X=", mwn, "x", mwnrhs
            call MATLAB_error(errmsg)
         endif
         nrhs = mwnrhs
         if(trans .or. .not.inv .or. nrhs.gt.1) then
            call matlab_warning("HSL_MA97: Converted sparse X vector to dense")
            ! Convert to dense
            allocate(xvec(n, nrhs))
            xvec(:,:) = 0
            do i = 1, nrhs
               do j = ptr(i), ptr(i+1)-1
                  xvec(row(j), i) = val(j)
               end do
            end do
            deallocate(ptr, row, val)
         endif
      else
         call matlab_to_fortran(prhs(x_in), xvec, mwn, mwnrhs, 'X')
         if(mwn .ne. n) then
            write(errmsg,"(4(a,i5))") "Dimensions of L and X inconsistent: L=",&
              n, "x", n, ", X=", mwn, "x", mwnrhs
            call MATLAB_error(errmsg)
         endif
         nrhs = mwnrhs
      endif

      ! Copy control in
      num_threads = 0
      matrix_type = HSL_MATRIX_CPLX_SYM
      if(nrhs_in.ge.control_in) then
         ! User supplied control structure
         if(MATLAB_is_structure(prhs(control_in))) then
            call copy_control_in(prhs(control_in), control, num_threads, &
               matrix_type)
         endif
      endif

      ! Set number of threads
      orig_num_threads = 0
      if(num_threads.ne.0) then
         ! User specified a specific number of threads
         if(num_threads.le.0) &
            call MATLAB_error("control.num_threads must be positive")
!$          orig_num_threads = omp_get_max_threads()
         if(orig_num_threads.eq.0 .and. num_threads.gt.1) then
            call MATLAB_warning( &
               "Parallel run requested, but compiled without OpenMP. Running in&
               & serial")
         endif
!$       if(num_threads .gt. orig_num_threads) &
!$          call MATLAB_warning( &
!$             "control.num_threads exceeds omp_get_max_threads(). Proceeding &
!$             &with control.num_threads anyway.")
!$       call omp_set_num_threads(num_threads)
      endif

      if(inv) then
         if(allocated(xvec)) then
            ! dense solve
            if(trans) then
               job = 3 ! bwd solve
            else
               job = 1 ! fwd solve
            endif
            ! Call solve
            call system_clock(t_start)
            call ma97_solve(nrhs, xvec, n, akeep, fkeep, control, info, job=job)
            call system_clock(t_stop, t_rate)
            select case (info%flag)
            case(0)
               ! Success. Do nothing.
            case(7)
               ! Original matrix was singular
               call MATLAB_warning("HSL_MA97: Matrix was singular and &
                  &consistency of rhs has not been checked.")
            case default
!$             if(orig_num_threads .ne. 0) &
!$                call omp_set_num_threads(orig_num_threads)
               call unknownError("ma97_solve", info%flag)
            end select
            stime = real(t_stop-t_start)/t_rate

            ! Copy solution out
            plhs(y_out) = fortran_to_matlab(xvec)
         else
            ! sparse forward solve

            ! First time only initialisation
            if(.not. associated(lflag)) then
               allocate(lflag(n), work(n))
               lflag(:) = .false.
               work(:) = 0
            endif

            ! allocate vectors
            allocate(xvec(n,1), xindex(n))

            ! Call solve
            do j = ptr(1), ptr(1+1)-1
               xvec(row(j), 1) = val(j)
            end do
            call system_clock(t_start)
            call ma97_sparse_fwd_solve(ptr(2)-ptr(1), &
               row(ptr(1):ptr(2)-1), xvec(:,1), order, lflag, nxi, xindex, &
               work, akeep, fkeep, control, info)
            call system_clock(t_stop, t_rate)
            select case (info%flag)
            case(0)
               ! Success. Do nothing.
            case(7)
               ! Original matrix was singular
               call MATLAB_warning("HSL_MA97: Matrix was singular and &
                  &consistency of rhs has not been checked.")
            case default
!$             if(orig_num_threads .ne. 0) &
!$                call omp_set_num_threads(orig_num_threads)
               call unknownError("ma97_solve", info%flag)
            end select
            ! extract soln, reset work arrays
            do j = 1, nxi
               lflag(xindex(j)) = .false.
               xvec(j,1) = work(xindex(j)) 
               work(xindex(j)) = 0
            end do
            stime = real(t_stop-t_start)/t_rate

            ! Copy solution out
            ptr(1) = 1; ptr(2) = nxi+1
            mwn = n
            plhs(y_out) = fortran_to_matlab(mwn, 1_mws_, ptr, xindex, xvec(:,1))
         endif
      else ! .not. inv
         allocate(yvec(n,nrhs))
         call system_clock(t_start)
         call ma97_lmultiply(trans, nrhs, xvec, n, yvec, n, akeep, fkeep, &
            control, info)
         select case (info%flag)
         case(0)
            ! Success. Do nothing.
         case(7)
            ! Original matrix was singular. Do nothing.
         case default
!$          if(orig_num_threads .ne. 0) &
!$             call omp_set_num_threads(orig_num_threads)
            call unknownError("ma97_lmultiply", info%flag)
         end select
         stime = real(t_stop-t_start)/t_rate
         call system_clock(t_stop, t_rate)

         ! Copy product out
         plhs(y_out) = fortran_to_matlab(yvec)
      endif

      ! Setup info structure if required
      if(nlhs_in.ge.info_out) then
         call copy_info_out(plhs(info_out), info, stime)
      endif

      ! Restore number of threads
!$    if(orig_num_threads .ne. 0) call omp_set_num_threads(orig_num_threads)

      return
      100 continue
      call MATLAB_error("Insufficient memory")
   contains
      subroutine copy_info_out(ml_info, info, stime)
         integer(mwp_), intent(inout) :: ml_info
         type(ma97_cinfo), intent(in) :: info
         real(wp), intent(in) :: stime

         integer(mwp_) :: comp
         integer*4 :: szf

         character(len=50), dimension(1) :: fields = (/ &
            "solve_time    " /)

         ml_info = MATLAB_create_structure(fields)

         call MATLAB_set_field(ml_info, 'solve_time', stime)
      end subroutine copy_info_out
   end subroutine ma97_matlab_multiply_complex


   ! ma97_expert('destroy', handle)
   ! arguments 'destroy' and handle already removed
   subroutine ma97_matlab_finalise_real(akeep, fkeep)
      type(ma97_rakeep), intent(inout) :: akeep
      type(ma97_rfkeep), intent(inout) :: fkeep

      type(ma97_rcontrol) :: control

      call ma97_finalise(akeep, fkeep)
   end subroutine ma97_matlab_finalise_real

   ! ma97_expert('destroy', handle)
   ! arguments 'destroy' and handle already removed
   subroutine ma97_matlab_finalise_complex(akeep, fkeep)
      type(ma97_cakeep), intent(inout) :: akeep
      type(ma97_cfkeep), intent(inout) :: fkeep

      type(ma97_ccontrol) :: control

      call ma97_finalise(akeep, fkeep)
   end subroutine ma97_matlab_finalise_complex

end module ma97_matlab_main

! This module looks after a SAVEd set of variables mapping integer handles
! to Fortran akeep and fkeep variables
module ma97_handles
   use hsl_ma97_double, only: ma97_rakeep => ma97_akeep, &
                              ma97_rfkeep => ma97_fkeep
   use hsl_ma97_double_complex, only: ma97_cakeep => ma97_akeep, &
                                      ma97_cfkeep => ma97_fkeep
   use ma97_matlab_main
   implicit none

   ! Data associated with the handle
   ! Considered to be empty if both associated(rakeep) and associated(cakeep)
   ! are .false.
   type ma97_hdl
      type(ma97_rakeep), pointer :: rakeep => null()
      type(ma97_rfkeep), pointer :: rfkeep => null()
      type(ma97_cakeep), pointer :: cakeep => null()
      type(ma97_cfkeep), pointer :: cfkeep => null()
      integer, pointer           :: order(:) => null()
      logical, pointer           :: lflag(:) => null()
      real(wp), pointer          :: rwork(:) => null()
      complex(wp), pointer       :: cwork(:) => null()
   end type ma97_hdl

   ! How many handles initally and how much increase once exhausted
   integer, parameter :: initial_handles = 5
   double precision, parameter :: multiplier = 2.0

   ! SAVEd data
   integer, save :: next_handle = 1
   integer, save :: total_handles = 0
   type(ma97_hdl), dimension(:), allocatable, save :: handles

contains

   integer function ma97_new_handle(complexFlag)
      logical, intent(in) :: complexFlag

      type(ma97_hdl), dimension(:), allocatable :: temp
      integer :: i

      ! Do we need to expand the number of available handles?
      if (next_handle .gt. total_handles) then
         if(total_handles.ne.0) then
            ! Need to expand existing handle selection
            allocate(temp(total_handles))
            do i = 1, total_handles
               temp(i)%rakeep => handles(i)%rakeep
               temp(i)%rfkeep => handles(i)%rfkeep
               temp(i)%cakeep => handles(i)%cakeep
               temp(i)%cfkeep => handles(i)%cfkeep
               temp(i)%order  => handles(i)%order
               temp(i)%lflag  => handles(i)%lflag
               temp(i)%rwork  => handles(i)%rwork
               temp(i)%cwork  => handles(i)%cwork
            end do
            deallocate(handles)
            total_handles = max(int(multiplier*total_handles), total_handles)
            allocate(handles(total_handles))
            do i = 1, size(temp)
               handles(i)%rakeep => temp(i)%rakeep
               handles(i)%rfkeep => temp(i)%rfkeep
               handles(i)%cakeep => temp(i)%cakeep
               handles(i)%cfkeep => temp(i)%cfkeep
               handles(i)%order  => temp(i)%order
               handles(i)%lflag  => temp(i)%lflag
               handles(i)%rwork  => temp(i)%rwork
               handles(i)%cwork  => temp(i)%cwork
            end do
            deallocate(temp)
         else
            ! First call since module loaded
            total_handles = initial_handles
            allocate(handles(total_handles))

            ! Register clean function
            call mexAtExit(cleanup_all_handles)
         endif
      endif

      ma97_new_handle = next_handle
      if(complexFlag) then
         allocate(handles(next_handle)%cakeep)
         allocate(handles(next_handle)%cfkeep)
      else
         allocate(handles(next_handle)%rakeep)
         allocate(handles(next_handle)%rfkeep)
      endif
      next_handle = next_handle + 1
   end function ma97_new_handle

   ! This routine is called at unload of this module from MATLAB.
   ! It shuold cleanup all SAVEd data
   subroutine cleanup_all_handles()
      integer :: i, st

      do i = 1, next_handle-1
         call cleanup_handle(i)
      end do
   end subroutine cleanup_all_handles

   ! Destroy the data associated with a handle.
   ! Recover all free pointers at end of handle list.
   subroutine cleanup_handle(handle)
      integer, intent(in) :: handle

      integer :: current
      integer :: st

      if(associated(handles(handle)%rakeep)) then
         call ma97_matlab_finalise(handles(handle)%rakeep, &
            handles(handle)%rfkeep)
         deallocate(handles(handle)%rakeep)
         deallocate(handles(handle)%rfkeep)
         if(associated(handles(handle)%rwork)) &
            deallocate(handles(handle)%rwork)
         nullify(handles(handle)%rakeep)
         nullify(handles(handle)%rfkeep)
         nullify(handles(handle)%rwork)
      endif
      if(associated(handles(handle)%cakeep)) then
         call ma97_matlab_finalise(handles(handle)%cakeep, &
            handles(handle)%cfkeep)
         deallocate(handles(handle)%cakeep)
         deallocate(handles(handle)%cfkeep)
         if(associated(handles(handle)%cwork)) &
            deallocate(handles(handle)%cwork)
         nullify(handles(handle)%cakeep)
         nullify(handles(handle)%cfkeep)
         nullify(handles(handle)%cwork)
      endif
      if(associated(handles(handle)%order)) &
         deallocate(handles(handle)%order)
      if(associated(handles(handle)%lflag)) &
         deallocate(handles(handle)%lflag)
      nullify(handles(handle)%order)
      nullify(handles(handle)%lflag)

      do current = handle, 1, -1
         if(current.ne.next_handle-1) exit
         if(.not.associated(handles(current)%rakeep) .and. .not. &
               associated(handles(current)%cakeep)) then
            ! Current "last" element is unallocated, make it next available
            ! element.
            next_handle = next_handle - 1
         endif
      end do
   end subroutine cleanup_handle

end module ma97_handles


! Gateway routine
! Strip first argument and converts any handles, then calls relevant routine
subroutine mexFunction(nlhs_in, plhs, nrhs_in, prhs)
   use hsl_matlab
   use hsl_mc68_double
   use hsl_mc69_double
   use hsl_ma97_double
   use ma97_handles
   use ma97_matlab_main
   implicit none

   integer*4 :: nlhs_in, nrhs_in, tlhs, trhs
   integer(mwp_) :: plhs(*), prhs(*)
   integer(mws_) :: mwstemp

   character(len=200) :: act
   integer :: handle
   integer, dimension(1) :: thandle
   integer :: st
   integer :: n

   character(len=200) :: errmsg

   if(nrhs_in.lt.1) call MATLAB_error("Insufficient input arguments")

   if(.not. MATLAB_is_character(prhs(1))) &
      call MATLAB_error("First argument must be string")
   mwstemp = len(act)
   call MATLAB_get_string(prhs(1), act, mwstemp)

   select case(trim(act))
   case('factor')
      ! [handle, info] = ma97_expert('factor', A[, P])
      ! At least handle is required for output

      ! Setup handle and store in first output argument
      if(nlhs_in.lt.1) call MATLAB_error("Insufficient output arguments")
      handle = ma97_new_handle(MATLAB_is_complex(prhs(2)))
      plhs(1) = fortran_to_matlab(handle)

      ! Call work routine
      if(associated(handles(handle)%rakeep)) then
         call ma97_matlab_analyse_factor_real(nlhs_in-1_int4_, plhs(2), &
            nrhs_in-1_int4_, prhs(2), handles(handle)%rakeep, &
            handles(handle)%rfkeep, handles(handle)%order)
      else
         call ma97_matlab_analyse_factor_complex(nlhs_in-1_int4_, plhs(2), &
            nrhs_in-1_int4_, prhs(2), handles(handle)%cakeep, &
            handles(handle)%cfkeep, handles(handle)%order)
      endif

   case('backslash')
      ! [x, info, handle] = ma97_expert('backslash', A, b[, control, P])
      ! At least soln is required for output

      ! Setup handle and store in third output argument, if present
      handle = ma97_new_handle(MATLAB_is_complex(prhs(2)))
      if(nlhs_in.ge.3) plhs(3) = fortran_to_matlab(handle)

      ! Call work routine
      if(associated(handles(handle)%rakeep)) then
         call ma97_matlab_backslash_real(nlhs_in, plhs, nrhs_in-1_int4_, &
            prhs(2), handles(handle)%rakeep, handles(handle)%rfkeep, &
            handles(handle)%order)
      else
         call ma97_matlab_backslash_complex(nlhs_in, plhs, nrhs_in-1_int4_, &
            prhs(2), handles(handle)%cakeep, handles(handle)%cfkeep, &
            handles(handle)%order)
      endif

      ! Destroy handle if it is not being returned to the user
      if(nlhs_in.lt.3) call cleanup_handle(handle)

   case('solve')
      ! [x, info, consist] = ma97_expert('solve', handle, b)

      ! Obtain and check handle
      if(nrhs_in.lt.2) call MATLAB_error("Insufficient input arguments")
      call matlab_to_fortran(prhs(2), handle, 'handle')
      if(handle.lt.1 .or. handle.gt.total_handles) &
         call MATLAB_error("Invalid handle")

      ! Call worker routine based on allocated parts of handle
      if(associated(handles(handle)%rakeep)) then
         call ma97_matlab_solve_real(nlhs_in, plhs(1), nrhs_in-2_int4_, &
            prhs(3), handles(handle)%rakeep, handles(handle)%rfkeep)
      elseif(associated(handles(handle)%cakeep)) then
         call ma97_matlab_solve_complex(nlhs_in, plhs(1), nrhs_in-2_int4_, &
            prhs(3), handles(handle)%cakeep, handles(handle)%cfkeep)
      else
         ! No parts allocated, probably a destroyed or unallocated handle
         call MATLAB_error("Invalid handle")
      endif

   case('multiply')
      ! [Y, info] = hsl_ma97_expert('multiply', handle, transpose, inverse,
      !                             X[, control])
      
      ! Obtain and check handle
      if(nrhs_in.lt.2) call MATLAB_error("Insufficient input arguments")
      call matlab_to_fortran(prhs(2), handle, 'handle')
      if(handle.lt.1 .or. handle.gt.total_handles) &
         call MATLAB_error("Invalid handle")

      ! Call worker routine based on allocated parts of handle
      if(associated(handles(handle)%rakeep)) then
         call ma97_matlab_multiply_real(nlhs_in, plhs(1), nrhs_in-2_int4_, &
            prhs(3), handles(handle)%rakeep, handles(handle)%rfkeep, &
            handles(handle)%order, handles(handle)%lflag, handles(handle)%rwork)
      elseif(associated(handles(handle)%cakeep)) then
         call ma97_matlab_multiply_complex(nlhs_in, plhs(1), nrhs_in-2_int4_, &
            prhs(3), handles(handle)%cakeep, handles(handle)%cfkeep, &
            handles(handle)%order, handles(handle)%lflag, handles(handle)%cwork)
      else
         ! No parts allocated, probably a destroyed or unallocated handle
         call MATLAB_error("Invalid handle")
      endif

   case('destroy')
      ! ma97_expert('destroy', handle)

      ! Obtain and check handle
      if(nrhs_in.lt.2) call MATLAB_error("Insufficient input arguments")
      call matlab_to_fortran(prhs(2), handle, 'handle')
      if(handle.lt.1 .or. handle.gt.total_handles) &
         call MATLAB_error("Invalid handle")

      ! Destroy anything that is there
      call cleanup_handle(handle)

   case default
      write(errmsg, "(3a)") "Unrecognised action: '", trim(act), "'"
      call MATLAB_error(errmsg)
   end select

   return
   100 continue
   call MATLAB_error("Insufficient memory")
end subroutine mexFunction

