! COPYRIGHT (c) 2010 Science and Technology Facilities Council
! Original date 18 January 2011. Version 1.0.0
!
! Written by: Jonathan Hogg, John Reid, and Sue Thorne
!
! Version 1.4.1
! For version history, see ChangeLog
!
module hsl_mc69_single
   implicit none

   private
   public :: HSL_MATRIX_UNDEFINED,                             &
      HSL_MATRIX_REAL_RECT, HSL_MATRIX_CPLX_RECT,              &
      HSL_MATRIX_REAL_UNSYM, HSL_MATRIX_CPLX_UNSYM,            &
      HSL_MATRIX_REAL_SYM_PSDEF, HSL_MATRIX_CPLX_HERM_PSDEF,   &
      HSL_MATRIX_REAL_SYM_INDEF, HSL_MATRIX_CPLX_HERM_INDEF,   &
      HSL_MATRIX_CPLX_SYM,                                     &
      HSL_MATRIX_REAL_SKEW, HSL_MATRIX_CPLX_SKEW
   public :: mc69_cscl_clean, mc69_verify, mc69_print, mc69_csclu_convert, &
      mc69_coord_convert, mc69_set_values, mc69_csrlu_convert, &
      mc69_cscl_convert, mc69_cscu_convert, mc69_csru_convert, &
      mc69_csrl_convert

   integer, parameter :: wp = kind(0.0)
   real(wp), parameter :: zero = 0.0_wp

   ! matrix types : real
   integer, parameter :: HSL_MATRIX_UNDEFINED      =  0 ! undefined/unknown
   integer, parameter :: HSL_MATRIX_REAL_RECT      =  1 ! real rectangular
   integer, parameter :: HSL_MATRIX_REAL_UNSYM     =  2 ! real unsymmetric
   integer, parameter :: HSL_MATRIX_REAL_SYM_PSDEF =  3 ! real symmetric pos def
   integer, parameter :: HSL_MATRIX_REAL_SYM_INDEF =  4 ! real symmetric indef
   integer, parameter :: HSL_MATRIX_REAL_SKEW      =  6 ! real skew symmetric

   ! matrix types : complex
   integer, parameter :: HSL_MATRIX_CPLX_RECT      = -1 ! complex rectangular
   integer, parameter :: HSL_MATRIX_CPLX_UNSYM     = -2 ! complex unsymmetric
   integer, parameter :: HSL_MATRIX_CPLX_HERM_PSDEF= -3 ! hermitian pos def
   integer, parameter :: HSL_MATRIX_CPLX_HERM_INDEF= -4 ! hermitian indef
   integer, parameter :: HSL_MATRIX_CPLX_SYM       = -5 ! complex symmetric
   integer, parameter :: HSL_MATRIX_CPLX_SKEW      = -6 ! complex skew symmetric

   ! Error flags
   integer, parameter :: MC69_SUCCESS                =  0
   integer, parameter :: MC69_ERROR_ALLOCATION       = -1
   integer, parameter :: MC69_ERROR_MATRIX_TYPE      = -2
   integer, parameter :: MC69_ERROR_N_OOR            = -3
   integer, parameter :: MC69_ERROR_M_NE_N           = -4
   integer, parameter :: MC69_ERROR_PTR_1            = -5
   integer, parameter :: MC69_ERROR_PTR_MONO         = -6
   integer, parameter :: MC69_ERROR_ROW_BAD_ORDER    = -7
   integer, parameter :: MC69_ERROR_ROW_OOR          = -8
   integer, parameter :: MC69_ERROR_ROW_DUP          = -9
   integer, parameter :: MC69_ERROR_ALL_OOR          = -10
   integer, parameter :: MC69_ERROR_MISSING_DIAGONAL = -11
   integer, parameter :: MC69_ERROR_IMAG_DIAGONAL    = -12
   integer, parameter :: MC69_ERROR_MISMATCH_LWRUPR  = -13
   integer, parameter :: MC69_ERROR_UPR_ENTRY        = -14
   integer, parameter :: MC69_ERROR_VAL_MISS         = -15
   integer, parameter :: MC69_ERROR_LMAP_MISS        = -16


   ! warning flags
   integer, parameter :: MC69_WARNING_IDX_OOR          = 1
   integer, parameter :: MC69_WARNING_DUP_IDX          = 2
   integer, parameter :: MC69_WARNING_DUP_AND_OOR      = 3
   integer, parameter :: MC69_WARNING_MISSING_DIAGONAL = 4
   integer, parameter :: MC69_WARNING_MISS_DIAG_OORDUP = 5

!            Possible error returns:

! MC69_ERROR_ALLOCATION         Allocation error
! MC69_ERROR_MATRIX_TYPE        Problem with matrix_type
! MC69_ERROR_N_OOR              n < 0 or m < 0
! MC69_ERROR_PTR_1              ptr(1) < 1
! MC69_ERROR_PTR_MONO           Error in ptr (not monotonic)
! MC69_ERROR_VAL_MISS           Only one of val_in and val_out is present
! MC69_ERROR_ALL_OOR            All the variables in a column are out-of-range
! MC69_ERROR_IMAG_DIAGONAL      Hermitian case and diagonal not real
! MC69_ERROR_ROW_BAD_ORDER      Entries within a column are not sorted by
!                               increasing row index 
! MC69_ERROR_MISMATCH_LWRUPR    Symmetric, skew symmetric or Hermitian: 
!                               entries in upper and lower
!                               triangles do not match
! MC69_ERROR_MISSING_DIAGONAL   Pos def and diagonal entries missing 
! MC69_ERROR_ROW_OOR            Row contains out-of-range entries      
! MC69_ERROR_ROW_DUP            Row contains duplicate entries         
! MC69_ERROR_M_NE_N             Square matrix and m .ne. n        

!           Possible warnings:

! MC69_WARNING_IDX_OOR          Out-of-range variable indices
! MC69_WARNING_DUP_IDX          Duplicated variable indices
! MC69_WARNING_DUP_AND_OOR      out of range and duplicated indices
! MC69_WARNING_MISSING_DIAGONAL Indefinite case and diagonal entries missing
! MC69_WARNING_MISS_DIAG_OORDUP As MC69_WARNING_MISSING_DIAGONAL, and 
!                               out-of-range and/or duplicates


!!!!!!!!!!!!!!!!!!!!!!!!
! Internal types
!!!!!!!!!!!!!!!!!!!!!!!!
type dup_list
   integer :: src
   integer :: dest
   type(dup_list), pointer :: next => null()
end type dup_list

interface mc69_verify
   module procedure mc69_verify_single
end interface mc69_verify
interface mc69_print
   module procedure mc69_print_single
end interface
interface mc69_cscl_clean
   module procedure mc69_cscl_clean_single
end interface mc69_cscl_clean
interface mc69_set_values
   module procedure mc69_set_values_single
end interface mc69_set_values
interface mc69_cscl_convert
   module procedure mc69_cscl_convert_single
end interface mc69_cscl_convert
interface mc69_cscu_convert
   module procedure mc69_cscu_convert_single
end interface mc69_cscu_convert
interface mc69_csclu_convert
   module procedure mc69_csclu_convert_single
end interface mc69_csclu_convert
interface mc69_csrl_convert
   module procedure mc69_csrl_convert_single
end interface mc69_csrl_convert
interface mc69_csru_convert
   module procedure mc69_csru_convert_single
end interface mc69_csru_convert
interface mc69_csrlu_convert
   module procedure mc69_csrlu_convert_single
end interface mc69_csrlu_convert
interface mc69_coord_convert
   module procedure mc69_coord_convert_single
end interface mc69_coord_convert

contains

!
! To verify that a matrix is in HSL format, or identify why it is not
!
subroutine mc69_verify_single(lp, matrix_type, m, n, ptr, row, flag, more, val)
   integer, intent(in) :: lp ! output unit
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: m ! number of rows
   integer, intent(in) :: n ! number of columns
   integer, dimension(*), intent(in) :: ptr ! column starts
   integer, dimension(*), intent(in) :: row ! row indices.
     ! Entries within each column must be sorted in order of 
     ! increasing row index. no duplicates and/or out of range entries allowed.
   integer, intent(out) :: flag ! return code
   integer, intent(out) :: more ! futher error information (or set to 0)
   real(wp), dimension(:), optional, intent(in) :: val ! matrix values,if any

   integer :: col ! current column
   character(50)  :: context  ! Procedure name (used when printing).
   logical :: diag ! flag for detection of diagonal
   integer :: i
   integer :: j
   integer :: k
   integer :: last ! last row index
   logical :: lwronly
   integer, dimension(:), allocatable :: ptr2
   integer :: st

   context = 'mc69_verify'
   flag = MC69_SUCCESS

   more = 0 ! ensure more is not undefined.

   ! Check matrix_type
   select case(matrix_type)
   case(0)
      ! Undefined matrix. OK, do nothing
   case(1:4,6)
      ! Real matrix. OK, do nothing
   case(-6:-1)
      ! Complex matrix. Issue error
      flag = MC69_ERROR_MATRIX_TYPE
      call mc69_print_flag(context,lp,flag)
      return
   case default
      ! Out of range value. Issue error
      flag = MC69_ERROR_MATRIX_TYPE
      call mc69_print_flag(context,lp,flag)
      return
   end select

   ! Check m and n are valid; skip remaining tests if n=0
   if(n < 0 .or. m < 0) then
      flag = MC69_ERROR_N_OOR
      call mc69_print_flag(context,lp,flag)
      return
   end if
   if(abs(matrix_type).ne.HSL_MATRIX_REAL_RECT .and. m.ne.n) then
      flag = MC69_ERROR_M_NE_N
      call mc69_print_flag(context,lp,flag)
      return
   endif
   if(n == 0 .or. m == 0) return

   ! Check ptr(1) is valid
   if(ptr(1) < 1) then
      more = ptr(1)
      flag = MC69_ERROR_PTR_1
      call mc69_print_flag(context,lp,flag)
      return
   end if

   ! Check ptr is monotonically increasing
   do i = 1, n
      if(ptr(i+1).ge.ptr(i)) cycle ! ptr is monotonically increasing, good.
      flag = MC69_ERROR_PTR_MONO
      more = i + 1
      call mc69_print_flag(context,lp,flag)
      return
   end do

   ! Count number of entries in each row. Also:
   ! * Check ordering of entries
   ! * Check for out-of-range entries
   ! * Check for duplicate entries
   ! * Lack of diagonal entries in real pos. def. case.

   ! ptr2(k+2) is set to hold the number of entries in row k
   allocate(ptr2(m+2), stat=st)
   if(st.ne.0) goto 100
   ptr2(:) = 0

   lwronly = (abs(matrix_type).ne.HSL_MATRIX_UNDEFINED) .and. &
             (abs(matrix_type).ne.HSL_MATRIX_REAL_RECT) .and. &
             (abs(matrix_type).ne.HSL_MATRIX_REAL_UNSYM)
   do col = 1, n
      last = -1
      diag = .false.
      do j = ptr(col), ptr(col+1)-1
         k = row(j)
         ! Check out-of-range
         if(k.lt.1 .or. k.gt.m) then
            flag = MC69_ERROR_ROW_OOR
            more = j
            call mc69_print_flag(context,lp,flag)
            return
         endif
         if(lwronly .and. k.lt.col) then
            flag = MC69_ERROR_UPR_ENTRY
            more = j
            call mc69_print_flag(context,lp,flag)
            return
         endif
         if(abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW .and. k.eq.col) then
            flag = MC69_ERROR_UPR_ENTRY
            more = j
            call mc69_print_flag(context,lp,flag)
            return
         endif
         ! Check duplicate
         if(k.eq.last) then
            flag = MC69_ERROR_ROW_DUP
            more = j-1
            call mc69_print_flag(context,lp,flag)
            return
         endif
         ! Check order
         if(k.lt.last) then
            flag = MC69_ERROR_ROW_BAD_ORDER
            more = j
            call mc69_print_flag(context,lp,flag)
            return
         endif
         ! Check for diagonal
         diag = diag .or. (k.eq.col)
         ! Increase count for row k
         ptr2(k+2) = ptr2(k+2) + 1
         ! Store value for next loop
         last = k
      end do
      ! If marked as positive definite, check if diagonal was present
      if(abs(matrix_type).eq.HSL_MATRIX_REAL_SYM_PSDEF .and. .not.diag) then
         flag = MC69_ERROR_MISSING_DIAGONAL
         more = col
         call mc69_print_flag(context,lp,flag)
         return
      endif
   end do

   if(present(val)) then
      select case(matrix_type)
      case(HSL_MATRIX_REAL_SYM_PSDEF)
         ! Check for positive diagonal entries
         do j = 1,n
            k = ptr(j)
            ! Note: column cannot be empty as previously checked pattern ok
            if(real(val(k))<=zero) then
               flag = MC69_ERROR_MISSING_DIAGONAL
               more = j
               call mc69_print_flag(context,lp,flag)
               return
            end if 
         end do
      end select
   endif

   return

   100 continue ! Allocation error
   flag = MC69_ERROR_ALLOCATION
   more = st
   call mc69_print_flag(context,lp,flag)
   return

end subroutine mc69_verify_single

!****************************************

!
! Pretty prints a matrix as best it can
!
subroutine mc69_print_single(lp, lines, matrix_type, m, n, ptr, row, val, cbase)
   integer, intent(in) :: lp ! unit to print on
   integer, intent(in) :: lines ! max number of lines to use (ignored if -ive)
   integer, intent(in) :: matrix_type ! type of matrix
   integer, intent(in) :: m ! number of rows in matrix
   integer, intent(in) :: n ! number of cols in matrix
   integer, dimension(n+1), intent(in) :: ptr ! column pointers
   integer, dimension(ptr(n+1)-1), intent(in) :: row ! row indices
   real(wp), dimension(ptr(n+1)-1), optional, intent(in) :: val ! matrix vals
   logical, optional, intent(in) :: cbase ! is true, rebase for C

   integer :: col, j, k
   integer :: llines
   integer, dimension(:,:), allocatable :: dmat
   character(len=5) :: mfrmt, nfrmt, nefrmt
   character(len=12) :: negfrmt, valfrmt, emptyfrmt
   integer ::  rebase

   if(lp.lt.0) return ! invalid unit

   ! Check if we need to rebase for C consistent output
   rebase = 0
   if(present(cbase)) then
      if(cbase) rebase = 1
   endif

   ! Calculate number of lines to play with
   llines = huge(llines)
   if(lines.gt.0) llines = lines

   ! Print a summary statement about the matrix
   mfrmt = digit_format(m)
   nfrmt = digit_format(n)
   nefrmt = digit_format(ptr(n+1)-1)

   select case(matrix_type)
   case(HSL_MATRIX_UNDEFINED)
      write(lp, "(a)", advance="no") &
         "Matrix of undefined type, dimension "
   case(HSL_MATRIX_REAL_RECT)
      write(lp, "(a)", advance="no") &
         "Real rectangular matrix, dimension "
   case(HSL_MATRIX_REAL_UNSYM)
      write(lp, "(a)", advance="no") &
         "Real unsymmetric matrix, dimension "
   case(HSL_MATRIX_REAL_SYM_PSDEF)
      write(lp, "(a)", advance="no") &
         "Real symmetric positive definite matrix, dimension "
   case(HSL_MATRIX_REAL_SYM_INDEF)
      write(lp, "(a)", advance="no") &
         "Real symmetric indefinite matrix, dimension "
   case(HSL_MATRIX_REAL_SKEW)
      write(lp, "(a)", advance="no") &
         "Real skew symmetric matrix, dimension "
   case default
      write(lp, "(a,i5)") &
         "Unrecognised matrix_type = ", matrix_type
      return
   end select
   write(lp, mfrmt, advance="no") m
   write(lp, "(a)", advance="no") "x"
   write(lp, nfrmt, advance="no") n
   write(lp, "(a)", advance="no") " with "
   write(lp, nefrmt, advance="no") ptr(n+1)-1
   write(lp, "(a)") " entries."

   ! immediate return if m = 0 or n = 0
   if (m == 0 .or. n == 0) return

   if(((present(val) .and. n.lt.10) .or. (.not.present(val) .and. n.lt.24)) &
         .and. m+1.le.llines) then
      ! Print the matrix as if dense
      allocate(dmat(m, n))
      dmat(:,:) = 0
      do col = 1, n
         do j = ptr(col), ptr(col+1) - 1
            k = row(j)
            if(abs(matrix_type).ge.HSL_MATRIX_REAL_SYM_PSDEF) then
               dmat(col, k) = -j
            endif
            dmat(k, col) = j
         end do
      end do

      select case(n)
      case(:6)
         valfrmt = "(1x,es12.4)"
         negfrmt = valfrmt
         emptyfrmt = "(1x,a12)"
      case(7)
         valfrmt = "(1x,es10.2)"
         negfrmt = valfrmt
         emptyfrmt = "(1x,a10)"
      case(8:)
         valfrmt = "(1x,es8.2)"
         negfrmt = "(1x,es8.1)"
         emptyfrmt = "(1x,a8)"
      end select

      do k = 1, m
         write(lp,mfrmt,advance="no") k-rebase
         write(lp,"(':')",advance="no")
         if(present(val)) then
            do j = 1, n
               if(dmat(k,j).eq.0) then 
                  ! nothing here
                  write(lp,emptyfrmt,advance="no") '                '
               elseif(dmat(k,j).gt.0) then
                  if(val(dmat(k,j)).gt.zero) then
                     write(lp,valfrmt,advance="no") val(dmat(k,j))
                  else
                     write(lp,negfrmt,advance="no") val(dmat(k,j))
                  endif
               else
                  ! in upper triangle
                  select case(matrix_type)
                  case(HSL_MATRIX_REAL_SYM_INDEF,HSL_MATRIX_REAL_SYM_PSDEF)
                     if(val(-dmat(k,j)).gt.zero) then
                        write(lp,valfrmt,advance="no") val(-dmat(k,j))
                     else
                        write(lp,negfrmt,advance="no") val(-dmat(k,j))
                     endif
                  case(HSL_MATRIX_REAL_SKEW)
                     if(-val(-dmat(k,j)).gt.zero) then
                        write(lp,valfrmt,advance="no") -val(-dmat(k,j))
                     else
                        write(lp,negfrmt,advance="no") -val(-dmat(k,j))
                     endif
                  end select
               endif
            end do
         else ! pattern only
            do j = 1, n
               if(dmat(k,j).eq.0) then 
                  ! nothing here
                  write(lp,"('  ')",advance="no")
               else
                  write(lp,"(1x,'x')",advance="no")
               endif
            end do
         end if
         write(lp,"()")
      end do
   else
      ! Output first x entries from each column
      llines = llines - 1 ! Discount first info line
      if(llines.le.2) return
      write(lp, "(a)") "First 4 entries in columns:"
      llines = llines - 1
      do col = 1, llines
         write(lp, "(a)", advance="no") "Col "
         write(lp, nfrmt, advance="no") col-rebase
         write(lp, "(':')", advance="no")
         do j = ptr(col), min(ptr(col+1)-1, ptr(col)+3)
            write(lp, "('  ')", advance="no")
            write(lp, mfrmt, advance="no") row(j)-rebase
            if(present(val)) then
               write(lp, "(' (',es12.4,')')", advance="no") val(j)
            endif
         end do
         write(lp, "()")
      end do
   endif
end subroutine mc69_print_single

!****************************************

character(len=5) function digit_format(x)
   integer, intent(in) :: x

   integer :: ndigit

   ndigit = int(log10(real(x))) + 1
   if(ndigit<10) then
      write(digit_format,"('(i',i1,')')") ndigit
   else
      write(digit_format,"('(i',i2,')')") ndigit
   endif
end function digit_format

!****************************************

! Convert a matrix held in lower compressed sparse column format to standard
! HSL format in-place. Duplicates are summed and out-of-range entries are
! remove.  For symmetric, skew-symmetric and Hermitian matrices, entries in 
! the upper triangle are removed and discarded.
!
! Note: this routine is in place! cscl_convert is the out of place version.
subroutine mc69_cscl_clean_single(matrix_type, m, n, ptr, row, flag, &
      val, lmap, map, lp, noor, ndup) 
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: m, n ! matrix dimensions
   integer, dimension(n+1), intent(inout) :: ptr ! column pointers
   integer, dimension(*), intent(inout) :: row ! row indices
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(inout) :: val ! values
   integer, optional, intent(out) :: lmap ! size of map
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives source: map(i) = j means val_out(i)=val(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed


   ! Local variables
   integer :: col ! current column
   character(50)  :: context  ! Procedure name (used when printing).
   integer :: i
   integer :: idiag ! number of columns with a diagonal entry
   integer :: idup ! number of duplicates summed
   integer :: ioor ! number of out-of-range entries
   integer :: j
   integer :: k
   integer :: k1
   integer :: k2
   integer :: l
   integer :: nout ! output unit (set to -1 if nout not present)
   integer :: st ! stat parameter
   integer, allocatable :: temp(:) ! work array

   context = 'mc69_cscl_clean'

   flag = MC69_SUCCESS

   nout = -1
   if(present(lp)) nout = lp

   ! ---------------------------------------------
   ! Check that restrictions are adhered to
   ! ---------------------------------------------

   ! Note: have to change this test for complex code
   if(matrix_type < 0 .or. matrix_type == 5 .or. matrix_type > 6) then
      flag = MC69_ERROR_MATRIX_TYPE
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(m < 0 .or. n < 0) then
      flag = MC69_ERROR_N_OOR
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(abs(matrix_type) > 1 .and. m /= n) then
      flag = MC69_ERROR_M_NE_N
      call mc69_print_flag(context,nout,flag)
      return
   end if
 
   if(ptr(1) < 1) then
      flag = MC69_ERROR_PTR_1
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(present(map).neqv.present(lmap)) then
      flag = MC69_ERROR_LMAP_MISS
      call mc69_print_flag(context,nout,flag)
      return
   end if

   ! ---------------------------------------------

   ! allocate work array

   allocate(temp(m),stat=st)
   if(st /= 0) then
      flag = MC69_ERROR_ALLOCATION
      call mc69_print_flag(context,lp,flag)
      return
   endif
   temp(:) = 0

   ! count duplicates and out-of-range indices
   idup = 0; ioor = 0; idiag = 0
   k = 0 ! last diagonal seen
   l = 1
   do col = 1, n
      if(ptr(col+1).lt.ptr(col)) then
         flag = MC69_ERROR_PTR_MONO
         call mc69_print_flag(context,nout,flag)
         return
      endif
      if(abs(matrix_type) > 2) l = col
      do i = ptr(col), ptr(col+1)-1
         j = row(i)
         if(j<l .or. j>m) then
            ioor = ioor + 1
            row(i) = m+1
         else
            if(temp(j)==col) idup = idup + 1
            temp(j) = col
         end if
         if(j.eq.col .and. k<col) then
            idiag = idiag + 1
            k = col
         endif
      end do
   end do
   if(present(ndup)) ndup = idup
   if(present(noor)) noor = ioor
   deallocate(temp,stat=st)
   

   k = 0
   l = ptr(n+1) 
   if(present(map)) then
      deallocate(map,stat=st)
      lmap = l - 1 - ioor + idup
      allocate(map(l-1+idup*2),stat=st)
      if(st /= 0) then
         flag = MC69_ERROR_ALLOCATION
         call mc69_print_flag(context,lp,flag)
         return
      endif
      do i = 1,ptr(n+1)-1
         map(i) = i
      end do
      if(present(val)) then
         do col = 1, n
            k1 = ptr(col)
            ptr(col) = k + 1
            k2 = ptr(col+1)-1
            if(k2-k1<0) cycle
            call sort (row(k1), k2-k1+1, map=map(k1:), val=val(k1))
            ! Squeeze out duplicates and out-of-range indices
            if(row(k1) /= m+1) then
               k = k + 1      
               row(k) = row(k1)
               map(k) = map(k1)
               val(k) = val(k1)
            end if
            do i = k1+1,k2
               if(row(i) == m+1) then
                  exit             
               else if(row(i)>row(i-1)) then
                  k = k + 1
                  row(k) = row(i)
                  map(k) = map(i)
                  val(k) = val(i)
               else
                  map(l) = k
                  map(l+1) = map(i)
                  val(k) = val(k) + val(i)
                  l = l + 2
               end if
            end do
         end do
      else 
         do col = 1, n
            k1 = ptr(col)
            ptr(col) = k + 1
            k2 = ptr(col+1)-1
            if(k2-k1<0) cycle
            call sort (row(k1), k2-k1+1, map=map(k1:))
            ! Squeeze out duplicates and out-of-range indices
            if(row(k1) /= m+1) then
               k = k + 1      
               row(k) = row(k1)
               map(k) = map(k1)
            end if
            do i = k1+1,k2
               if(row(i) == m+1) then
                  exit             
               else if(row(i)>row(i-1)) then
                  k = k + 1
                  row(k) = row(i)
                  map(k) = map(i)
               else
                  map(l) = k
                  map(l+1) = map(i)
                  l = l + 2
               end if
            end do
         end do
      end if
      l = ptr(n+1) - 1
      ptr(n+1) = k + 1
      ! Move duplicate pair in map forward
      do i = 1, idup*2
         map(k+i) = map(l+i)
      end do
   else if(present(val)) then
      do col = 1, n
         k1 = ptr(col)
         ptr(col) = k + 1
         k2 = ptr(col+1)-1
         if(k2-k1<0) cycle
         call sort (row(k1), k2-k1+1, val=val(k1))
         ! Squeeze out duplicates and out-of-range indices
         if(row(k1) /= m+1) then
            k = k + 1      
            row(k) = row(k1)
            val(k) = val(k1)
         end if
         do i = k1+1,k2
            if(row(i) == m+1) then
               exit             
            else if(row(i)>row(i-1)) then
               k = k + 1
               row(k) = row(i)
               val(k) = val(i)
            else
               val(k) = val(k) + val(i)
               l = l + 2
            end if
         end do
      end do
      ptr(n+1) = k + 1
   else 
      do col = 1, n
         k1 = ptr(col)
         ptr(col) = k + 1
         k2 = ptr(col+1)-1
         if(k2-k1<0) cycle
         call sort (row(k1), k2-k1+1)
         ! Squeeze out duplicates and out-of-range indices
         if(row(k1) /= m+1) then
            k = k + 1      
            row(k) = row(k1)
         end if
         do i = k1+1,k2
            if(row(i) == m+1) then
               exit             
            else if(row(i)>row(i-1)) then
               k = k + 1
               row(k) = row(i)
            end if
         end do
      end do
      ptr(n+1) = k + 1
   end if

   select case(matrix_type)
   case(HSL_MATRIX_REAL_SYM_PSDEF)
      ! Check for positive diagonal entries
      do j = 1,n
         k = ptr(j)
         ! Note: we're OK if the column is empty - row(k) is still a diagonal
         ! entry, however we must be careful that we've not gone past the
         ! end of the matrix
         if(k.ge.ptr(n+1)) exit
         if(row(k)/=j)then
            flag = MC69_ERROR_MISSING_DIAGONAL
            call mc69_print_flag(context,nout,flag)
            return
         end if
         if(present(val))then
            if(val(k)<=zero)then
               flag = MC69_ERROR_MISSING_DIAGONAL
               call mc69_print_flag(context,nout,flag)
               return
            end if 
         end if 
      end do
   end select

   ! Check whether a warning needs to be raised
   if(ioor > 0 .or. idup > 0 .or. idiag < n) then
     if(ioor > 0) flag = MC69_WARNING_IDX_OOR
     if(idup > 0) flag = MC69_WARNING_DUP_IDX
     if(idup > 0 .and. ioor > 0) flag = MC69_WARNING_DUP_AND_OOR
     if(abs(matrix_type) .ne. HSL_MATRIX_REAL_SKEW) then
         if(idiag < n .and. ioor > 0) then
            flag = MC69_WARNING_MISS_DIAG_OORDUP
         else if(idiag < n .and. idup > 0) then
            flag = MC69_WARNING_MISS_DIAG_OORDUP
         else if(idiag < n) then
            flag = MC69_WARNING_MISSING_DIAGONAL
         end if
      endif
      call mc69_print_flag(context,nout,flag)
   end if

end subroutine mc69_cscl_clean_single

!****************************************

!
! Converts CSC with lower entries only for symmetric, skew-symmetric and
! Hermitian matrices to HSL standard format
!
subroutine mc69_cscl_convert_single(matrix_type, m, n, ptr_in, row_in, ptr_out,&
      row_out, flag, val_in, val_out, lmap, map, lp, noor, ndup) 
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: m ! number of rows
   integer, intent(in) :: n ! number of columns
   integer, dimension(*), intent(in) :: ptr_in ! column pointers on input
   integer, dimension(*), intent(in) :: row_in ! row indices on input.
      ! These may be unordered within each column and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(*), intent(out) :: ptr_out ! col ptr output
   integer, allocatable, dimension(:), intent(out) :: row_out ! row indices out
      ! Duplicates and out-of-range entries are dealt with and
      ! the entries within each column are ordered by increasing row index.
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(in) :: val_in ! values input
   real(wp), optional, allocatable, dimension(:) :: val_out
      ! values on output
   integer, optional, intent(out) :: lmap
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives src: map(i) = j means val_out(i)=val_in(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val_in(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed

   ! Local variables
   character(50)  :: context  ! Procedure name (used when printing).
   integer :: nout ! output unit (set to -1 if nout not present)

   context = 'mc69_cscl_convert'

   nout = -1
   if(present(lp)) nout = lp

   ! Note: have to change this test for complex code
   if(matrix_type < 0 .or. matrix_type == 5 .or. matrix_type > 6) then
      flag = MC69_ERROR_MATRIX_TYPE
      call mc69_print_flag(context,nout,flag)
      return
   end if

   call mc69_cscl_convert_main(context, 1, matrix_type, m, n, &
      ptr_in, row_in, ptr_out, row_out, flag, val_in, val_out, lmap, map, &
      lp, noor, ndup) 
end subroutine mc69_cscl_convert_single

!****************************************

!
! Converts CSC (with lower entries only for symmetric, skew-symmetric and
! Hermitian matrices) to HSL standard format.
! Also used for symmetric, skew-symmetric and Hermitian matrices in upper
! CSR format
!
subroutine mc69_cscl_convert_main(context, multiplier, matrix_type, m, n, &
      ptr_in, row_in, ptr_out, row_out, flag, val_in, val_out, lmap, map, &
      lp, noor, ndup) 
   character(50), intent(in) :: context ! Procedure name (used when printing).
   integer, intent(in) :: multiplier ! -1 or 1, differs for csc/csr
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: m ! number of rows
   integer, intent(in) :: n ! number of columns
   integer, dimension(*), intent(in) :: ptr_in ! column pointers on input
   integer, dimension(*), intent(in) :: row_in ! row indices on input.
      ! These may be unordered within each column and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(*), intent(out) :: ptr_out ! col ptr output
   integer, allocatable, dimension(:), intent(out) :: row_out ! row indices out
      ! Duplicates and out-of-range entries are dealt with and
      ! the entries within each column are ordered by increasing row index.
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(in) :: val_in ! values input
   real(wp), optional, allocatable, dimension(:) :: val_out
      ! values on output
   integer, optional, intent(out) :: lmap
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives src: map(i) = j means val_out(i)=val_in(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val_in(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed


   ! Local variables
   integer :: col ! current column
   integer :: i
   integer :: idiag
   integer :: idup
   integer :: ioor
   integer :: j
   integer :: k
   integer :: nout ! output unit (set to -1 if nout not present)
   integer :: st ! stat parameter
   integer :: minidx

   type(dup_list), pointer :: dup
   type(dup_list), pointer :: duphead

   ! ---------------------------------------------
   ! Check that restrictions are adhered to
   ! ---------------------------------------------

   nullify(dup, duphead)

   flag = MC69_SUCCESS

   nout = -1
   if(present(lp)) nout = lp

   if(n < 0 .or. m < 0) then
      flag = MC69_ERROR_N_OOR
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(ptr_in(1) < 1) then
      flag = MC69_ERROR_PTR_1
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(present(val_in).neqv.present(val_out)) then
      flag = MC69_ERROR_VAL_MISS
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(present(map).neqv.present(lmap)) then
      flag = MC69_ERROR_LMAP_MISS
      call mc69_print_flag(context,nout,flag)
      return
   end if

   ! ---------------------------------------------

   ! ensure output arrays are not allocated

   deallocate(row_out,stat=st)
   if(present(val_out)) deallocate(val_out,stat=st)
   if(present(map)) deallocate(map,stat=st)

   idup = 0; ioor = 0; idiag = 0

   allocate(row_out(ptr_in(n+1)-1),stat=st)
   if(st.ne.0) goto 100
   if(present(map)) then
      ! Allocate map for worst case where all bar one are duplicates
      allocate(map(2*ptr_in(n+1)-2),stat=st)
      k = 1 ! insert location
      do col = 1, n
         ptr_out(col) = k
         if(ptr_in(col+1).lt.ptr_in(col)) then
            flag = MC69_ERROR_PTR_MONO
            call mc69_print_flag(context,nout,flag)
            call cleanup_dup(duphead)
            return
         endif
         minidx = 1
         if(abs(matrix_type).ge.HSL_MATRIX_REAL_SYM_PSDEF) minidx = col
         if(abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW) minidx = col + 1
         ! Loop over column, copy across while dropping any out of range entries
         do i = ptr_in(col), ptr_in(col+1)-1
            j = row_in(i)
            if(j.lt.minidx .or. j.gt.m) then
               ! out of range, ignore
               ioor = ioor + 1
               cycle
            endif
            row_out(k) = row_in(i)
            map(k) = multiplier*i
            k = k + 1
         end do
         ! Sort entries into order
         i = k - ptr_out(col)
         if(i.eq.0 .and. ptr_in(col+1)-ptr_in(col).ne.0) then
            flag = MC69_ERROR_ALL_OOR
            call mc69_print_flag(context,nout,flag)
            call cleanup_dup(duphead)
            return
         endif
         if(i.ne.0) then
            call sort(row_out(ptr_out(col):k-1), i, map=map(ptr_out(col):k-1))
            ! Loop over sorted list and drop duplicates
            i = k-1 ! last entry in column
            k = ptr_out(col)+1 ! insert position
            ! Note: we are skipping the first entry as it cannot be a duplicate
            if(row_out(ptr_out(col)).eq.col) then
               idiag = idiag + 1
            elseif(abs(matrix_type).eq.HSL_MATRIX_REAL_SYM_PSDEF) then
               flag = MC69_ERROR_MISSING_DIAGONAL
               call mc69_print_flag(context,nout,flag)
               call cleanup_dup(duphead)
               return
            endif
            do i = ptr_out(col)+1, i
               if(row_out(i).eq.row_out(i-1)) then
                  ! duplicate, drop
                  idup = idup + 1
                  allocate(dup,stat=st)
                  if(st.ne.0) goto 100
                  dup%next => duphead
                  duphead => dup
                  dup%src = map(i)
                  dup%dest = k-1
                  cycle
               endif
               if(row_out(i).eq.col) idiag = idiag + 1
               row_out(k) = row_out(i)
               map(k) = map(i)
               k = k + 1
            end do
         elseif(abs(matrix_type).eq.HSL_MATRIX_REAL_SYM_PSDEF) then
            flag = MC69_ERROR_MISSING_DIAGONAL
            call mc69_print_flag(context,nout,flag)
            call cleanup_dup(duphead)
            return
         endif
      end do
      ptr_out(n+1) = k
      lmap = k-1
   elseif(present(val_out)) then
      allocate(val_out(ptr_in(n+1)-1),stat=st)
      if(st.ne.0) goto 100
      k = 1 ! insert location
      do col = 1, n
         ptr_out(col) = k
         if(ptr_in(col+1).lt.ptr_in(col)) then
            flag = MC69_ERROR_PTR_MONO
            call mc69_print_flag(context,nout,flag)
            return
         endif
         minidx = 1
         if(abs(matrix_type).ge.HSL_MATRIX_REAL_SYM_PSDEF) minidx = col
         if(abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW) minidx = col + 1
         ! Loop over column, copy across while dropping any out of range entries
         select case(matrix_type)
         case(HSL_MATRIX_REAL_SKEW)
            do i = ptr_in(col), ptr_in(col+1)-1
               j = row_in(i)
               if(j.lt.minidx .or. j.gt.m) then
                  ! out of range, ignore
                  ioor = ioor + 1
                  cycle
               endif
               row_out(k) = row_in(i)
               val_out(k) = multiplier*val_in(i)
               k = k + 1
            end do
         case default
            do i = ptr_in(col), ptr_in(col+1)-1
               j = row_in(i)
               if(j.lt.minidx .or. j.gt.m) then
                  ! out of range, ignore
                  ioor = ioor + 1
                  cycle
               endif
               row_out(k) = row_in(i)
               val_out(k) = val_in(i)
               k = k + 1
            end do
         end select
         ! Sort entries into order
         i = k - ptr_out(col)
         if(i.eq.0 .and. ptr_in(col+1)-ptr_in(col).ne.0) then
            flag = MC69_ERROR_ALL_OOR
            call mc69_print_flag(context,nout,flag)
            return
         endif
         if(i.ne.0) then
            call sort(row_out(ptr_out(col):k-1), i, &
               val=val_out(ptr_out(col):k-1))
            ! Loop over sorted list and drop duplicates
            i = k-1 ! last entry in column
            k = ptr_out(col)+1 ! insert position
            ! Note: we are skipping the first entry as it cannot be a duplicate
            if(row_out(ptr_out(col)).eq.col) then
               idiag = idiag + 1
            elseif(abs(matrix_type).eq.HSL_MATRIX_REAL_SYM_PSDEF) then
               flag = MC69_ERROR_MISSING_DIAGONAL
               call mc69_print_flag(context,nout,flag)
               return
            endif
            do i = ptr_out(col)+1, i
               if(row_out(i).eq.row_out(i-1)) then
                  ! duplicate, sum then drop from pattern
                  idup = idup + 1
                  val_out(i-1) = val_out(i-1) + val_out(i)
                  cycle
               endif
               if(row_out(i).eq.col) idiag = idiag + 1
               row_out(k) = row_out(i)
               val_out(k) = val_out(i)
               k = k + 1
            end do
         elseif(abs(matrix_type).eq.HSL_MATRIX_REAL_SYM_PSDEF) then
            flag = MC69_ERROR_MISSING_DIAGONAL
            call mc69_print_flag(context,nout,flag)
            return
         endif
      end do
      ptr_out(n+1) = k
   else ! pattern only
      k = 1 ! insert location
      do col = 1, n
         ptr_out(col) = k
         if(ptr_in(col+1).lt.ptr_in(col)) then
            flag = MC69_ERROR_PTR_MONO
            call mc69_print_flag(context,nout,flag)
            return
         endif
         minidx = 1
         if(abs(matrix_type).ge.HSL_MATRIX_REAL_SYM_PSDEF) minidx = col
         if(abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW) minidx = col + 1
         ! Loop over column, copy across while dropping any out of range entries
         do i = ptr_in(col), ptr_in(col+1)-1
            j = row_in(i)
            if(j.lt.minidx .or. j.gt.m) then
               ! out of range, ignore
               ioor = ioor + 1
               cycle
            endif
            row_out(k) = row_in(i)
            k = k + 1
         end do
         ! Sort entries into order
         i = k - ptr_out(col)
         if(i.eq.0 .and. ptr_in(col+1)-ptr_in(col).ne.0) then
            flag = MC69_ERROR_ALL_OOR
            call mc69_print_flag(context,nout,flag)
            return
         endif
         if(i.ne.0) then
            call sort(row_out(ptr_out(col):k-1), i)
            ! Loop over sorted list and drop duplicates
            i = k-1 ! last entry in column
            k = ptr_out(col)+1 ! insert position
            ! Note: we are skipping the first entry as it cannot be a duplicate
            if(row_out(ptr_out(col)).eq.col) then
               idiag = idiag + 1
            elseif(abs(matrix_type).eq.HSL_MATRIX_REAL_SYM_PSDEF) then
               flag = MC69_ERROR_MISSING_DIAGONAL
               call mc69_print_flag(context,nout,flag)
               return
            endif
            do i = ptr_out(col)+1, i
               if(row_out(i).eq.row_out(i-1)) then
                  ! duplicate, drop
                  idup = idup + 1
                  cycle
               endif
               if(row_out(i).eq.col) idiag = idiag + 1
               row_out(k) = row_out(i)
               k = k + 1
            end do
         elseif(abs(matrix_type).eq.HSL_MATRIX_REAL_SYM_PSDEF) then
            flag = MC69_ERROR_MISSING_DIAGONAL
            call mc69_print_flag(context,nout,flag)
            return
         end if
      end do
      ptr_out(n+1) = k
   endif

   if(present(map)) then
      ! Append duplicates to map
      do while(associated(duphead))
         idup = idup + 1
         map(lmap+1) = duphead%dest
         map(lmap+2) = duphead%src
         lmap = lmap + 2
         dup => duphead%next
         deallocate(duphead)
         duphead => dup
      end do
      if(present(val_out)) then
         allocate(val_out(ptr_out(n+1)-1), stat=st)
         if(st.ne.0) goto 100
         call mc69_set_values(matrix_type, lmap, map, val_in, ptr_out(n+1)-1, &
            val_out)
      endif
   endif

   if(present(val_out)) then
      select case(matrix_type)
      case(HSL_MATRIX_REAL_SYM_PSDEF)
         ! Check for positive diagonal entries
         do j = 1,n
            k = ptr_out(j)
            ! ps def case - can't reach here unless all entries have diagonal
            if(real(val_out(k))<=zero) then
               flag = MC69_ERROR_MISSING_DIAGONAL
               call mc69_print_flag(context,nout,flag)
               return
            end if 
         end do
      end select
   endif

   ! Check whether a warning needs to be raised
   if(ioor > 0 .or. idup > 0 .or. idiag < n) then
     if(ioor > 0) flag = MC69_WARNING_IDX_OOR
     if(idup > 0) flag = MC69_WARNING_DUP_IDX
     if(idup > 0 .and. ioor > 0) flag = MC69_WARNING_DUP_AND_OOR
     if(abs(matrix_type) .ne. HSL_MATRIX_REAL_SKEW) then
         if(idiag < n .and. ioor > 0) then
            flag = MC69_WARNING_MISS_DIAG_OORDUP
         else if(idiag < n .and. idup > 0) then
            flag = MC69_WARNING_MISS_DIAG_OORDUP
         else if(idiag < n) then
            flag = MC69_WARNING_MISSING_DIAGONAL
         end if
      endif
      call mc69_print_flag(context,nout,flag)
   end if

   if(present(noor)) noor = ioor
   if(present(ndup)) ndup = idup

   return

100 if(st /= 0) then
      flag = MC69_ERROR_ALLOCATION
      call mc69_print_flag(context,nout,flag)
    end if

end subroutine mc69_cscl_convert_main

!****************************************

!
! Converts CSC with only upper entries present to CSC with 
! only lower entries present and entries
! within each column ordered by increasing row index (HSL standard format)
!
subroutine mc69_cscu_convert_single(matrix_type, n, ptr_in, row_in, ptr_out, &
      row_out, flag, val_in, val_out, lmap, map, lp, noor, ndup) 
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: n ! matrix dimension
   integer, dimension(*), intent(in) :: ptr_in ! row pointers on input
   integer, dimension(*), intent(in) :: row_in ! col indices on input.
      ! These may be unordered within each row and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(*), intent(out) :: ptr_out ! col ptr output
   integer, allocatable, dimension(:), intent(out) :: row_out ! row indices out
      ! Duplicates and out-of-range entries are dealt with and
      ! the entries within each column are ordered by increasing row index.
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(in) :: val_in ! values input
   real(wp), optional, allocatable, dimension(:) :: val_out
      ! values on output
   integer, optional, intent(out) :: lmap
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives src: map(i) = j means val_out(i)=val_in(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val_in(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed

   ! Local variables
   character(50)  :: context  ! Procedure name (used when printing).
   integer :: nout ! output unit (set to -1 if nout not present)

   context = 'mc69_cscu_convert'

   nout = -1
   if(present(lp)) nout = lp

   ! Note: have to change this test for complex code
   if(matrix_type < 3 .or. matrix_type == 5 .or. matrix_type > 6) then
      flag = MC69_ERROR_MATRIX_TYPE
      call mc69_print_flag(context,nout,flag)
      return
   end if

   call mc69_csrl_convert_main(context, -1, matrix_type, n, n, &
      ptr_in, row_in, ptr_out, row_out, flag, val_in, val_out, lmap, map, &
      lp, noor, ndup) 
end subroutine mc69_cscu_convert_single

!****************************************

!
! Converts CSC with both lower and upper entries present to CSC with 
! only lower entries present and entries
! within each column ordered by increasing row index (HSL standard format)
!
subroutine mc69_csclu_convert_single(matrix_type, n, ptr_in, row_in, ptr_out, &
      row_out, flag, val_in, val_out, lmap, map, lp, noor, ndup) 
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: n ! matrix dimension
   integer, dimension(*), intent(in) :: ptr_in ! column pointers on input
   integer, dimension(*), intent(in) :: row_in ! row indices on input.
      ! These may be unordered within each column and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(*), intent(out) :: ptr_out ! col ptr output
   integer, allocatable, dimension(:), intent(out) :: row_out ! row indices out
      ! Duplicates and out-of-range entries are dealt with and
      ! the entries within each column are ordered by increasing row index.
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(in) :: val_in ! values input
   real(wp), optional, allocatable, dimension(:) :: val_out
      ! values on output
   integer, optional, intent(out) :: lmap
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives src: map(i) = j means val_out(i)=val_in(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val_in(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed

   ! Local variables
   character(50)  :: context  ! Procedure name (used when printing).
   integer :: nout ! output unit (set to -1 if nout not present)

   context = 'mc69_csclu_convert'

   nout = -1
   if(present(lp)) nout = lp

   ! Note: have to change this test for complex code
   if(matrix_type < 3 .or. matrix_type == 5 .or. matrix_type > 6) then
      flag = MC69_ERROR_MATRIX_TYPE
      call mc69_print_flag(context,nout,flag)
      return
   end if

   call mc69_csclu_convert_main(context, -1, matrix_type, n, &
      ptr_in, row_in, ptr_out, row_out, flag, val_in, val_out, lmap, map, &
      lp, noor, ndup) 
end subroutine mc69_csclu_convert_single

!****************************************

subroutine mc69_csclu_convert_main(context, multiplier, matrix_type, n, &
      ptr_in, row_in, ptr_out, row_out, flag, val_in, val_out, lmap, map, &
      lp, noor, ndup) 
   character(50), intent(in) :: context  ! Procedure name (used when printing).
   integer, intent(in) :: multiplier ! -1 for upr source, +1 for lwr
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: n ! matrix dimension
   integer, dimension(*), intent(in) :: ptr_in ! column pointers on input
   integer, dimension(*), intent(in) :: row_in ! row indices on input.
      ! These may be unordered within each column and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(*), intent(out) :: ptr_out ! col ptr output
   integer, allocatable, dimension(:), intent(out) :: row_out ! row indices out
      ! Duplicates and out-of-range entries are dealt with and
      ! the entries within each column are ordered by increasing row index.
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(in) :: val_in ! values input
   real(wp), optional, allocatable, dimension(:) :: val_out
      ! values on output
   integer, optional, intent(out) :: lmap
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives src: map(i) = j means val_out(i)=val_in(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val_in(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed


   ! Local variables
   integer :: col ! current column
   integer :: i
   integer :: idiag
   integer :: idup
   integer :: ioor
   integer :: j
   integer :: k
   integer :: nlwr ! number of strictly lower triangular entries
   integer :: npre
   integer :: nupr ! number of strictly upper triangular entries
   integer :: ne    
   integer :: nout ! output unit (set to -1 if nout not present)
   integer :: st ! stat parameter

   type(dup_list), pointer :: dup
   type(dup_list), pointer :: duphead

   nullify(dup, duphead)

   ! ---------------------------------------------
   ! Check that restrictions are adhered to
   ! ---------------------------------------------

   flag = MC69_SUCCESS

   nout = -1
   if(present(lp)) nout = lp

   if(n < 0) then
      flag = MC69_ERROR_N_OOR
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(ptr_in(1) < 1) then
      flag = MC69_ERROR_PTR_1
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(present(val_in).neqv.present(val_out)) then
      flag = MC69_ERROR_VAL_MISS
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(present(map).neqv.present(lmap)) then
      flag = MC69_ERROR_LMAP_MISS
      call mc69_print_flag(context,nout,flag)
      return
   end if

   ! ---------------------------------------------

   ! allocate output and work arrays

   deallocate(row_out,stat=st)
   if(present(val_out)) deallocate(val_out,stat=st)
   if(present(map)) deallocate(map,stat=st)


   ! ---------------------------------------------
   ! check for duplicates and/or out-of-range. check diagonal present.
   ! first do the case where values not present

   idup = 0; ioor = 0; idiag = 0
   nlwr = 0; nupr = 0

   !
   ! First pass, count number of entries in each row of the matrix.
   ! Count is at an offset of 1 to allow us to play tricks
   ! (ie. ptr_out(i+1) set to hold number of entries in column i
   ! of expanded matrix).
   ! Excludes out of range entries. Includes duplicates.
   !
   ptr_out(1:n+1) = 0
   k = 0 ! last diagonal found
   do col = 1, n
      if(ptr_in(col+1).lt.ptr_in(col)) then
         flag = MC69_ERROR_PTR_MONO
         call mc69_print_flag(context,nout,flag)
         return
      endif
      npre = nlwr + nupr + idiag
      if(ptr_in(col+1).eq.ptr_in(col)) cycle ! no entries in column
      do i = ptr_in(col), ptr_in(col+1)-1
         j = row_in(i)
         if(j.lt.1 .or. j.gt.n) then
            ioor = ioor + 1
            cycle
         endif
         if(j.eq.col .and. abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW) then
            ioor = ioor + 1
            cycle
         endif
         if(j.gt.col) then
            nlwr = nlwr + 1
            cycle
         endif
         if(j.lt.col) then
            nupr = nupr + 1
         elseif(j.ne.k) then ! diagonal entry (not second in column)
            idiag = idiag + 1
            k = col
         endif
         ptr_out(j+1) = ptr_out(j+1) + 1
      end do
      if(nlwr+nupr+idiag.eq.npre) then
         ! Column contains only out of range entries
         flag = MC69_ERROR_ALL_OOR
         call mc69_print_flag(context,nout,flag)
         return
      endif
   end do

   ! Check for missing diagonals in pos def case
   ! Note: change this test for complex case
   if(abs(matrix_type) == HSL_MATRIX_REAL_SYM_PSDEF) then
     if(idiag < n) then
         flag = MC69_ERROR_MISSING_DIAGONAL
         call mc69_print_flag(context,nout,flag)
         return
      end if
   end if

   ! Check number in lower triangle = number in upper triangle
   if(nlwr .ne. nupr) then
      flag = MC69_ERROR_MISMATCH_LWRUPR
      call mc69_print_flag(context,nout,flag)
      return
   endif

   ! Determine column starts for transposed matrix such 
   ! that column i starts at ptr_out(i+1)
   ne = 0
   ptr_out(1) = 1
   do i = 1, n
      ne = ne + ptr_out(i+1)
      ptr_out(i+1) = ptr_out(i) + ptr_out(i+1)
   end do
   do i = n,1,-1
      ptr_out(i+1) = ptr_out(i)
   end do

   !
   ! Second pass, drop entries into place for conjugate of transposed
   ! matrix
   !
   allocate(row_out(ne), stat=st)
   if(st.ne.0) goto 100
   if(present(map)) then
      ! Needs to be this big for worst case: all entries are repeat of same
      allocate(map(2*(ptr_in(n+1)-1)), stat=st)
      if(st.ne.0) goto 100
      do col = 1, n
         do i = ptr_in(col), ptr_in(col+1)-1
            j = row_in(i)
            if(j.lt.1 .or. j.gt.col) cycle ! ignore oor and lwr triangle entries
            if(j.eq.col .and. abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW) cycle
            k = ptr_out(j+1)
            ptr_out(j+1) = k + 1
            row_out(k) = col
            map(k) = multiplier*i
         end do
      end do
   elseif(present(val_out)) then
      allocate(val_out(2*(ptr_in(n+1)-1)), stat=st)
      if(st.ne.0) goto 100
      select case(matrix_type)
      case(HSL_MATRIX_REAL_SKEW)
         do col = 1, n
            do i = ptr_in(col), ptr_in(col+1)-1
               j = row_in(i)
               if(j.lt.1 .or. j.ge.col) cycle
               k = ptr_out(j+1)
               ptr_out(j+1) = k + 1
               row_out(k) = col
               val_out(k) = multiplier*val_in(i)
            end do
         end do
      case(HSL_MATRIX_REAL_SYM_PSDEF, HSL_MATRIX_REAL_SYM_INDEF)
         do col = 1, n
            do i = ptr_in(col), ptr_in(col+1)-1
               j = row_in(i)
               if(j.lt.1 .or. j.gt.col) cycle
               k = ptr_out(j+1)
               ptr_out(j+1) = k + 1
               row_out(k) = col
               val_out(k) = val_in(i)
            end do
         end do
      end select
   else ! neither val_out nor map present
      do col = 1, n
         do i = ptr_in(col), ptr_in(col+1)-1
            j = row_in(i)
            if(j.lt.1 .or. j.gt.col) cycle
            if(j.eq.col .and. abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW) cycle
            k = ptr_out(j+1)
            ptr_out(j+1) = k + 1
            row_out(k) = col
         end do
      end do
   endif

   !
   ! Third pass, removal of duplicates (no sort required by construction)
   !
   if(present(map)) then
      k = 1 ! insert position
      do col = 1, n
         i = ptr_out(col)
         ptr_out(col) = k
         if(ptr_out(col+1).eq.i) cycle ! no entries
         ! Move first entry of column forward
         row_out(k) = row_out(i)
         map(k) = map(i)
         k = k + 1
         ! Loop over remaining entries
         do i = i+1, ptr_out(col+1)-1
            if(row_out(i).eq.row_out(k-1)) then
               ! Duplicate
               allocate(dup,stat=st)
               if(st.ne.0) goto 100
               dup%next => duphead
               duphead => dup
               dup%src = map(i)
               dup%dest = k-1
               cycle
            endif
            ! Pull entry forwards
            row_out(k) = row_out(i)
            map(k) = map(i)
            k = k + 1
         end do
      end do
      ptr_out(n+1) = k
      lmap = k-1
   elseif(present(val_out)) then
      k = 1 ! insert position
      do col = 1, n
         i = ptr_out(col)
         ptr_out(col) = k
         if(ptr_out(col+1).eq.i) cycle ! no entries
         ! Move first entry of column forward
         row_out(k) = row_out(i)
         val_out(k) = val_out(i)
         k = k + 1
         ! Loop over remaining entries
         do i = i+1, ptr_out(col+1)-1
            if(row_out(i).eq.row_out(k-1)) then
               ! Duplicate
               idup = idup + 1
               val_out(k-1) = val_out(k-1) + val_out(i)
               cycle
            endif
            ! Pull entry forwards
            row_out(k) = row_out(i)
            val_out(k) = val_out(i)
            k = k + 1
         end do
      end do
      ptr_out(n+1) = k
   else ! neither val_out nor map are present
      k = 1 ! insert position
      do col = 1, n
         i = ptr_out(col)
         ptr_out(col) = k
         if(ptr_out(col+1).eq.i) cycle ! no entries
         ! Move first entry of column forward
         row_out(k) = row_out(i)
         k = k + 1
         ! Loop over remaining entries
         do i = i+1, ptr_out(col+1)-1
            if(row_out(i).eq.row_out(k-1)) then
               ! Duplicate
               idup = idup + 1
               cycle
            endif
            ! Pull entry forwards
            row_out(k) = row_out(i)
            k = k + 1
         end do
      end do
      ptr_out(n+1) = k
   endif

   if(present(map)) then
      ! Append duplicates to map
      do while(associated(duphead))
         idup = idup + 1
         map(lmap+1) = duphead%dest
         map(lmap+2) = duphead%src
         lmap = lmap + 2
         dup => duphead%next
         deallocate(duphead)
         duphead => dup
      end do
      if(present(val_out)) then
         allocate(val_out(ptr_out(n+1)-1), stat=st)
         if(st.ne.0) goto 100
         call mc69_set_values(matrix_type, lmap, map, val_in, ptr_out(n+1)-1, &
            val_out)
      endif
   endif

   if(present(val_out)) then
      select case(matrix_type)
      case(HSL_MATRIX_REAL_SYM_PSDEF)
         ! Check for positive diagonal entries
         do j = 1,n
            k = ptr_out(j)
            ! ps def case - can't reach here unless all entries have diagonal
            if(real(val_out(k))<=zero) then
               flag = MC69_ERROR_MISSING_DIAGONAL
               call mc69_print_flag(context,nout,flag)
               return
            end if 
         end do
      end select
   endif

   ! Check whether a warning needs to be raised
   if(ioor > 0 .or. idup > 0 .or. idiag < n) then
     if(ioor > 0) flag = MC69_WARNING_IDX_OOR
     if(idup > 0) flag = MC69_WARNING_DUP_IDX
     if(idup > 0 .and. ioor > 0) flag = MC69_WARNING_DUP_AND_OOR
     if(abs(matrix_type) .ne. HSL_MATRIX_REAL_SKEW) then
         if(idiag < n .and. ioor > 0) then
            flag = MC69_WARNING_MISS_DIAG_OORDUP
         else if(idiag < n .and. idup > 0) then
            flag = MC69_WARNING_MISS_DIAG_OORDUP
         else if(idiag < n) then
            flag = MC69_WARNING_MISSING_DIAGONAL
         end if
      endif
      call mc69_print_flag(context,nout,flag)
   end if

   if(present(noor)) noor = ioor
   if(present(ndup)) ndup = idup

   return

100 if(st /= 0) then
      flag = MC69_ERROR_ALLOCATION
      call mc69_print_flag(context,nout,flag)
    end if

end subroutine mc69_csclu_convert_main

!****************************************

!
! Converts CSR with only lower entries present to CSC with 
! only lower entries present and entries
! within each column ordered by increasing row index (HSL standard format)
!
subroutine mc69_csrl_convert_single(matrix_type, m, n, ptr_in, col_in, ptr_out,&
      row_out, flag, val_in, val_out, lmap, map, lp, noor, ndup)
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: m ! matrix dimension
   integer, intent(in) :: n ! matrix dimension
   integer, dimension(*), intent(in) :: ptr_in ! row pointers on input
   integer, dimension(*), intent(in) :: col_in ! col indices on input.
      ! These may be unordered within each row and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(*), intent(out) :: ptr_out ! col ptr output
   integer, allocatable, dimension(:), intent(out) :: row_out ! row indices out
      ! Duplicates and out-of-range entries are dealt with and
      ! the entries within each column are ordered by increasing row index.
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(in) :: val_in ! values input
   real(wp), optional, allocatable, dimension(:) :: val_out
      ! values on output
   integer, optional, intent(out) :: lmap
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives src: map(i) = j means val_out(i)=val_in(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val_in(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed

   ! Local variables
   character(50)  :: context  ! Procedure name (used when printing).
   integer :: nout ! output unit (set to -1 if nout not present)

   context = 'mc69_csrl_convert'

   nout = -1
   if(present(lp)) nout = lp

   ! Note: have to change this test for complex code
   if(matrix_type < 1 .or. matrix_type == 5 .or. matrix_type > 6) then
      flag = MC69_ERROR_MATRIX_TYPE
      call mc69_print_flag(context,nout,flag)
      return
   end if

   call mc69_csrl_convert_main(context, 1, matrix_type, m, n, &
      ptr_in, col_in, ptr_out, row_out, flag, val_in, val_out, lmap, map, &
      lp, noor, ndup) 
end subroutine mc69_csrl_convert_single

!****************************************

subroutine mc69_csrl_convert_main(context, multiplier, matrix_type, m, n, &
      ptr_in, col_in, ptr_out, row_out, flag, val_in, val_out, lmap, map, &
      lp, noor, ndup) 
   character(50), intent(in) :: context  ! Procedure name (used when printing).
   integer, intent(in) :: multiplier ! 1 for lwr triangle source, -1 for upr
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: m ! matrix dimension
   integer, intent(in) :: n ! matrix dimension
   integer, dimension(*), intent(in) :: ptr_in ! row pointers on input
   integer, dimension(*), intent(in) :: col_in ! col indices on input.
      ! These may be unordered within each row and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(*), intent(out) :: ptr_out ! col ptr output
   integer, allocatable, dimension(:), intent(out) :: row_out ! row indices out
      ! Duplicates and out-of-range entries are dealt with and
      ! the entries within each column are ordered by increasing row index.
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(in) :: val_in ! values input
   real(wp), optional, allocatable, dimension(:) :: val_out
      ! values on output
   integer, optional, intent(out) :: lmap
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives src: map(i) = j means val_out(i)=val_in(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val_in(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed


   ! Local variables
   integer :: row ! current row
   integer :: col ! current col
   integer :: i
   integer :: idiag
   integer :: idup
   integer :: ioor
   integer :: j
   integer :: k
   integer :: ne    
   integer :: nout ! output unit (set to -1 if nout not present)
   integer :: npre ! number of out-of-range entries prior to this column
   integer :: st ! stat parameter
   integer :: maxv ! maximum value before out of range

   type(dup_list), pointer :: dup
   type(dup_list), pointer :: duphead

   nullify(dup, duphead)

   ! ---------------------------------------------
   ! Check that restrictions are adhered to
   ! ---------------------------------------------

   flag = MC69_SUCCESS

   nout = -1
   if(present(lp)) nout = lp

   if(n < 0 .or. m < 0) then
      flag = MC69_ERROR_N_OOR
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(ptr_in(1) < 1) then
      flag = MC69_ERROR_PTR_1
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(present(val_in).neqv.present(val_out)) then
      flag = MC69_ERROR_VAL_MISS
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(present(map).neqv.present(lmap)) then
      flag = MC69_ERROR_LMAP_MISS
      call mc69_print_flag(context,nout,flag)
      return
   end if

   ! ---------------------------------------------

   ! allocate output and work arrays

   deallocate(row_out,stat=st)
   if(present(val_out)) deallocate(val_out,stat=st)
   if(present(map)) deallocate(map,stat=st)


   ! ---------------------------------------------
   ! check for duplicates and/or out-of-range. check diagonal present.
   ! first do the case where values not present

   idup = 0; ioor = 0; idiag = 0

   !
   ! First pass, count number of entries in each col of the matrix.
   ! Count is at an offset of 1 to allow us to play tricks
   ! (ie. ptr_out(i+1) set to hold number of entries in column i
   ! of expanded matrix).
   ! Excludes out of range entries. Includes duplicates.
   !
   ptr_out(1:n+1) = 0
   k = 0 ! last diagonal found
   maxv = n
   do row = 1, m
      if(ptr_in(row+1).lt.ptr_in(row)) then
         flag = MC69_ERROR_PTR_MONO
         call mc69_print_flag(context,nout,flag)
         return
      endif
      npre = ioor
      if(abs(matrix_type).ge.HSL_MATRIX_REAL_SYM_PSDEF) maxv = row
      if(abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW) maxv = row-1
      do i = ptr_in(row), ptr_in(row+1)-1
         j = col_in(i)
         if(j.lt.1 .or. j.gt.maxv) then
            ioor = ioor + 1
            cycle
         endif
         if(j.eq.row .and. j.ne.k) then ! diagonal entry (not second in column)
            idiag = idiag + 1
            k = row
         endif
         ptr_out(j+1) = ptr_out(j+1) + 1
      end do
      if(ioor-npre.ne.0 .and. ptr_in(row+1)-ptr_in(row).eq.ioor-npre) then
         ! Column contains only out of range entries
         flag = MC69_ERROR_ALL_OOR
         call mc69_print_flag(context,nout,flag)
         return
      endif
   end do

   ! Check for missing diagonals in pos def case
   ! Note: change this test for complex case
   if(abs(matrix_type) == HSL_MATRIX_REAL_SYM_PSDEF) then
     if(idiag < n) then
         flag = MC69_ERROR_MISSING_DIAGONAL
         call mc69_print_flag(context,nout,flag)
         return
      end if
   end if

   ! Determine column starts for matrix such 
   ! that column i starts at ptr_out(i+1)
   ne = 0
   ptr_out(1) = 1
   do i = 1, n
      ne = ne + ptr_out(i+1)
      ptr_out(i+1) = ptr_out(i) + ptr_out(i+1)
   end do
   do i = n,1,-1
      ptr_out(i+1) = ptr_out(i)
   end do

   !
   ! Second pass, drop entries into place for matrix
   !
   allocate(row_out(ne), stat=st)
   if(st.ne.0) goto 100
   if(present(map)) then
      ! Needs to be this big for worst case: all entries are repeat of same
      allocate(map(2*(ptr_in(m+1)-1)), stat=st)
      if(st.ne.0) goto 100
      maxv = n
      do row = 1, m
         if(abs(matrix_type).ge.HSL_MATRIX_REAL_SYM_PSDEF) maxv = row
         if(abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW) maxv = row-1
         do i = ptr_in(row), ptr_in(row+1)-1
            j = col_in(i)
            if(j.lt.1 .or. j.gt.maxv) cycle ! ignore oor and upr triangle
            k = ptr_out(j+1)
            ptr_out(j+1) = k + 1
            row_out(k) = row
            map(k) = multiplier*i
         end do
      end do
   elseif(present(val_out)) then
      allocate(val_out(2*(ptr_in(m+1)-1)), stat=st)
      if(st.ne.0) goto 100
      select case(matrix_type)
      case(HSL_MATRIX_REAL_SKEW)
         do row = 1, m
            do i = ptr_in(row), ptr_in(row+1)-1
               j = col_in(i)
               if(j.lt.1 .or. j.ge.row) cycle
               k = ptr_out(j+1)
               ptr_out(j+1) = k + 1
               row_out(k) = row
               val_out(k) = multiplier*val_in(i)
            end do
         end do
      case default
         maxv= n
         do row = 1, m
            if(abs(matrix_type).ge.HSL_MATRIX_REAL_SYM_PSDEF) maxv = row
            if(abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW) maxv = row-1
            do i = ptr_in(row), ptr_in(row+1)-1
               j = col_in(i)
               if(j.lt.1 .or. j.gt.maxv) cycle
               k = ptr_out(j+1)
               ptr_out(j+1) = k + 1
               row_out(k) = row
               val_out(k) = val_in(i)
            end do
         end do
      end select
   else ! neither val_out nor map present
      maxv = n
      do row = 1, m
         if(abs(matrix_type).ge.HSL_MATRIX_REAL_SYM_PSDEF) maxv = row
         if(abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW) maxv = row-1
         do i = ptr_in(row), ptr_in(row+1)-1
            j = col_in(i)
            if(j.lt.1 .or. j.gt.maxv) cycle
            k = ptr_out(j+1)
            ptr_out(j+1) = k + 1
            row_out(k) = row
         end do
      end do
   endif

   !
   ! Third pass, removal of duplicates (no sort required by construction)
   !
   if(present(map)) then
      k = 1 ! insert position
      do col = 1, n
         i = ptr_out(col)
         ptr_out(col) = k
         if(ptr_out(col+1).eq.i) cycle ! no entries
         ! Move first entry of column forward
         row_out(k) = row_out(i)
         map(k) = map(i)
         k = k + 1
         ! Loop over remaining entries
         do i = i+1, ptr_out(col+1)-1
            if(row_out(i).eq.row_out(k-1)) then
               ! Duplicate
               allocate(dup,stat=st)
               if(st.ne.0) goto 100
               dup%next => duphead
               duphead => dup
               dup%src = map(i)
               dup%dest = k-1
               cycle
            endif
            ! Pull entry forwards
            row_out(k) = row_out(i)
            map(k) = map(i)
            k = k + 1
         end do
      end do
      ptr_out(n+1) = k
      lmap = k-1
   elseif(present(val_out)) then
      k = 1 ! insert position
      do col = 1, n
         i = ptr_out(col)
         ptr_out(col) = k
         if(ptr_out(col+1).eq.i) cycle ! no entries
         ! Move first entry of column forward
         row_out(k) = row_out(i)
         val_out(k) = val_out(i)
         k = k + 1
         ! Loop over remaining entries
         do i = i+1, ptr_out(col+1)-1
            if(row_out(i).eq.row_out(k-1)) then
               ! Duplicate
               idup = idup + 1
               val_out(k-1) = val_out(k-1) + val_out(i)
               cycle
            endif
            ! Pull entry forwards
            row_out(k) = row_out(i)
            val_out(k) = val_out(i)
            k = k + 1
         end do
      end do
      ptr_out(n+1) = k
   else ! neither val_out nor map are present
      k = 1 ! insert position
      do col = 1, n
         i = ptr_out(col)
         ptr_out(col) = k
         if(ptr_out(col+1).eq.i) cycle ! no entries
         ! Move first entry of column forward
         row_out(k) = row_out(i)
         k = k + 1
         ! Loop over remaining entries
         do i = i+1, ptr_out(col+1)-1
            if(row_out(i).eq.row_out(k-1)) then
               ! Duplicate
               idup = idup + 1
               cycle
            endif
            ! Pull entry forwards
            row_out(k) = row_out(i)
            k = k + 1
         end do
      end do
      ptr_out(n+1) = k
   endif

   if(present(map)) then
      ! Append duplicates to map
      do while(associated(duphead))
         idup = idup + 1
         map(lmap+1) = duphead%dest
         map(lmap+2) = duphead%src
         lmap = lmap + 2
         dup => duphead%next
         deallocate(duphead)
         duphead => dup
      end do
      if(present(val_out)) then
         allocate(val_out(ptr_out(n+1)-1), stat=st)
         if(st.ne.0) goto 100
         call mc69_set_values(matrix_type, lmap, map, val_in, ptr_out(n+1)-1, &
            val_out)
      endif
   endif

   if(present(val_out)) then
      select case(matrix_type)
      case(HSL_MATRIX_REAL_SYM_PSDEF)
         ! Check for positive diagonal entries
         do j = 1,n
            k = ptr_out(j)
            ! ps def case - can't reach here unless all entries have diagonal
            if(real(val_out(k))<=zero) then
               flag = MC69_ERROR_MISSING_DIAGONAL
               call mc69_print_flag(context,nout,flag)
               return
            end if 
         end do
      end select
   endif

   ! Check whether a warning needs to be raised
   if(ioor > 0 .or. idup > 0 .or. idiag < n) then
     if(ioor > 0) flag = MC69_WARNING_IDX_OOR
     if(idup > 0) flag = MC69_WARNING_DUP_IDX
     if(idup > 0 .and. ioor > 0) flag = MC69_WARNING_DUP_AND_OOR
     if(abs(matrix_type) .ne. HSL_MATRIX_REAL_SKEW) then
         if(idiag < n .and. ioor > 0) then
            flag = MC69_WARNING_MISS_DIAG_OORDUP
         else if(idiag < n .and. idup > 0) then
            flag = MC69_WARNING_MISS_DIAG_OORDUP
         else if(idiag < n) then
            flag = MC69_WARNING_MISSING_DIAGONAL
         end if
      endif
      call mc69_print_flag(context,nout,flag)
   end if

   if(present(noor)) noor = ioor
   if(present(ndup)) ndup = idup

   return

100 if(st /= 0) then
      flag = MC69_ERROR_ALLOCATION
      call mc69_print_flag(context,nout,flag)
    end if

end subroutine mc69_csrl_convert_main

!****************************************

!
! Converts CSR with upper entries only for symmetric, skew-symmetric and
! Hermitian matrices to HSL standard format
!
subroutine mc69_csru_convert_single(matrix_type, n, ptr_in, col_in, ptr_out, &
      row_out, flag, val_in, val_out, lmap, map, lp, noor, ndup) 
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: n ! number of columns
   integer, dimension(*), intent(in) :: ptr_in ! column pointers on input
   integer, dimension(*), intent(in) :: col_in ! row indices on input.
      ! These may be unordered within each column and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(*), intent(out) :: ptr_out ! col ptr output
   integer, allocatable, dimension(:), intent(out) :: row_out ! row indices out
      ! Duplicates and out-of-range entries are dealt with and
      ! the entries within each column are ordered by increasing row index.
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(in) :: val_in ! values input
   real(wp), optional, allocatable, dimension(:) :: val_out
      ! values on output
   integer, optional, intent(out) :: lmap
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives src: map(i) = j means val_out(i)=val_in(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val_in(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed

   ! Local variables
   character(50)  :: context  ! Procedure name (used when printing).
   integer :: nout ! output unit (set to -1 if lp not present)

   context = 'mc69_csru_convert'

   nout = -1
   if(present(lp)) nout = lp

   ! Note: have to change this test for complex code
   if(matrix_type < 3 .or. matrix_type == 5 .or. matrix_type > 6) then
      flag = MC69_ERROR_MATRIX_TYPE
      call mc69_print_flag(context,nout,flag)
      return
   end if

   call mc69_cscl_convert_main(context, -1, matrix_type, n, n, &
      ptr_in, col_in, ptr_out, row_out, flag, val_in, val_out, lmap, map, &
      lp, noor, ndup) 
end subroutine mc69_csru_convert_single

!****************************************

!
! Converts CSR with both lower and upper entries present to CSC with 
! only lower entries present and entries
! within each column ordered by increasing row index (HSL standard format)
!
subroutine mc69_csrlu_convert_single(matrix_type, n, ptr_in, col_in, ptr_out, &
      row_out, flag, val_in, val_out, lmap, map, lp, noor, ndup) 
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: n ! matrix dimension
   integer, dimension(*), intent(in) :: ptr_in ! row pointers on input
   integer, dimension(*), intent(in) :: col_in ! col indices on input.
      ! These may be unordered within each row and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(*), intent(out) :: ptr_out ! col ptr output
   integer, allocatable, dimension(:), intent(out) :: row_out ! row indices out
      ! Duplicates and out-of-range entries are dealt with and
      ! the entries within each column are ordered by increasing row index.
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(in) :: val_in ! values input
   real(wp), optional, allocatable, dimension(:) :: val_out
      ! values on output
   integer, optional, intent(out) :: lmap
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives src: map(i) = j means val_out(i)=val_in(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val_in(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed

   ! Local variables
   character(50)  :: context  ! Procedure name (used when printing).
   integer :: nout ! output unit (set to -1 if nout not present)

   context = 'mc69_csrlu_convert'

   nout = -1
   if(present(lp)) nout = lp

   ! Note: have to change this test for complex code
   if(matrix_type < 3 .or. matrix_type == 5 .or. matrix_type > 6) then
      flag = MC69_ERROR_MATRIX_TYPE
      call mc69_print_flag(context,nout,flag)
      return
   end if

   call mc69_csclu_convert_main(context, 1, matrix_type, n, &
      ptr_in, col_in, ptr_out, row_out, flag, val_in, val_out, lmap, map, &
      lp, noor, ndup) 
end subroutine mc69_csrlu_convert_single

!****************************************

!
! Converts COOR format to CSC with only lower entries present for
! (skew-)symmetric problems. Entries within each column ordered by increasing
! row index (HSL standard format)
!
subroutine mc69_coord_convert_single(matrix_type, m, n, ne, row, col, ptr_out, &
      row_out, flag, val_in, val_out, lmap, map, lp, noor, ndup) 
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: m ! number of rows in matrix
   integer, intent(in) :: n ! number of columns in matrix
   integer, intent(in) :: ne ! number of input nonzero entries
   integer, dimension(:), intent(in) :: row(ne) ! row indices on input.
      ! These may be unordered within each column and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(:), intent(in) :: col(ne) ! column indices on input.
      ! These may be unordered within each column and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(:), intent(out) :: ptr_out(n+1) ! col ptr output
   integer, allocatable, dimension(:), intent(out) :: row_out ! row indices out
      ! Duplicates and out-of-range entries are dealt with and
      ! the entries within each column are ordered by increasing row index.
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(in) :: val_in ! values input
   real(wp), optional, allocatable, dimension(:) :: val_out
      ! values on output
   integer, optional, intent(out) :: lmap
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives source: map(i) = j means
      ! val_out(i)=val_in(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val_in(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed


   ! Local variables
   character(50)  :: context  ! Procedure name (used when printing).
   integer :: i
   integer :: idiag
   integer :: idup
   integer :: ioor
   integer :: j
   integer :: k, l1, l2
   integer :: l
   integer :: ne_new
   integer :: nout ! output unit (set to -1 if lp not present)
   integer :: st ! stat parameter

   type(dup_list), pointer :: dup
   type(dup_list), pointer :: duphead

   nullify(dup, duphead)

   context = 'mc69_coord_convert'

   flag = MC69_SUCCESS

   nout = -1
   if(present(lp)) nout = lp

   ! ---------------------------------------------
   ! Check that restrictions are adhered to
   ! ---------------------------------------------

   ! Note: have to change this test for complex code
   if(matrix_type < 0 .or. matrix_type == 5 .or. matrix_type > 6) then
      flag = MC69_ERROR_MATRIX_TYPE
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(m < 0 .or. n < 0) then
      flag = MC69_ERROR_N_OOR
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(abs(matrix_type).ge.HSL_MATRIX_REAL_UNSYM .and. m.ne.n) then
      flag = MC69_ERROR_M_NE_N
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(present(val_in).neqv.present(val_out)) then
      flag = MC69_ERROR_VAL_MISS
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(present(map).neqv.present(lmap)) then
      flag = MC69_ERROR_LMAP_MISS
      call mc69_print_flag(context,nout,flag)
      return
   end if

   ! ---------------------------------------------

   ! allocate output and work arrays

   deallocate(row_out,stat=st)
   if(present(val_out)) deallocate(val_out,stat=st)
   if(present(map)) deallocate(map,stat=st)


   ! ---------------------------------------------
   ! check for duplicates and/or out-of-range. check diagonal present.
   ! first do the case where values not present

   idup = 0; ioor = 0; idiag = 0

   !
   ! First pass, count number of entries in each col of the matrix
   ! matrix. Count is at an offset of 1 to allow us to play tricks
   ! (ie. ptr_out(i+1) set to hold number of entries in column i
   ! of expanded matrix).
   ! Excludes out of range entries. Includes duplicates.
   !
   ptr_out(1:n+1) = 0
   do l = 1, ne
      i = row(l)
      j = col(l)
      if(j.lt.1 .or. j.gt.n .or. i.lt.1 .or. i.gt.m) then
         ioor = ioor + 1
         cycle
      endif

      if(abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW .and. i.eq.j) then
         ioor = ioor + 1
         cycle
      endif
   
      select case (abs(matrix_type))
      case (HSL_MATRIX_REAL_SYM_PSDEF:)
         if(i.ge.j) then
            ptr_out(j+1) = ptr_out(j+1) + 1
         else
            ptr_out(i+1) = ptr_out(i+1) + 1
         end if
      case default
          ptr_out(j+1) = ptr_out(j+1) + 1
      end select
   end do


   ! Determine column starts for transposed expanded matrix such 
   ! that column i starts at ptr_out(i)
   ne_new = 0
   ptr_out(1) = 1
   do i = 2, n+1
      ne_new = ne_new + ptr_out(i)
      ptr_out(i) = ptr_out(i) + ptr_out(i-1)
   end do

   ! Check whether all entries out of range
   if(ne.gt.0 .and. ne_new.eq.0) then
      flag = MC69_ERROR_ALL_OOR
      call mc69_print_flag(context,nout,flag)
      return
   end if

   !
   ! Second pass, drop entries into place for conjugate of transposed
   ! expanded matrix
   !
   allocate(row_out(ne_new), stat=st)
   if(st.ne.0) goto 100
   if(present(map)) then
      if(allocated(map)) deallocate(map,stat=st)
      allocate(map(2*ne), stat=st)
      if(st.ne.0) goto 100
      map(:) = 0
      select case(matrix_type)
      case(HSL_MATRIX_REAL_SKEW)
         do l = 1, ne
            i = row(l)
            j = col(l)
            if(j.lt.1 .or. j.gt.n .or. i.lt.1 .or. i.gt.m .or. i.eq.j) cycle
            if(i.ge.j) then
               k=ptr_out(j)
               ptr_out(j) = k+1
               row_out(k) = i
               map(k) = l
            else
               k=ptr_out(i)
               ptr_out(i) = k+1
               row_out(k) = j
               map(k) = -l
            end if
         end do

      case(HSL_MATRIX_REAL_SYM_PSDEF, HSL_MATRIX_REAL_SYM_INDEF)
         do l = 1, ne
            i = row(l)
            j = col(l)
            if(j.lt.1 .or. j.gt.n .or. i.lt.1 .or. i.gt.m) cycle
            if(i.ge.j) then
               k=ptr_out(j)
               ptr_out(j) = k+1
               row_out(k) = i
               map(k) = l
            else
               k=ptr_out(i)
               ptr_out(i) = k+1
               row_out(k) = j
               map(k) = l
            end if
         end do
      case default
         do l = 1, ne
            i = row(l)
            j = col(l)
            if(j.lt.1 .or. j.gt.n .or. i.lt.1 .or. i.gt.m) cycle
            k=ptr_out(j)
            ptr_out(j) = k+1
            row_out(k) = i
            map(k) = l
         end do
      end select
   elseif(present(val_out)) then
      allocate(val_out(ne_new), stat=st)
      if(st.ne.0) goto 100
      select case(matrix_type)
      case(HSL_MATRIX_REAL_SKEW)
         do l = 1, ne
            i = row(l)
            j = col(l)
            if(j.lt.1 .or. j.gt.n .or. i.lt.1 .or. i.gt.m .or. i.eq.j) cycle
            if(i.ge.j) then
               k=ptr_out(j)
               ptr_out(j) = k+1
               row_out(k) = i
               val_out(k) = val_in(l)
            else
               k=ptr_out(i)
               ptr_out(i) = k+1
               row_out(k) = j
               val_out(k) = -val_in(l)
            end if
         end do

      case(HSL_MATRIX_REAL_SYM_PSDEF, HSL_MATRIX_REAL_SYM_INDEF)
         do l = 1, ne
            i = row(l)
            j = col(l)
            if(j.lt.1 .or. j.gt.n .or. i.lt.1 .or. i.gt.m) cycle
            if(i.ge.j) then
               k=ptr_out(j)
               ptr_out(j) = k+1
               row_out(k) = i
               val_out(k) = val_in(l)
            else
               k=ptr_out(i)
               ptr_out(i) = k+1
               row_out(k) = j
               val_out(k) = val_in(l)
            end if
         end do
      case default
         do l = 1, ne
            i = row(l)
            j = col(l)
            if(j.lt.1 .or. j.gt.n .or. i.lt.1 .or. i.gt.m) cycle
            k=ptr_out(j)
            ptr_out(j) = k+1
            row_out(k) = i
            val_out(k) = val_in(l)
         end do
      end select


   else
      ! neither val_out or map present
      select case(matrix_type)
      case(HSL_MATRIX_REAL_SKEW)
         do l = 1, ne
            i = row(l)
            j = col(l)
            if(j.lt.1 .or. j.gt.n .or. i.lt.1 .or. i.gt.m .or. i.eq.j) cycle
            if(i.ge.j) then
               k=ptr_out(j)
               ptr_out(j) = k+1
               row_out(k) = i
            else
               k=ptr_out(i)
               ptr_out(i) = k+1
               row_out(k) = j
            end if
         end do

      case(HSL_MATRIX_REAL_SYM_PSDEF, HSL_MATRIX_REAL_SYM_INDEF)
          do l = 1, ne
            i = row(l)
            j = col(l)
            if(j.lt.1 .or. j.gt.n .or. i.lt.1 .or. i.gt.m) cycle
            if(i.ge.j) then
               k=ptr_out(j)
               ptr_out(j) = k+1
               row_out(k) = i
            else
               k=ptr_out(i)
               ptr_out(i) = k+1
               row_out(k) = j
            end if
         end do
      case default
         do l = 1, ne
            i = row(l)
            j = col(l)
            if(j.lt.1 .or. j.gt.n .or. i.lt.1 .or. i.gt.m) cycle
            k=ptr_out(j)
            ptr_out(j) = k+1
            row_out(k) = i
         end do
      end select
   endif

   do j=n,2,-1
      ptr_out(j) = ptr_out(j-1)
   end do
   ptr_out(1) = 1

   !
   ! Third pass, in place sort and removal of duplicates
   ! Also checks for diagonal entries in pos. def. case.
   ! Uses a modified insertion sort for best speed on already ordered data
   !
   idup=0
   if(present(map)) then
      do j = 1, n
         l1 = ptr_out(j)
         l2 = ptr_out(j+1)-1
         ! sort(row_out(l1:l2)) and permute map(l1:l2) accordingly
         l = l2-l1+1
         if(l.gt.1) call sort ( row_out(l1:l2),l, map=map(l1:l2) )
      end do

      ! work through removing duplicates
      k = 1 ! insert position      
      do j = 1, n
         l1 = ptr_out(j)
         l2 = ptr_out(j+1)-1
         ptr_out(j) = k
         ! sort(row_out(l1:l2)) and permute map(l1:l2) accordingly
         l = l2-l1+1
         if(l.eq.0) cycle ! no entries

         ! Move first entry of column forward
         if(row_out(l1).eq.j) idiag = idiag + 1
         row_out(k) = row_out(l1)
         map(k) = map(l1)
         k = k + 1
         ! Loop over remaining entries
         do i = l1+1, l2
            if(row_out(i).eq.row_out(k-1)) then
               ! Duplicate
               idup = idup + 1
               allocate(dup,stat=st)
               if(st.ne.0) goto 100
               dup%next => duphead
               duphead => dup
               dup%src = map(i)
               dup%dest = k-1
               cycle
            endif
            ! Pull entry forwards
            if(row_out(i).eq.j) idiag = idiag + 1
            row_out(k) = row_out(i)
            map(k) = map(i)
            k = k + 1
         end do
      end do
      ptr_out(n+1) = k
      lmap = k-1
   else if(present(val_out)) then
      ! ADD
      do j = 1, n
         l1 = ptr_out(j)
         l2 = ptr_out(j+1)-1
         l = l2-l1+1
         if(l.gt.1) call sort( row_out(l1:l2),l,val=val_out(l1:l2) )
         ! ADD
      end do 

      ! work through removing duplicates
      k = 1 ! insert position      
      do j = 1, n
         l1 = ptr_out(j)
         l2 = ptr_out(j+1)-1
         ptr_out(j) = k
         ! sort(row_out(l1:l2)) and permute map(l1:l2) accordingly
         l = l2-l1+1
         if(l.eq.0) cycle ! no entries

         ! Move first entry of column forward
         if(row_out(l1).eq.j) idiag = idiag + 1
         row_out(k) = row_out(l1)
         val_out(k) = val_out(l1)
         k = k + 1
         ! Loop over remaining entries
         do i = l1+1, l2
            if(row_out(i).eq.row_out(k-1)) then
              idup = idup + 1
              val_out(k-1) = val_out(k-1)+val_out(i)
              cycle
            end if
            ! Pull entry forwards
            if(row_out(i).eq.j) idiag = idiag + 1
            row_out(k) = row_out(i)
            val_out(k) = val_out(i)
            k = k + 1
         end do
      end do
      ptr_out(n+1) = k
   else
      do j = 1, n
         l1 = ptr_out(j)
         l2 = ptr_out(j+1)-1
         l = l2-l1+1
         if(l.gt.1)call sort ( row_out(l1:l2),l)
         ! ADD
      end do

      ! work through removing duplicates
      k = 1 ! insert position      
      do j = 1, n
         l1 = ptr_out(j)
         l2 = ptr_out(j+1)-1
         ptr_out(j) = k
         ! sort(row_out(l1:l2)) and permute map(l1:l2) accordingly
         l = l2-l1+1
         if(l.eq.0) cycle ! no entries

         ! Move first entry of column forward
         if(row_out(l1).eq.j) idiag = idiag + 1
         row_out(k) = row_out(l1)
         k = k + 1
         ! Loop over remaining entries
         do i = l1+1, l2
            if(row_out(i).eq.row_out(k-1)) then
              idup = idup + 1
              cycle
            end if
            ! Pull entry forwards
            if(row_out(i).eq.j) idiag = idiag + 1
            row_out(k) = row_out(i)
            k = k + 1
         end do
      end do
      ptr_out(n+1) = k
   
   endif


   ! Check for missing diagonals in pos def and indef cases
   ! Note: change this test for complex case
   if(abs(matrix_type) == HSL_MATRIX_REAL_SYM_PSDEF) then
      do j = 1,n
         if(ptr_out(j).lt.ptr_out(n+1)) then
            if(row_out(ptr_out(j)) .ne. j) then
               flag = MC69_ERROR_MISSING_DIAGONAL
               call mc69_print_flag(context,nout,flag)
               return
            end if
         end if
      end do
   end if

   if(present(map)) then
      ! Append duplicates to map
      do while(associated(duphead))
         map(lmap+1) = duphead%dest
         map(lmap+2) = duphead%src
         lmap = lmap + 2
         dup => duphead%next
         deallocate(duphead)
         duphead => dup
      end do
      if(present(val_out)) then
         allocate(val_out(ptr_out(n+1)-1), stat=st)
         if(st.ne.0) goto 100
         call mc69_set_values(matrix_type, lmap, map, val_in, ptr_out(n+1)-1, &
            val_out)
      endif
   endif

   if(present(val_out)) then
      select case(matrix_type)
      case(HSL_MATRIX_REAL_SYM_PSDEF)
         ! Check for positive diagonal entries
         do j = 1,n
            k = ptr_out(j)
            ! ps def case - can't reach here unless all entries have diagonal
            if(real(val_out(k))<=zero) then
               flag = MC69_ERROR_MISSING_DIAGONAL
               call mc69_print_flag(context,nout,flag)
               return
            end if 
         end do
      end select
   endif

   ! Check whether a warning needs to be raised
   if(ioor > 0 .or. idup > 0 .or. idiag < n) then
     if(ioor > 0) flag = MC69_WARNING_IDX_OOR
     if(idup > 0) flag = MC69_WARNING_DUP_IDX
     if(idup > 0 .and. ioor > 0) flag = MC69_WARNING_DUP_AND_OOR
     if(abs(matrix_type) .ne. HSL_MATRIX_REAL_SKEW) then
         if(idiag < n .and. ioor > 0) then
            flag = MC69_WARNING_MISS_DIAG_OORDUP
         else if(idiag < n .and. idup > 0) then
            flag = MC69_WARNING_MISS_DIAG_OORDUP
         else if(idiag < n) then
            flag = MC69_WARNING_MISSING_DIAGONAL
         end if
      endif
      call mc69_print_flag(context,nout,flag)
   end if

   ! Check whether a warning needs to be raised
   if(ioor > 0 .or. idup > 0 .or. idiag < n) then
      if(ioor > 0) flag = MC69_WARNING_IDX_OOR
      if(idup > 0) flag = MC69_WARNING_DUP_IDX
      if(idup > 0 .and. ioor > 0) flag = MC69_WARNING_DUP_AND_OOR
      if(abs(matrix_type).eq.HSL_MATRIX_REAL_SYM_INDEF .and. &
            idiag<n .and. ioor>0) then
         flag = MC69_WARNING_MISS_DIAG_OORDUP
      else if(abs(matrix_type).eq.HSL_MATRIX_REAL_SYM_INDEF .and. idiag<n .and.&
            idup>0) then
         flag = MC69_WARNING_MISS_DIAG_OORDUP
      else if(abs(matrix_type).eq.HSL_MATRIX_REAL_SYM_INDEF .and. idiag<n) then
         flag = MC69_WARNING_MISSING_DIAGONAL
      end if
      call mc69_print_flag(context,nout,flag)
   end if

   if(present(noor)) noor = ioor
   if(present(ndup)) ndup = idup
   return

   100 continue
   if(st /= 0) then
      flag = MC69_ERROR_ALLOCATION
      call mc69_print_flag(context,nout,flag)
   end if
   return

end subroutine mc69_coord_convert_single

!*************************************************

!
! This subroutine will use map to translate the values of val to val_out
!
subroutine mc69_set_values_single(matrix_type, lmap, map, val, ne, val_out)
   integer, intent(in) :: matrix_type
   integer, intent(in) :: lmap
   integer, dimension(lmap), intent(in) :: map
   real(wp), dimension(*), intent(in) :: val
   integer, intent(in) :: ne
   real(wp), dimension(ne), intent(out) :: val_out

   integer :: i, j, k

   select case(matrix_type)
   case default
      !
      ! Rectangular, Unsymmetric or Symmetric Matrix
      !

      ! First set val_out using first part of map
      do i = 1, ne
         j = abs(map(i))
         val_out(i) = val(j)
      end do

      ! Second examine list of duplicates
      do i = ne+1, lmap, 2
         j = abs(map(i))
         k = abs(map(i+1))
         val_out(j) = val_out(j) + val(k)
      end do
   case(HSL_MATRIX_REAL_SKEW)
      !
      ! Skew symmetric Matrix
      !

      ! First set val_out using first part of map
      do i = 1, ne
         j = abs(map(i))
         val_out(i) = sign(1.0,real(map(i)))*val(j)
      end do

      ! Second examine list of duplicates
      do i = ne+1, lmap, 2
         j = abs(map(i))
         k = abs(map(i+1))
         val_out(j) = val_out(j) + sign(1.0,real(map(i+1)))*val(k)
      end do
   end select
end subroutine mc69_set_values_single

!*************************************************

subroutine mc69_print_flag(context,nout,flag)
   integer, intent (in) :: flag, nout
   character (len=*), optional, intent(in) :: context

   if(nout < 0) return
   if(flag < 0) then
      write (nout,'(/3a,i3)') ' Error return from ',trim(context), &
         '. Error flag = ', flag
   else
      write (nout,'(/3a,i3)') ' Warning from ',trim(context), &
         '. Warning flag = ', flag
   end if

   select case(flag)
   !
   ! Errors
   !
   case(MC69_ERROR_ALLOCATION)
      write (nout,'(a)') ' Allocation error'
   case(MC69_ERROR_MATRIX_TYPE)
       write (nout,'(a)') ' matrix_type has invalid value'
   case(MC69_ERROR_N_OOR)
      write (nout,'(a)') ' m or n is out-of-range'
   case(MC69_ERROR_ALL_OOR)
      write (nout,'(a)') ' All entries in a column out-of-range'
   case(MC69_ERROR_PTR_MONO)
      write (nout,'(a)') ' ptr not monotonic'
   case(MC69_ERROR_PTR_1)
      write (nout,'(a)') ' ptr(1) < 1'
   case(MC69_ERROR_IMAG_DIAGONAL)
      write (nout,'(a)') ' one or more diagonal entries is not real'
   case(MC69_ERROR_MISSING_DIAGONAL)
      write (nout,'(a)') ' one or more diagonal entries are not positive'
   case(MC69_ERROR_VAL_MISS)
      write (nout,'(a)') ' Only one of val and val_out is present'
   case(MC69_ERROR_LMAP_MISS)
      write (nout,'(a)') ' Only one of lmap and map is present'
   case(MC69_ERROR_UPR_ENTRY)
      write (nout,'(a)') ' Entry in upper triangle'
   case(MC69_ERROR_M_NE_N)
      write (nout,'(a)') ' m is not equal to n'
   !
   ! Warnings
   !
   case(MC69_WARNING_IDX_OOR)
      write (nout,'(a)') ' out-of-range indices detected'
   case(MC69_WARNING_DUP_IDX)
      write (nout,'(a)') ' duplicate entries detected'
   case(MC69_WARNING_DUP_AND_OOR)
      write (nout,'(a)') &
         ' out-of-range indices detected and duplicate entries detected'
   case(MC69_WARNING_MISSING_DIAGONAL)
      write (nout,'(a)') ' one or more diagonal entries is missing'
   case(MC69_WARNING_MISS_DIAG_OORDUP)
      write (nout,'(a)') ' one or more diagonal entries is missing and'
      write (nout,'(a)') ' out-of-range and/or duplicate entries detected'
   end select

end subroutine mc69_print_flag

!************************************************************************

!
!   Sort an integer array by heapsort into ascending order.
!
subroutine sort( array, n, map, val )
   integer, intent(in) :: n       ! Size of array to be sorted
   integer, dimension(n), intent(inout) :: array ! Array to be sorted
   integer, dimension(n), optional, intent(inout) :: map
   real(wp), dimension(n), optional, intent(inout) :: val ! Apply same
      ! permutation to val

   integer :: i
   integer :: temp
   real(wp) :: vtemp
   integer :: root

   if(n.le.1) return ! nothing to do

   !
   ! Turn array into a heap with largest element on top (as this will be pushed
   ! on the end of the array in the next phase)
   !
   ! We start at the bottom of the heap (i.e. 3 above) and work our way
   ! upwards ensuring the new "root" of each subtree is in the correct
   ! location
   root = n / 2
   do root = root, 1, -1
      call pushdown(root, n, array, val=val, map=map)
   end do

   !
   ! Repeatedly take the largest value and swap it to the back of the array
   ! then repeat above code to sort the array
   !
   do i = n, 2, -1
      ! swap a(i) and head of heap a(1)
      temp = array(1)
      array(1) = array(i)
      array(i) = temp
      if(present(val)) then
         vtemp = val(1)
         val(1) = val(i)
         val(i) = vtemp
      endif
      if(present(map)) then
         temp = map(1)
         map(1) = map(i)
         map(i) = temp
      endif
      call pushdown(1,i-1, array, val=val, map=map)
   end do
end subroutine sort

!****************************************

! This subroutine will assume everything below head is a heap and will
! push head downwards into the correct location for it
subroutine pushdown(root, last, array, val, map)
   integer, intent(in) :: root
   integer, intent(in) :: last
   integer, dimension(last), intent(inout) :: array
   real(wp), dimension(last), optional, intent(inout) :: val
   integer, dimension(last), optional, intent(inout) :: map

   integer :: insert ! current insert position
   integer :: test ! current position to test
   integer :: root_idx ! value of array(root) at start of iteration
   real(wp) :: root_val ! value of val(root) at start of iteration
   integer :: root_map ! value of map(root) at start of iteration

   ! NB a heap is a (partial) binary tree with the property that given a
   ! parent and a child, array(child)>=array(parent).
   ! If we label as
   !                      1
   !                    /   \
   !                   2     3
   !                  / \   / \
   !                 4   5 6   7
   ! Then node i has nodes 2*i and 2*i+1 as its children

   if(present(val) .and. present(map)) then ! both val and map
      root_idx = array(root)
      root_val = val(root)
      root_map = map(root)
      insert = root
      test = 2*insert
      do while(test.le.last)
         ! First check for largest child branch to descend
         if(test.ne.last) then
            if(array(test+1).gt.array(test)) test = test + 1
         endif
         if(array(test).le.root_idx) exit ! root gets tested here
         ! Otherwise, move on to next level down, percolating value up
         array(insert) = array(test);
         val(insert) = val(test);
         map(insert) = map(test)
         insert = test
         test = 2*insert
      end do
      ! Finally drop root value into location found
      array(insert) = root_idx
      val(insert) = root_val
      map(insert) = root_map
   elseif(present(val)) then ! val only, not map
      root_idx = array(root)
      root_val = val(root)
      insert = root
      test = 2*insert
      do while(test.le.last)
         ! First check for largest child branch to descend
         if(test.ne.last) then
            if(array(test+1).gt.array(test)) test = test + 1
         endif
         if(array(test).le.root_idx) exit ! root gets tested here
         ! Otherwise, move on to next level down, percolating value up
         array(insert) = array(test)
         val(insert) = val(test)
         insert = test
         test = 2*insert
      end do
      ! Finally drop root value into location found
      array(insert) = root_idx
      val(insert) = root_val
   elseif(present(map)) then ! map only, not val
      root_idx = array(root)
      root_map = map(root)
      insert = root
      test = 2*insert
      do while(test.le.last)
         ! First check for largest child branch to descend
         if(test.ne.last) then
            if(array(test+1).gt.array(test)) test = test + 1
         endif
         if(array(test).le.root_idx) exit ! root gets tested here
         ! Otherwise, move on to next level down, percolating mapue up
         array(insert) = array(test)
         map(insert) = map(test)
         insert = test
         test = 2*insert
      end do
      ! Finally drop root mapue into location found
      array(insert) = root_idx
      map(insert) = root_map
   else ! neither map nor val
      root_idx = array(root)
      insert = root
      test = 2*insert
      do while(test.le.last)
         ! First check for largest child branch to descend
         if(test.ne.last) then
            if(array(test+1).gt.array(test)) test = test + 1
         endif
         if(array(test).le.root_idx) exit ! root gets tested here
         ! Otherwise, move on to next level down, percolating value up
         array(insert) = array(test)
         insert = test
         test = 2*insert
      end do
      ! Finally drop root value into location found
      array(insert) = root_idx
   endif

end subroutine pushdown

!****************************************

subroutine cleanup_dup(duphead)
   type(dup_list), pointer :: duphead ! NB: can't have both intent() and pointer

   type(dup_list), pointer :: dup

   do while(associated(duphead))
      dup => duphead%next
      deallocate(duphead)
      duphead => dup
   end do
end subroutine cleanup_dup

end module hsl_mc69_single
! COPYRIGHT (c) 2006 Council for the Central Laboratory
!                    of the Research Councils
! This package may be copied and used in any application, provided no
! changes are made to these or any other lines.
! Original date 21 February 2006. Version 1.0.0.
! 6 March 2007 Version 1.1.0. Argument stat made non-optional

MODULE HSL_ZD11_single

!  ==========================
!  Sparse matrix derived type
!  ==========================

  TYPE, PUBLIC :: ZD11_type
    INTEGER :: m, n, ne
    CHARACTER, ALLOCATABLE, DIMENSION(:) :: id, type
    INTEGER, ALLOCATABLE, DIMENSION(:) :: row, col, ptr
    REAL ( KIND( 1.0E+0 ) ), ALLOCATABLE, DIMENSION(:) :: val
  END TYPE

CONTAINS

   SUBROUTINE ZD11_put(array,string,stat)
     CHARACTER, allocatable :: array(:)
     CHARACTER(*), intent(in) ::  string
     INTEGER, intent(OUT) ::  stat

     INTEGER :: i,l

     l = len_trim(string)
     if (allocated(array)) then
        deallocate(array,stat=stat)
        if (stat/=0) return
     end if
     allocate(array(l),stat=stat)
     if (stat/=0) return
     do i = 1, l
       array(i) = string(i:i)
     end do

   END SUBROUTINE ZD11_put

   FUNCTION ZD11_get(array)
     CHARACTER, intent(in):: array(:)
     CHARACTER(size(array)) ::  ZD11_get
! Give the value of array to string.

     integer :: i
     do i = 1, size(array)
        ZD11_get(i:i) = array(i)
     end do

   END FUNCTION ZD11_get

END MODULE HSL_ZD11_single


! COPYRIGHT (c) 2007 Science and Technology Facilities Council
!           and Jacko Koster (Trondheim, Norway)
!
! Version 2.3.1
! History: See ChangeLog
!

! To convert from double to single:
! * s/_single/_single/g
! * Change wp
! * Change: MC64A, MC64D, MC64E, MC64F, MC64R, MC64Q, mc64a, mc64w
!

module hsl_mc64_single
   use hsl_mc34_single
   use hsl_mc69_single
   use hsl_zd11_single
   implicit none

   private
   public :: mc64_control, mc64_info
   public :: mc64_initialize, mc64_matching

   integer, parameter :: wp = kind(0.0)

   type mc64_control
!     real(wp) :: relax = 0d0 ! Relaxes matching
      integer :: lp = 6       ! Unit for error messages
      integer :: wp = 6       ! Unit for warning messages
      integer :: sp = -1      ! Unit for statistical output
      integer :: ldiag = 2    ! Controls level of diagnostic output
      integer :: checking = 0 ! Control for checking input
   end type mc64_control

   type mc64_info
      integer :: flag   ! Flags success or failure case
      integer :: more    ! More information on failure
      integer :: strucrank ! Structural rank
      integer :: stat    ! STAT value after allocate failure
   end type mc64_info

   interface mc64_matching
      module procedure mc64_matching_zd11_single
      module procedure mc64_matching_hslstd_single
   end interface mc64_matching

contains

   subroutine mc64_initialize(control)
      type(mc64_control), intent(out) :: control

      type(mc64_control) :: default

!       control%relax = default%relax
        control%lp = default%lp
        control%wp = default%wp
        control%sp = default%sp
        control%ldiag = default%ldiag
        control%checking = default%checking
    end subroutine mc64_initialize

   subroutine mc64_matching_zd11_single(job,matrix,control,info,perm, &
         scale)
      integer, intent(in) :: job ! Control parameter for algorithm choice
      type(zd11_type), intent(in) :: matrix
      type(mc64_control), intent(in) :: control
      type(mc64_info), intent(out) :: info
      integer, intent(out) :: perm(matrix%m + matrix%n)
      real(wp), optional, intent(out) :: scale(matrix%m + matrix%n)

      integer :: matrix_type

      matrix_type = HSL_MATRIX_UNDEFINED

      if (allocated(matrix%id)) then
        if (matrix%id(1) == 'S' .or. matrix%id(1) == 's') then
          matrix_type = HSL_MATRIX_REAL_SYM_INDEF
        endif
      endif

      if (control%checking == 0) then
        if (matrix%n .lt. 1) then
          info%flag = -2
          info%more = matrix%n
          if (control%ldiag>0 .and. control%lp>0 ) &
              write (control%lp,'(/a,i5/a,i5)') &
              'Error return from MC64 with info%flag =', &
              info%flag, &
              'Value of matrix%n out-of-range', info%more
          return
        endif

        if(matrix%ne .ne. matrix%ptr(matrix%n+1)-1) then
          info%flag = -3
          info%more = matrix%ptr(matrix%n+1)-1
          if (control%ldiag>0 .and. control%lp>0 ) &
              write (control%lp,'(/a,i5/a,i5)') &
              'Error return from MC64 with info%flag =', &
              info%flag, &
              'Value of ptr(n+1)-1!=ne', info%more
          return
        endif
      endif

      call mc64_matching_hslstd_single(job,matrix_type,matrix%m, &
         matrix%n,matrix%ptr,matrix%row,matrix%val,control,info,perm,scale)
   end subroutine mc64_matching_zd11_single

   subroutine mc64_matching_hslstd_single(job,matrix_type,m,n,ptr,row, &
         val,control,info,perm,scale)
      integer, intent(in) :: job ! Control parameter for algorithm choice
      integer, intent(in) :: matrix_type
      integer, intent(in) :: m
      integer, intent(in) :: n
      integer, dimension(n+1), intent(in) :: ptr
      integer, dimension(*), intent(in) :: row
      real(wp), dimension(*), intent(in) :: val
      type(mc64_control), intent(in) :: control
      type(mc64_info), intent(out) :: info
      integer, intent(out) :: perm(m + n)
      real(wp), optional, intent(out) :: scale(m + n)

      integer, allocatable :: iw(:)
      real(wp), allocatable :: dw(:)
      real(wp) :: cntl(10)
      integer :: i,j,k,ne,lliw,liw,lldw,ldw,ndiag,num,stat &
                 ,icntl(10),info64(10)
      real(wp), parameter :: zero = 0.0
      logical sym

      stat = 0

! Check input data if checking is requested
      if (control%checking == 0) then
        if (job .lt. 1 .or. job .gt. 5) then
          info%flag = -1
          info%more = job
          if (control%ldiag>0 .and. control%lp>0 ) &
              write (control%lp,'(/a,i5/a,i5)') &
              'Error return from MC64 with info%flag =', &
              info%flag, &
              'Value of job out-of-range', info%more
          return
        endif

        if (n .lt. 1) then
          info%flag = -2
          info%more = n
          if (control%ldiag>0 .and. control%lp>0 ) &
              write (control%lp,'(/a,i5/a,i5)') &
              'Error return from MC64 with info%flag =', &
              info%flag, &
              'Value of matrix%n out-of-range', info%more
          return
        endif

        if (ptr(n+1)-1 .lt. 1) then
          info%flag = -3
          info%more = ptr(n+1)-1
          if (control%ldiag>0 .and. control%lp>0 ) &
              write (control%lp,'(/a,i5/a,i5)') &
              'Error return from MC64 with info%flag =', &
              info%flag, &
              'Value of ptr(n+1)-1 out-of-range', info%more
          return
        endif

        if (m .lt. n) then
          info%flag = -4
          info%more = m
          if (control%ldiag>0 .and. control%lp>0 ) &
              write (control%lp,'(/a,i5/a,i5)') &
              'Error return from MC64 with info%flag =', &
              info%flag, &
              'Value of matrix%m less than matrix%n', info%more
          return
        endif

! end of checking
      endif

      ne = ptr(n+1)-1

      icntl(1) = control%lp
      icntl(2) = control%wp
      icntl(3) = control%sp

! Action if matrix is symmetric
      sym = (matrix_type.ge.HSL_MATRIX_REAL_SYM_PSDEF)
      if(sym .and. job.eq.5) then
         call sym_maxprod_match(n,ptr,row,val,control,info,perm,scale)
         return
      endif
      if (sym) then
        ndiag = 0
! Matrix is symmetric ... check only lower triangle supplied
        if (control%checking == 0) then
          do j = 1,n
            do k = ptr(j),ptr(j+1)-1
              if (row(k) .lt. j) then
                info%flag = -8
                info%more = j
                if (control%ldiag>0 .and. control%lp>0 ) &
                  write (control%lp,'(/a,i5/a,i5)') &
                  'Error return from MC64 with info%flag =', &
                   info%flag, &
       'Input symmetric matrix has entries in upper triangle in column', &
                   info%more
                return
              endif
              if (row(k) == j) ndiag = ndiag + 1
            enddo
          enddo
        else
          do j = 1,n
            do k = ptr(j),ptr(j+1)-1
              if (row(k) == j) ndiag = ndiag + 1
            enddo
          enddo
        endif
! Reset ne for the expanded symmetric matrix
        ne = 2*ne - ndiag
      endif


      if (job == 1) liw = 4*n + m
      if (job == 2) liw = 2*n + 2*m
      if (job == 3) liw = 8*n+2*m+ne
      if (job == 4) liw = 3*n+2*m
      if (job == 5) liw = 3*n+2*m
!     if (job == 6) liw = 3*n+2*m+ne

      if (job == 1) ldw = 1
      if (job == 2) ldw = m
      if (job == 3) ldw = ne
      if (job == 4) ldw = 2*m+ne
      if (job == 5) ldw = 2*m+n+ne
!     if (job == 6) ldw = 3*m+n+ne

! Expand matrix if symmetric
      if (sym) then
        lliw = liw+n+1+2*ne
        lldw = ldw+2*ne
        allocate (iw(lliw),stat=stat)
        if (stat/=0) go to 100
        allocate (dw(lldw),stat=stat)
        if (stat/=0) go to 100

        iw(n+2:n+1+ptr(n+1)-1) = row(1:ptr(n+1)-1)
        dw(1:ptr(n+1)-1) = val(1:ptr(n+1)-1)
        iw(1:n+1) = ptr(1:n+1)
        !call mc34ad(n,iw(n+2),iw,.true.,dw,iw(n+2+2*ne))
        call mc34_expand(n,iw(n+2:),iw,iw(n+2+2*ne:),a=dw)
      else
        allocate (iw(liw),stat=stat)
        if (stat/=0) go to 100
        allocate (dw(ldw),stat=stat)
        if (stat/=0) go to 100
      endif

      icntl(4) = control%checking
      icntl(5:10) = 0

!     if (m==n .and. control%relax == zero .and. job .le. 5) then
      if (m==n) then
! Needed because no ldiag in mc64
        if (control%ldiag .lt. 1) then
          icntl(1) = -1
          icntl(2) = -1
          icntl(3) = -1
        endif
        if (control%ldiag .eq. 1) then
          icntl(2) = -1
          icntl(3) = -1
        endif
        if (control%ldiag .eq. 2) then
          icntl(3) = -1
        endif
        if (sym)  then
! Call HSL F77 code
          call mc64a(job,n,ne,iw,iw(n+2),dw,num,perm, &
                      liw,iw(n+2+2*ne),ldw,dw(2*ne+1),icntl,info64)
          info%flag   = info64(1)
          info%more   = info64(2)
          if (info64(1) == 2) info%flag = -9
          if (info%flag .lt. 0) go to 80
! Note that in this case m=n
! Set column permutation and row permutation to be the same
          perm(m+1:m+n) = perm(1:m)
        else
          call mc64a(job,n,ne,ptr,row,val,num,  &
                      perm(m+1),liw,iw,ldw,dw,icntl,info64)
          info%flag   = info64(1)
          info%more   = info64(2)
          if (info%flag .lt. 0) go to 80
! Set row permutation to identity
          do i = 1,m
            perm(i) = i
          enddo
! Invert column permutation
          do i = 1,n
            iw(abs(perm(m+i))) = sign(i,perm(m+i))
          enddo
          do i = 1,n
            perm(m+i) = iw(i)
          enddo
        endif
      else
        icntl(5) = control%ldiag
        icntl(6:10) = 0
!       cntl(1) = control%relax
        cntl(1) = zero
        call mc64_hsl_a(job,m,n,ne,ptr,row,val, &
                  num,perm,liw,iw,ldw,dw,icntl,cntl,info64)
          info%flag   = info64(1)
          info%more   = info64(2)
          if (info%flag .lt. 0) go to 80
! Set column permutation to identity
        do i = 1,n
          perm(m+i) = i
        enddo
      endif

      info%strucrank  = num

      if (present(scale) .and. job == 5) then
      !!! Commented as handled by seperate subroutine now
      !!!  if (sym) then
      !!!    scale(1:m) = (dw(2*ne+1:2*ne+n)+dw(2*ne+n+1:2*ne+2*n))/2
      !!!    scale(m+1:m+n) = scale(1:m)
      !!!  else
          scale(1:m) = dw(1:m)
          scale(m+1:m+n) = dw(m+1:m+n)
      !!!  endif
      endif

  80  deallocate (iw, stat=stat)
      if (stat/=0) go to 100
      deallocate (dw, stat=stat)
      if (stat/=0) go to 100

      return

  100 info%flag = -5
      info%stat = stat
      if (control%ldiag>0 .and. control%lp>0 ) &
          write (control%lp,'(/a,i5/a,i5)') &
         'Error return from MC64 with info%flag =', info%flag, &
         'Allocate or deallocate failed with STAT=',stat

   end subroutine mc64_matching_hslstd_single

   subroutine sym_maxprod_match(n,ptr,row,val,control,info,perm,scale)
      integer, intent(in) :: n
      integer, dimension(n+1), intent(in) :: ptr
      integer, dimension(*), intent(in) :: row
      real(wp), dimension(*), intent(in) :: val
      type(mc64_control), intent(in) :: control
      type(mc64_info), intent(out) :: info
      integer, intent(out) :: perm(2*n)
      real(wp), optional, intent(out) :: scale(2*n)

      integer, allocatable :: ptr2(:), row2(:), iw(:), new_to_old(:), &
         old_to_new(:), cperm(:)
      real(wp), allocatable :: val2(:), dw(:), cmax(:), cscale(:)
      real(wp) :: colmax
      integer :: i,j,k,ne,num,stat,nn,j1,j2,jj
      real(wp), parameter :: zero = 0.0

      stat = 0
      ne = ptr(n+1)-1

! Matrix is symmetric ... check only lower triangle supplied
      if (control%checking == 0) then
        do j = 1,n
          do k = ptr(j),ptr(j+1)-1
            if (row(k) .lt. 1 .or. row(k).gt.n) then
              info%flag = -6
              info%more = j
              if (control%ldiag>0 .and. control%lp>0 ) &
                write (control%lp,'(/a,i5/a,i5)') &
                  'Error return from MC64 with info%flag =', &
                  info%flag, &
                  'Input symmetric matrix has entries with invald row index &
                  &in position', &
                  info%more
              return
            endif
            if (row(k) .lt. j) then
              info%flag = -8
              info%more = j
              if (control%ldiag>0 .and. control%lp>0 ) &
                write (control%lp,'(/a,i5/a,i5)') &
                  'Error return from MC64 with info%flag =', &
                   info%flag, &
                   'Input symmetric matrix has entries in upper triangle in &
                   &position', &
                   info%more
              return
            endif
          enddo
        enddo
      endif
! Reset ne for the expanded symmetric matrix
      ne = 2*ne

! Expand matrix, drop explicit zeroes and take log absolute values
      allocate (ptr2(n+1), row2(ne), val2(ne), &
                iw(5*n), dw(2*n), cmax(n), stat=stat)
      if (stat/=0) then
         info%flag = -5
         info%stat = stat
         if (control%ldiag>0 .and. control%lp>0 ) &
             write (control%lp,'(/a,i5/a,i5)') &
            'Error return from MC64 with info%flag =', info%flag, &
            'Allocate or deallocate failed with STAT=',stat
         return
      endif

      if(control%checking==0) then
         ! checking enabled, look for duplicates
         iw(1:n) = 0
         do i = 1, n
           do j = ptr(i), ptr(i+1)-1
             if(iw(row(j)).ge.i) then
               ! duplicate detected
               info%flag = -7
               info%more = j
               if (control%ldiag>0 .and. control%lp>0 ) &
                 write (control%lp,'(/a,i5/a,i5)') &
                   'Error return from MC64 with info%flag =', &
                    info%flag, &
                    'Input symmetric matrix has duplicate entries in &
                    &position', &
                    info%more
               return
             endif
             iw(row(j)) = i
           end do
         end do
      endif

      k = 1
      do i = 1, n
         ptr2(i) = k
         do j = ptr(i), ptr(i+1)-1
            if(val(j).eq.zero) cycle
            row2(k) = row(j)
            val2(k) = abs(val(j))
            k = k + 1
         end do
         ! Following log is seperated from above loop to expose expensive
         ! log operation to vectorization.
         val2(ptr2(i):k-1) = log(val2(ptr2(i):k-1))
      end do
      ptr2(n+1) = k
      call mc34_expand(n, row2, ptr2, iw, a=val2)

! Compute column maximums
      do i = 1, n
         colmax = maxval(val2(ptr2(i):ptr2(i+1)-1))
         cmax(i) = colmax
         val2(ptr2(i):ptr2(i+1)-1) = colmax - val2(ptr2(i):ptr2(i+1)-1)
      end do

      call mc64w(n,ne,ptr2,row2,val2,perm,num,iw(1),iw(n+1), &
                  iw(2*n+1),iw(3*n+1),iw(4*n+1),dw(1),dw(n+1))

      info%flag   = 0
      info%strucrank  = num

      if(num.eq.n) then ! Full rank
! Note that in this case m=n
! Set column permutation and row permutation to be the same
         perm(n+1:n+n) = perm(1:n)


         if (present(scale)) then
             scale(1:n) = (dw(1:n)+dw(n+1:2*n)-cmax(1:n))/2
             !scale(n+1:n+n) = scale(1:n)
             ! above line causes segfault in ifort, workaround:
             do i = 1, n
               scale(n+i) = scale(i)
             end do
         endif
         return
      endif

      ! If we reach this point then structually rank deficient:
      ! Build a full rank submatrix and call mc64w on it.
      info%flag = 1 ! structually singular warning

      allocate(old_to_new(n),new_to_old(n),cperm(n),stat=info%stat)
      if (info%stat.ne.0) then
         info%flag = -5
         info%stat = stat
         if (control%ldiag>0 .and. control%lp>0 ) &
             write (control%lp,'(/a,i5/a,i5)') &
            'Error return from MC64 with info%flag =', info%flag, &
            'Allocate or deallocate failed with STAT=',stat
         return
      end if

      j = num + 1
      k = 0
      do i = 1,n
         if (perm(i) < 0) then
            ! row i is not part of the matching
            old_to_new(i) = -j
            j = j + 1
         else
            k = k + 1
            ! old_to_new(i) holds the new index for variable i after
            ! removal of singular part and new_to_old(k) is the
            ! original index for k
            old_to_new(i) = k
            new_to_old(k) = i
         end if
      end do

      ! Overwrite ptr2, row2 and val2
      ne = 0
      k = 0
      ptr2(1) = 1
      j2 = 1
      do i = 1,n
         j1 = j2
         j2 = ptr2(i+1)
         ! skip over unmatched entries
         if (perm(i) < 0) cycle
         k = k + 1
         do j = j1,j2-1
            jj = row2(j)
            if (perm(jj) < 0) cycle
            ne = ne + 1
            row2(ne) = old_to_new(jj)
            val2(ne) = val2(j)
         end do
         ptr2(k+1) = ne + 1
      end do
      ! nn is order of non-singular part.
      nn = k
      call mc64w(nn,ne,ptr2,row2,val2,cperm,num, &
                  iw(1),iw(nn+1),iw(2*nn+1),iw(3*nn+1),iw(4*nn+1), &
                  dw(1),dw(nn+1))

      if(present(scale)) then
          do i = 1,n
             j = old_to_new(i)
             if (j < 0) then
                scale(i) = -huge(scale)
             else
               ! Note: we need to subtract col max using old matrix numbering
               scale(i) = (dw(j)+dw(nn+j)-cmax(i))/2
            end if
         end do
      endif

      perm(1:n) = -1
      do i = 1,nn
         j = cperm(i)
         perm(new_to_old(i)) = j
      end do

      do i = 1, n
         if(perm(i).eq.-1) then
            perm(i) = old_to_new(i)
         endif
      end do

      perm(n+1:n+n) = perm(1:n)

      ! Apply Duff and Pralet correction to unmatched row scalings
      if(present(scale)) then
         allocate(cscale(n), stat=stat)
         if(stat/=0) then
            info%flag = -5
            info%stat = stat
            if (control%ldiag>0 .and. control%lp>0 ) &
                write (control%lp,'(/a,i5/a,i5)') &
               'Error return from MC64 with info%flag =', info%flag, &
               'Allocate or deallocate failed with STAT=',stat
            return
         endif
         ! For columns i not in the matched set I, set
         !     s_i = 1 / (max_{k in I} | a_ik s_k |)
         ! with convention that 1/0 = 1
         cscale(1:n) = scale(1:n)
         do i = 1,n
            do j = ptr(i), ptr(i+1)-1
               k = row(j)
               if(cscale(i).eq.-huge(scale).and.cscale(k).ne.-huge(scale)) then
                  ! i not in I, k in I
                  scale(i) = max(scale(i), log(abs(val(j)))+scale(k))
               endif
               if(cscale(k).eq.-huge(scale).and.cscale(i).ne.-huge(scale)) then
                  ! k not in I, i in I
                  scale(k) = max(scale(k), log(abs(val(j)))+scale(i))
               endif
            end do
         end do
         do i = 1,n
            if(cscale(i).ne.-huge(scale)) cycle ! matched part
            if(scale(i).eq.-huge(scale)) then
               scale(i) = zero
            else
               scale(i) = -scale(i)
            endif
         end do
      endif

   end subroutine sym_maxprod_match

!**********************************************************************
!CCCC LAST UPDATE Tue Nov 26 03:20:26 MET 2002
SUBROUTINE MC64_HSL_A (JOB, M, N, NE, IP, IRN, A, NUM, PERM, LIW,&
      IW, LDW, DW, ICNTL, CNTL, INFO)
      IMPLICIT NONE
!
! *** Copyright (c) 2002  I.S. Duff and J. Koster                   ***
! *** Although every effort has been made to ensure robustness and  ***
! *** reliability of the subroutines in this MC64 suite, we         ***
! *** disclaim any liability arising through the use or misuse of   ***
! *** any of the subroutines.                                       ***
!
! Purpose
! =======
!
! This subroutine attempts to find a permutation for an MxN, M>=N,
! sparse matrix A = {a_ij} that makes the permuted matrix have N
! entries on its diagonal.
! If the matrix is structurally nonsingular, the subroutine optionally
! returns a permutation that maximizes the smallest element on the
! diagonal, maximizes the sum of the diagonal entries, or maximizes
! the product of the diagonal entries of the permuted matrix.
! For the latter option, the subroutine also finds scaling factors that
! may be used to scale the matrix so that the nonzero diagonal entries
! of the permuted matrix are one in absolute value and all the
! off-diagonal entries are less than or equal to one in absolute value.
! The natural logarithms of the scaling factors u(i), i=1..M, for the
! rows and v(j), j=1..N, for the columns are returned so that the
! scaled matrix B = {b_ij} has entries b_ij = a_ij * EXP(u_i + v_j).
! The scaling factors are returned by this subroutine, but the actual
! scaling of the matrix has to be performed by the calling program.
!
! Parameters
! ==========
!
      INTEGER NICNTL, NCNTL, NINFO
      PARAMETER (NICNTL = 10, NCNTL = 10, NINFO = 10)

      INTEGER JOB, M, N, NE, NUM, LIW, LDW
      INTEGER IP (N + 1), IRN (NE), PERM (M), IW (LIW)
      INTEGER ICNTL(NICNTL), INFO (NINFO)
      REAL(WP) A (NE)
!
! JOB is an INTEGER variable which must be set by the user to control
! the action. It is not altered by the subroutine.
! Possible values for JOB are:
!   1 Compute a column permutation of the matrix so that the
!     permuted matrix has as many entries on its diagonal as possible.
!     The values on the diagonal are of arbitrary size. HSL subroutine
!     MC21A/MC64Z is used for this. See [1].
!   2 Compute a column permutation of the matrix so that the smallest
!     value on the diagonal of the permuted matrix is maximized.
!     See [3].
!   3 Compute a column permutation of the matrix so that the smallest
!     value on the diagonal of the permuted matrix is maximized.
!     The algorithm differs from the one used for JOB = 2 and may
!     have quite a different performance. See [2].
!   4 Compute a column permutation of the matrix so that the sum
!     of the diagonal entries of the permuted matrix is maximized.
!     See [3].
!   5 Compute a column permutation of the matrix so that the product
!     of the diagonal entries of the permuted matrix is maximized
!     and vectors to scale the matrix so that the nonzero diagonal
!     entries of the permuted matrix are one in absolute value and
!     all the off-diagonal entries are less than or equal to one in
!     absolute value. See [3].
!  Restriction: 1 <= JOB <= 5.
!
! M is an INTEGER variable which must be set by the user to the
!   number of rows of the matrix A. It is not altered by the
!   subroutine. Restriction: M >= N.
!
! N is an INTEGER variable which must be set by the user to the
!   number of columns of the matrix A. It is not altered by the
!   subroutine. Restriction: N >= 1.
!
! NE is an INTEGER variable which must be set by the user to the
!   number of entries in the matrix. It is not altered by the
!   subroutine. Restriction: NE >= 1.
!
! IP is an INTEGER array of length N+1.
!   IP(J), J=1..N, must be set by the user to the position in array IRN
!   of the first row index of an entry in column J. IP(N+1) must be set
!   to NE+1. It is not altered by the subroutine.
!
! IRN is an INTEGER array of length NE.
!   IRN(K), K=1..NE, must be set by the user to hold the row indices of
!   the entries of the matrix. Those belonging to column J must be
!   stored contiguously in the positions IP(J)..IP(J+1)-1. The ordering
!   of the row indices within each column is unimportant. Repeated
!   entries are not allowed. The array IRN is not altered by the
!   subroutine.
!
! A is a REAL array of length NE.
!   The user must set A(K), K=1..NE, to the numerical value of the
!   entry that corresponds to IRN(K).
!   It is not used by the subroutine when JOB = 1.
!   It is not altered by the subroutine.
!
! NUM is an INTEGER variable that need not be set by the user.
!   On successful exit, NUM will be the number of entries on the
!   diagonal of the permuted matrix.
!   If NUM < N, the matrix is structurally singular.
!
! PERM is an INTEGER array of length M that need not be set by the
!   user. On successful exit, PERM can be interpreted in any of the
!   following ways:
!
!   1. If M=N, PERM contains the column permutation.
!      Column PERM(I) of the original matrix is column I in the
!      permuted matrix, I=1..N.
!      (This was the definition of parameter CPERM in versions of
!      MC64A before version 1.2b)
!
!   2. If M>=N, PERM contains the row permutation.
!      Row I of the original matrix is row ABS(PERM(I)) in the
!      permuted matrix, I=1..M.
!      The rows where PERM(I) is positive constitute an N by N matrix
!      the scaled version of which has ones on the diagonal.
!
! LIW is an INTEGER variable that must be set by the user to
!   the dimension of array IW. It is not altered by the subroutine.
!   Restriction:
!     JOB = 1 :  LIW >=  4N +  M
!     JOB = 2 :  LIW >=  2N + 2M
!     JOB = 3 :  LIW >=  8N + 2M + NE
!     JOB = 4 :  LIW >=  3N + 2M
!     JOB = 5 :  LIW >=  3N + 2M
!
! IW is an INTEGER array of length LIW that is used for workspace.
!
! LDW is an INTEGER variable that must be set by the user to the
!   dimension of array DW. It is not altered by the subroutine.
!   Restriction:
!     JOB = 1 :  LDW not used
!     JOB = 2 :  LDW >=      M
!     JOB = 3 :  LDW >=          NE
!     JOB = 4 :  LDW >=     2M + NE
!     JOB = 5 :  LDW >= N + 2M + NE
!
! DW is a REAL array of length LDW used for workspace.
!   If JOB = 5, on return, DW(i) contains u_i, i=1..M, and
!   DW(M+j) contains v_j, j=1..N.
!
! ICNTL is an INTEGER array of length NICNTL.
!   Its components control the output of MC64A and must be set by the
!   user before calling MC64A. They are not altered by the subroutine.
!
!   ICNTL(1) must be set to specify the output stream for
!   error messages. If ICNTL(1) < 0, messages are suppressed.
!
!   ICNTL(2) must be set by the user to specify the output stream for
!   warning messages. If ICNTL(2) < 0, messages are suppressed.
!
!   ICNTL(3) must be set by the user to specify the output stream for
!   diagnostic messages. If ICNTL(3) < 0, messages are suppressed.
!
!   ICNTL(4) must be set by the user to a value other than 0 to avoid
!   checking of the input data.  Setting ICNTL(4) to any
!   other will avoid the checks but is likely to cause problems
!   later if out-of-range indices or duplicates are present.
!   The user should set ICNTL(4) nonzero, if the data is known not
!   to contain such problems. The code will exhibit undefined
!   behaviour in case data checking is not done and the
!   input data does not satisfy the restrictions as listed
!   elsewhere.
!
!   ICNTL(5) must be set by the user to control the printing of
!   diagnostic messages.
!   If ICNTL(5) <= 0, no messages are output.
!   If ICNTL(5) = 1, only error messages are output.
!   If ICNTL(5) = 2, error and warning messages output.
!   If ICNTL(5) = 3, as for 2 plus scalar parameters, the first 
!   ten entries of array parameters, and the control parameters on 
!   the first entry.
!   If ICNTL(5) > 3, full data will be printed on entry and exit.
!
! CNTL is a REAL array of length NCNTL.
!   Its components control the output of MC64A and must be set by the
!   user before calling MC64A. They are not altered by the subroutine.
!
!   CNTL(1) must be set to specify the relaxation parameter.
!   It is used by MC64 only if JOB = 3,4,5.
!   It must be set to a non-negative value (usually close to zero).
!   If CNTL(1) < 0.0, it is treated as 0.0.
!
!   CNTL(1) is a relaxation parameter. A positive value will lead to
!   matchings computed by MC64A that are not optimal/maximal in some
!   sense but only nearly so. However, these non-optimal matchings are
!   often computed more quickly. Appropriate values for CNTL(1) are
!   problem dependent but usually slightly larger than 0.0.
!
!
! INFO is an INTEGER array of length NINFO which need not be set by the
!   user. INFO(1) is set non-negative to indicate success. A negative
!   value is returned if an error occurred, a positive value if a
!   warning occurred. INFO(2) holds further information on the error.
!   On exit from the subroutine, INFO(1) will take one of the
!   following values:
!    0 : successful entry (for structurally nonsingular matrix).
!   +1 : successful entry (for structurally singular matrix).
!   +2 : the returned scaling factors are large and may cause
!        overflow when used to scale the matrix.
!        (For JOB = 4,5 entries only.)
!   +4 : CNTL(1) is negative and treated as zero.
!   -1 : JOB < 1 or JOB > 5.  Value of JOB held in INFO(2).
!   -2 : N < 1.  Value of invalid N held in INFO(2).
!   -3 : NE < 1.  Value of NE held in INFO(2).
!   -4 : M < N. Value of M held in INFO(2).
!   -6 : entries are found whose row indices are out of range. INFO(2)
!        contains the position in arrays A/IRN in which first entry is found.
!        (This value can be returned only if ICNTL(4) was set to zero.)
!   -7 : repeated entries are found. INFO(2) contains the position in arrays
!        A/IRN in which first entry is found.
!        (This value can be returned only if ICNTL(4) was set to zero.)
!
!   A return with one of the values INFO(1)=+3,+5,+6,+7 is also possible
!   These values are combinations of the above warnings (+1,+2,+4) and
!   correspond to the sum of the constituent warnings.
!
!   INFO(3) to INFO(NINFO) are not currently used and are set to zero
!        by the routine.
!
      REAL(WP) DW (LDW), CNTL (NCNTL)
! References:
!  [1] I. S. Duff, (1981),
!      "Algorithm 575. Permutations for a zero-free diagonal",
!      ACM Trans. Math. Software 7(3), 387-390.
!  [2] I. S. Duff and J. Koster, (1998),
!      "The design and use of algorithms for permuting large
!      entries to the diagonal of sparse matrices",
!      SIAM J. Matrix Anal. Appl., vol. 20, no. 4, pp. 889-901.
!  [3] I. S. Duff and J. Koster, (2001),
!      "On algorithms for permuting large entries to the diagonal
!      of sparse matrices",
!      SIAM J. Matrix Anal. Appl., vol. 22, no. 4, pp. 973-996.

! Local variables and parameters
      INTEGER I, J, K, WARN1, WARN2, WARN4
      REAL(WP) FACT, RINF
      REAL(WP), PARAMETER :: ZERO = 0.0_WP
      REAL(WP), PARAMETER :: ONE = 1.0_WP
! Intrinsic functions
      INTRINSIC ABS, LOG

! Set RINF to largest positive real number (infinity)
      RINF = HUGE (RINF)
      RINF = RINF / N
      WARN1 = 0
      WARN2 = 0
      WARN4 = 0
! Check value of JOB
      IF (ICNTL(4).EQ.0) THEN
! Check input data
        IF (JOB.LT.1.OR.JOB.GT.5) THEN
           INFO (1) = - 1
           INFO (2) = JOB
           IF (ICNTL(1).GE.0 .AND. ICNTL(5).GT.0)   &
             WRITE (ICNTL(1) , 9001) INFO (1) , 'JOB', JOB
           GOTO 99
        ENDIF
! Check value of N
        IF (N.LT.1) THEN
           INFO (1) = - 2
           INFO (2) = N
           IF (ICNTL(1).GE.0 .AND. ICNTL(5).GT.0)   &
             WRITE (ICNTL(1) , 9001) INFO (1) , 'N', N
           GOTO 99
        ENDIF
! Check value of NE
        IF (NE.LT.1) THEN
           INFO (1) = - 3
           INFO (2) = NE
           IF (ICNTL(1) .GE.0 .AND. ICNTL(5).GT.0)   &
             WRITE (ICNTL(1) , 9001) INFO (1) , 'NE', NE
           GOTO 99
        ENDIF
! Check value of M
        IF (M.LT.N) THEN
           INFO (1) = - 4
           INFO (2) = M
           IF (ICNTL(1) .GE.0 .AND. ICNTL(5).GT.0)   &
             WRITE (ICNTL(1) , 9001) INFO (1) , 'M', M
           GOTO 99
        ENDIF
! Check LIW
!     IF (JOB.EQ.1) K = 4*N +   M
!     IF (JOB.EQ.2) K = 2*N + 2*M
!     IF (JOB.EQ.3) K = 8*N + 2*M + NE
!     IF (JOB.EQ.4) K = 3*N + 2*M
!     IF (JOB.EQ.5) K = 3*N + 2*M
!     IF (JOB.EQ.6) K = 3*N + 2*M + NE
!     IF (LIW.LT.K) THEN
!       INFO(1) = -5
!       INFO(2) = K
!       IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9004) INFO(1),K
!       GO TO 99
!     ENDIF
! Check LDW; If JOB = 1, do not check
!     IF (JOB.GT.1) THEN
!       IF (JOB.EQ.2) K =       M
!       IF (JOB.EQ.3) K =           NE
!       IF (JOB.EQ.4) K =     2*M + NE
!       IF (JOB.EQ.5) K = N + 2*M + NE
!       IF (JOB.EQ.6) K = N + 3*M + NE
!       IF (LDW.LT.K) THEN
!         INFO(1) = -5
!         INFO(2) = K
!         IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9005) INFO(1),K
!         GO TO 99
!       ENDIF
!     ENDIF
      ENDIF
      IF (ICNTL(4) .EQ.0) THEN
! Check row indices. Use IW(1:M) as workspace
         DO 3 I = 1, M
            IW (I) = 0
    3    END DO
         DO 6 J = 1, N
            DO 4 K = IP (J), IP (J + 1) - 1
               I = IRN (K)
! Check for row indices that are out of range
               IF (I.LT.1.OR.I.GT.M) THEN
                  INFO (1) = - 6
                  INFO (2) = K
                  IF (ICNTL(1).GE.0 .AND. ICNTL(5).GT.0)   &
                    WRITE (ICNTL(1), 9006) INFO (1), J, I
                  GOTO 99
               ENDIF
! Check for repeated row indices within a column
               IF (IW (I) .EQ.J) THEN
                  INFO (1) = - 7
                  INFO (2) = K
                  IF (ICNTL(1).GE.0 .AND. ICNTL(5).GT.0)   &
                    WRITE (ICNTL(1), 9007) INFO (1), J, I
                  GOTO 99
               ELSE
                  IW (I) = J
               ENDIF
    4       END DO
    6    END DO
      ENDIF

! Print diagnostics on input
      IF (ICNTL(5).GT.2) THEN
        IF (ICNTL(3) .GE.0) THEN
          WRITE (ICNTL(3), 9020) JOB, M, N, NE
          IF (ICNTL(5).EQ.3) THEN
            WRITE (ICNTL(3), 9021) (IP (J), J = 1, MIN (10, N + 1) )
            WRITE (ICNTL(3), 9022) (IRN (J), J = 1, MIN (10, NE) )
            IF (JOB.GT.1) WRITE (ICNTL(3), 9023) (A (J), J = 1, MIN &
            (10, NE) )
          ELSE
            WRITE (ICNTL(3), 9021) (IP (J), J = 1, N + 1)
            WRITE (ICNTL(3), 9022) (IRN (J), J = 1, NE)
            IF (JOB.GT.1) WRITE (ICNTL(3), 9023) (A (J), J = 1, NE)
          ENDIF
          WRITE (ICNTL(3), 9024) (ICNTL(J), J = 1, NICNTL)
          WRITE (ICNTL(3), 9025) CNTL(1)
        ENDIF
      ENDIF

! Set components of INFO to zero
      DO I = 1, NINFO
         INFO (I) = 0
      END DO

! Compute maximum matching
      IF (JOB.EQ.1) THEN
! Put length of column J in IW(J)
         DO J = 1, N
            IW (J) = IP (J + 1) - IP (J)
         END DO
! IW(N+1:3N+M+N) is workspace
         CALL MC64_HSL_Z (M, N, IRN, NE, IP, IW (1), PERM, NUM, IW (N +&
         1), IW (2 * N + 1), IW (3 * N + 1), IW (3 * N + M + 1) )
         GOTO 90
      ENDIF

! Compute bottleneck matching
      IF (JOB.EQ.2) THEN
! Pass CNTL(1) to MC64B through DW(1)
         DW (1) = MAX (ZERO, CNTL (1) )
! IW(1:2N+2M), DW(1:M) are workspaces
         CALL MC64_HSL_B (M, N, NE, IP, IRN, A, PERM, NUM, IW (1),     &
         IW (N + 1), IW (2 * N + 1), IW (2 * N + M + 1), DW, RINF)
         GOTO 90
      ENDIF

! Compute bottleneck matching
      IF (JOB.EQ.3) THEN
! Copy IRN(K) into IW(K), ABS(A(K)) into DW(K), K=1..NE
         DO K = 1, NE
            IW (K) = IRN (K)
            DW (K) = ABS (A (K) )
         END DO
! Sort entries in each column by decreasing value.
         CALL MC64R (N, NE, IP, IW, DW)
! Pass CNTL(1) to MC64S through FACT
         FACT = MAX (ZERO, CNTL (1) )
! IW(NE+1:NE+5N+M+(3N+M)) is workspace
         CALL MC64_HSL_S (M, N, NE, IP, IW (1), DW, PERM, NUM, IW (    &
         NE+1), IW (NE+N + 1), IW (NE+2 * N + 1), IW (NE+3 * N + 1),    &
         IW (NE+4 * N + 1), IW (NE+5 * N + 1), IW (NE+5 * N + M + 1),   &
         FACT, RINF)
         GOTO 90
      ENDIF

      IF (JOB.EQ.4) THEN
         DO J = 1, N
            FACT = ZERO
            DO K = IP (J), IP (J + 1) - 1
               IF (ABS (A (K) ) .GT.FACT) FACT = ABS (A (K) )
            END DO
            DO K = IP (J), IP (J + 1) - 1
               DW (2 * M + K) = FACT - ABS (A (K) )
            END DO
         END DO
! B = DW(2M+1:2M+NE); IW(1:3N+2M) and DW(1:2M) are workspaces
! Pass CNTL(1) to MC64W through DW(1)
! Pass JOB to MC64W through IW(1)
         DW (1) = MAX (ZERO, CNTL (1) )
         IW (1) = JOB
! Call MC64W
         CALL MC64_HSL_W (M, N, NE, IP, IRN, DW (2 * M + 1), PERM, NUM,&
         IW (1), IW (N + 1), IW (2 * N + 1), IW (3 * N + 1), IW (3 * N +&
         M + 1), DW (1), DW (M + 1), RINF)
         GOTO 90
      ENDIF

      IF (JOB.EQ.5.or.JOB.EQ.6) THEN
         IF (JOB.EQ.5) THEN
            DO 75 J = 1, N
               FACT = ZERO
               DO K = IP (J), IP (J + 1) - 1
                  DW (2 * M + N + K) = ABS (A (K) )
                  IF (DW (2 * M + N + K) .GT.FACT) FACT = DW (2 * M + N &
                  + K)
               END DO
               DW (2 * M + J) = FACT
!CC Significant change made here so that column with single
!   zero gets set to RINF and not 1.0
               IF (FACT.NE.ZERO) THEN
                  FACT = LOG (FACT)
               ELSE
                  FACT = RINF
               ENDIF
               DO K = IP (J), IP (J + 1) - 1
                  IF (DW (2 * M + N + K) .NE.ZERO) THEN
                     DW (2 * M + N + K) = FACT - LOG (DW (2 * M + N + K)&
                     )
                  ELSE
!                  write(*,*) 'set diag to ',RINF
                     DW (2 * M + N + K) = RINF
!5.0D+14
!*RINF/(N+1)
                  ENDIF
               END DO
!           ELSE
!             DO 71 K = IP(J),IP(J+1)-1
!               DW(2*M+N+K) = ONE
!  71         CONTINUE
!           ENDIF
   75       END DO
         ENDIF

!       IF (JOB.EQ.6) THEN
!         DO 175 K = 1,NE
!           IW(3*N+2*M+K) = IRN(K)
!           DW(2*M+N+K) = ABS(A(K))
! 175     CONTINUE
!         DO 61 I = 1,M
!           DW(2*M+N+NE+I) = ZERO
!  61     CONTINUE
!         DO 63 J = 1,N
!           DO 62 K = IP(J),IP(J+1)-1
!             I = IRN(K)
!             IF (DW(2*M+N+K).GT.DW(2*M+N+NE+I)) THEN
!               DW(2*M+N+NE+I) = DW(2*M+N+K)
!             ENDIF
!  62       CONTINUE
!  63     CONTINUE
!         DO 64 I = 1,M
!           IF (DW(2*M+N+NE+I).NE.ZERO) THEN
!             DW(2*M+N+NE+I) = 1/DW(2*M+N+NE+I)
!           ENDIF
!  64     CONTINUE
!         DO 66 J = 1,N
!           DO 65 K = IP(J),IP(J+1)-1
!             I = IRN(K)
!             DW(2*M+N+K) = DW(2*M+N+NE+I) * DW(2*M+N+K)
!  65       CONTINUE
!  66     CONTINUE
!         CALL MC64R(N,NE,IP,IW(3*N+2*M+1),DW(2*M+N+1))
!         DO 176 J = 1,N
!           IF (IP(J).NE.IP(J+1)) THEN
!             FACT = DW(2*M+N+IP(J))
!           ELSE
!             FACT = ZERO
!           ENDIF
!           DW(2*M+J) = FACT
!           IF (FACT.NE.ZERO) THEN
!             FACT = LOG(FACT)
!             DO 170 K = IP(J),IP(J+1)-1
!               IF (DW(2*M+N+K).NE.ZERO) THEN
!                 DW(2*M+N+K) = FACT - LOG(DW(2*M+N+K))
!               ELSE
!                  write(*,*) 'set diag to ',RINF
!                 DW(2*M+N+K) = RINF
!5.0D+14
!0.5* RINF/(N+1)
!               ENDIF
! 170         CONTINUE
!           ELSE
!             DO 171 K = IP(J),IP(J+1)-1
!               DW(2*M+N+K) = ONE
! 171         CONTINUE
!           ENDIF
! 176     CONTINUE
!       ENDIF
! Pass CNTL(1) to MC64W through DW(1)
! Pass JOB to MC64W through IW(1)
         DW (1) = MAX (ZERO, CNTL (1) )
         IW (1) = JOB
! Call MC64W
         IF (JOB.EQ.5) THEN
            CALL MC64_HSL_W (M, N, NE, IP, IRN, DW (2 * M + N + 1),    &
            PERM, NUM, IW (1), IW (N + 1), IW (2 * N + 1), IW (3 * N +  &
            1), IW (3 * N + M + 1), DW (1), DW (M + 1), RINF)
         ENDIF

!       IF (JOB.EQ.6) THEN
!         CALL MC64_HSL_W(M,N,NE,IP,IW(3*N+2*M+1),DW(2*M+N+1),PERM,NUM,
!    &         IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(3*N+M+1),
!    &         DW(1),DW(M+1),RINF)
!       ENDIF
!       IF (JOB.EQ.6) THEN
!         DO 79 I = 1,M
!           IF (DW(2*M+N+NE+I).NE.0) THEN
!             DW(I) = DW(I) + LOG(DW(2*M+N+NE+I))
!           ENDIF
!  79     CONTINUE
!       ENDIF
         IF (NUM.EQ.N) THEN
            DO    J = 1, N
               IF (DW (2 * M + J) .NE.ZERO) THEN
                  DW (M + J) = DW (M + J) - LOG (DW (2 * M + J) )
               ELSE
                  DW (M + J) = ZERO
               ENDIF
            END DO
         ENDIF
! Check size of row and column scaling factors
         FACT = 0.5 * LOG (RINF)
         DO J = 1, N
            IF (DW (M + J) .LT.FACT) CYCLE
            WARN2 = 2
! Scaling factor is large, return with warning
            GOTO 90
         END DO
         DO I = 1, M
            IF (DW (I) .LT.FACT) CYCLE
            WARN2 = 2
! Scaling factor is large, return with warning
            GOTO 90
         END DO
!       GO TO 90
      ENDIF

! If matrix is structurally singular, return with warning
   90 IF (NUM.LT.N) WARN1 = 1

! If CNTL(1) is negative and treated as zero, return with warning
      IF (JOB.EQ.4.OR.JOB.EQ.5.OR.JOB.EQ.6) THEN
         IF (CNTL (1) .LT.ZERO) WARN4 = 4
      ENDIF

! Set warning flag and print warnings (only if no errors were found)
      IF (INFO (1) .EQ.0) THEN
         INFO (1) = WARN1 + WARN2 + WARN4
         IF (INFO (1).GT.0 .AND. ICNTL(2).GE.0 .AND. ICNTL(5).GT.1) THEN
           WRITE (ICNTL(2), 9010) INFO (1)
           IF (WARN1.EQ.1) WRITE (ICNTL(2), 9011)
           IF (WARN2.EQ.2) WRITE (ICNTL(2), 9012)
           IF (WARN4.EQ.4) WRITE (ICNTL(2), 9014)
         ENDIF
      ENDIF
! Print diagnostics on output
      IF (ICNTL(5).GT.2) THEN
        IF (ICNTL(3).GE.0) THEN
          WRITE (ICNTL(3), 9030) (INFO (J), J = 1, 2)
          WRITE (ICNTL(3), 9031) NUM
          IF (ICNTL(5) .EQ.3) THEN
            WRITE (ICNTL(3), 9032) (PERM (J), J = 1, MIN (10, M) )
            IF (JOB.EQ.5.OR.JOB.EQ.6) THEN
              WRITE (ICNTL(3), 9033) (DW (J), J = 1, MIN (10, M) )
              WRITE (ICNTL(3), 9034) (DW (M + J), J = 1, MIN (10, N))
            ENDIF
          ELSE
            WRITE (ICNTL(3), 9032) (PERM (J), J = 1, M)
            IF (JOB.EQ.5.OR.JOB.EQ.6) THEN
              WRITE (ICNTL(3), 9033) (DW (J), J = 1, M)
              WRITE (ICNTL(3), 9034) (DW (M + J), J = 1, N)
            ENDIF
          ENDIF
        ENDIF
      ENDIF

! Return from subroutine.
   99 RETURN

 9001 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2,               &
     &        ' because ',(A),' = ',I10)
!9004 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/               &
!    &        '        LIW too small, must be at least ',I8)
!9005 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/               &
!    &        '        LDW too small, must be at least ',I8)
 9006 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/               &
     &        '        Column ',I8,                                     &
     &        ' contains an entry with invalid row index ',I8)
 9007 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/               &
     &        '        Column ',I8,                                     &
     &        ' contains two or more entries with row index ',I8)

 9010 FORMAT (' ****** Warning from MC64A/AD. INFO(1) = ',I2)
 9011 FORMAT ('        - The matrix is structurally singular.')
 9012 FORMAT ('        - Some scaling factors may be too large.')
 9014 FORMAT ('        - CNTL(1) is negative and was treated as zero.')

 9020 FORMAT (' ****** Input parameters for MC64A:'/                   &
     &        ' JOB =',I10/' M   =',I10/' N   =',I10/' NE  =',I10)
 9021 FORMAT (' IP(1:N+1)  = ',8I8/(14X,8I8))
 9022 FORMAT (' IRN(1:NE)  = ',8I8/(14X,8I8))
 9023 FORMAT (' A(1:NE)    = ',4(1PD14.4)/(14X,4(1PD14.4)))
 9024 FORMAT (' ICNTL(1:10)= ',8I8/(14X,2I8))
 9025 FORMAT (' CNTL(1)    = ',1PD14.4)
 9030 FORMAT (' ****** Output parameters for MC64A:'/                  &
     &        ' INFO(1:2)  = ',2I8)
 9031 FORMAT (' NUM        = ',I8)
 9032 FORMAT (' PERM(1:M)  = ',8I8/(14X,8I8))
 9033 FORMAT (' DW(1:M)    = ',5(F11.3)/(14X,5(F11.3)))
 9034 FORMAT (' DW(M+1:M+N)= ',5(F11.3)/(14X,5(F11.3)))
END SUBROUTINE MC64_HSL_A

!**********************************************************************
!CCCC LAST UPDATE Tue Nov 26 03:20:26 MET 2002
SUBROUTINE MC64_HSL_B (M, N, NE, IP, IRN, A, IPERM, NUM, JPERM,  &
      PR, Q, L, D, RINF)
      IMPLICIT NONE
!
! *** Copyright (c) 2002  I.S. Duff and J. Koster                   ***
! *** Although every effort has been made to ensure robustness and  ***
! *** reliability of the subroutines in this MC64 suite, we         ***
! *** disclaim any liability arising through the use or misuse of   ***
! *** any of the subroutines.                                       ***
!
      INTEGER M, N, NE, NUM
      INTEGER IP (N + 1), IRN (NE), IPERM (M), JPERM (N), PR (N),       &
      Q (M), L (M)
      REAL(WP) A (NE)
      REAL(WP) D (M), RINF

! N, NE, IP, IRN are described in MC64A.
! A is a REAL array of length NE.
!   A(K), K=1..NE, must be set to the value of the entry
!   that corresponds to IRN(K). It is not altered.
! IPERM is an INTEGER array of length M. On exit, it contains the
!    matching: IPERM(I) = 0 or row I is matched to column IPERM(I).
! NUM is INTEGER variable. On exit, it contains the cardinality of the
!    matching stored in IPERM.
! D is a REAL work array of length M.
!    On entry, D(1) contains the relaxation parameter RLX.
! RINF is the largest positive real number

! Local variables
      INTEGER I, II, J, JJ, JORD, Q0, QLEN, IDUM, JDUM, ISP, JSP, K, KK,&
      KK1, KK2, I0, UP, LOW, LPOS
      REAL(WP) CSP, DI, DNEW, DQ0, AI, A0, BV, TBV, RLX
! Local parameters
      REAL(WP), PARAMETER :: ZERO=0.0_WP
      REAL(WP), PARAMETER :: MINONE=-1.0_WP
      INTRINSIC ABS, MIN
! External subroutines and/or functions
      EXTERNAL MC64D, MC64E, MC64F

! Initialize variables to eliminate copmiler warnings
      I0 = -HUGE(I0); ISP = -HUGE(ISP); JSP = -HUGE(JSP)

      RLX = D (1)
! Initialization
      NUM = 0
      BV = RINF
      DO K = 1, N
         JPERM (K) = 0
         PR (K) = IP (K)
      END DO
      DO K = 1, M
         IPERM (K) = 0
         D (K) = ZERO
      END DO

      DO J = 1, N
         A0 = MINONE
         DO K = IP (J), IP (J + 1) - 1
            I = IRN (K)
            AI = ABS (A (K) )
            IF (AI.GT.D (I) ) D (I) = AI
            IF (JPERM (J) .NE.0) CYCLE
            IF (AI.GE.BV) THEN
               A0 = BV
               IF (IPERM (I) .NE.0) CYCLE
               JPERM (J) = I
               IPERM (I) = J
               NUM = NUM + 1
            ELSE
               IF (AI.LE.A0) CYCLE
               A0 = AI
               I0 = I
            ENDIF
         END DO
         IF (A0.NE.MINONE.AND.A0.LT.BV) THEN
            BV = A0
            IF (IPERM (I0) .NE.0) CYCLE
            IPERM (I0) = J
            JPERM (J) = I0
            NUM = NUM + 1
         ENDIF
      END DO

      IF (M.EQ.N) THEN
! Update BV with smallest of all the largest maximum absolute values
! of the rows. D(I) contains the largest absolute value in row I.
         DO I = 1, M
            BV = MIN (BV, D (I) )
         END DO
      ENDIF

! Shortcut if all columns are matched at this stage.
      IF (NUM.EQ.N) GOTO 1000

! Rescan unassigned columns; improve initial assignment
      DO J = 1, N
         IF (JPERM (J) .NE.0) CYCLE
         DO K = IP (J), IP (J + 1) - 1
            I = IRN (K)
            AI = ABS (A (K) )
            IF (AI.LT.BV) CYCLE
            IF (IPERM (I) .EQ.0) GOTO 90
            JJ = IPERM (I)
            KK1 = PR (JJ)
            KK2 = IP (JJ + 1) - 1
            IF (KK1.GT.KK2) CYCLE
            DO KK = KK1, KK2
               II = IRN (KK)
               IF (IPERM (II) .NE.0) CYCLE
               IF (ABS (A (KK) ) .GE.BV) GOTO 80
            END DO
            PR (JJ) = KK2 + 1
         END DO
         CYCLE
   80    JPERM (JJ) = II
         IPERM (II) = JJ
         PR (JJ) = KK + 1
   90    NUM = NUM + 1
         JPERM (J) = I
         IPERM (I) = J
         PR (J) = K + 1
      END DO

! Shortcut if all columns are matched at this stage.
      IF (NUM.EQ.N) GOTO 1000

! Prepare for main loop
      DO I = 1, M
         D (I) = MINONE
         L (I) = 0
      END DO
! TBV is a relaxed value of BV (ie TBV is slightly smaller than BV).
      TBV = BV * (1 - RLX)

! Main loop ... each pass round this loop is similar to Dijkstra's
! algorithm for solving the single source shortest path problem

      DO JORD = 1, N

         IF (JPERM (JORD) .NE.0) CYCLE
         QLEN = 0
         LOW = M + 1
         UP = M + 1
! CSP is cost of shortest path to any unassigned row
! ISP is matrix position of unassigned row element in shortest path
! JSP is column index of unassigned row element in shortest path
         CSP = MINONE
! Build shortest path tree starting from unassigned column JORD
         J = JORD
         PR (J) = - 1

! Scan column J
         DO K = IP (J), IP (J + 1) - 1
            I = IRN (K)
            DNEW = ABS (A (K) )
            IF (CSP.GE.DNEW) CYCLE
            IF (IPERM (I) .EQ.0) THEN
! Row I is unassigned; update shortest path info
               CSP = DNEW
               ISP = I
               JSP = J
               IF (CSP.GE.TBV) GOTO 160
            ELSE
               D (I) = DNEW
               IF (DNEW.GE.TBV) THEN
! Add row I to Q2
                  LOW = LOW - 1
                  Q (LOW) = I
               ELSE
! Add row I to Q, and push it
                  QLEN = QLEN + 1
                  L (I) = QLEN
                  CALL MC64D (I, M, Q, D, L, 1)
               ENDIF
               JJ = IPERM (I)
               PR (JJ) = J
            ENDIF
         END DO

         DO JDUM = 1, NUM
! If Q2 is empty, extract new rows from Q
            IF (LOW.EQ.UP) THEN
               IF (QLEN.EQ.0) EXIT
               I = Q (1)
               IF (CSP.GE.D (I) ) EXIT
               BV = D (I)
               TBV = BV * (1 - RLX)
               DO IDUM = 1, M
                  CALL MC64E (QLEN, M, Q, D, L, 1)
                  L (I) = 0
                  LOW = LOW - 1
                  Q (LOW) = I
                  IF (QLEN.EQ.0) EXIT
                  I = Q (1)
                  IF (D (I) .LT.TBV) EXIT
               END DO
! End of dummy loop; this point is never reached
            ENDIF
! Move row Q0
            UP = UP - 1
            Q0 = Q (UP)
            DQ0 = D (Q0)
            L (Q0) = UP
! Scan column that matches with row Q0
            J = IPERM (Q0)
            DO K = IP (J), IP (J + 1) - 1
               I = IRN (K)
! Update D(I); only if row I is not marked
               IF (L (I) .GE.UP) CYCLE
               DNEW = MIN (DQ0, ABS (A (K) ) )
               IF (CSP.GE.DNEW) CYCLE
               IF (IPERM (I) .EQ.0) THEN
! Row I is unassigned; update shortest path info
                  CSP = DNEW
                  ISP = I
                  JSP = J
                  IF (CSP.GE.TBV) GOTO 160
               ELSE
                  DI = D (I)
                  IF (DI.GE.TBV.OR.DI.GE.DNEW) CYCLE
                  D (I) = DNEW
                  IF (DNEW.GE.TBV) THEN
! Delete row I from Q (if necessary); add row I to Q2
                     IF (DI.NE.MINONE) THEN
                        LPOS = L (I)
                        CALL MC64F (LPOS, QLEN, M, Q, D, L, 1)
                     ENDIF
                     L (I) = 0
                     LOW = LOW - 1
                     Q (LOW) = I
                  ELSE
! Add row I to Q (if necessary); push row I up Q
                     IF (DI.EQ.MINONE) THEN
                        QLEN = QLEN + 1
                        L (I) = QLEN
                     ENDIF
                     CALL MC64D (I, M, Q, D, L, 1)
                  ENDIF
! Update tree
                  JJ = IPERM (I)
                  PR (JJ) = J
               ENDIF
            END DO
         END DO

! If CSP = MINONE, no augmenting path is found
  160    IF (CSP.EQ.MINONE) GOTO 190
! Update bottleneck value
         BV = MIN (BV, CSP)
         TBV = BV * (1 - RLX)
! Find augmenting path by tracing backward in PR; update IPERM,JPERM
         NUM = NUM + 1
         I = ISP
         J = JSP
         DO JDUM = 1, NUM + 1
            I0 = JPERM (J)
            JPERM (J) = I
            IPERM (I) = J
            J = PR (J)
            IF (J.EQ. - 1) EXIT
            I = I0
         END DO
! End of dummy loop; this point is never reached
  190    DO 191 KK = UP, M
            I = Q (KK)
            D (I) = MINONE
            L (I) = 0
  191    END DO
         DO 192 KK = LOW, UP - 1
            I = Q (KK)
            D (I) = MINONE
  192    END DO
         DO 193 KK = 1, QLEN
            I = Q (KK)
            D (I) = MINONE
            L (I) = 0
  193    END DO

      END DO
! End of main loop

! BV is now bottleneck value of final matching

! IPERM is complete if M = N and NUM = N
 1000 IF (M.EQ.N.and.NUM.EQ.N) GOTO 2000

! Complete IPERM; L, JPERM are work arrays
      CALL MC64_HSL_X (M, N, IPERM, L, JPERM)

 2000 RETURN
END SUBROUTINE MC64_HSL_B

!**********************************************************************
!CCCC LAST UPDATE Tue Nov 26 03:20:26 MET 2002
SUBROUTINE MC64_HSL_S (M, N, NE, IP, IRN, A, IPERM, NUMX, W, LEN,&
      LENL, LENH, FC, IW, IW4, RLX, RINF)
      IMPLICIT NONE
!
! *** Copyright (c) 2002  I.S. Duff and J. Koster                   ***
! *** Although every effort has been made to ensure robustness and  ***
! *** reliability of the subroutines in this MC64 suite, we         ***
! *** disclaim any liability arising through the use or misuse of   ***
! *** any of the subroutines.                                       ***
!
      INTEGER M, N, NE, NUMX
      INTEGER IP (N + 1), IRN (NE), IPERM (M), W (N), LEN (N), LENL (N),&
      LENH (N), FC (N), IW (M), IW4 (3 * N + M)
      REAL(WP) A (NE), RLX, RINF

! M, N, NE, IP, IRN, are described in MC64A.
! A is a REAL array of length NE.
!   A(K), K=1..NE, must be set to the value of the entry that
!   corresponds to IRN(k). The entries in each column must be
!   non-negative and ordered by decreasing value.
! IPERM is an INTEGER array of length M. On exit, it contains the
!   bottleneck matching: IPERM(I) - 0 or row I is matched to column
!   IPERM(I).
! NUMX is an INTEGER variable. On exit, it contains the cardinality
!   of the matching stored in IPERM.

! FC is an integer array of length N that contains the list of
!   unmatched columns.
! LEN(J), LENL(J), LENH(J) are integer arrays of length N that point
!   to entries in matrix column J.
!   In the matrix defined by the column parts IP(J)+LENL(J) we know
!   a matching does not exist; in the matrix defined by the column
!   parts IP(J)+LENH(J) we know one exists.
!   LEN(J) lies between LENL(J) and LENH(J) and determines the matrix
!   that is tested for a maximum matching.
! W is an integer array of length N and contains the indices of the
!   columns for which LENL /= LENH.
! WLEN is number of indices stored in array W.
! IW is integer work array of length M.
! IW4 is integer work array of length 3N+M used by MC64U.
!
! RLX is a REAL variable. It is a relaxation
!   parameter for finding the optimal matching.
!
! RINF is the largest positive real number

      INTEGER NUM, NVAL, WLEN, II, I, J, K, L, CNT, MOD, IDUM1, IDUM2,  &
      IDUM3
      REAL(WP) BVAL, BMIN, BMAX
! External subroutines and/or functions
      EXTERNAL MC64Q
! Intrinsic functions

! BMIN and BMAX are such that a maximum matching exists for the input
!   matrix in which all entries smaller than BMIN are dropped.
!   For BMAX, a maximum matching does not exist.
! BVAL is a value between BMIN and BMAX.
! CNT is the number of calls made to MC64U so far.
! NUM is the cardinality of last matching found.

! Compute a first maximum matching from scratch on whole matrix.
      DO J = 1, N
         FC (J) = J
         LEN (J) = IP (J + 1) - IP (J)
      END DO
      DO I = 1, M
         IW (I) = 0
      END DO
! The first call to MC64U
      CNT = 1
      MOD = 1
      NUMX = 0
      CALL MC64_HSL_U (CNT, MOD, M, N, IRN, NE, IP, LEN, FC, IW, NUMX, &
      N, IW4 (1), IW4 (N + 1), IW4 (2 * N + 1), IW4 (2 * N + M + 1) )

! IW contains a maximum matching of length NUMX.
      NUM = NUMX

      IF (NUM.NE.N) THEN
! Matrix is structurally singular
         BMAX = RINF
      ELSE
! Matrix is structurally nonsingular, NUM=NUMX=N;
! Set BMAX just above the smallest of all the maximum absolute
! values of the columns
         BMAX = RINF
         DO J = 1, N
            BVAL = 0.0
            DO K = IP (J), IP (J + 1) - 1
               IF (A (K) .GT.BVAL) BVAL = A (K)
            END DO
            IF (BVAL.LT.BMAX) BMAX = BVAL
         END DO
! ... should print warning if BMAX == RINF
         BMAX = 1.001 * BMAX
      ENDIF

! Initialize BVAL,BMIN
      BVAL = 0.0
      BMIN = 0.0
! Initialize LENL,LEN,LENH,W,WLEN according to BMAX.
! Set LEN(J), LENH(J) just after last entry in column J.
! Set LENL(J) just after last entry in column J with value >= BMAX.
      WLEN = 0
      DO J = 1, N
         L = IP (J + 1) - IP (J)
         LENH (J) = L
         LEN (J) = L
         DO K = IP (J), IP (J + 1) - 1
            IF (A (K) .LT.BMAX) GOTO 46
         END DO
! Column J is empty or all entries are >= BMAX
         K = IP (J + 1)
   46    LENL (J) = K - IP (J)
! Add J to W if LENL(J) /= LENH(J)
         IF (LENL (J) .EQ.L) CYCLE
         WLEN = WLEN + 1
         W (WLEN) = J
      END DO

! Main loop
      DO IDUM1 = 1, NE
         IF (NUM.EQ.NUMX) THEN
! We have a maximum matching in IW; store IW in IPERM
            DO  I = 1, M
               IPERM (I) = IW (I)
            END DO
! Keep going round this loop until matching IW is no longer maximum.
            DO 80 IDUM2 = 1, NE
               BMIN = BVAL
               IF (BMAX - BMIN.LE.RLX) GOTO 1000
! Find splitting value BVAL
               CALL MC64Q (IP, LENL, LEN, W, WLEN, A, NVAL, BVAL)
               IF (NVAL.LE.1) GOTO 1000
! Set LEN such that all matrix entries with value < BVAL are
! discarded. Store old LEN in LENH. Do this for all columns W(K).
! Each step, either K is incremented or WLEN is decremented.
               K = 1
               DO IDUM3 = 1, N
                  IF (K.GT.WLEN) EXIT
                  J = W (K)
                  DO II = IP (J) + LEN (J) - 1, IP (J) + LENL (J),   &
                  - 1
                     IF (A (II) .GE.BVAL) EXIT
                     I = IRN (II)
                     IF (IW (I) .NE.J) CYCLE
! Remove entry from matching
                     IW (I) = 0
                     NUM = NUM - 1
                     FC (N - NUM) = J
                  END DO
                  LENH (J) = LEN (J)
! IP(J)+LEN(J)-1 is last entry in column >= BVAL
                  LEN (J) = II - IP (J) + 1
! If LENH(J) = LENL(J), remove J from W
                  IF (LENL (J) .EQ.LENH (J) ) THEN
                     W (K) = W (WLEN)
                     WLEN = WLEN - 1
                  ELSE
                     K = K + 1
                  ENDIF
               END DO
               IF (NUM.LT.NUMX) EXIT
   80       END DO
! End of dummy loop; this point is never reached
! Set mode for next call to MC64U
            MOD = 1
         ELSE
! We do not have a maximum matching in IW.
            BMAX = BVAL
! BMIN is the bottleneck value of a maximum matching;
! for BMAX the matching is not maximum, so BMAX>BMIN
! and following condition is always false if RLX = 0.0
            IF (BMAX - BMIN.LE.RLX) GOTO 1000
! Find splitting value BVAL
            CALL MC64Q (IP, LEN, LENH, W, WLEN, A, NVAL, BVAL)
            IF (NVAL.EQ.0.OR.BVAL.EQ.BMIN) GOTO 1000
! Set LEN such that all matrix entries with value >= BVAL are
! inside matrix. Store old LEN in LENL. Do this for all columns W(K).
! Each step, either K is incremented or WLEN is decremented.
            K = 1
            DO 87 IDUM3 = 1, N
               IF (K.GT.WLEN) GOTO 88
               J = W (K)
               DO 85 II = IP (J) + LEN (J), IP (J) + LENH (J) - 1
                  IF (A (II) .LT.BVAL) GOTO 86
   85          END DO
   86          LENL (J) = LEN (J)
               LEN (J) = II - IP (J)
               IF (LENL (J) .EQ.LENH (J) ) THEN
                  W (K) = W (WLEN)
                  WLEN = WLEN - 1
               ELSE
                  K = K + 1
               ENDIF
   87       END DO
! End of dummy loop; this point is never reached
! Set mode for next call to MC64U
   88       MOD = 0
         ENDIF
         CNT = CNT + 1
         CALL MC64_HSL_U (CNT, MOD, M, N, IRN, NE, IP, LEN, FC, IW,    &
         NUM, NUMX, IW4 (1), IW4 (N + 1), IW4 (2 * N + 1), IW4 (2 * N + &
         M + 1) )

! IW contains maximum matching of length NUM
      END DO
! End of dummy loop; this point is never reached

! BMIN is now bottleneck value of final matching

! IPERM is complete if M = N and NUMX = N
 1000 IF (M.EQ.N.and.NUMX.EQ.N) GOTO 2000

! Complete IPERM; IW, W are work arrays
      CALL MC64_HSL_X (M, N, IPERM, IW, W)

 2000 RETURN
END SUBROUTINE MC64_HSL_S

!**********************************************************************
!CCCC LAST UPDATE Tue Nov 26 03:20:26 MET 2002
SUBROUTINE MC64_HSL_U (ID, MOD, M, N, IRN, LIRN, IP, LENC, FC,   &
      IPERM, NUM, NUMX, PR, ARP, CV, OUT)
      IMPLICIT NONE
!
! *** Copyright (c) 2002  I.S. Duff and J. Koster                   ***
! *** Although every effort has been made to ensure robustness and  ***
! *** reliability of the subroutines in this MC64 suite, we         ***
! *** disclaim any liability arising through the use or misuse of   ***
! *** any of the subroutines.                                       ***
!
      INTEGER ID, MOD, M, N, LIRN, NUM, NUMX
      INTEGER ARP (N), CV (M), IRN (LIRN), IP (N), FC (N), IPERM (M),   &
      LENC (N), OUT (N), PR (N)

! PR(J) is the previous column to J in the depth first search.
!   Array PR is used as workspace in the sorting algorithm.
! Elements (I,IPERM(I)) I=1,..,M are entries at the end of the
!   algorithm unless N assignments have not been made in which case
!   N-NUM pairs (I,IPERM(I)) will not be entries in the matrix.
! CV(I) is the most recent loop number (ID+JORD) at which row I
!   was visited.
! ARP(J) is the number of entries in column J which have been scanned
!   when looking for a cheap assignment.
! OUT(J) is one less than the number of entries in column J which have
!   not been scanned during one pass through the main loop.
! NUMX is maximum possible size of matching.

      INTEGER I, II, IN1, IN2, J, J1, JORD, K, KK, LAST, NFC, NUM0,     &
      NUM1, NUM2, ID0, ID1

! Initialize variables to eliminate copmiler warnings
      I = -1; II = -1

      IF (ID.EQ.1) THEN
! The first call to MC64U.
! Initialize CV and ARP; parameters MOD, NUMX are not accessed
         DO I = 1, M
            CV (I) = 0
         END DO
         DO J = 1, N
            ARP (J) = 0
         END DO
         NUM1 = N
         NUM2 = N
      ELSE
! Not the first call to MC64U.
! Re-initialize ARP if entries were deleted since last call to MC64U
         IF (MOD.EQ.1) THEN
            DO J = 1, N
               ARP (J) = 0
            END DO
         ENDIF
         NUM1 = NUMX
         NUM2 = N - NUMX
      ENDIF
      NUM0 = NUM

! NUM0 is size of input matching
! NUM1 is maximum possible size of matching
! NUM2 is maximum allowed number of unassigned rows/columns
! NUM is size of current matching

! Quick return if possible
!      IF (NUM.EQ.N) GO TO 199
! NFC is number of rows/columns that could not be assigned
      NFC = 0
! Integers ID0+1 to ID0+N are unique numbers for call ID to MC64U,
! so 1st call uses 1..N, 2nd call uses N+1..2N, etc
      ID0 = (ID-1) * N

! Main loop. Each pass round this loop either results in a new
! assignment or gives a column with no assignment

      OUTER: DO JORD = NUM0 + 1, N

! Each pass uses unique number ID1
         ID1 = ID0 + JORD
! J is unmatched column
         J = FC (JORD-NUM0)
         PR (J) = - 1
         SCAN: DO K = 1, JORD
! Look for a cheap assignment
            IF (ARP (J) .GE.LENC (J) ) GOTO 30
            IN1 = IP (J) + ARP (J)
            IN2 = IP (J) + LENC (J) - 1
            DO II = IN1, IN2
               I = IRN (II)
               IF (IPERM (I) .EQ.0) GOTO 80
            END DO
! No cheap assignment in row
            ARP (J) = LENC (J)
! Begin looking for assignment chain starting with row J
   30       OUT (J) = LENC (J) - 1
! Inner loop.  Extends chain by one or backtracks
            DO KK = 1, JORD
               IN1 = OUT (J)
               IF (IN1.LT.0) GOTO 50
               IN2 = IP (J) + LENC (J) - 1
               IN1 = IN2 - IN1
! Forward scan
               DO II = IN1, IN2
                  I = IRN (II)
                  IF (CV (I) .EQ.ID1) CYCLE
! Column J has not yet been accessed during this pass
                  J1 = J
                  J = IPERM (I)
                  CV (I) = ID1
                  PR (J) = J1
                  OUT (J1) = IN2 - II - 1
                  CYCLE SCAN
               END DO
! Backtracking step.
   50          J1 = PR (J)
               IF (J1.EQ. - 1) THEN
! No augmenting path exists for column J.
                  NFC = NFC + 1
                  FC (NFC) = J
                  IF (NFC.GT.NUM2) THEN
! A matching of maximum size NUM1 is not possible
                     LAST = JORD
                     GOTO 101
                  ENDIF
                  CYCLE OUTER
               ENDIF
               J = J1
            END DO
! End of dummy loop; this point is never reached
         END DO SCAN
! End of dummy loop; this point is never reached

! New assignment is made.
   80    IPERM (I) = J
         ARP (J) = II - IP (J) + 1
         NUM = NUM + 1
         DO 90 K = 1, JORD
            J = PR (J)
            IF (J.EQ. - 1) GOTO 95
            II = IP (J) + LENC (J) - OUT (J) - 2
            I = IRN (II)
            IPERM (I) = J
   90    END DO
! End of dummy loop; this point is never reached

   95    IF (NUM.EQ.NUM1) THEN
! A matching of maximum size NUM1 is found
            LAST = JORD
            GOTO 101
         ENDIF
!
      END DO OUTER

! All unassigned columns have been considered
      LAST = N

! Now, a transversal is computed or is not possible.
! Complete FC before returning.
  101 DO 110 JORD = LAST + 1, N
         NFC = NFC + 1
         FC (NFC) = FC (JORD-NUM0)
  110 END DO

!  199 RETURN
      RETURN
END SUBROUTINE MC64_HSL_U

!**********************************************************************
!CCCC LAST UPDATE Tue Nov 26 03:20:26 MET 2002
SUBROUTINE MC64_HSL_W (M, N, NE, IP, IRN, A, IPERM, NUM, JPERM,  &
      OUT, PR, Q, L, U, D, RINF)
      IMPLICIT NONE
!
! *** Copyright (c) 2002  I.S. Duff and J. Koster                   ***
! *** Although every effort has been made to ensure robustness and  ***
! *** reliability of the subroutines in this MC64 suite, we         ***
! *** disclaim any liability arising through the use or misuse of   ***
! *** any of the subroutines.                                       ***
!
      INTEGER M, N, NE, NUM
      INTEGER IP (N + 1), IRN (NE), IPERM (M), JPERM (N), OUT (N),      &
      PR (N), Q (M), L (M)
      REAL(WP) A (NE), U (M), D (M), RINF

! M, N, NE, IP, IRN are described in MC64A.
! A is a REAL array of length NE.
!   A(K), K=1..NE, must be set to the value of the entry that
!   corresponds to IRN(K). It is not altered.
!   All values A(K) must be non-negative.
! IPERM is an INTEGER array of length M. On exit, it contains the
!   weighted matching: IPERM(I) = 0 or row I is matched to column
!   IPERM(I).
! NUM is an INTEGER variable. On exit, it contains the cardinality of
!   the matching stored in IPERM.
! D is a REAL array of length M.
!   On exit, V = D(1:N) contains the dual column variable.
!   If U(1:M) denotes the dual row variable and if the matrix
!   is structurally nonsingular (NUM = N), the following holds:
!      U(I)+V(J) <= A(I,J)  if IPERM(I) /= J
!      U(I)+V(J)  = A(I,J)  if IPERM(I)  = J
!      U(I) = 0  if IPERM(I) = 0
!      V(J) = 0  if there is no I for which IPERM(I) = J
!   On entry, U(1) contains the relaxation parameter RLX.
! RINF is the largest positive real number

! Local variables
      INTEGER I, I0, II, J, JJ, JORD, Q0, QLEN, JDUM, ISP, JSP, K, K0,  &
      K1, K2, KK, KK1, KK2, UP, LOW, LPOS
      REAL(WP) CSP, DI, DMIN, DNEW, DQ0, VJ
! Local parameters
      REAL(WP),PARAMETER :: ZERO=0.0_WP
! External subroutines and/or functions
      EXTERNAL MC64D, MC64E, MC64F

! Initialize variables to eliminate copmiler warnings
      ISP = -HUGE(ISP); JSP = -HUGE(JSP)

! Set RINF to largest positive real number
      RINF = HUGE(RINF)

! Initialization
      NUM = 0
      DO K = 1, N
         D (K) = ZERO
         JPERM (K) = 0
         PR (K) = IP (K)
      END DO
      DO K = 1, M
         U (K) = RINF
         IPERM (K) = 0
         L (K) = 0
      END DO
! Initialize U(I)
      DO J = 1, N
         DO K = IP (J), IP (J + 1) - 1
            I = IRN (K)
            IF (A (K) .GT.U (I) ) CYCLE
            U (I) = A (K)
            IPERM (I) = J
            L (I) = K
         END DO
      END DO
      DO I = 1, M
         J = IPERM (I)
         IF (J.EQ.0) CYCLE
! Row I is not empty
         IPERM (I) = 0
         IF (JPERM (J) .NE.0) CYCLE
! Don't choose cheap assignment from dense columns
         IF (IP (J + 1) - IP (J) .GT.N / 10.AND.N.GT.50) CYCLE
! Assignment of column J to row I
         NUM = NUM + 1
         IPERM (I) = J
         JPERM (J) = L (I)
      END DO

!      write(14,*) 'Number of cheap assignments ',NUM

      IF (NUM.EQ.N) GOTO 1000
! Scan unassigned columns; improve assignment
      DO J = 1, N
! JPERM(J) ne 0 iff column J is already assigned
         IF (JPERM (J) .NE.0) CYCLE
         K1 = IP (J)
         K2 = IP (J + 1) - 1
! Continue only if column J is not empty
         IF (K1.GT.K2) CYCLE
!       VJ = RINF
! Changes made to allow for NaNs
         I0 = IRN (K1)
         VJ = A (K1) - U (I0)
         K0 = K1
         DO K = K1 + 1, K2
            I = IRN (K)
            DI = A (K) - U (I)
            IF (DI.GT.VJ) CYCLE
            IF (.NOT.(DI.LT.VJ.OR.DI.EQ.RINF)) THEN
               IF (IPERM (I) .NE.0.OR.IPERM (I0) .EQ.0) CYCLE
            ENDIF
            VJ = DI
            I0 = I
            K0 = K
         END DO
         D (J) = VJ
         K = K0
         I = I0
         IF (IPERM (I) .EQ.0) GOTO 90
         DO K = K0, K2
            I = IRN (K)
            IF (A (K) - U (I) .GT.VJ) CYCLE
            JJ = IPERM (I)
! Scan remaining part of assigned column JJ
            KK1 = PR (JJ)
            KK2 = IP (JJ + 1) - 1
            IF (KK1.GT.KK2) CYCLE
            DO KK = KK1, KK2
               II = IRN (KK)
               IF (IPERM (II) .GT.0) CYCLE
               IF (A (KK) - U (II) .LE.D (JJ) ) GOTO 80
            END DO
            PR (JJ) = KK2 + 1
         END DO
         CYCLE
   80    JPERM (JJ) = KK
         IPERM (II) = JJ
         PR (JJ) = KK + 1
   90    NUM = NUM + 1
         JPERM (J) = K
         IPERM (I) = J
         PR (J) = K + 1
      END DO

!     write(14,*) 'Number of improved  assignments ',NUM

      IF (NUM.EQ.N) GOTO 1000

! Prepare for main loop
      DO I = 1, M
         D (I) = RINF
         L (I) = 0
      END DO

! Main loop ... each pass round this loop is similar to Dijkstra's
! algorithm for solving the single source shortest path problem

      DO JORD = 1, N

         IF (JPERM (JORD) .NE.0) CYCLE
! JORD is next unmatched column
! DMIN is the length of shortest path in the tree
         DMIN = RINF
         QLEN = 0
         LOW = N + 1
         UP = N + 1
! CSP is the cost of the shortest augmenting path to unassigned row
! IRN(ISP). The corresponding column index is JSP.
         CSP = RINF
! Build shortest path tree starting from unassigned column (root) JORD
         J = JORD
         PR (J) = - 1

! Scan column J
         DO K = IP (J), IP (J + 1) - 1
!         IF (N.EQ.3) THEN
!           write(14,*) 'Scanning column ',J
!           write(14,*) 'IP  ',IP(1:4)
!           write(14,*) 'IRN ',IRN(1:6)
!           write(14,*) 'A   ',A(1:6)
!           write(14,*) 'U ',U(1:3)
!           write(14,*) 'IPERM ',IPERM(1:3)
!         ENDIF
            I = IRN (K)
            DNEW = A (K) - U (I)
            IF (DNEW.GE.CSP) CYCLE
            IF (IPERM (I) .EQ.0) THEN
               CSP = DNEW
               ISP = K
               JSP = J
            ELSE
               IF (DNEW.LT.DMIN) DMIN = DNEW
               D (I) = DNEW
               QLEN = QLEN + 1
               Q (QLEN) = K
            ENDIF
         END DO
! Initialize heap Q and Q2 with rows held in Q(1:QLEN)
         Q0 = QLEN
         QLEN = 0
         DO KK = 1, Q0
            K = Q (KK)
            I = IRN (K)
            IF (CSP.LE.D (I) ) THEN
               D (I) = RINF
               CYCLE
            ENDIF
            IF (D (I) .LE.DMIN) THEN
               LOW = LOW - 1
               Q (LOW) = I
               L (I) = LOW
            ELSE
               QLEN = QLEN + 1
               L (I) = QLEN
               CALL MC64D (I, M, Q, D, L, 2)
            ENDIF
! Update tree
            JJ = IPERM (I)
            OUT (JJ) = K
            PR (JJ) = J
         END DO

         DO JDUM = 1, NUM

! If Q2 is empty, extract rows from Q
            IF (LOW.EQ.UP) THEN
               IF (QLEN.EQ.0) GOTO 160
               I = Q (1)
               IF (D (I) .GE.CSP) GOTO 160
               DMIN = D (I)
  152          CALL MC64E (QLEN, M, Q, D, L, 2)
               LOW = LOW - 1
               Q (LOW) = I
               L (I) = LOW
               IF (QLEN.EQ.0) GOTO 153
               I = Q (1)
               IF (D (I) .GT.DMIN) GOTO 153
               GOTO 152
            ENDIF
! Q0 is row whose distance D(Q0) to the root is smallest
  153       Q0 = Q (UP - 1)
            DQ0 = D (Q0)
! Exit loop if path to Q0 is longer than the shortest augmenting path
            IF (DQ0.GE.CSP) GOTO 160
            UP = UP - 1

! Scan column that matches with row Q0
            J = IPERM (Q0)
            VJ = DQ0 - A (JPERM (J) ) + U (Q0)
            DO K = IP (J), IP (J + 1) - 1
               I = IRN (K)
               IF (L (I) .GE.UP) CYCLE
! DNEW is new cost
               DNEW = VJ + A (K) - U (I)
! Do not update D(I) if DNEW ge cost of shortest path
               IF (DNEW.GE.CSP) CYCLE
               IF (IPERM (I) .EQ.0) THEN
! Row I is unmatched; update shortest path info
                  CSP = DNEW
                  ISP = K
                  JSP = J
               ELSE
! Row I is matched; do not update D(I) if DNEW is larger
                  DI = D (I)
                  IF (DI.LE.DNEW) CYCLE
                  IF (L (I) .GE.LOW) CYCLE
                  D (I) = DNEW
                  IF (DNEW.LE.DMIN) THEN
                     LPOS = L (I)
                     IF (LPOS.NE.0) CALL MC64F (LPOS, QLEN, M, Q, D, L,&
                     2)
                     LOW = LOW - 1
                     Q (LOW) = I
                     L (I) = LOW
                  ELSE
                     IF (L (I) .EQ.0) THEN
                        QLEN = QLEN + 1
                        L (I) = QLEN
                     ENDIF
                     CALL MC64D (I, M, Q, D, L, 2)
                  ENDIF
! Update tree
                  JJ = IPERM (I)
                  OUT (JJ) = K
                  PR (JJ) = J
               ENDIF
            END DO
         END DO

! If CSP = RINF, no augmenting path is found
  160    IF (CSP.EQ.RINF) GOTO 190
! Find augmenting path by tracing backward in PR; update IPERM,JPERM
         NUM = NUM + 1

!       write(14,*) 'NUM = ',NUM
!       write(14,*) 'Augmenting path found from unmatched col ',JORD

         I = IRN (ISP)
         IPERM (I) = JSP
         JPERM (JSP) = ISP
         J = JSP
         DO 170 JDUM = 1, NUM
            JJ = PR (J)
            IF (JJ.EQ. - 1) GOTO 180
            K = OUT (J)
            I = IRN (K)
            IPERM (I) = JJ
            JPERM (JJ) = K
            J = JJ
  170    END DO
! End of dummy loop; this point is never reached

! Update U for rows in Q(UP:N)
  180    DO 185 KK = UP, N
            I = Q (KK)
            U (I) = U (I) + D (I) - CSP
  185    END DO
  190    DO 191 KK = LOW, N
            I = Q (KK)
            D (I) = RINF
            L (I) = 0
  191    END DO
         DO 193 KK = 1, QLEN
            I = Q (KK)
            D (I) = RINF
            L (I) = 0
  193    END DO

      END DO
! End of main loop


! Set dual column variable in D(1:N)
 1000 DO 200 J = 1, N
         K = JPERM (J)
         IF (K.NE.0) THEN
            D (J) = A (K) - U (IRN (K) )
         ELSE
            D (J) = ZERO
         ENDIF
  200 END DO
      DO 1201 I = 1, M
         IF (IPERM (I) .EQ.0) U (I) = ZERO
 1201 END DO

      IF (NUM.EQ.N.AND.M.EQ.N) GOTO 1100
! Complete IPERM; L, JPERM are work arrays
      CALL MC64_HSL_X (M, N, IPERM, L, JPERM)

! The matrix is structurally singular, complete IPERM.
! JPERM, OUT are work arrays
!     DO 300 J = 1,N
!       JPERM(J) = 0
! 300 CONTINUE
!     K = 0
!     DO 310 I = 1,N
!       IF (IPERM(I).EQ.0) THEN
!         K = K + 1
!         OUT(K) = I
!       ELSE
!         J = IPERM(I)
!         JPERM(J) = I
!       ENDIF
! 310 CONTINUE
!     K = 0
!     DO 320 J = 1,N
!       IF (JPERM(J).NE.0) GO TO 320
!       K = K + 1
!       JDUM = OUT(K)
!       IPERM(JDUM) = - J
! 320 CONTINUE
 1100 RETURN
END SUBROUTINE MC64_HSL_W


!**********************************************************************
!CCCC LAST UPDATE Tue Nov 26 03:20:26 MET 2002
SUBROUTINE MC64_HSL_Z (M, N, IRN, LIRN, IP, LENC, IPERM, NUM, PR,&
      ARP, CV, OUT)
      IMPLICIT NONE
!
! *** Copyright (c) 2002  I.S. Duff and J. Koster                   ***
! *** Although every effort has been made to ensure robustness and  ***
! *** reliability of the subroutines in this MC64 suite, we         ***
! *** disclaim any liability arising through the use or misuse of   ***
! *** any of the subroutines.                                       ***
!
! PR(I) is the previous row to I in the depth first search.
!   It is used as a work array in the sorting algorithm.
! Elements (IPERM(I),I) I=1,...M  are non-zero at the end of the
!   algorithm unless N assignments have not been made.  In which case
!   (IPERM(I),I) will be zero for N-NUM entries.
! CV(I) is the most recent row extension at which column I was visited.
! ARP(I) is one less than the number of non-zeros in row I
!   which have not been scanned when looking for a cheap assignment.
! OUT(I) is one less than the number of non-zeros in row I
!   which have not been scanned during one pass through the main loop.
!
      INTEGER LIRN, M, N, NUM
      INTEGER ARP (N), CV (M), IRN (LIRN), IP (N), IPERM (M), LENC (N), &
      OUT (N), PR (N)

      INTEGER I, II, IN1, IN2, J, J1, JORD, K, KK

! Initialize variables to eliminate copmiler warnings
      II = -HUGE(II); IN2 = -HUGE(IN2)

      DO I = 1, M
         CV (I) = 0
         IPERM (I) = 0
      END DO
      DO J = 1, N
         ARP (J) = LENC (J) - 1
      END DO
      NUM = 0
!
! Main loop. Each pass round this loop either results in a new
! assignment or gives a row with no assignment.
!
      OUTER: DO JORD = 1, N
!
         J = JORD
         PR (J) = - 1
         SCAN: DO K = 1, JORD
! Look for a cheap assignment
            IN1 = ARP (J)
            IF (IN1.GE.0) THEN
               IN2 = IP (J) + LENC (J) - 1
               IN1 = IN2 - IN1
               DO II = IN1, IN2
                  I = IRN (II)
                  IF (IPERM (I) .EQ.0) EXIT scan
               END DO
! No cheap assignment in row.
               ARP (J) = - 1
            END IF
! Begin looking for assignment chain starting with row J.
            OUT (J) = LENC (J) - 1
! Inner loop.  Extends chain by one or backtracks.
            DO KK = 1, JORD
               IN1 = OUT (J)
               IF (IN1.GE.0) THEN
                  IN2 = IP (J) + LENC (J) - 1
                  IN1 = IN2 - IN1
! Forward scan.
                  DO II = IN1, IN2
                     I = IRN (II)
                     IF (CV (I) .EQ.JORD) CYCLE
! Column I has not yet been accessed during this pass.
                     J1 = J
                     J = IPERM (I)
                     CV (I) = JORD
                     PR (J) = J1
                     OUT (J1) = IN2 - II - 1
                     CYCLE SCAN
                  END DO
               ENDIF
! Backtracking step.
               J = PR (J)
               IF (J.EQ. - 1) CYCLE OUTER
            END DO
         END DO SCAN
!
! New assignment is made.
         IPERM (I) = J
         ARP (J) = IN2 - II - 1
         NUM = NUM + 1
         DO K = 1, JORD
            J = PR (J)
            IF (J.EQ. - 1) EXIT
            II = IP (J) + LENC (J) - OUT (J) - 2
            I = IRN (II)
            IPERM (I) = J
         END DO
!
      END DO OUTER


! IPERM is complete if M = N and NUM = N
      IF (M.EQ.N.and.NUM.EQ.N) RETURN

! Complete IPERM; CV, ARP are work arrays
      CALL MC64_HSL_X (M, N, IPERM, CV, ARP)

END SUBROUTINE MC64_HSL_Z

!**********************************************************************
!CCCC LAST UPDATE Tue Nov 26 03:20:26 MET 2002
SUBROUTINE MC64_HSL_X (M, N, IPERM, RW, CW)
      IMPLICIT NONE
!
! *** Copyright (c) 2002  I.S. Duff and J. Koster                   ***
! *** Although every effort has been made to ensure robustness and  ***
! *** reliability of the subroutines in this MC64 suite, we         ***
! *** disclaim any liability arising through the use or misuse of   ***
! *** any of the subroutines.                                       ***
!
! Complete the (incomplete) row permutation in IPERM.
!
      INTEGER M, N
      INTEGER RW (M), CW (N), IPERM (M)

      INTEGER I, J, K

! If M=N, the matrix is structurally singular, complete IPERM
! If M>N, the matrix is rectangular, complete IPERM

! RW, CW are work arrays;
! Store indices of unmatched rows in RW
! Mark matched columns in CW

      DO J = 1, N
         CW (J) = 0
      END DO
      K = 0
      DO I = 1, M
         IF (IPERM (I) .EQ.0) THEN
            K = K + 1
            RW (K) = I
         ELSE
            J = IPERM (I)
            CW (J) = I
         ENDIF
      END DO
      K = 0
      DO J = 1, N
         IF (CW (J) .NE.0) CYCLE
         K = K + 1
         I = RW (K)
         IPERM (I) = - J
      END DO
      DO J = N + 1, M
         K = K + 1
         I = RW (K)
         IPERM (I) = - J
      END DO
END SUBROUTINE MC64_HSL_X

end module hsl_mc64_single
! COPYRIGHT (c) 2012 Science and Technology Facilities Council
! Original date 6 June 2012 as version 1.0.0
!
! Version 1.1.0
! See ChangeLog for version history.
!
! Written by: Jonathan Hogg and Jennifer Scott

! Given a sparse symmetric  matrix A, this package 
! uses a matching algorithm to compute an elimination
! order that is suitable for use with a sparse direct solver. 
! It optionally computes scaling factors.

! Note: this version uses mc64 to compute the matching
!       and follows approach of Duff and Pralet (2005).
!       It does not compute the optimal matching in the
!       structurally singular case as we found the extra work
!       required to do this did not result in a better
!       pivot order.

! To convert from double to single:
! * Change wp
! * Replace _single by _single
! * Replace call to mc64w to mc64w

module hsl_mc80_single
   use hsl_mc34_single
   use hsl_mc68_integer
   ! Also calls mc64
   implicit none

   private
   public :: mc80_order, mc80_order_full, mc80_control, mc80_info

   integer, parameter :: wp = kind(0.0)
   real(wp), parameter :: one = 1.0_wp
   real(wp), parameter :: zero = 0.0_wp
   real(wp), parameter :: rinf = huge(0.0_wp)

   ! Error flags
   integer, parameter :: MC80_SUCCESS               = 0
   integer, parameter :: MC80_ERROR_ALLOCATION      = -1
   integer, parameter :: MC80_ERROR_A_N_OOR         = -2
   integer, parameter :: MC80_ERROR_SINGULAR        = -3
   integer, parameter :: MC80_ERROR_MC68            = -4
   integer, parameter :: MC80_ERROR_ORD_OOR         = -5
   integer, parameter :: MC80_ERROR_NO_METIS        = -6

   ! warning flags
   integer, parameter :: MC80_WARNING_SINGULAR      = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Data type for information returned by code
   !
   type mc80_info
      integer :: compress_rank = 0 ! order of compressed matrix (ordering
         ! applied to this matrix)
      integer :: flag = 0 ! Takes one of the enumerated flag values:
         ! Possible  following values:
         !    0 : successful entry (for structurally nonsingular matrix).
         !   -1 : allocation error
         !   -2 : Metis not available when required
         !   -3 : either n or ord is out of range (immediate return)
         !   -4 : unexpected error returned by hsl_mc68
         !   +1 : successful entry (for structurally singular matrix).
      integer :: flag68 = 0 ! error flag from hsl_mc68
      integer :: max_cycle = 0 ! maximum cycle length
      integer :: struct_rank = 0 ! structural rank 
      integer :: stat = 0 ! stat parameter

   end type mc80_info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   type mc80_control
      logical :: action = .true. ! if set to .true. and the matrix
         ! is found to be structurally singular, the code exits immediately
         ! with an error flag set.
      logical :: unmatched_scale_zero = .false. ! If true, then in the singular
         ! case, the scaling factors associated with the unmatched part of the
         ! matrix are set to zero. This will breakdown if the matched submatrix
         ! is numerically singular, but if no any factorization will likely have
         ! far fewer delayed pivots. If false, then Duff and Pralet are followed
         ! and the value is set such that the corresponding scaled entries are
         ! <= 1 in absolute value.
      logical :: unmatched_last = .false. ! If true, then in singular case, rows
         ! and columns associated with the unmatched part are ordered last in
         ! the elimination order. If false, they are placed in a fill-minimising
         ! position.
   end type mc80_control

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   interface mc80_order
      ! to be used if lower triangle of A available
      module procedure mc80_order_single
   end interface

   interface mc80_order_full
      ! to be used if lower and upper triangles of A available
      module procedure mc80_order_full_single
   end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! This routine computes the matching, splits long cycles,
! compresses the matrix, applies AMD, Min Deg or Metis to the compressed
! matrix and then returns an ordering for the original
! matrix. This ordering flags 2x2 pivots using negative signs
! (as used on input to hsl_ma77).
! If scale is present, the mc64 scaling factors are returned.
! The matrix may be singular.
!
! Input (ptr, row , val) is the ** lower triangular part ** of the matrix.
! Diagonal entries need not be present.
! There is no matrix data checking.
!
subroutine mc80_order_single(ord, n, ptr, row, val, order, &
      control, info, scale)
   integer, intent(in) :: ord ! controls ordering on compressed matrix
         ! 1 An approximate minimum degree ordering is used.
         ! 2 A minimum degree ordering is used (as in MA27).
         ! 3 METIS ordering with default settings is used.
   integer, intent(in) :: n
   integer, dimension(:), intent(in) :: ptr
   integer, dimension(:), intent(in) :: row
   real(wp), dimension(:), intent(in) :: val
   integer, dimension(:), intent(out) :: order ! |order(i)|  holds the position
      ! of variable i in the elimination order (pivot sequence). If a
      ! 1x1 pivot i is obtained,  order(i)>0. If a
      ! 2x2 pivot involving  i and j is obtained, 
      ! order(i)<0, order(j)<0 and |order(j)|=|order(i)|+1. 
   type(mc80_control), intent(in) :: control
   type (mc80_info), intent(out) :: info ! used to hold information
   real(wp), dimension(n), intent(out), optional :: scale ! if present,
      ! returns the mc64 symmetric scaling

   integer, dimension(:), allocatable :: cperm ! used to hold matching
   integer, dimension(:), allocatable :: ptr2 ! column pointers for expanded
      ! matrix. 
   integer, dimension(:), allocatable :: row2 ! row indices for expanded matrix.
   real(wp), dimension(:), allocatable :: val2 ! entries of expanded matrix.
   real(wp), dimension(:), allocatable :: scale2 ! holds scaling factors
      ! if scale is not present.

   integer :: i, j, k, ne

   info%compress_rank = 0
   info%flag = 0
   info%flag68 = 0
   info%stat = 0
   info%max_cycle = 0 
   info%struct_rank = n

   ! check n has valid value
   if (n < 0) then
      info%flag = MC80_ERROR_A_N_OOR
      return
   end if

   ! check ord has valid value
   if (ord < 1 .or. ord > 3) then
      info%flag = MC80_ERROR_ORD_OOR
      return
   end if

   ! just return with no action if n = 0
   if (n.eq. 0) return

   !
   ! Take absolute values, expand out, removing any explicit zeroes
   !
   ne = ptr(n+1) - 1
   allocate(ptr2(n+1), row2(2*ne), val2(2*ne), cperm(n), stat=info%stat)
   if (info%stat.ne.0) then
      info%flag = MC80_ERROR_ALLOCATION
      return
   end if

   k = 1
   do i = 1, n
      ptr2(i) = k
      do j = ptr(i), ptr(i+1)-1
         if (val(j).eq.zero) cycle
         row2(k) = row(j)
         val2(k) = abs(val(j))
         k = k + 1
      end do
   end do
   ptr2(n+1) = k

   call mc34_expand(n, row2, ptr2, cperm, a=val2)
   
   ! Compute matching and scaling

   if (present(scale)) then 
      call mc80_scale(n,ptr2,row2,val2,scale,control,info, &
           perm=cperm)
   else
      allocate(scale2(n), stat=info%stat)
      if (info%stat.ne.0) then
         info%flag = MC80_ERROR_ALLOCATION
         return
      end if
      call mc80_scale(n,ptr2,row2,val2,scale2,control,info, &
           perm=cperm)
      deallocate(scale2, stat=info%stat)
   end if
   deallocate(val2, stat=info%stat)

   if (info%flag.lt.0) return

   ! Note: row j is matched with column cperm(j)
   !
   ! Split matching into 1- and 2-cycles only and then
   ! compress matrix and order.

   call mc80_split(ord,n,row2,ptr2,order,cperm,control,info)

   if(present(scale)) scale(1:n) = exp( scale(1:n) )

end subroutine mc80_order_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! This routine is the same as mc80_order EXCEPT on
! input ptr, row , val hold the ** lower AND upper ** triangular 
! parts of the matrix.
! this reduces amount of copies of matrix required (so slightly
! more efficient on memory and does not need to expand supplied matrix)
!  
subroutine mc80_order_full_single(ord, n, ptr, row, val, order,  &
      control, info, scale)
   integer, intent(in) :: ord ! controls ordering on compressed matrix
         ! 1 An approximate minimum degree ordering is used.
         ! 2 A minimum degree ordering is used (as in MA27).
         ! 3 METIS ordering with default settings is used.
   integer, intent(in) :: n
   integer, dimension(:), intent(in) :: ptr
   integer, dimension(:), intent(in) :: row
   real(wp), dimension(:), intent(in) :: val
   integer, dimension(:), intent(out) :: order ! |order(i)|  holds the position
      ! of variable i in the elimination order (pivot sequence). If a
      ! 1x1 pivot i is obtained,  order(i)>0. If a
      ! 2x2 pivot involving  i and j is obtained, 
      ! order(i)<0, order(j)<0 and |order(j)|=|order(i)|+1. 
   type(mc80_control), intent(in) :: control
   type (mc80_info), intent(out) :: info ! used to hold information
   real(wp), dimension(n), intent(out), optional :: scale ! if present,
      ! returns the mc64 symmetric scaling

   integer, dimension(:), allocatable :: cperm ! used to hold matching
   integer, dimension(:), allocatable :: ptr2 ! column pointers for expanded 
      ! matrix. 
   integer, dimension(:), allocatable :: row2 ! row indices for expanded matrix 
   real(wp), dimension(:), allocatable :: val2 ! entries of expanded matrix.
   real(wp), dimension(:), allocatable :: scale2 ! holds scaling factors
      ! if scale is not present.

   integer :: i, j, k, ne

   info%compress_rank = 0
   info%flag = 0
   info%flag68 = 0
   info%stat = 0
   info%max_cycle = 0 
   info%struct_rank = n

   ! check n has valid value
   if (n < 0) then
     info%flag = MC80_ERROR_A_N_OOR
     return
   end if

   ! check ord has valid value
   if (ord < 1 .or. ord > 3) then
     info%flag = MC80_ERROR_ORD_OOR
     return
   end if

   ! just return with no action if n = 0
   if (n.eq.0) return

   !
   ! Remove any explicit zeroes and take absolute values
   !
   ne = ptr(n+1) - 1
   allocate(ptr2(n+1), row2(ne), val2(ne), cperm(n), stat=info%stat)
   if (info%stat.ne.0) then
      info%flag = MC80_ERROR_ALLOCATION
      return
   end if

   k = 1
   do i = 1, n
      ptr2(i) = k
      do j = ptr(i), ptr(i+1)-1
         if (val(j).eq.zero) cycle
         row2(k) = row(j)
         val2(k) = abs(val(j))
         k = k + 1
      end do
   end do
   ptr2(n+1) = k

   ! Compute matching and scaling

   if (present(scale)) then 
      call mc80_scale(n,ptr2,row2,val2,scale,control,info, &
         perm=cperm)
   else
      allocate(scale2(n), stat=info%stat)
      if (info%stat.ne.0) then
         info%flag = MC80_ERROR_ALLOCATION
         return
      end if
      call mc80_scale(n,ptr2,row2,val2,scale2,control,info, &
         perm=cperm)
      deallocate(scale2, stat=info%stat)
   end if
   deallocate(val2, stat=info%stat)

   if (info%flag.lt.0) return

   ! Note: row j is matched with column cperm(j)
   ! write (*,'(a,15i4)') 'cperm',cperm(1:min(15,n))
   !
   ! Split matching into 1- and 2-cycles only and then
   ! compress matrix and order.

   call mc80_split(ord,n,row2,ptr2,order,cperm,control,info)

   if(present(scale)) scale(1:n) = exp( scale(1:n) )

end subroutine mc80_order_full_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Split matching into 1- and 2-cycles only and then
! compress matrix and order using mc68.
!
! Input (ptr2, row2 , val2) holds the ** lower and upper triangles ** 
! of the matrix (with explicit zeros removed).
! Overwritten in the singular case
!
subroutine mc80_split(ord,n,row2,ptr2,order,cperm,control,info)
   integer, intent(in) :: ord ! controls choice of ordering algorithm 
      ! on compressed matrix
   integer, intent(in) :: n
   integer, dimension(:), intent(in) :: ptr2 
   integer, dimension(:), intent(in) :: row2 
   integer, dimension(n), intent(out) :: order ! used to hold ordering
   integer, dimension(n), intent(inout) :: cperm ! used to hold matching 
   type (mc80_control), intent(in) :: control
   type (mc80_info), intent(inout) :: info ! used to hold information

   type (mc68_control) :: control68
   type (mc68_info) :: info68

   integer, dimension(:), allocatable :: iwork ! work array
   integer, dimension(:), allocatable :: old_to_new, new_to_old
      ! holds mapping between original matrix indices and those in condensed
      ! matrix.
   integer, dimension(:), allocatable :: ptr3 ! column pointers for condensed 
      ! matrix.
   integer, dimension(:), allocatable :: row3 ! row indices for condensed 
      ! matrix.


   integer :: csz ! current cycle length
   integer :: i, j, j1, j2, jj, k, krow
   integer :: max_csz ! maximum cycle length
   integer :: ncomp ! order of compressed matrix
   integer :: ncomp_matched ! order of compressed matrix (matched entries only)
   integer :: ne ! number of non zeros

   ! Use iwork to track what has been matched:
   ! -2 unmatched
   ! -1 matched as singleton
   !  0 not yet seen
   ! >0 matched with specified node

   ne = ptr2(n+1) - 1
   allocate(ptr3(n+1), row3(ne), old_to_new(n), new_to_old(n), iwork(n), &
      stat=info%stat)
   if (info%stat.ne.0) return

   iwork(1:n) = 0
   max_csz = 0
   do i = 1, n
      if (iwork(i).ne.0) cycle
      j = i
      csz = 0
      do
         if (cperm(j).eq.-1) then
            ! unmatched by MC64
            iwork(j) = -2
            csz = csz + 1
            exit
         else if (cperm(j).eq.i) then
            ! match as singleton, unmatched or finished
            iwork(j) = -1
            csz = csz + 1
            exit
         end if
         ! match j and cperm(j)
         jj = cperm(j)
         iwork(j) = jj
         iwork(jj) = j
         csz = csz + 2
         ! move onto next start of pair
         j = cperm(jj)
         if (j.eq.i) exit
      end do
      max_csz = max(max_csz, csz)
   end do
   info%max_cycle = max_csz

   ! Overwrite cperm with new matching
   cperm(1:n) = iwork(1:n)

   !
   ! Build maps for new numbering schemes
   !
   k = 1
   do i = 1,n
      j = cperm(i)
      if (control%unmatched_last .and. j.eq.-2) cycle
      if (j<i .and. j.gt.0) cycle
      old_to_new(i) = k
      new_to_old(k) = i ! note: new_to_old only maps to first of a pair
      if (j.gt.0) old_to_new(j) = k   
      k = k + 1
   end do
   ncomp_matched = k-1
   
   ! Place unmatched columns to rear
   if(control%unmatched_last) then
      do i = 1, n
         if(cperm(i).eq.-2) then
            old_to_new(i) = k
            new_to_old(k) = i
            k = k + 1
         endif
      end do
   end if

   !
   ! Produce a condensed version of the matrix for ordering.
   ! Hold pattern using ptr3 and row3.
   !
   ptr3(1) = 1
   iwork(:) = 0 ! Use to indicate if entry is in a paired column
   ncomp = 1
   jj = 1
   do i = 1, n
      j = cperm(i)
      if (j<i .and. j.gt.0) cycle ! already seen
      if (control%unmatched_last .and. j.eq.-2) cycle ! column not participating
      do k = ptr2(i), ptr2(i+1)-1
         krow = old_to_new(row2(k))
         if (iwork(krow).eq.i) cycle ! already added to column
         if (krow>ncomp_matched) cycle ! unmatched row not participating
         row3(jj) = krow
         jj = jj + 1
         iwork(krow) = i
      end do
      if (j.gt.0) then
         ! Also check column cperm(i)
         do k = ptr2(j), ptr2(j+1)-1
            krow = old_to_new(row2(k))
            if (iwork(krow).eq.i) cycle ! already added to column
            if (krow>ncomp_matched) cycle ! unmatched row not participating
            row3(jj) = krow
            jj = jj + 1
            iwork(krow) = i
         end do
      end if
      ptr3(ncomp+1) = jj
      ncomp = ncomp + 1
   end do
   ncomp = ncomp - 1
   info%compress_rank = ncomp

   if(control%unmatched_last) ncomp = ncomp_matched

   ! store just lower triangular part for input to hsl_mc68
   ptr3(1) = 1
   jj = 1
   j1 = 1
   do i = 1, ncomp
      j2 = ptr3(i+1)
      do k = j1, j2-1
         krow = row3(k)
         if ( krow.lt.i ) cycle ! already added to column
         row3(jj) = krow
         jj = jj + 1
      end do
      ptr3(i+1) = jj
      j1 = j2
   end do

   ! reorder the compressed matrix using hsl_mc68.
   ! switch off hsl_mc68 printing
   control68%lp = -1
   control68%wp = -1
   control68%mp = -1
   control68%print_level = -1
   call mc68_order(ord,ncomp,ptr3,row3,order,control68,info68)

   if (info68%flag < 0) then
      info%flag68 = info68%flag
      select case(info68%flag)
      case(-1)
         info%flag = MC80_ERROR_ALLOCATION
         info%stat = info68%stat
      case(-5)
         info%flag = MC80_ERROR_NO_METIS
      case default
         info%flag = MC80_ERROR_MC68
      end select
      return
   end if

   do i = 1, ncomp
      j = order(i)
      iwork(j) = i
   end do

   !
   ! Translate inverse permutation in iwork back to 
   ! permutation for original variables.
   ! Set negative signs for 2x2 pivots (exploited by hsl_ma77).
   !
   k = 1
   do i = 1, ncomp
      j = new_to_old( iwork(i) )
      order(j) = k
      k = k + 1
      if (cperm(j).gt.0) then
         order(j) = -order(j)
         j = cperm(j)
         order(j) = -k
         k = k + 1
      end if
   end do

   ! Place unmatched columns last
   if(control%unmatched_last) then
      do i = ncomp+1, ncomp + (n-info%struct_rank)
         j = new_to_old( i )
         order(j) = k
         k = k + 1
      end do
   endif

end subroutine mc80_split

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Scale the matrix using MC64, accounting for singular matrices using the
! approach of Duff and Pralet
!
! Expects a full matrix as input
!
subroutine mc80_scale(n, ptr, row, val, scale, control, info, perm)
   integer, intent(in) :: n
   integer, dimension(:), intent(in) :: ptr
   integer, dimension(:), intent(in) :: row
   real(wp), dimension(:), intent(in) :: val
   real(wp), dimension(n), intent(out) :: scale ! returns the symmetric scaling
   type(mc80_control), intent(in) :: control
   type(mc80_info), intent(inout) :: info ! used to hold information
   integer, dimension(n), intent(out), optional :: perm ! if present, returns
      ! the matching

   integer, dimension(:), allocatable :: ptr2 ! column pointers after 
      ! zeros removed.
   integer, dimension(:), allocatable :: row2 ! row indices after zeros
      !  removed. 
   real(wp), dimension(:), allocatable :: val2 ! matrix of absolute values
      ! (zeros removed).
   real(wp), dimension(:), allocatable :: cscale ! temporary copy of scaling
      ! factors. Only needed if A rank deficient. allocated to have size n.

   integer :: i, j, k, ne

   info%struct_rank = n

   !
   ! Remove any explicit zeroes and take absolute values
   !
   ne = ptr(n+1) - 1
   allocate(ptr2(n+1), row2(ne), val2(ne),stat=info%stat)
   if (info%stat.ne.0) then
      info%flag = MC80_ERROR_ALLOCATION
      return
   end if

   k = 1
   do i = 1, n
      ptr2(i) = k
      do j = ptr(i), ptr(i+1)-1
         if (val(j).eq.zero) cycle
         row2(k) = row(j)
         val2(k) = abs(val(j))
         k = k + 1
      end do
   end do
   ptr2(n+1) = k

   call mc80_match(n,row2,ptr2,val2,scale,control,info,perm=perm)
   if (info%flag.lt.0) return

   if (info%struct_rank.ne.n) then
      ! structurally singular case. At this point, scaling factors
      ! for rows in corresponding to rank deficient part are set to 
      ! zero. The following is to set them according to Duff and Pralet.
      deallocate(ptr2, stat=info%stat)
      deallocate(row2, stat=info%stat)
      deallocate(val2, stat=info%stat)
      if(.not.control%unmatched_scale_zero) then
         allocate(cscale(n),stat=info%stat)
         if (info%stat.ne.0) then
            info%flag = MC80_ERROR_ALLOCATION
            return
         end if
         cscale(1:n) = scale(1:n)
         do i = 1,n
            if (cscale(i).ne.-huge(scale)) cycle
            do j = ptr(i), ptr(i+1)-1
               k = row(j)
               if (cscale(k).eq.-huge(scale)) cycle
               scale(i) = max(scale(i), val(j)+scale(k))
            end do
            if(scale(i).eq.-huge(scale)) then
               scale(i) = zero
            else
               scale(i) = -scale(i)
            endif
         end do
      end if
   end if

end subroutine mc80_scale

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Input (ptr2, row2 , val2) holds the ** lower and upper triangles ** 
! of the matrix (with explicit zeros removed).
! val2 holds absolute values of matrix entries.
! Overwritten in the singular case
!
subroutine mc80_match(n,row2,ptr2,val2,scale,control,info,perm)
   integer, intent(in) :: n
   integer, dimension(:), intent(inout) :: ptr2 ! In singular case, overwritten
      ! by column pointers for non singular part of matrix.
   integer, dimension(:), intent(inout) :: row2 ! In singular case, overwritten
      ! by row indices for non singular part of matrix.
   real(wp), dimension(:), intent(inout) :: val2 ! In singular case, overwritten
      ! by entries for non singular part of matrix.
   real(wp), dimension(n), intent(out) :: scale ! returns the symmetric scaling
   type(mc80_control), intent(in) :: control
   type (mc80_info), intent(inout) :: info ! used to hold information
   integer, dimension(n), intent(out), optional :: perm ! if present, returns
      ! the matching

   integer, dimension(:), allocatable :: iwork ! work array
   integer, dimension(:), allocatable :: cperm ! used to hold matching
   integer, dimension(:), allocatable :: old_to_new, new_to_old
      ! holds mapping between original matrix indices and those in reduced
      ! non singular matrix. 
   real(wp), dimension(:), allocatable :: cmax ! (log) column maximum
   real(wp), dimension(:), allocatable :: dw ! array used by mc64

   integer :: i, j, j1, j2, jj, k
   integer :: ne ! number of non zeros
   integer :: nn ! Holds number of rows/cols in non singular part of matrix
   integer :: nne ! Only used in singular case. Holds number of non zeros
     ! in non-singular part of matrix.
   integer :: rank ! returned by mc64
   real(wp) :: colmax ! max. entry in col. of expanded matrix

   allocate(iwork(5*n), cperm(n), dw(2*n), cmax(n), stat=info%stat)
   if (info%stat.ne.0) then
      info%flag = MC80_ERROR_ALLOCATION
      return
   end if
    
   ! Compute column maximums    
   do i = 1,n
      colmax = max(zero,maxval(val2(ptr2(i):ptr2(i+1)-1)))
      if (colmax.ne.zero) colmax = log(colmax)
      cmax(i) = colmax
   end do

   do i = 1,n
      val2(ptr2(i):ptr2(i+1)-1) = cmax(i) - log(val2(ptr2(i):ptr2(i+1)-1))
   end do

   ne = ptr2(n+1)-1
   call mc64w(n,ne,ptr2,row2,val2,cperm,rank, &
      iwork(1),iwork(n+1),iwork(2*n+1),iwork(3*n+1),iwork(4*n+1), &
      dw(1),dw(n+1))

   if (rank.eq.n) then
      do i = 1,n
         scale(i) = (dw(i)+dw(n+i)-cmax(i))/2
      end do
      if (present(perm)) perm(1:n) = cperm(1:n)
      return
   end if

   !!!! we have to handle the singular case. Either immediate exit
   ! or set warning, squeeze out the unmatched entries and recall mc64w.

   info%struct_rank = rank
   if (.not.control%action) then
      info%flag = MC80_ERROR_SINGULAR
      return
   end if

   info%flag = MC80_WARNING_SINGULAR

   allocate(old_to_new(n), new_to_old(n),stat=info%stat)
   if (info%stat.ne.0) then
      info%flag = MC80_ERROR_ALLOCATION
      return
   end if

   k = 0
   do i = 1,n
      if (cperm(i) < 0) then
         ! row i and col j are not part of the matching
         old_to_new(i) = -1
      else
         k = k + 1
         ! old_to_new(i) holds the new index for variable i after
         ! removal of singular part and new_to_old(k) is the
         ! original index for k
         old_to_new(i) = k
         new_to_old(k) = i
      end if
   end do

   ! Overwrite ptr2, row2 and val2
   nne = 0
   k = 0
   ptr2(1) = 1
   j2 = 1
   do i = 1,n
      j1 = j2
      j2 = ptr2(i+1)
      ! skip over unmatched entries
      if (cperm(i) < 0) cycle
      k = k + 1
      do j = j1,j2-1
         jj = row2(j)
         if (cperm(jj) < 0) cycle
         nne = nne + 1
         row2(nne) = old_to_new(jj)
         val2(nne) = val2(j)
      end do
      ptr2(k+1) = nne + 1
    end do
    ! nn is order of non-singular part.
    nn = k
    call mc64w(nn,nne,ptr2,row2,val2,cperm,rank, &
       iwork(1),iwork(nn+1),iwork(2*nn+1),iwork(3*nn+1),iwork(4*nn+1), &
       dw(1),dw(nn+1))

    do i = 1,n
       j = old_to_new(i)
       if (j < 0) then
          scale(i) = -huge(scale)
       else
         ! Note: we need to subtract col max using old matrix numbering
         scale(i) = (dw(j)+dw(nn+j)-cmax(i))/2
      end if
   end do

   if (present(perm)) then
      perm(1:n) = -1
      do i = 1,nn
         j = cperm(i)
         perm(new_to_old(i)) = new_to_old(j)
      end do
   end if

end subroutine mc80_match

!**********************************************************************
end module hsl_mc80_single
