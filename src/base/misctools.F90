!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

!---------------------------------------------------------------
!
! Miscallaneous tools
!
!  subroutine dist
!     compute distance between points on the globe
!  subroutine dist_corr_fun
!     calculate a distance-dependent correlation function
!  subroutine eigen_symm
!     calculate eigenvectors/-values of symmetric matrix
!  function gasdev
!     Calculate normally distributed deviate with zero mean
!     and unit variance (from numerical recipes)
!  function ran_n
!     random number generator
!  subroutine check_dir
!     check if the directory for a file exists, and if not, create it
!  subroutine collapse_loop
!     collapse nested iterators into a single iterator
!
!---------------------------------------------------------------

module misctools

  use GO, only : gol, goPr, goErr

  ! --- in/out -----------------------------

  implicit none

  private

  public :: dist
  public :: dist_corr_fun
  public :: gasdev
  public :: ran_n
  public :: check_dir
  public :: collapse_loop, T_region_iter
  public :: find_unique_strings

  interface collapse_loop
    module procedure make_loop_2d
    module procedure make_loop_3d
    module procedure collapse_loop_2d
    module procedure collapse_loop_3d
  end interface collapse_loop

  type T_region_iter
    integer, allocatable :: counter_2d(:,:)
    integer, allocatable :: counter_3d(:,:)
    integer              :: n_iter
  end type T_region_iter

  ! --- const ------------------------------

  character(len=*), parameter      :: mname = 'misctools'

contains

  subroutine find_unique_strings(input_str_array, num_unique, unique_string_indices)
    ! Given an array of strings, find the unique strings
    ! I/O
    character(len=*), intent(in)                :: input_str_array(:)
    integer, intent(out)                        :: num_unique
    integer, intent(out), allocatable           :: unique_string_indices(:)
    ! Local
    character(len=*), parameter                 :: rname = mname//'/find_unique_string'
    integer                                     :: len_arr, i, j
    integer, parameter                          :: MAX_STRING_LEN = 512
    character(len=MAX_STRING_LEN), allocatable  :: temp_arr_str(:)
    integer, allocatable                        :: temp_arr_idx(:)
    logical                                     :: already_counted

    len_arr = size(input_str_array, 1)
    allocate(temp_arr_str(len_arr), temp_arr_idx(len_arr))

    num_unique = 0
    do i = 1, len_arr
        already_counted = .false.

        do j = 1, num_unique
            if (trim(input_str_array(i)) == trim(temp_arr_str(j))) already_counted = .true.
        end do

        if (.not. already_counted) then
            ! input_str_array(i) is a new string
            num_unique = num_unique + 1
            temp_arr_str(num_unique) = input_str_array(i)
            temp_arr_idx(num_unique) = i
        end if
    end do

    allocate(unique_string_indices(num_unique))
    unique_string_indices(:) = temp_arr_idx(1:num_unique)

    deallocate(temp_arr_str, temp_arr_idx)

  end subroutine find_unique_strings

  subroutine make_loop_2d(i1, i2, j1, j2, output_iter)
    !---------------------------------------------------------------
    ! Like collapse_loop_2d, except that the inputs are not iterators, but simply the limits
    !---------------------------------------------------------------
    ! I/O
    integer, intent(in)                     :: i1, i2, j1, j2
    integer, dimension(:,:), intent(out)    :: output_iter
    ! Local
    character(len=*), parameter             :: rname = mname//'/make_loop_2d'
    integer, allocatable, dimension(:)      :: iter_1, iter_2
    integer                                 :: N1, N2, i, status

    N1 = i2-i1+1
    N2 = j2-j1+1

    allocate(iter_1(N1), iter_2(N2))

    do i = 1, N1
        iter_1(i) = i1+(i-1)
    end do

    do i = 1, N2
        iter_2(i) = j1+(i-1)
    end do

    call collapse_loop_2d(iter_1, iter_2, output_iter)

    deallocate(iter_1, iter_2)

  end subroutine make_loop_2d

  subroutine make_loop_3d(i1, i2, j1, j2, k1, k2, output_iter)
    ! I/O
    integer, intent(in)                 :: i1, i2, j1, j2, k1, k2
    integer, intent(out)                :: output_iter(:,:)
    ! Local
    character(len=*), parameter         :: rname = mname//'/make_loop_3d'
    integer, allocatable, dimension(:)  :: iter_1, iter_2, iter_3
    integer                             :: N1, N2, N3, i, status

    N1 = i2-i1+1
    N2 = j2-j1+1
    N3 = k2-k1+1

    allocate(iter_1(N1), iter_2(N2), iter_3(N3))

    do i = 1, N1
        iter_1(i) = i1+(i-1)
    end do

    do i = 1, N2
        iter_2(i) = j1+(i-1)
    end do

    do i = 1, N3
        iter_3(i) = k1+(i-1)
    end do

    call collapse_loop_3d(iter_1, iter_2, iter_3, output_iter)

    deallocate(iter_1, iter_2, iter_3)

  end subroutine make_loop_3d

  subroutine collapse_loop_2d(iter_1, iter_2, output_iter)
    !---------------------------------------------------------------
    ! Take two iterators iter_1 and iter_2, and create a single iterator of
    ! dimension (n_iter_1 x n_iter_2, 2) to loop over all the indices
    !---------------------------------------------------------------
    ! I/O
    integer, dimension(:), intent(in)       :: iter_1, iter_2
    integer, dimension(:,:), intent(out)    :: output_iter ! assume it has the right shape
    ! Local
    integer                                 :: i, j, Ni, Nj, counter, status
    character(len=*), parameter             :: rname = mname//'/collapse_loop_2d'

    Ni = size(iter_1, 1)
    Nj = size(iter_2, 1)

    if (size(output_iter, 1) /= Ni*Nj) then
        write(0, '("Length of output array ", i7, " is not the product of the two input lengths ", i7, " and ", i7)') size(output_iter, 1), Ni, Nj
        IF_NOTOK_RETURN(status=1)
    end if

    counter = 1
    do i = lbound(iter_1, 1), ubound(iter_1, 1)
        do j = lbound(iter_2, 1), ubound(iter_2, 1)
            output_iter(counter, 1) = iter_1(i)
            output_iter(counter, 2) = iter_2(j)
            counter = counter + 1
        end do
    end do

  end subroutine collapse_loop_2d

  subroutine collapse_loop_3d(iter_1, iter_2, iter_3, output_iter)
    !---------------------------------------------------------------
    ! Take three iterators iter_1, iter_2 and iter_3, and create a single iterator
    ! of dimension (n_iter_1 x n_iter_2 x n_iter_3, 3) to loop over all the indices
    !---------------------------------------------------------------
    ! I/O
    integer, dimension(:), intent(in)       :: iter_1, iter_2, iter_3
    integer, dimension(:,:), intent(out)    :: output_iter ! assume it has the right shape
    ! Local
    integer                                 :: i, j, k, Ni, Nj, Nk, counter, status
    character(len=*), parameter             :: rname = mname//'/collapse_loop_3d'

    Ni = size(iter_1, 1)
    Nj = size(iter_2, 1)
    Nk = size(iter_3, 1)

    if (size(output_iter, 1) /= Ni*Nj*Nk) then
        write(0, '("Length of output array ", i7, " is not the product of the three input lengths ", i7, ", ", i7, " and ", i7)') size(output_iter, 1), Ni, Nj, Nk
        IF_NOTOK_RETURN(status=1)
    end if

    counter = 1
    do i = lbound(iter_1, 1), ubound(iter_1, 1)
        do j = lbound(iter_2, 1), ubound(iter_2, 1)
            do k = lbound(iter_3, 1), ubound(iter_3, 1)
                output_iter(counter, 1) = iter_1(i)
                output_iter(counter, 2) = iter_2(j)
                output_iter(counter, 3) = iter_3(k)
                counter = counter + 1
            end do
        end do
    end do

  end subroutine collapse_loop_3d

  subroutine check_dir(file_name)
    !---------------------------------------------------------------
    ! Check if the directory for a file exists or not, and if not, create it
    !---------------------------------------------------------------
    use os_specs, only : MAX_FILENAME_LEN
    implicit none

    ! I/O
    character(len=*), intent(in)    :: file_name
    ! Local
    character(len=MAX_FILENAME_LEN) :: dir_name
    logical                         :: dir_exists
    integer                         :: sep_pos
    ! Begin

    sep_pos = index(file_name, '/', .true.)
    ! If the file refers to the current folder or the root folder, no need to do anything more
    if (sep_pos .gt. 1) then
        dir_name = file_name(1:sep_pos-1)
        inquire(file=trim(dir_name)//'/.', exist=dir_exists)
        if (.not. dir_exists) call system('mkdir -p '//trim(dir_name))
    end if

  end subroutine check_dir

  subroutine dist( x1g, y1g, x2g, y2g, ddg )
    !---------------------------------------------------------------
    !  Compute distance of two points on the globe.
    !
    !  This formula is correct, also for points close together.
    !
    !  Input:
    !   x1g, y1g : longitude and latitude of first point  (degrees)
    !   x2g, y2g : longitude and latitude of second point (degrees)
    !
    !  Output:
    !   ddg      : distance between these points          (km)
    !----------------------------------------------------------------

    ! --- modules ---------------------------------------------

    use binas, only                 : pi, ae

    ! --- in/out ----------------------------------------------

    real, intent(in)               :: x1g, x2g, y1g, y2g
    real, intent(out)              :: ddg

    ! --- const -----------------------------------------------

    character(len=*), parameter    :: name = mname//', dist'

    ! --- local -----------------------------------------------

    real                           :: x1, x2, y1, y2, dd
    real                           :: cc, dx, dy, dx2, dy2, dd2

    ! --- begin -----------------------------------------------

    cc = pi / 180.0
    !
    x1 = x1g * cc
    x2 = x2g * cc
    y1 = y1g * cc
    y2 = y2g * cc
    !
    !      print *,'x1,x2,y1,y2', x1,x2,y1,y2
    dy = 0.5 * ( y2 - y1 )
    dy = sin(dy)
    dy2 = dy * dy
    !
    dx = 0.5 * ( x2 - x1 )
    dx = sin(dx)
    dx2 = dx * dx * cos(y1) * cos(y2)
    !
    !     print *,'dx,dy,dx2,dy2',dx,dy,dx2,dy2
    dd2 = dx2 + dy2
    dd = sqrt(dd2)
    dd = 2.0 * asin(dd)
    !
    ddg = dd / cc
    ddg = (ddg * 2 * pi * 0.001 * ae) / 360
    !     print *,'dd2,dd,ddg', dd2,dd,ddg
    !
  end subroutine dist


  subroutine dist_corr_fun( choice, corlen, rho, status )

    !-----------------------------------------------------------------------
    !
    ! Calculate a distance-dependent correlation function
    !
    ! choice [in]  - single character describing correlation function:
    !    'e' : exponential distance dependence, with linear tail
    !    'g' : Gaussian distance dependence
    !    't' : Thiebaux autoregressive correlation
    ! corlen [in]  - correlation length in km (1/e length)
    ! rho   [out] - correlation function for 0, 10, ... 21000 km
    !
    !-----------------------------------------------------------------------

    ! --- in/out ----------------------------------------------

    character(1), intent(in)              :: choice
    real, intent(in)                      :: corlen
    real, dimension(0:2100), intent(out)  :: rho
    integer, intent(out)                  :: status

    ! --- const -----------------------------------------------

    character(len=*), parameter    :: name = mname//', dist_corr_fun'

    ! --- local -----------------------------------------------

    integer                        :: i, icorlen

    ! --- begin -----------------------------------------------

    if( corlen <= 0.0 )then
       write (*,'("ERROR - negative correlation length")')
       write (*,'("ERROR in ",a)') name; status=1; return
    end if
    if( corlen < 100.0 )then
       write (*,'("ERROR - correlation length < 100 km")')
       write (*,'("ERROR in ",a)') name; status=1; return
    end if

    rho(0:2100) = 0.0

    icorlen = nint( corlen/10.0 )
    if( 5*icorlen > 2100 )then
       write (*,'("ERROR - correlation length too large")')
       write (*,'("ERROR in ",a)') name; status=1; return
    end if

    select case ( choice )
    case ( 'e' )   ! exponential correlation
       if( icorlen .gt. 0 )then
          do i=0,4*icorlen-1
             rho(i) = exp(-float(i)/float(icorlen))
          end do
          do i=4*icorlen,5*icorlen-1
             rho(i) = exp(-4.0)*(1.0-float(i-4*icorlen)/float(icorlen))
          end do
       end if
    case ( 'g' )   ! Gaussian correlation
       do i=0,5*icorlen-1
          rho(i) = exp(-(float(i)/float(icorlen))**2)
       end do
    case ( 't' )   ! Thiebaux autoregressive correlation
       do i=0,5*icorlen-1
          rho(i) = (1.0+2.0*float(i)/float(icorlen))* &
               exp( -2.0*float(i)/float(icorlen) )
       end do
    case default
       write (*,'("ERROR - choice `",a,"` invalid")') choice
       write (*,'("ERROR in ",a)') name; status=1; return
    end select

    status = 0

  end subroutine dist_corr_fun

  real function gasdev()
    !
    ! Returns a normally distributed deviate with zero mean and unit variance
    ! From Numerical Recipes in Fortran, 2nd edition, p.280
    !
    implicit none
    !
    integer,save  :: iset=0
    real,save     :: gset
    real  :: fac,rsq,v1,v2,rnd1,rnd2
    !
    if (iset.eq.0) then
1      call random_number(rnd1)
       call random_number(rnd2)
       v1=2.*rnd1-1.
       v2=2.*rnd2-1.
       rsq=v1**2+v2**2
       if (rsq.ge.1..or.rsq.eq.0.) goto 1
       fac=sqrt(-2.*log(rsq)/rsq)
       gset=v1*fac
       gasdev=v2*fac
       iset=1
    else
       gasdev=gset
       iset=0
    endif
    !
  end function gasdev

! <mgv>
  real FUNCTION ran_n(idum)
    IMPLICIT NONE
    INTEGER, PARAMETER :: K4B=selected_int_kind(9)
    INTEGER(K4B), INTENT(INOUT) :: idum

    INTEGER(K4B), PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836
    REAL, SAVE :: am
    INTEGER(K4B), SAVE :: ix=-1,iy=-1,k
    if (idum <= 0 .or. iy < 0) then            ! Initialize.
      am=nearest(1.0,-1.0)/IM
      iy=ior(ieor(888889999,abs(idum)),1)
      ix=ieor(777755555,abs(idum))
      idum=abs(idum)+1                         ! Set idum positive.
    end if
    ix=ieor(ix,ishft(ix,13))                   ! Marsaglia shift sequence with period 232 - 1.
    ix=ieor(ix,ishft(ix,-17))
    ix=ieor(ix,ishft(ix,5))
    k=iy/IQ                                    ! Park-Miller sequence by Schrage’s method,
    iy=IA*(iy-k*IQ)-IR*k                       ! period 231 - 2.
    if (iy < 0) iy=iy+IM
    ran_n=am*ior(iand(IM,ieor(ix,iy)),1)         ! Combine the two generators with masking to
  END FUNCTION ran_n                           ! ensure nonzero value

! </mgv>

end module misctools
