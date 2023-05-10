!###############################################################################
!
! 4D-var IO tools.
!
! Routines:
!
!   check_gridconsistency
!     Compares the grid attributes in a hdf file with the actual values
!     to see if they match.
!
!
!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module Var4D_IO_Tools

  use GO, only : gol, goPr, goErr

  implicit none


  ! --- in/out -----------------------------------

  private

  ! Check consistency of grid definitions
  public :: check_gridconsistency

  ! add time profile description to HDF file:
  public  ::  Time_Profile_Write_to_HDF


  ! --- const ------------------------------------

  character(len=*), parameter   ::  mname = 'Var4D_IO_Tools'


  ! --- types ------------------------------------


  ! --- var --------------------------------------


contains


  !===========================================================================================================


  subroutine check_gridconsistency(nc_id, status)

    use file_hdf,            only: THdfFile, ReadAttribute
    use file_hdf,            only: Done
    use file_netcdf

    use toolbox,             only: escape_tm

    use dims,                only: im, jm, lm
    use dims,                only: nregions, region_name, len_region_name
    use dims,                only: dx,dy
    use dims,                only: xref, yref, tref
    use dims,                only: xbeg, xend, ybeg, yend


    implicit none

    !__IO___________________________________________________________________

    integer, intent(in)                         :: nc_id
    integer, intent(out)                        :: status

    !__LOCAL_VARIABLES______________________________________________________

    character(len=*), parameter        :: rname = mname//'/check_gridconsistency'

    integer                                     :: region

    integer                                     :: nregions_in
    character(len=len_region_name), allocatable :: region_name_in(:)
    integer, dimension(:), allocatable          :: im_in, jm_in, lm_in
    integer                                     :: dx_in, dy_in
    integer, dimension(:), allocatable          :: xref_in, yref_in, tref_in
    integer, dimension(:), allocatable          :: xbeg_in, xend_in, ybeg_in, yend_in


    !__START_SUBROUTINE______________________________________________________

    !write (gol,'(a,": begin")') rname; call goPr

    ! somthing with debug info ...
    status=0

    nregions_in = nc_get_attr(nc_id, 'nregions', status)

    if(nregions_in /= nregions) then
       status=-1
       print *, 'nregions', nregions_in, nregions
    else
       allocate(region_name_in(nregions))
       allocate(im_in(nregions))
       allocate(jm_in(nregions))
       allocate(lm_in(nregions))
       allocate(xref_in(0:nregions))
       allocate(yref_in(0:nregions))
       allocate(tref_in(0:nregions))
       allocate(xbeg_in(nregions))
       allocate(xend_in(nregions))
       allocate(ybeg_in(nregions))
       allocate(yend_in(nregions))

       im_in   = nc_get_attr(nc_id, 'im', status)
       jm_in   = nc_get_attr(nc_id, 'jm', status)
       lm_in   = nc_get_attr(nc_id, 'lm', status)
       xref_in = nc_get_attr(nc_id, 'xref', status)
       yref_in = nc_get_attr(nc_id, 'yref', status)
       tref_in = nc_get_attr(nc_id, 'tref', status)
       xbeg_in = nc_get_attr(nc_id, 'xbeg', status)
       xend_in = nc_get_attr(nc_id, 'xend', status)
       ybeg_in = nc_get_attr(nc_id, 'ybeg', status)
       yend_in = nc_get_attr(nc_id, 'yend', status)


       do region=1, nregions

         ! if(region_name_in(region) /= region_name(region)) then
         !    status=-1
         !    print *, 'region ', region, ' region_name ', region_name_in(:), region_name(:)
         ! endif

          if(im_in(region) /= im(region)) then
             status=-1
             print *, 'region ', region, ' im ', im_in(region), im(region)
          endif

          if(jm_in(region) /= jm(region)) then
             status=-1
             print *, 'region ', region, ' jm ', jm_in(region), jm(region)
          endif

          if(lm_in(region) /= lm(region)) then
             status=-1
             print *, 'region ', region, ' lm ', lm_in(region), lm(region)
          endif

          if(xref_in(region) /= xref(region)) then
             status=-1
             print *, 'region ', region, ' xref ', xref_in(region), xref(region)
          endif

          if(yref_in(region) /= yref(region)) then
             status=-1
             print *, 'region ', region, ' yref ', yref_in(region), yref(region)
          endif

          if(tref_in(region) /= tref(region)) then
             status=-1
             print *, 'region ', region, ' tref ', tref_in(region), tref(region)
          endif

          if(xbeg_in(region) /= xbeg(region)) then
             status=-1
             print *, 'region ', region, ' xbeg ', xbeg_in(region), xbeg(region)
          endif

          if(xend_in(region) /= xend(region)) then
             status=-1
             print *, 'region ', region, ' xend ', xend_in(region), xend(region)
          endif

          if(ybeg_in(region) /= ybeg(region)) then
             status=-1
             print *, 'region ', region, ' ybeg ', ybeg_in(region), ybeg(region)
          endif

          if(yend_in(region) /= yend(region)) then
             status=-1
             print *, 'region ', region, ' yend ', yend_in(region), yend(region)
          endif


       enddo

       deallocate(region_name_in)
       deallocate(im_in)
       deallocate(jm_in)
       deallocate(lm_in)
       deallocate(xref_in)
       deallocate(yref_in)
       deallocate(tref_in)
       deallocate(xbeg_in)
       deallocate(xend_in)
       deallocate(ybeg_in)
       deallocate(yend_in)


    endif

    if(status /= 0) then
       call escape_tm('check_gridconsistency: definitions not consistent')
    endif

    !write (gol,'(a,": end")') rname; call goPr

    ! ok
    status = 0

  end subroutine check_gridconsistency


  ! ***


  !
  ! add time profile info to hdf file
  !
  ! variables written are:
  !   <prefix>_time1(6,nt)
  !   <prefix>_time2(6,nt)
  !

  subroutine Time_Profile_Write_to_HDF( tp, hdf, prefix, status )

    use GO      , only : Get
    use GO      , only : T_Time_Profile
    use file_hdf, only : THdfFile, TSds, Init, Done, WriteData

    ! --- in/out ---------------------------------

    type(T_Time_Profile), intent(in)          ::  tp
    type(THdfFile), intent(inout)             ::  hdf
    character(len=*), intent(in)              ::  prefix
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/Time_Profile_Write_to_HDF'

    ! --- local ----------------------------------

    integer                  ::  i_period
    integer, allocatable     ::  itimes(:,:)
    type(TSds)               ::  Sds

    ! --- begin ----------------------------------

    ! storage for 6 time values per period:
    allocate( itimes(6,tp%n_period) )

    ! fill:
    do i_period = 1, tp%n_period
      call Get( tp%period(i_period)%t1, time6=itimes(:,i_period) )
    end do

    ! write:
    call Init( Sds, hdf, trim(prefix)//'_time1', shape(itimes), 'int', status )
    IF_NOTOK_RETURN(status=1)
    call WriteData( Sds, itimes, status )
    IF_NOTOK_RETURN(status=1)
    call Done( Sds, status )
    IF_NOTOK_RETURN(status=1)

    ! fill:
    do i_period = 1, tp%n_period
      call Get( tp%period(i_period)%t2, time6=itimes(:,i_period) )
    end do

    ! write:
    call Init( Sds, hdf, trim(prefix)//'_time2', shape(itimes), 'int', status )
    IF_NOTOK_RETURN(status=1)
    call WriteData( Sds, itimes, status )
    IF_NOTOK_RETURN(status=1)
    call Done( Sds, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( itimes )

    ! ok:
    status = 0

  end subroutine Time_Profile_Write_to_HDF


end module Var4D_IO_Tools


