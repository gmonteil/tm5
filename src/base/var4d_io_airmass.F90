!###############################################################################
!
! Routines to read/write air mass and pressure arrays.
! Might be useful for debugging.
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

module Var4D_IO_Mass

  use GO, only : gol, goPr, goErr
  use os_specs, only : MAX_FILENAME_LEN

  implicit none


  ! --- in/out -----------------------------------

  private

  ! public routines:
  public :: Var4D_IO_Mass_Init, Var4D_IO_Mass_Done
  public :: save_masses, restore_masses
!  public :: save_pressure, restore_pressure
!  public :: WriteAirMass

  ! public variables:


  ! --- const ------------------------------------

  character(len=*), parameter   ::  mname = ' Var4D_IO_Mass'


  ! --- types ------------------------------------


  ! --- var --------------------------------------

  character(len=MAX_FILENAME_LEN)   :: outdir_mass  ! output dir for air mass

contains


  ! ============================================================================


  subroutine Var4D_IO_Mass_Init( rcF, status )

    use GO,                        only : goPathJoin
    use GO,                        only : TrcFile, ReadRc
    use global_data,               only : outdir

    ! --- in/out ----------------------------------------------

    type(TrcFile), intent(in)           ::  rcF
    integer, intent(out)                ::  status

    ! --- const -----------------------------------------------

    character(len=*), parameter         ::  rname = mname//'/Var4D_IO_Mass_Init'

    ! --- local -----------------------------------------------

    character(len=MAX_FILENAME_LEN)     :: subdir

    ! --- begin -----------------------------------------------

    ! output directory for air mass
    call ReadRc( rcF, 'mass.output.subdir', subdir, status )
    IF_NOTOK_RETURN(status=1)
    call goPathJoin( outdir, subdir, outdir_mass, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Var4D_IO_Mass_Init


  ! ***


  subroutine Var4D_IO_Mass_Done( status )

    ! --- in/out ----------------------------------------------

    integer, intent(out)                ::  status

    ! --- const -----------------------------------------------

    character(len=*), parameter         ::  rname = mname//'/Var4D_IO_Mass_Done'

    ! --- local -----------------------------------------------

    ! --- begin -----------------------------------------------

    ! ok
    status = 0

  end subroutine Var4D_IO_Mass_Done


  !==============================================================================================
  !==============================================================================================

  subroutine save_masses( status )

    use dims,                       only : nregions, idate, region_name
    use MeteoData                 , only : m_dat
    use misctools,                  only : check_dir
    use toolbox,                    only : escape_tm
    use file_netcdf

    ! --- in/out ----------------------------------------------

    integer, intent(out)                        ::  status

    ! --- const -----------------------------------------------

    character(len=*), parameter        :: rname = mname//'/save_masses'

    ! --- local -----------------------------------------------

    integer                         :: region, dimx, dimy, dimz
    integer, dimension(3)           :: lowbound, upbound
    real, dimension(:,:,:), pointer :: m
!    integer                         :: sfstart, io, sfend, sfcreate
!    integer                         :: sds_id, istat, sfsnatt, sfendacc, sfwdata
!    integer                         :: dimid0, dimid1, dimid2, sfsdmname, sfdimid
    integer                         :: nc_id, grp_id

    ! --- begin -----------------------------------------------

    call check_dir(trim(outdir_mass)//'/m_final_forward.nc4')

    nc_id = nc_open(trim(outdir_mass)//'/m_final_forward.nc4', 'c', status)
    IF_NOTOK_RETURN(status=1)

    do region = 1, nregions
        m => m_dat(region)%data
        lowbound = lbound(m)
        upbound = ubound(m)
        dimx = upbound(1)-lowbound(1) + 1
        dimy = upbound(2)-lowbound(2) + 1
        dimz = upbound(3)-lowbound(3) + 1
        grp_id = nc_create_group(nc_id, trim(region_name(region)))
        call nc_create_dim(grp_id, 'dimx', dimx)
        call nc_create_dim(grp_id, 'dimy', dimy)
        call nc_create_dim(grp_id, 'dimz', dimz)
        call nc_dump_var(grp_id, 'm', (/'dimx','dimy','dimz'/), m(:,:,:))
        nullify(m)
    end do ! region
    call nc_close(nc_id)

!    do region = 1, nregions
!       m=> mass_dat(region)%m_t
!       lowbound = lbound(m)
!       upbound = ubound(m)
!       dimx = upbound(1)-lowbound(1) + 1
!       dimy = upbound(2)-lowbound(2) + 1
!       dimz = upbound(3)-lowbound(3) + 1
!       io = sfstart(trim(outdir_mass)//'/m_final_forward'//trim(region_name(region))//'.hdf', DFACC_CREATE)
!       sds_id = sfcreate(io,'m', DFNT_FLOAT64, 3, (/dimx,dimy,dimz/))
!       istat =  sfsnatt(sds_id,'idate',  DFNT_INT32, 6, idate)
!       dimid2 = sfdimid(sds_id, 0)
!       istat = sfsdmname(dimid2,'dimx')
!       dimid1 = sfdimid(sds_id, 1)
!       istat = sfsdmname(dimid1,'dimy')
!       dimid0 = sfdimid(sds_id, 2)
!       istat = sfsdmname(dimid0,'dimz')

!       istat  = sfwdata( sds_id, (/0,0,0/), (/1,1,1/), (/dimx,dimy,dimz/), m)
!       istat = sfendacc(sds_id)

!       istat = sfend(io)
!       if (istat == SUCCEED) then
!          print *, 'Created file m_final_forward for region ', region
!       else
!          call escape_tm('Error in creating m_final_forward')
!       endif
!       nullify(m)
!    enddo

    ! ok
    status = 0

  end subroutine save_masses

  !==============================================================================================
  !==============================================================================================

  subroutine restore_masses( status, only_region )

    use dims,                       only : nregions, idate, region_name
    use MeteoData                 , only : m_dat
!    use io_hdf,                     only : DFACC_READ, DFNT_FLOAT64, DFNT_INT32, SUCCEED
    use toolbox,                    only : escape_tm
    use file_netcdf

    ! --- in/out ----------------------------------------------

    integer, intent(out)             ::  status
    integer, intent(in), optional    ::  only_region

    ! --- const -----------------------------------------------

    character(len=*), parameter        :: rname = mname//'/restore_masses'

    ! --- local -----------------------------------------------

    integer                         :: region, dimx, dimy, dimz, i
    integer, dimension(3)           :: lowbound, upbound
    real, dimension(:,:,:), pointer :: m
    real, allocatable, dimension(:,:,:) :: dummy_arr
    integer                         :: nc_id, grp_id

    ! --- begin -----------------------------------------------

    ! open file:
    nc_id = nc_open(trim(outdir_mass)//'/m_final_forward.nc4', 'r', status)
    IF_NOTOK_RETURN(status=1)
    
    ! loop over regions:
    do region = 1, nregions

      ! skip?
      if ( present(only_region) ) then
        if ( region /= only_region ) exit
      end if

      ! mass array:
      m => m_dat(region)%data
      ! dims:
        lowbound = lbound(m)
        upbound = ubound(m)
        dimx = upbound(1)-lowbound(1) + 1
        dimy = upbound(2)-lowbound(2) + 1
        dimz = upbound(3)-lowbound(3) + 1
      ! open group:
        grp_id = nc_get_group(nc_id, trim(region_name(region)))
      ! check:
        if (dimx /= nc_get_dim(grp_id, 'dimx') .or. dimy /= nc_get_dim(grp_id, 'dimy') .or. dimz /= nc_get_dim(grp_id, 'dimz')) then
            write(0, '(a, " :: shape mismatch trying to read mass file for region ", a)') rname, trim(region_name(region))
            write(0, '("Expected array of shape ", i3, " x ", i3, " x ", i3)') dimx, dimy, dimz
            write(0, '("Encountered array of shape ", i3, " x ", i3, " x ", i3)') nc_get_dim(grp_id, 'dimx'), nc_get_dim(grp_id, 'dimy'), nc_get_dim(grp_id, 'dimz')
            call nc_close(nc_id)
            TRACEBACK; status=1; return
        end if
      ! read:
        dummy_arr = nc_read_var(grp_id, 'm')
      ! store:
        m(:,:,:) = dummy_arr(:,:,:)
      ! clear:
        if (allocated(dummy_arr)) deallocate(dummy_arr)
        nullify(m)

    end do ! region
    ! close:
    call nc_close(nc_id)

!    do region = 1, nregions
!       m=> mass_dat(region)%m_t
!       lowbound = lbound(m)
!       upbound = ubound(m)
!       dimx = upbound(1)-lowbound(1) + 1
!       dimy = upbound(2)-lowbound(2) + 1
!       dimz = upbound(3)-lowbound(3) + 1
!       write (fname,'(a,"/m_final_forward",a,".hdf")') trim(indir_mass), trim(region_name(region))
!       io = sfstart( trim(fname), DFACC_READ )
!       if ( io < 0 ) then
!         write (gol,'("could not open hdf file : ",a)') trim(fname); call goErr
!         TRACEBACK; status=1; return
!       end if
!       idx = sfn2index(io,'m')
!       if ( idx < 0 ) then
!         write (gol,'("could not find dataset `m` in hdf file : ",a)') trim(fname); call goErr
!         TRACEBACK; status=1; return
!       end if
!       sds_id = sfselect(io, idx)
!       istat = sfginfo(sds_id, xname, rank, dim_sizes, num_type, attributes)
!       if (rank /= 3 .or. xname /= 'm' .or. num_type /= DFNT_FLOAT64) then
!         write (gol,'("properties of `m` from `",a,"` do not match")') trim(fname); call goErr
!         TRACEBACK; status=1; return
!       endif
!       if (dim_sizes(1) /= dimx .or. dim_sizes(2) /= dimy .or. dim_sizes(3) /= dimz) then
!         write (gol,'("dimension of `m` from `",a,"` do not match")') trim(fname); call goErr
!         TRACEBACK; status=1; return
!       endif
!       attr_index = sffattr(sds_id, 'idate')
!       istat = sfgainfo(sds_id, attr_index, attr_name, data_type,  n_values)
!       if(data_type == DFNT_INT32 .and. n_values == 6) then
!          istat = sfrnatt(sds_id, attr_index, idate_file)
!       else
!          call escape_tm('Failed to read date when restoring m')
!       endif


!       if ( any((idate-idate_file) /= 0)) then
!          write (gol,'("Date in restore mass file is wrong:")'); call goErr
!          write (gol,'("  idate          : ",i4,5i3)') idate; call goErr
!          write (gol,'("  idate in file  : ",i4,5i3)') idate_file; call goErr
!          TRACEBACK; status=1; return
!       endif

!       istat  = sfrdata( sds_id, (/0,0,0/), (/1,1,1/), (/dimx,dimy,dimz/), m)

!       if (istat == SUCCEED) then
!          print *, 'Restored m from file m_final_forward for region ', region
!       else
!          call escape_tm('Error in reading the m_final_forward')
!       endif
!       istat = sfendacc(sds_id)
!       istat = sfend(io)

!       nullify(m)
!    enddo

    ! ok
    status = 0

  end subroutine restore_masses

!  !==============================================================================================
!  !==============================================================================================
!
!  subroutine save_pressure(region)
!    use dims,                       only : idate, region_name
!    use global_data,                only : mass_dat
!    use io_hdf,                     only : DFACC_CREATE, DFNT_FLOAT64, DFNT_INT32, SUCCEED
!    use toolbox,                    only : escape_tm
!    implicit none
!    integer, intent(IN)             :: region
!    integer                         :: dimx, dimy
!    integer, dimension(2)           :: lowbound, upbound
!    real, dimension(:,:), pointer   :: ps
!    integer                         :: sfstart, io, sfend, sfcreate
!    integer                         :: sds_id, istat, sfsnatt, sfendacc, sfwdata
!    integer                         :: dimid2, dimid1, sfsdmname, sfdimid
!
!
!
!    ps=> mass_dat(region)%p
!    lowbound = lbound(ps)
!    upbound = ubound(ps)
!    dimx = upbound(1)-lowbound(1) + 1
!    dimy = upbound(2)-lowbound(2) + 1
!    io = sfstart('p_final_forward'//trim(region_name(region))//'.hdf', DFACC_CREATE)
!    sds_id = sfcreate(io,'ps', DFNT_FLOAT64, 2, (/dimx,dimy/))
!    istat =  sfsnatt(sds_id,'idate',  DFNT_INT32, 6, idate)
!    dimid2 = sfdimid(sds_id, 0)
!    istat = sfsdmname(dimid2,'dimx')
!    dimid1 = sfdimid(sds_id, 1)
!    istat = sfsdmname(dimid1,'dimy')
!
!    istat  = sfwdata( sds_id, (/0,0/), (/1,1/), (/dimx,dimy/), ps)
!    istat = sfendacc(sds_id)
!
!    istat = sfend(io)
!    if (istat == SUCCEED) then
!       print *, 'Created file p_final_forward for region ', region
!    else
!       call escape_tm('Error in creating p_final_forward')
!    endif
!    nullify(ps)
!
!
!  end subroutine save_pressure
!
!  !==============================================================================================
!  !==============================================================================================
!
!  subroutine restore_pressure(region)
!    use dims,                       only : idate, region_name
!    use global_data,                only : mass_dat
!    use io_hdf,                     only : DFACC_READ, DFNT_FLOAT64, DFNT_INT32, SUCCEED
!    use toolbox,                    only : escape_tm
!    implicit none
!    integer, intent(in)             :: region
!
!    integer                         :: dimx, dimy, i
!    integer, dimension(2)           :: lowbound, upbound
!    real, dimension(:,:), pointer   :: ps
!    integer                         :: sfstart, io, sfend
!    integer                         :: sds_id, istat, sfsnatt, sfendacc, sfwdata
!    integer                         :: idx
!    integer, dimension(6)           :: idate_file
!
!    integer, parameter   :: MAX_VAR_DIMS = 32
!    character(len=64)    :: xname, attr_name
!    integer              :: rank
!    integer              :: attributes, num_type
!    integer              :: sffinfo, sfselect, sfginfo
!    integer              :: sfrnatt
!    integer              :: sfrdata, sfn2index
!    integer,dimension(MAX_VAR_DIMS) :: dim_sizes
!    integer              :: attr_index, sffattr, data_type, n_values
!    integer              :: sfgainfo
!
!
!
!    ps => mass_dat(region)%p
!    lowbound = lbound(ps)
!    upbound = ubound(ps)
!    dimx = upbound(1)-lowbound(1) + 1
!    dimy = upbound(2)-lowbound(2) + 1
!    io = sfstart('p_final_forward'//trim(region_name(region))//'.hdf', DFACC_READ)
!    idx = sfn2index(io,'ps')
!    if(idx == -1) then
!       call escape_tm('Failed to restore ps from file')
!    endif
!    sds_id = sfselect(io, idx)
!    istat = sfginfo(sds_id, xname, rank, dim_sizes, num_type, attributes)
!    if (rank /= 2 .or. xname /= 'ps' .or. num_type /= DFNT_FLOAT64) then
!       call escape_tm('Failed to restore ps from file')
!    endif
!    if (dim_sizes(1) /= dimx .or. dim_sizes(2) /= dimy ) then
!       call escape_tm('Failed to restore ps from file')
!    endif
!    attr_index = sffattr(sds_id, 'idate')
!    istat = sfgainfo(sds_id, attr_index, attr_name, data_type,  n_values)
!    if(data_type == DFNT_INT32 .and. n_values == 6) then
!       istat = sfrnatt(sds_id, attr_index, idate_file)
!    else
!       call escape_tm('Failed to read date when restoring ps')
!    endif
!
!
!    if(any((idate-idate_file) /= 0)) then
!       call escape_tm('Date in restore pressure file is wrong')
!    endif
!
!    istat  = sfrdata( sds_id, (/0,0/), (/1,1/), (/dimx,dimy/), ps)
!
!    if (istat == SUCCEED) then
!       print *, 'Restored ps from file p_final_forward for region ', region
!    else
!       call escape_tm('Error in reading the p_final_forward')
!    endif
!    istat = sfendacc(sds_id)
!    istat = sfend(io)
!
!    nullify(ps)
!
!  end subroutine restore_pressure
!
!  !==============================================================================================
!  !==============================================================================================
!
!  subroutine WriteAirMass( fname_base )
!
!    ! Writes air mass field to ascii file (one file per region).
!    ! Useful for testing whether air mass in forward and adjoint
!    ! mode is equivalent.
!
!    use dims,                       only : nregions, im, jm, lm
!    use global_data,                only : mass_dat
!
!    character(len=*), intent(in)    :: fname_base
!
!    integer, parameter              :: kout = 777
!    integer                         :: region, i, j, l, imr, jmr, lmr
!    integer, dimension(3)           :: lowbound, upbound
!    real, dimension(:,:,:), pointer :: m
!    character(len=80)               :: fname
!
!    do region = 1, nregions
!       imr = im(region)
!       jmr = jm(region)
!       lmr = lm(region)
!       write( fname, '(a,"_",i2.2,".dat")' ) trim(fname_base), region
!       open( unit=kout, form='formatted', file=fname )
!       m=> mass_dat(region)%m_t
!       lowbound = lbound(m)
!       upbound = ubound(m)
!!       do i = lowbound(1), upbound(1)
!!          do j = lowbound(2), upbound(2)
!!             do l = lowbound(3), upbound(3)
!       do i = 1, imr
!          do j = 1, jmr
!             do l = 1, lmr
!                write( kout, '(3i3.3,e20.10)' ) i, j, l, m(i,j,l)
!             end do
!          end do
!       end do
!       nullify(m)
!       close( kout )
!    enddo
!
!  end subroutine WriteAirMass
!
!  !==============================================================================================
!  !==============================================================================================


end module Var4D_IO_Mass

