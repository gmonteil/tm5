!###############################################################################
!
! contains routines to read emissions and
! to read and write the main model state from/to file
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

module io_save_netcdf

  use GO, only : gol, goPr, goErr

  implicit none

  private

  public :: savenetcdf
  public :: readnetcdf
  public :: save_filename

  character(len=*), parameter :: mname = 'io_save_netcdf'

contains

  subroutine readnetcdf(file_name, status)
    ! Read in files written out by savenetcdf
    use dims,        only : region_name, im, jm, lm, nregions, kmain, itau
    use global_data, only : mass_dat
    use MeteoData  , only : sp_dat, m_dat
    use chem_param,  only : ntrace, ntracet
    use toolbox,     only : escape_tm
    use global_data, only : rcF
    use Go,          only : ReadRc
    use datetime,    only : tstamp
    use file_netcdf
    use os_specs,    only : WRITE_STR_LEN

    implicit none

    ! in/out
    character(len=*), intent(in)    :: file_name
    integer, intent(out)            :: status

    ! local
    character(len=*), parameter :: rname = mname//'/readnetcdf'

    real,dimension(:,:,:,:),pointer :: rm,rxm,rym,rzm
#ifdef secmom
    real,dimension(:,:,:,:),pointer :: rxxm,rxym,rxzm,ryym,ryzm,rzzm
#endif
    real,dimension(:,:,:),pointer   :: m
    real, dimension(:,:,:), pointer :: sp
    integer :: region, nc_id, grp_id, imr, jmr, lmr, n, dimx, dimy, dimz, xub, xlb, yub, ylb
    real, allocatable :: field4d(:,:,:,:), saved_mass(:,:,:), field2d(:,:)
    logical :: read_full ! whether to read the full rm, m etc arrays or only 1:imr, 1:jmr, etc.
    logical :: correct_mix ! whether to make mixing ratios continuous or not
    logical :: overwrite_mass ! whether to overwrite mass and surface pressure
    character(len=WRITE_STR_LEN)  :: write_string

    call ReadRc( rcF, 'restart.correct.mixing.ratio', correct_mix, status, default=.true.)
    IF_ERROR_RETURN(status=1)
    call ReadRc( rcF, 'restart.overwrite.mass', overwrite_mass, status, default=.false.)
    IF_ERROR_RETURN(status=1)

    nc_id = nc_open(file_name, 'r', status)
    IF_NOTOK_RETURN(status=1)

    do region = 1, nregions
        imr = im(region)
        jmr = jm(region)
        lmr = lm(region)

        m  =>  m_dat(region)%data
        sp => sp_dat(region)%data
        rm => mass_dat(region)%rm_t
        rxm => mass_dat(region)%rxm_t
        rym => mass_dat(region)%rym_t
        rzm => mass_dat(region)%rzm_t

        grp_id = nc_get_group(nc_id, region_name(region))
        saved_mass = nc_read_var(grp_id, 'm')

        ! Is this with or without padding?
        if ((size(saved_mass, 1) == imr) .and. (size(saved_mass, 2) == jmr)) then
            read_full = .false.
        else if ((size(saved_mass, 1) == imr+4) .and. (size(saved_mass, 2) == jmr+4)) then
            read_full = .true.
        else
            write(0, '("For region ", i1, ", lateral size of mass in save file is ", i3, " x ", i3)') &
                region, size(saved_mass, 1), size(saved_mass, 2)
            write(0, '("    This is neither ", i3, " x ", i3, ", nor ", i3, " x ", i3)') imr, jmr, imr+4, jmr+4
            status = 1
            IF_NOTOK_RETURN(status=1)
        end if

        if (read_full) then
            xlb = -1
            xub = imr+2
            ylb = -1
            yub = jmr+2
        else
            xlb = 1
            xub = imr
            ylb = 1
            yub = jmr
        end if

        if (overwrite_mass) then
            m(xlb:xub, ylb:yub, 1:lmr) = saved_mass
            field2d = nc_read_var(grp_id, 'sp')
            sp(xlb:xub, ylb:yub,1) = field2d
            deallocate(field2d)
        end if

        field4d = nc_read_var(grp_id, 'rm')
        rm(xlb:xub,ylb:yub,1:lmr,1:ntracet) = field4d(:,:,:,:)
        deallocate(field4d)

        field4d = nc_read_var(grp_id, 'rxm')
        rxm(xlb:xub,ylb:yub,1:lmr,1:ntracet) = field4d(:,:,:,:)
        deallocate(field4d)

        field4d = nc_read_var(grp_id, 'rym')
        rym(xlb:xub,ylb:yub,1:lmr,1:ntracet) = field4d(:,:,:,:)
        deallocate(field4d)

        field4d = nc_read_var(grp_id, 'rzm')
        rzm(xlb:xub,ylb:yub,1:lmr,1:ntracet) = field4d(:,:,:,:)
        deallocate(field4d)

        print*, 'init rm:', sum(rm), sum(rxm), sum(rzm), sum(rym), sum(m), sum(saved_mass)

        ! Now scale the rm, rxm, rym and rzm by the ratio of masses
        if (correct_mix) then
            do n = 1, ntracet
                rm(xlb:xub,ylb:yub,1:lmr,n) = rm(xlb:xub,ylb:yub,1:lmr,n) * m(xlb:xub,ylb:yub,1:lmr) / saved_mass(xlb:xub,ylb:yub,1:lmr)
                rxm(xlb:xub,ylb:yub,1:lmr,n) = rxm(xlb:xub,ylb:yub,1:lmr,n) * m(xlb:xub,ylb:yub,1:lmr) / saved_mass(xlb:xub,ylb:yub,1:lmr)
                rym(xlb:xub,ylb:yub,1:lmr,n) = rym(xlb:xub,ylb:yub,1:lmr,n) * m(xlb:xub,ylb:yub,1:lmr) / saved_mass(xlb:xub,ylb:yub,1:lmr)
                rzm(xlb:xub,ylb:yub,1:lmr,n) = rzm(xlb:xub,ylb:yub,1:lmr,n) * m(xlb:xub,ylb:yub,1:lmr) / saved_mass(xlb:xub,ylb:yub,1:lmr)
            end do
        end if

        print*, 'init rm:', sum(rm)

!        ! SBi
!        write(*, '("Maximum mismatch between saved and model airmass = ", es20.12)') &
!            maxval(m(1:imr,1:jmr,1:lmr) - saved_mass(1:imr,1:jmr,1:lmr))
!        write(*, '("Minimum mismatch between saved and model airmass = ", es20.12)') &
!            minval(m(1:imr,1:jmr,1:lmr) - saved_mass(1:imr,1:jmr,1:lmr))
!        write(*, '("Total saved airmass = ", es20.12)') sum(saved_mass(1:imr,1:jmr,1:lmr))
!        write(*, '("Total model airmass = ", es20.12)') sum(m(1:imr,1:jmr,1:lmr))
!        write(*, '("Maximum deviation from unit ratio = ", es20.12)') maxval(m(1:imr,1:jmr,1:lmr) / saved_mass(1:imr,1:jmr,1:lmr) - 1.0)
!        write(*, '("Minimum deviation from unit ratio = ", es20.12)') minval(m(1:imr,1:jmr,1:lmr) / saved_mass(1:imr,1:jmr,1:lmr) - 1.0)
!        ! SBf

        nullify(rm,rxm,rym,rzm,sp)
#ifdef secmom
        rxxm => mass_dat(region)%rxxm_t
        rxym => mass_dat(region)%rxym_t
        rxzm => mass_dat(region)%rxzm_t
        ryym => mass_dat(region)%ryym_t
        ryzm => mass_dat(region)%ryzm_t
        rzzm => mass_dat(region)%rzzm_t

        field4d = nc_read_var(grp_id, 'rxxm')
        rxxm(xlb:xub,ylb:yub,1:lmr,1:ntracet) = field4d(:,:,:,:)
        deallocate(field4d)

        field4d = nc_read_var(grp_id, 'rxym')
        rxym(xlb:xub,ylb:yub,1:lmr,1:ntracet) = field4d(:,:,:,:)
        deallocate(field4d)

        field4d = nc_read_var(grp_id, 'rxzm')
        rxzm(xlb:xub,ylb:yub,1:lmr,1:ntracet) = field4d(:,:,:,:)
        deallocate(field4d)

        field4d = nc_read_var(grp_id, 'ryym')
        ryym(xlb:xub,ylb:yub,1:lmr,1:ntracet) = field4d(:,:,:,:)
        deallocate(field4d)

        field4d = nc_read_var(grp_id, 'ryzm')
        ryzm(xlb:xub,ylb:yub,1:lmr,1:ntracet) = field4d(:,:,:,:)
        deallocate(field4d)

        field4d = nc_read_var(grp_id, 'rzzm')
        rzzm(xlb:xub,ylb:yub,1:lmr,1:ntracet) = field4d(:,:,:,:)
        deallocate(field4d)

        if (correct_mix) then
            do n = 1, ntracet
                rxxm(xlb:xub,ylb:yub,1:lmr,n) = rxxm(xlb:xub,ylb:yub,1:lmr,n) * m(xlb:xub,ylb:yub,1:lmr) / saved_mass(xlb:xub,ylb:yub,1:lmr)
                rxym(xlb:xub,ylb:yub,1:lmr,n) = rxym(xlb:xub,ylb:yub,1:lmr,n) * m(xlb:xub,ylb:yub,1:lmr) / saved_mass(xlb:xub,ylb:yub,1:lmr)
                rxzm(xlb:xub,ylb:yub,1:lmr,n) = rxzm(xlb:xub,ylb:yub,1:lmr,n) * m(xlb:xub,ylb:yub,1:lmr) / saved_mass(xlb:xub,ylb:yub,1:lmr)
                ryym(xlb:xub,ylb:yub,1:lmr,n) = ryym(xlb:xub,ylb:yub,1:lmr,n) * m(xlb:xub,ylb:yub,1:lmr) / saved_mass(xlb:xub,ylb:yub,1:lmr)
                ryzm(xlb:xub,ylb:yub,1:lmr,n) = ryzm(xlb:xub,ylb:yub,1:lmr,n) * m(xlb:xub,ylb:yub,1:lmr) / saved_mass(xlb:xub,ylb:yub,1:lmr)
                rzzm(xlb:xub,ylb:yub,1:lmr,n) = rzzm(xlb:xub,ylb:yub,1:lmr,n) * m(xlb:xub,ylb:yub,1:lmr) / saved_mass(xlb:xub,ylb:yub,1:lmr)
            end do
        end if
        nullify(rxxm,rxym,rxzm,ryzm,ryym,rzzm)
#endif
        deallocate(saved_mass)
        nullify(m)
    end do
    call nc_close(nc_id)

    write(write_string, '("Starting masses read in from ", a)') trim(file_name)
    call tstamp(kmain, itau, trim(write_string))
    status = 0

  end subroutine readnetcdf

  subroutine save_filename(file_name, status)
    use global_data,    only : RcF, outdir
    use go,             only : pathsep, ReadRc
    use dims,           only : idate, revert
    use os_specs,       only : MAX_FILENAME_LEN, MAX_RCKEY_LEN

    implicit none

    character(len=MAX_FILENAME_LEN), intent(out) :: file_name
    character(len=MAX_RCKEY_LEN)    :: save_subdir
    integer, intent(out)            :: status

    character(len=*), parameter     :: rname = mname//'/save_filename'

    call ReadRc(rcF, 'save.output.subdir', save_subdir, status, default='save')
    IF_ERROR_RETURN(status=1)

    if (revert == 1) then
        write( file_name, '(4a,"save_",i4,4i2.2,".nc4")' ) trim(outdir), pathsep, trim(save_subdir), pathsep, idate(1:5)
    else
        write( file_name, '(4a,"adj_save_",i4,4i2.2,".nc4")' ) trim(outdir), pathsep, trim(save_subdir), pathsep, idate(1:5)
    end if

  end subroutine save_filename

  subroutine savenetcdf(file_name, status)
    ! Save masses or adjoint masses at the end of a run
    ! Identical to savehdf, except that this writes a netcdf-4 file
    use dims,        only : nregions, im, jm, lm, region_name, kmain, itau
    use global_data, only : mass_dat, rcF
    use MeteoData  , only : sp_dat, m_dat
    use Go, only          : ReadRc
    use chem_param
    use file_netcdf
    use datetime,    only : tstamp
    use toolbox,     only : escape_tm
    use netcdf,      only : nf90_put_att, NF90_GLOBAL
    use os_specs,    only : WRITE_STR_LEN

    implicit none
    ! in/out
    character(len=*),intent(in) :: file_name
    integer, intent(out)        :: status

    ! local
    character(len=*), parameter :: rname = mname//'/savenetcdf'

    real,dimension(:,:,:,:),pointer :: rm,rxm,rym,rzm
#ifdef secmom
    real,dimension(:,:,:,:),pointer :: rxxm,rxym,rxzm,ryym,ryzm,rzzm
#endif
    real,dimension(:,:,:),  pointer :: m
    real, dimension(:,:,:), pointer :: sp ! surface pressure
    integer :: io_status, region, imr, jmr, lmr, nc_id, grp_id
    integer :: dimx, dimy, dimz, xlb, xub, ylb, yub
    logical :: store_full
    character(len=WRITE_STR_LEN) :: write_string

    call ReadRc(rcF, 'restart.file.contains.padded.mass', store_full, status, default=.false.)
    IF_ERROR_RETURN(status=1)

    nc_id = nc_open(file_name, 'c', status)
    IF_NOTOK_RETURN(status=1)
    
    ! Write the global attributes first
    io_status = nf90_put_att(nc_id, NF90_GLOBAL, 'ntrace', ntrace)
    io_status = nf90_put_att(nc_id, NF90_GLOBAL, 'ra', ra)
    io_status = nf90_put_att(nc_id, NF90_GLOBAL, 'fscale', fscale)
    ! Create the global dimensions
    call nc_create_dim(nc_id, 'ntracet', ntracet)
    ! Now loop over the regions
    do region = 1, nregions
        grp_id = nc_create_group(nc_id, region_name(region))

        m  =>  m_dat(region)%data
        sp => sp_dat(region)%data
        rm => mass_dat(region)%rm_t
        rxm => mass_dat(region)%rxm_t
        rym => mass_dat(region)%rym_t
        rzm => mass_dat(region)%rzm_t

        imr = im(region)
        jmr = jm(region)
        lmr = lm(region)

        ! The dimensions
        if (store_full) then
            dimx = imr+4
            dimy = jmr+4
        else
            dimx = imr
            dimy = jmr
        end if
        dimz = lmr
        call nc_create_dim(grp_id, 'dimx', dimx)
        call nc_create_dim(grp_id, 'dimy', dimy)
        call nc_create_dim(grp_id, 'dimz', dimz)

        if (store_full) then
            xlb = -1
            xub = imr+2
            ylb = -1
            yub = jmr+2
        else
            xlb = 1
            xub = imr
            ylb = 1
            yub = jmr
        end if

        ! The variables
        call nc_dump_var(grp_id, 'm', (/'dimx','dimy','dimz'/), m(xlb:xub, ylb:yub, 1:lmr))
        call nc_dump_var(grp_id,  'rm', (/'dimx   ','dimy   ','dimz   ','ntracet'/),  rm(xlb:xub,ylb:yub,1:lmr,1:ntracet))
        call nc_dump_var(grp_id, 'rxm', (/'dimx   ','dimy   ','dimz   ','ntracet'/), rxm(xlb:xub,ylb:yub,1:lmr,1:ntracet))
        call nc_dump_var(grp_id, 'rym', (/'dimx   ','dimy   ','dimz   ','ntracet'/), rym(xlb:xub,ylb:yub,1:lmr,1:ntracet))
        call nc_dump_var(grp_id, 'rzm', (/'dimx   ','dimy   ','dimz   ','ntracet'/), rzm(xlb:xub,ylb:yub,1:lmr,1:ntracet))
        call nc_dump_var(grp_id, 'sp', (/'dimx', 'dimy'/), sp(xlb:xub, ylb:yub,1) )

        nullify(m, rm, rxm, rym, rzm, sp)
#ifdef secmom
        rxxm => mass_dat(region)%rxxm_t
        rxym => mass_dat(region)%rxym_t
        rxzm => mass_dat(region)%rxzm_t
        ryym => mass_dat(region)%ryym_t
        ryzm => mass_dat(region)%ryzm_t
        rzzm => mass_dat(region)%rzzm_t
        call nc_dump_var(grp_id, 'rxxm', (/'dimx   ','dimy   ','dimz   ','ntracet'/), rxxm(xlb:xub,ylb:yub,1:lmr,1:ntracet))
        call nc_dump_var(grp_id, 'rxym', (/'dimx   ','dimy   ','dimz   ','ntracet'/), rxym(xlb:xub,ylb:yub,1:lmr,1:ntracet))
        call nc_dump_var(grp_id, 'rxzm', (/'dimx   ','dimy   ','dimz   ','ntracet'/), rxzm(xlb:xub,ylb:yub,1:lmr,1:ntracet))
        call nc_dump_var(grp_id, 'ryym', (/'dimx   ','dimy   ','dimz   ','ntracet'/), ryym(xlb:xub,ylb:yub,1:lmr,1:ntracet))
        call nc_dump_var(grp_id, 'ryzm', (/'dimx   ','dimy   ','dimz   ','ntracet'/), ryzm(xlb:xub,ylb:yub,1:lmr,1:ntracet))
        call nc_dump_var(grp_id, 'rzzm', (/'dimx   ','dimy   ','dimz   ','ntracet'/), rzzm(xlb:xub,ylb:yub,1:lmr,1:ntracet))
        nullify(rxxm, rxym, rxzm, ryym, ryzm, rzzm)
#endif
    end do ! region
    call nc_close(nc_id)

    write(write_string, '("Ending masses saved to ", a)') trim(file_name)
    call tstamp(kmain, itau, trim(write_string))

    status = 0

  end subroutine savenetcdf

end module io_save_netcdf
