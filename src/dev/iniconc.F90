#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if

module iniconc_module

    use dims,           only : nregions, region_name, im, jm, lm, adv_scheme, itaur
    use GO,             only : readrc, gol, goerr
    use global_data,    only : rcf, mass_dat
    use file_netcdf
    use zoom_tools,     only : update_parent
    use grid_type_ll,   only : init_grid => init, tllgridinfo
    use grid_type_hyb,  only : init_levels => init, tlevelinfo
    use grid_3d,        only : regrid_3d_mix
    use chem_param,     only : ntracet, names, mixrat_unit
    use tm5_geometry,   only : lli, levi
    use datetime,       only : tau2date
    use MeteoData,      only : m_dat
    use chem_param,     only : fscale

    implicit none

    public :: read_iniconc_fitic
    private

    character(len=*), parameter :: mname = 'iniconc_module'

    contains

        subroutine read_iniconc_fitic(status)

            integer, intent(out)            :: status
            integer                         :: itrac
            character(len=200)              :: fname
            character(len=10)               :: ftype
            integer                         :: region
            character(len=*), parameter     :: rname = mname//'/read_iniconc_fitic'
            
            status = 0

            do itrac = 1, ntracet
                call readrc(rcf, 'start.' // trim(names(itrac)) // '.type', ftype, status)
                select case (trim(ftype))
                    case ('cams')
                        call read_iniconc_cams(itrac, trim(names(itrac)), status)
                    case ('constant')
                        call set_constant_iniconc(itrac, trim(names(itrac)), status)
                end select
                do region = 1, nregions
                    mass_dat(region)%rm_t = mass_dat(region)%rm_t / mixrat_unit(itrac)
                    if ( adv_scheme == 'slope' ) then
                        mass_dat(region)%rxm_t(:, :, :, itrac) = 0.0
                        mass_dat(region)%rym_t(:, :, :, itrac) = 0.0
                        mass_dat(region)%rzm_t(:, :, :, itrac) = 0.0
                    end if
                enddo
            enddo

            do region = nregions,2,-1
                call update_parent(region)
            end do

        end subroutine read_iniconc_fitic


        subroutine read_iniconc_cams(itrac, tracname, status)

            use ISO_FORTRAN_ENV

            character(len=*), intent(in)    :: tracname
            integer, intent(in)             :: itrac
            integer, intent(out)            :: status

            ! netCDF variables
            character(len=200)      :: filename
            integer                 :: ncf
            real(kind=4), dimension(:), allocatable             :: hyai, hybi
            integer(int16), dimension(:, :, :, :), allocatable  :: mixglo1x1_packed
            real, dimension(:, :, :, :), allocatable            :: mixglo1x1
            real, dimension(:), allocatable                     :: time 
            integer(int16), dimension(:, :, :), allocatable     :: ps_packed
            real, dimension(:,:,:), allocatable                 :: ps
            integer, dimension(6)   :: idate
            integer                 :: nhours_since_start_of_month
            integer                 :: itime
            integer                 :: region, nlev
            type(tllgridinfo)       :: hor_grid_in
            type(tlevelinfo)        :: ver_grid_in
            character(len=*), parameter    :: rname = mname//'/read_iniconc_cams'

            !!MVO, 2024-12-17:
            ! CAMS initial concentration file uses packed format.
            ! so the physical data needs to be determined as     
            ! unpacked_value = packed_value * scale_factor + add_offset
!            type(nc_uni_att) :: attr
            real(kind=4) :: add_offset
            real(kind=4) :: scale_factor
            integer :: mixglo1x1_shape(4) !lon/lat/level/time
            integer :: ps_shape(3)        !lon/lat/time

            status = 0
            
            call readrc(rcf, 'start.' // trim(tracname) // '.filename', filename, status)
            print*, rname//"::reading initial condition from "//trim(filename)

            ! Read the relevant info from the netCDF file
            ncf = nc_open(trim(filename), 'r', status)

            !-- convert packed mixing ratio to physical value
            mixglo1x1_packed = nc_read_var(ncf, trim(names(itrac)))
            add_offset = nc_get_attr(ncf, trim(names(itrac)), 'add_offset', status)
            scale_factor = nc_get_attr(ncf, trim(names(itrac)), 'scale_factor', status)
            print*, rname//"::converted packed mixing ratio to physical values with "//&
                 "scale_factor=", scale_factor, " and add_offset=", add_offset, 'for tracer ' // trim(names(itrac))

            mixglo1x1_shape = shape(mixglo1x1_packed)
            allocate(mixglo1x1(mixglo1x1_shape(1), mixglo1x1_shape(2), mixglo1x1_shape(3), mixglo1x1_shape(4)))
            mixglo1x1 = real(mixglo1x1_packed) * scale_factor + add_offset
            
            time = nc_read_var(ncf, 'time')
            hyai = nc_read_var(ncf, 'hyai')
            hybi = nc_read_var(ncf, 'hybi')
            ps_packed = nc_read_var(ncf, 'ps')

            !- convert packed surface pressure to physical value
            add_offset = nc_get_attr(ncf, 'ps', 'add_offset', status)
            scale_factor = nc_get_attr(ncf, 'ps', 'scale_factor', status)
            ps_shape = shape(ps_packed)
            ps = real(ps_packed * scale_factor + add_offset)
            print*, rname//"::converted packed surface pressure to physical values with "//&
                 "scale_factor=",scale_factor," and add_offset=",add_offset
            call nc_close(ncf)

            ! Select the proper time index:
            call tau2date(itaur(1), idate)
            nhours_since_start_of_month = (idate(3) - 1) * 24 + idate(4)
            itime = minloc(time, 1, time - time(1) >= nhours_since_start_of_month)

            ! MVO-DEBUG
            print*, rname//'::mixglo1x1@itime=',itime,&
                 'min/mean/max=', &
                 minval(mixglo1x1(:,:,:,itime)), &
                 sum(mixglo1x1(:,:,:,itime))/size(mixglo1x1(:,:,:,itime)),&
                 maxval(mixglo1x1(:,:,:,itime))

                 ! Get the number of levels of the input data
            nlev = size(hyai) - 1

            ! Create a set of coordinates for the input field:
            ! hyai and hybi are real(4) variables, so add real(0, 8) for on-the-spot conversion
            call init_grid(hor_grid_in, -179.5, 1.0, 360, -89.5, 1.0, 180, status)
            IF_NOTOK_RETURN(status=1)
            call init_levels(ver_grid_in, nlev, hyai + real(0, 8), hybi + real(0, 8), status)
            IF_NOTOK_RETURN(status=1)

            ! Propagate the global 1x1 field to the regions
            do region = 1, nregions
                call regrid_3d_mix( &
                    hor_grid_in, lli(region), ver_grid_in, levi, &
                    ps(:, :, itime), mixglo1x1(:, :, :, itime), &
                    mass_dat(region)%rm_t(1:im(region), 1:jm(region), 1:lm(region), itrac), &
                    status, .true. &
                )
                mass_dat(region)%rm_t(:, :, :, itrac) = mass_dat(region)%rm_t(:, :, :, itrac) * m_dat(region)%data / fscale(itrac)
                print*, mass_dat(region)%rm_t(3, 3, 3, itrac)
            enddo


            status = 0

        end subroutine read_iniconc_cams

        subroutine set_constant_iniconc(itrac, tracname, status)

            integer, intent(in)     :: itrac
            integer, intent(out)    :: status
            character(len=*), intent(in)    :: tracname
            integer :: region
            real    :: mix
            character(len=*), parameter    :: rname = mname//'/sec_constant_iniconc'

            status = 0

            call readrc(rcf, 'start.' // trim(tracname) // '.mix_ratio', mix, status)
            IF_NOTOK_RETURN(status=1)

            do region = 1, nregions
                mass_dat(region)%rm_t(1:im(region), 1:jm(region), 1:lm(region), itrac) = mix
            enddo

        end subroutine set_constant_iniconc

end module iniconc_module
