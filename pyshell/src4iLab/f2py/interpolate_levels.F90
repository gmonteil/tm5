#define TRACEBACK write (0,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__
#define PRINT_ERROR(string) TRACEBACK; write(0,'(a)') trim(string)

module indexing

implicit none

public :: findIndex, findUnique

character(len=*), parameter :: mname = 'indexing'

contains

subroutine findUnique(len_source, src_lat, src_lon, tolerance, uniq_lat, uniq_lon, frequency, n_uniq)

    ! Given a sequence of lat/lon positions, find the unique tuples, which are
    ! unique up to some tolerance value, which makes it easier to plot all
    ! the locations
    integer, intent(in)     :: len_source
    real, intent(in)        :: src_lat(len_source), src_lon(len_source)
    real, intent(in)        :: tolerance
    integer, intent(out)    :: n_uniq
    real, intent(out)       :: uniq_lat(len_source), uniq_lon(len_source)
    integer, intent(out)   :: frequency(len_source)
    ! local variables
    integer                 :: i, j
    logical                 :: al_exists
    character(len=*), parameter :: rname = mname//'/findUnique'

    uniq_lat = 0.0
    uniq_lon = 0.0
    n_uniq = 0
    frequency = 0

    ! n_uniq always contains the number of unique locations stored so far
    do i = 1, len_source
        al_exists = .false.
        do j = 1, n_uniq
            if ((abs(src_lat(i)-uniq_lat(j)) .le. tolerance) .and. (abs(src_lon(i)-uniq_lon(j)) .le. tolerance)) then
                ! the i-th point in (src_lat, src_lon) corresponds to the j-th (uniq_lat, uniq_lon)
                frequency(j) = frequency(j) + 1
                al_exists = .true.
                exit
            end if
        end do
        if (.not. al_exists) then
            n_uniq = n_uniq + 1
            uniq_lat(n_uniq) = src_lat(i)
            uniq_lon(n_uniq) = src_lon(i)
            frequency(n_uniq) = frequency(n_uniq) + 1
        end if
    end do

    ! Have all points been accounted for?
!    write(*,'("Out of ", i8, " input points, ", i8, " have been accounted for")') len_source, sum(frequency)

end subroutine findUnique

subroutine findIndex(len_source, len_seek, source_array, seek_array, tolerance, indices_array)

    ! Given a source_array(len_source) and a seek_array(len_seek), return a list of
    ! indices, indices_array(len_seek), that contain the indices of source_array
    ! containing the elements of seek_array. The first occurrence is returned. If an
    ! element of seek_array is not present in source_array, the corresponding index
    ! is set to -1. Since we will use this to match ACOS IDs, which have fourteen
    ! digits, we need to use long integer arrays. We will match to within a tolerance,
    ! which can be zero.
    integer, intent(in)     :: len_source, len_seek, tolerance
    integer(8), intent(in)  :: source_array(len_source), seek_array(len_seek)
    integer, intent(out)    :: indices_array(len_seek)
    ! local variables
    integer                 :: i, j
    integer(8)              :: tol_long
    character(len=*), parameter :: rname = mname//'/findIndex'

    tol_long = abs(tolerance)

    ! Initially, set all elements of indices_array to -1
    indices_array(:) = -1

    do i = 1, len_seek
        do j = 1, len_source
            if ((seek_array(i) .ge. (source_array(j)-tol_long)) .and. (seek_array(i) .le. (source_array(j)+tol_long))) then
                indices_array(i) = j-1 ! switch from Fortran to Python indices
                exit
            end if
        end do
    end do

end subroutine findIndex

end module indexing

module redistrib_flux

implicit none

public :: extrapolate3D, regrid_fluxes, regrid_flux_series, regrid_flux_submul

character(len=*), parameter :: mname = 'redistrib_flux'

double precision, parameter :: R_e = 6371000.0 ! radius of earth in meters
double precision, parameter :: pi = 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214 ! http://oeis.org/A000796/constant
double precision, parameter :: deg_to_rad = pi/180.0
integer, parameter          :: ksp = selected_real_kind(2)
integer, parameter          :: kdp = selected_real_kind(15)

interface searchSorted
    module procedure D_searchSorted, R_searchSorted, I_searchSorted
end interface searchSorted

contains

subroutine extrapolate3D(num_yrs_exist, num_yrs_predict, time_dim, lat_dim, lon_dim, ip_array, yrs_predict, op_array)

    ! This is for exrapolating CO2 prior fluxes to years that don't have prior fluxes. The existing prior fluxes are fed in
    ! through ip_array or dimension num_yrs_exist x time_dim x lat_dim x lon_dim, and the op_array is num_yrs_predict x time_dim x
    ! lat_dim x lon_dim, containing the predicted fluxes. The years to be predicted are specified as follows. Suppose
    ! fluxes exist for 2000 through 2010, so num_yrs_exist = 11. 2000 is then year 1, 2001 is year 2, etc. So if yrs_predict
    ! is (13,14,16), then fluxes are required for years 2012, 2013 and 2015. This allows us to produce fluxes for multiple
    ! arbitrary years in the future, instead of sticking to an integer sequence of years immediately after existing years.
    integer, intent(in)             :: time_dim, lat_dim, lon_dim, num_yrs_exist, num_yrs_predict
    double precision, intent(in)    :: ip_array(num_yrs_exist,time_dim,lat_dim,lon_dim)
    integer, intent(in)             :: yrs_predict(num_yrs_predict)
    double precision, intent(out)   :: op_array(num_yrs_predict,time_dim,lat_dim,lon_dim)
    ! local variables
    integer                         :: i, j, t
    double precision                :: y(num_yrs_exist), x(num_yrs_exist) ! x for year index, y for flux
    double precision                :: x_ext(num_yrs_predict) ! year indices for the years to predict
    double precision                :: a0, a1, n, xbar, ybar
    character(len=*), parameter     :: rname = mname//'/extrapolate3D'

    do t = 1, num_yrs_exist
        x(t) = t
    end do
    n = num_yrs_exist

    do t = 1, num_yrs_predict
        x_ext(t) = yrs_predict(t)
    end do

    do t = 1, time_dim
        do j = 1, lat_dim
            do i = 1, lon_dim
                y = ip_array(:,t,j,i)
                ! Fit a straight line, y vs x
                xbar = sum(x)/n
                ybar = sum(y)/n
                a1 = (sum(x*y)/n - xbar*ybar)/(sum(x*x)/n - xbar**2)
                a0 = ybar - a1*xbar
                ! Extrapolate to x_ext
                op_array(:,t,j,i) = a0 + a1*x_ext
            end do
        end do
    end do

end subroutine extrapolate3D

subroutine regrid_fluxes_from_grid(jm_ip, im_ip, jm_op, im_op, ip_lats, ip_lons, op_lats, op_lons, ip_flux, per_area, op_flux)

    ! Given a TM5 flux array, in mass/gridbox/time, redistribute it over the output grid in
    ! units of mass/m^2/time or mass/gridbox/time. The grids are defined by ip_lats/ip_lons
    ! and op_lats/op_lons, which define the boundary lats/lons, not the centers.
    integer, intent(in)             :: jm_ip, im_ip, jm_op, im_op ! jm corresponds to lats, im corresponds to lons
    double precision, intent(in)    :: ip_lats(jm_ip+1), ip_lons(im_ip+1), op_lats(jm_op+1), op_lons(im_op+1)
    double precision, intent(in)    :: ip_flux(jm_ip,im_ip)
    logical, intent(in)             :: per_area ! is the output per area or per gridbox?
    double precision, intent(out)   :: op_flux(jm_op, im_op)
    ! local variables
    double precision                :: ip_surf_area(jm_ip,im_ip), dummy_flux(jm_ip,im_ip)
    double precision                :: dlon_ip(im_ip), dlat_ip(jm_ip)
    integer                         :: i,j
    character(len=*), parameter     :: rname = mname//'/regrid_fluxes_from_grid'

    ! easiest way to do this is to convert the input flux to per area, then pass it to regrid_fluxes
    ! calculate the surface area per gridbox for the input grid

    dlon_ip = deg_to_rad * (ip_lons(2:im_ip+1) - ip_lons(1:im_ip))
    dlat_ip = sin(deg_to_rad*ip_lats(2:jm_ip+1)) - sin(deg_to_rad*ip_lats(1:jm_ip))
    do j = 1, jm_ip
        do i = 1, im_ip
            ip_surf_area(j,i) = dlon_ip(i) * dlat_ip(j)
        end do
    end do
    ip_surf_area = R_e * R_e * ip_surf_area

    dummy_flux = ip_flux/ip_surf_area
    ! now redistribute
    call regrid_fluxes(jm_ip, im_ip, jm_op, im_op, ip_lats, ip_lons, op_lats, op_lons, dummy_flux, per_area, op_flux)

end subroutine regrid_fluxes_from_grid

subroutine regrid_fluxes(jm_ip, im_ip, jm_op, im_op, ip_lats, ip_lons, op_lats, op_lons, ip_flux, per_area, op_flux)

    ! Given a flux in a 2D array ip_flux, in units of mass/m^2/time, redistribute it over the output
    ! grid in units of mass/gridbox/time and mass/m^2/time. The grids are defined by ip_lats/ip_lons
    ! and op_lats/op_lons, which define the boundary lats/lons, not the centers.
    integer, intent(in)             :: jm_ip, im_ip, jm_op, im_op ! jm corresponds to lats, im corresponds to lons
    double precision, intent(in)    :: ip_lats(jm_ip+1), ip_lons(im_ip+1), op_lats(jm_op+1), op_lons(im_op+1)
    double precision, intent(in)    :: ip_flux(jm_ip,im_ip)
    logical, intent(in)             :: per_area
    double precision, intent(out)   :: op_flux(jm_op, im_op)
    ! local variables
    double precision                :: op_surf_area(jm_op,im_op)
    double precision                :: dlon_ip(im_ip), dlat_ip(jm_ip), dlon_op(im_op), dlat_op(jm_op)
    integer                         :: i,j,k,l
    double precision                :: ovlap_lat_up, ovlap_lat_dn, ovlap_lon_up, ovlap_lon_dn, ovlap_area
    character(len=*), parameter     :: rname = mname//'/regrid_fluxes'

    ! calculate the surface area per gridbox for the output grid

    if (per_area) then
        dlon_op = deg_to_rad * (op_lons(2:im_op+1) - op_lons(1:im_op))
        dlat_op = sin(deg_to_rad*op_lats(2:jm_op+1)) - sin(deg_to_rad*op_lats(1:jm_op))
        do j = 1, jm_op
            do i = 1, im_op
                op_surf_area(j,i) = dlon_op(i) * dlat_op(j)
            end do
        end do
        op_surf_area = R_e * R_e * op_surf_area
    end if

    ! time for redistribution
    op_flux = 0.0
    do i = 1, im_ip
        do j = 1, jm_ip
            ! the box with input flux is between latitudes ip_lats(j) and ip_lats(j+1), and between longitudes ip_lons(i) and ip_lons(i+1)
            do k = 1, im_op ! lon limits are op_lons(k) and op_lons(k+1)
                if ((op_lons(k) .ge. ip_lons(i+1)) .or. (op_lons(k+1) .le. ip_lons(i))) cycle
                do l = 1, jm_op ! lat limits are op_lats(l) and op_lats(l+1)
                    if ((op_lats(l) .ge. ip_lats(j+1)) .or. (op_lats(l+1) .le. ip_lats(j))) cycle
                    ! what area of the box ((ip_lons(i),ip_lons(i+1)), (ip_lats(j),ip_lats(j+1))) lies inside
                    ! the box ((op_lons(k),op_lons(k+1)), (op_lats(l),op_lats(l+1)))?
                    ovlap_lat_up = min(ip_lats(j+1),op_lats(l+1))
                    ovlap_lat_dn = max(ip_lats(j),op_lats(l))
                    ovlap_lon_up = min(ip_lons(i+1),op_lons(k+1))
                    ovlap_lon_dn = max(ip_lons(i),op_lons(k))
                    ovlap_area = R_e * R_e * deg_to_rad * (ovlap_lon_up-ovlap_lon_dn) * (sin(deg_to_rad*ovlap_lat_up) - sin(deg_to_rad*ovlap_lat_dn))
                    op_flux(l,k) = op_flux(l,k) + ip_flux(j,i) * ovlap_area
                end do
            end do
        end do
    end do

    ! now create a per-area output flux map as well
    if (per_area) op_flux = op_flux/op_surf_area

end subroutine regrid_fluxes

subroutine regrid_flux_submul(jm_ip, im_ip, jm_op, im_op, ip_lats, ip_lons, op_lats, op_lons, ip_flux, per_area, op_flux)

    ! Do the same as regrid_fluxes, but assume that the output grid is a submultiple of the input grid,
    ! and both grids are regular. For example, the input grid could be global 1x1, whereas the output
    ! grid could be regional 3x2. This speeds up the regridding considerably, but does not work for general
    ! cases such as nonuniform grids or where the input grid is (say) 1x1 and the output grid is (say) 3.25x1.75.
    integer, intent(in)             :: jm_ip, im_ip, jm_op, im_op ! jm corresponds to lats, im corresponds to lons
    double precision, intent(in)    :: ip_lats(jm_ip+1), ip_lons(im_ip+1), op_lats(jm_op+1), op_lons(im_op+1)
    double precision, intent(in)    :: ip_flux(jm_ip,im_ip)
    logical, intent(in)             :: per_area
    double precision, intent(out)   :: op_flux(jm_op, im_op)
    ! local variables
    double precision                :: op_surf_area(jm_op,im_op), ip_surf(jm_ip,im_ip)
    double precision                :: dlon_ip(im_ip), dlat_ip(jm_ip), dlon_op(im_op), dlat_op(jm_op)
    integer                         :: lat_mul, lon_mul, i_beg, j_beg, i, j
    character(len=*), parameter     :: rname = mname//'/regrid_flux_submul'

    ! calculate the surface area per gridbox for the output grid

    if (per_area) then
        dlon_op = deg_to_rad * (op_lons(2:im_op+1) - op_lons(1:im_op))
        dlat_op = sin(deg_to_rad*op_lats(2:jm_op+1)) - sin(deg_to_rad*op_lats(1:jm_op))
        do j = 1, jm_op
            do i = 1, im_op
                op_surf_area(j,i) = dlon_op(i) * dlat_op(j)
            end do
        end do
        op_surf_area = R_e * R_e * op_surf_area
    end if

    ! calculate the surface area of the input grid
    dlon_ip = deg_to_rad * (ip_lons(2:im_ip+1) - ip_lons(1:im_ip))
    dlat_ip = sin(deg_to_rad*ip_lats(2:jm_ip+1)) - sin(deg_to_rad*ip_lats(1:jm_ip))
    do j = 1, jm_ip
        do i = 1, im_ip
            ip_surf(j,i) = dlon_ip(i) * dlat_ip(j)
        end do
    end do
    ip_surf = R_e * R_e * ip_surf

    ! convert input flux to mass/gridcell/time, store it in the variable ip_surf since ip_flux is intent(in)
    ip_surf = ip_flux * ip_surf

    ! what submultiple is the o/p grid of the i/p grid?
    lat_mul = int(((op_lats(jm_op+1)-op_lats(1))/jm_op) / ((ip_lats(jm_ip+1)-ip_lats(1))/jm_ip))
    lon_mul = int(((op_lons(im_op+1)-op_lons(1))/im_op) / ((ip_lons(im_ip+1)-ip_lons(1))/im_ip))

    ! at what indices of the i/p array does the o/p array begin?
    i_beg = searchSorted(ip_lons, op_lons(1))
    j_beg = searchSorted(ip_lats, op_lats(1))
    op_flux = 0.0
    do i = 1, im_op
        do j = 1, jm_op
            op_flux(j,i) = op_flux(j,i) + sum(ip_surf(j_beg+(j-1)*lat_mul:j_beg+j*lat_mul-1, i_beg+(i-1)*lon_mul:i_beg+i*lon_mul-1))
        end do
    end do

    ! now create a per-area output flux map
    if (per_area) op_flux = op_flux/op_surf_area

end subroutine regrid_flux_submul

subroutine regrid_flux_series(jm_ip, im_ip, jm_op, im_op, timesteps, ip_lats, ip_lons, op_lats, op_lons, ip_flux, per_area, op_flux, flux_totals)

    ! Do the same thing as regrid_fluxes, but then for a series of 2D flux maps
    integer, intent(in)               :: jm_ip, im_ip, jm_op, im_op, timesteps ! jm corresponds to lats, im corresponds to lons
    logical, intent(in)               :: per_area ! True means that the flux output should be per area and not per gridbox
    double precision, intent(in)    :: ip_lats(jm_ip+1), ip_lons(im_ip+1), op_lats(jm_op+1), op_lons(im_op+1)
    double precision, intent(in)    :: ip_flux(timesteps,jm_ip,im_ip)
    double precision, intent(out)   :: op_flux(timesteps,jm_op, im_op), flux_totals(timesteps,2)
    ! local variables
    double precision                :: op_surf_area(jm_op,im_op), ip_surf_area(jm_ip,im_ip)
    double precision                :: dlon_ip(im_ip), dlat_ip(jm_ip), dlon_op(im_op), dlat_op(jm_op)
    integer                           :: i,j,k,l,t
    double precision                :: ovlap_lat_up, ovlap_lat_dn, ovlap_lon_up, ovlap_lon_dn, ovlap_area
    character(len=*), parameter     :: rname = mname//'/regrid_flux_series'

    if (per_area) then
        dlon_op = deg_to_rad * (op_lons(2:im_op+1) - op_lons(1:im_op))
        dlat_op = sin(deg_to_rad*op_lats(2:jm_op+1)) - sin(deg_to_rad*op_lats(1:jm_op))
        do j = 1, jm_op
            do i = 1, im_op
                op_surf_area(j,i) = dlon_op(i) * dlat_op(j)
            end do
        end do
        op_surf_area = R_e * R_e * op_surf_area
    end if

    dlon_ip = deg_to_rad * (ip_lons(2:im_ip+1) - ip_lons(1:im_ip))
    dlat_ip = sin(deg_to_rad*ip_lats(2:jm_ip+1)) - sin(deg_to_rad*ip_lats(1:jm_ip))
    do j = 1, jm_ip
        do i = 1, im_ip
            ip_surf_area(j,i) = dlon_ip(i) * dlat_ip(j)
        end do
    end do
    ip_surf_area = R_e * R_e * ip_surf_area

    ! time for redistribution
    op_flux = 0.0
    do i = 1, im_ip
        do j = 1, jm_ip
            ! the box with input flux is between latitudes ip_lats(j) and ip_lats(j+1), and between longitudes ip_lons(i) and ip_lons(i+1)
            do k = 1, im_op ! lon limits are op_lons(k) and op_lons(k+1)
                if ((op_lons(k) .ge. ip_lons(i+1)) .or. (op_lons(k+1) .le. ip_lons(i))) cycle
                do l = 1, jm_op ! lat limits are op_lats(l) and op_lats(l+1)
                    if ((op_lats(l) .ge. ip_lats(j+1)) .or. (op_lats(l+1) .le. ip_lats(j))) cycle
                    ! what area of the box ((ip_lons(i),ip_lons(i+1)), (ip_lats(j),ip_lats(j+1))) lies inside
                    ! the box ((op_lons(k),op_lons(k+1)), (op_lats(l),op_lats(l+1)))?
                    ovlap_lat_up = min(ip_lats(j+1),op_lats(l+1))
                    ovlap_lat_dn = max(ip_lats(j),op_lats(l))
                    ovlap_lon_up = min(ip_lons(i+1),op_lons(k+1))
                    ovlap_lon_dn = max(ip_lons(i),op_lons(k))
                    ovlap_area = R_e * R_e * deg_to_rad * (ovlap_lon_up-ovlap_lon_dn) * (sin(deg_to_rad*ovlap_lat_up) - sin(deg_to_rad*ovlap_lat_dn))
                    op_flux(:,l,k) = op_flux(:,l,k) + ip_flux(:,j,i) * ovlap_area
                end do
            end do
        end do
    end do

    ! op_choice = 0 corresponds to per gridbox fluxes, op_choice = 1 corresponds to per area fluxes
    if (per_area) then
        do t = 1, timesteps
            op_flux(t,:,:) = op_flux(t,:,:)/op_surf_area
        end do
    end if

    ! Debugging code to compare flux totals before and after regridding
    do t = 1, timesteps
        flux_totals(t,1) = sum(ip_flux(t,:,:) * ip_surf_area)
        if (per_area) then
            flux_totals(t,2) = sum(op_flux(t,:,:) * op_surf_area)
        else
            flux_totals(t,2) = sum(op_flux(t,:,:))
        end if
    end do
    ! End debugging code

end subroutine regrid_flux_series

function D_searchSorted(in_array, value) result(idx)
    ! __________________________________________________________
    !   D_searchSorted = return the index in array in_array such that
    !   in_array(index) <= value < in_array(index+1). If value is outside
    !   in_array then either lbound(in_array)-1 or ubound(in_array) is returned.
    ! __________________________________________________________
    ! __________________________________________________________
    real (kind=kdp), dimension(:), intent(in) :: in_array
    real (kind=kdp), intent(in)               :: value
    integer                                   :: idx
    ! __________________________________________________________
    integer :: ub, lb, i
    character(len=*), parameter     :: rname = mname//'/D_searchSorted'

    ub = ubound(in_array,1)
    lb = lbound(in_array,1)
    idx = lb-1
    do i=lb,ub
        if (value .lt. in_array(i)) exit
        idx = i
    end do
end function D_searchSorted

function R_searchSorted(in_array, value) result(idx)
    ! __________________________________________________________
    !   R_searchSorted = return the index in array in_array such that
    !   in_array(index) <= value < in_array(index+1). If value is outside
    !   in_array then either lbound(in_array)-1 or ubound(in_array) is returned.
    ! __________________________________________________________
    ! __________________________________________________________
    real (kind=ksp), dimension(:), intent(in) :: in_array
    real (kind=ksp), intent(in)               :: value
    integer                                   :: idx
    ! __________________________________________________________
    integer :: ub, lb, i
    character(len=*), parameter     :: rname = mname//'/R_searchSorted'

    ub = ubound(in_array,1)
    lb = lbound(in_array,1)
    idx = lb-1
    do i=lb,ub
        if (value .lt. in_array(i)) exit
        idx = i
    end do
end function R_searchSorted

function I_searchSorted(in_array, value) result(idx)
    ! __________________________________________________________
    !   I_searchSorted = return the index in array in_array such that
    !   in_array(index) <= value < in_array(index+1). If value is outside
    !   in_array then either lbound(in_array)-1 or ubound(in_array) is returned.
    ! __________________________________________________________
    ! __________________________________________________________
    integer, dimension(:), intent(in) :: in_array
    integer, intent(in)               :: value
    integer                           :: idx
    ! __________________________________________________________
    integer :: ub, lb, i
    character(len=*), parameter     :: rname = mname//'/I_searchSorted'

    ub = ubound(in_array,1)
    lb = lbound(in_array,1)
    idx = lb-1
    do i=lb,ub
        if (value .lt. in_array(i)) exit
        idx = i
    end do
end function I_searchSorted

end module redistrib_flux

module precon

use omp_lib

implicit none

public :: xc_to_x, g_to_gc
public :: mult_M, mult_MT

!interface matmul_fast
!    module procedure matmul_fast_matmat
!    module procedure matmul_fast_matvec
!end interface matmul_fast

character(len=*), parameter :: mname = 'precon'

contains

subroutine xc_to_x(n_state, nt, n_hor, G_state, Temp_L, Hor_L, x_c, ipos, num_threads, x)

    use omp_lib
    ! Converting x_c to x, especially for larger resolutions such as global 3x2
    ! or global 1x1, takes a lot of time in python. So here we do the matrix
    ! manipulations in fortran. Turns out, it's just as slow, so we're going
    ! to use OpenMP.
    integer, intent(in)             :: n_state, nt, n_hor, ipos, num_threads
    double precision, intent(in)    :: G_state(n_state), Temp_L(nt,nt), Hor_L(n_hor,n_hor), x_c(n_state)
    double precision, intent(out)   :: x(n_state)
    ! local variables
    integer                         :: i, j, k
    integer, dimension(nt*nt,2)     :: idx
    character(len=*), parameter     :: rname = mname//'/xc_to_x'

    x = 0.0
    k = 1
    do i = 1, nt
        do j = 1, nt
            idx(k,:) = (/i,j/)
            k = k + 1
        end do
    end do

    call omp_set_num_threads(num_threads)

    !$omp parallel private(i,j)
    !write(*,'("xc_to_x operating on ", i2, " threads")') omp_get_num_threads()
    !$omp do schedule(guided)
    do i = 1, nt
        do j = 1, nt
            x(ipos+(i-1)*n_hor+1:ipos+i*n_hor) = x(ipos+(i-1)*n_hor+1:ipos+i*n_hor) + G_state(ipos+(i-1)*n_hor+1:ipos+i*n_hor) * matmul(Temp_L(i,j)*Hor_L, x_c(ipos+(j-1)*n_hor+1:ipos+j*n_hor))
        end do
    end do
    !$omp end do
    !$omp end parallel

end subroutine xc_to_x

subroutine g_to_gc(n_state, nt, n_hor, G_state, Temp_Lt, Hor_Lt, g, ipos, num_threads, g_c)

    use omp_lib
    ! Convering g to g_c, especially for larger resolutions such as global 3x2
    ! or global 1x1, takes a lot of time in python. So here we do the matrix
    ! manipulations in fortran. Turns out, it's just as slow, so we're going
    ! to use OpenMP.
    integer, intent(in)             :: n_state, nt, n_hor, ipos, num_threads
    double precision, intent(in)    :: G_state(n_state), Temp_Lt(nt,nt), Hor_Lt(n_hor,n_hor), g(n_state)
    double precision, intent(out)   :: g_c(n_state)
    ! local variables
    integer                         :: i, j, k
    integer, dimension(nt*nt,2)     :: idx
    character(len=*), parameter     :: rname = mname//'/g_to_gc'

    g_c = 0.0
    k = 1
    do i = 1, nt
        do j = 1, nt
            idx(k,:) = (/i,j/)
            k = k + 1
        end do
    end do

    call omp_set_num_threads(num_threads)

    !$omp parallel private(i,j)
    !write(*,'("g_to_gc operating on ", i2, " threads")') omp_get_num_threads()
    !$omp do schedule(guided)
    do i = 1, nt
        do j = 1, nt
            g_c(ipos+(i-1)*n_hor+1:ipos+i*n_hor) = g_c(ipos+(i-1)*n_hor+1:ipos+i*n_hor) + matmul(Temp_Lt(i,j)*Hor_Lt, G_state(ipos+(j-1)*n_hor+1:ipos+j*n_hor) * g(ipos+(j-1)*n_hor+1:ipos+j*n_hor))
        end do
    end do
    !$omp end do
    !$omp end parallel

end subroutine g_to_gc

subroutine mult_MT(n_block, nt, n_hor, G_state, G_mask, Temp_L, Hor_L, v_in, num_threads, v_out)

    integer, intent(in)             :: n_block, nt, n_hor, num_threads
    double precision, intent(in)    :: G_state(n_block), Temp_L(nt,nt), Hor_L(n_hor,n_hor)
    logical, intent(in)             :: G_mask(n_block)
    double precision, intent(in), target    :: v_in(n_block)
    double precision, intent(out), target   :: v_out(n_block)
    ! local variables
    integer                         :: i, j, k, nvalid
    integer, allocatable            :: idx(:)
    logical, allocatable            :: sub_mask(:)
    double precision, pointer       :: vip(:), vop(:)
    character(len=*), parameter     :: rname = mname//'/mult_MT'

    v_out = 0.0
    call omp_set_num_threads(num_threads)

    !$omp parallel do private(i,j,nvalid,idx,k,sub_mask,vip,vop) schedule(static)
    do i = 1, nt
        allocate(sub_mask(n_hor))
        vop => v_out((i-1)*n_hor+1:i*n_hor)
        do j = 1, nt
            ! What are the valid G_state elements?
            sub_mask = G_mask((j-1)*n_hor+1:j*n_hor)
            nvalid = count(sub_mask)
            allocate(idx(nvalid))
            idx = pack([(k, k=1,n_hor)], sub_mask)
            vip => v_in((j-1)*n_hor+1:j*n_hor)

            vop = vop + matmul(Temp_L(i,j)*Hor_L(:,idx), 1.0/pack(G_state((j-1)*n_hor+1:j*n_hor), sub_mask) * vip(idx))

            deallocate(idx)
            nullify(vip)
        end do
        deallocate(sub_mask)
        nullify(vop)
    end do
    !$omp end parallel do

end subroutine mult_MT

subroutine mult_M(n_block, nt, n_hor, G_state, G_mask, Temp_L, Hor_L, v_in, num_threads, v_out)

    integer, intent(in)             :: n_block, nt, n_hor, num_threads
    double precision, intent(in)    :: G_state(n_block), Temp_L(nt,nt), Hor_L(n_hor,n_hor)
    logical, intent(in)             :: G_mask(n_block)
    double precision, intent(in), target    :: v_in(n_block)
    double precision, intent(out), target   :: v_out(n_block)
    ! local variables
    integer                         :: i, j, k, nvalid
    integer, allocatable            :: idx(:)
    logical, allocatable            :: sub_mask(:)
    double precision, pointer       :: vip(:), vop(:)
    character(len=*), parameter     :: rname = mname//'/mult_M'

    v_out = 0.0
    call omp_set_num_threads(num_threads)

    !$omp parallel do private(i,j,k,nvalid,idx,sub_mask,vip,vop) schedule(static)
    do i = 1, nt
        allocate(sub_mask(n_hor))
        ! What are the valid G_state elements?
        sub_mask = G_mask((i-1)*n_hor+1:i*n_hor)
        nvalid = count(sub_mask)
        allocate(idx(nvalid))
        idx = pack([(k, k=1,n_hor)], sub_mask)
        vop => v_out((i-1)*n_hor+1:i*n_hor)

        do j = 1, nt
            vip => v_in((j-1)*n_hor+1:j*n_hor)
            vop(idx) = vop(idx) + 1.0/pack(G_state((i-1)*n_hor+1:i*n_hor), sub_mask) * matmul(Temp_L(i,j)*Hor_L(idx,:), vip)
            nullify(vip)
        end do

        deallocate(idx,sub_mask)
        nullify(vop)
    end do
    !$omp end parallel do

end subroutine mult_M

end module precon

module sat_utils

implicit none

double precision, parameter   :: earth_radius = 6371.0 ! km, to be consistent with what's in TM5
double precision, parameter   :: pi = 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214 ! http://oeis.org/A000796/constant
double precision, parameter   :: deg_to_rad = pi/180.0

public :: errorInflation

character(len=*), parameter :: mname = 'sat_utils'

contains

subroutine errorInflation(num_records, times, in_errors, locations, corr_time, corr_length, bin_choice, err_floor, out_errors, num_samples)

    ! We need to inflate satellite measurement errors based on how many other observations are near a certain observation
    ! This will be used in lieu of using the full observation covariance matrix
    ! num_records -- total number of samples
    ! times       -- array of times, in ascending order
    ! in_errors   -- array of standard deviations
    ! locations   -- (lat, lon) for each observation
    ! corr_length -- floating point number, radius of circle inside which to average
    ! corr_time   -- floating point number, half the window width
    ! bin_choice  -- 'b' for the algorithm of Basu et al (2013), 'k' for Susan Kulawik's estimate of random and systematic error
    ! err_floor   -- estimate of systematic error, 0.7 ppm according to Kulawik, not applicable for Basu
    ! out_errors  -- array of inflated errors
    ! num_samples -- number of samples from the original time array used to inflate each error
    integer, intent(in)           :: num_records
    double precision, intent(in)  :: times(num_records), in_errors(num_records), locations(num_records, 2)
    double precision, intent(in)  :: corr_time, corr_length, err_floor
    character(len=1), intent(in)  :: bin_choice
    double precision, intent(out) :: out_errors(num_records)
    integer, intent(out)          :: num_samples(num_records)
    ! local variables
    logical                       :: input_ok
    integer                       :: i, j, w_beg, w_end
    double precision              :: distance, sigma_sum, var_sum, lat1, lat2, lon1, lon2, random_error
    character(len=*), parameter   :: rname = mname//'/errorInflation'

    input_ok = .true.
    ! first check whether the time array is sorted or not
    do i = 1, num_records-1
        if (times(i) > times(i+1)) input_ok = .false.
    end do
    if ( .not. input_OK ) then
        write(0,*) 'Times are not strictly increasing'
        stop
    end if
    ! then check if all the errors are positive or not
    do i = 1, num_records
        if (in_errors(i) .lt. 0.0) input_ok = .false.
    end do
    if ( .not. input_OK ) then
        write(0,*) 'Some errors are negative'
        stop
    end if
    ! now check if the time_window is positive or not
    if (corr_length .le. 0.0) input_ok = .false.
    if (corr_time .le. 0.0) input_ok = .false.
    if ( .not. input_OK ) then
        write(0,*) 'Correlation length and time are not positive'
        stop
    end if

    !$omp parallel do schedule(guided) &
    !$omp private(i, w_beg, w_end, sigma_sum, var_sum, lat1, lon1, j, distance, lat2, lon2, random_error)
    do i = 1, num_records
        w_beg = 1
        do while (times(w_beg) .lt. times(i) - corr_time)
            w_beg = w_beg + 1
        end do
        w_end = num_records
        do while (times(w_end) .ge. times(i) + corr_time)
            w_end = w_end - 1
        end do

        ! select by spatial distance
        sigma_sum = 0.0
        var_sum = 0.0
        lat1 = locations(i,1) * deg_to_rad
        lon1 = locations(i,2) * deg_to_rad
        num_samples(i) = 0

        select case (bin_choice)

        case ('B', 'b')

            do j = w_beg, w_end
                if (i == j) then
                    distance = 0.0
                else
                    lat2 = locations(j,1) * deg_to_rad
                    lon2 = locations(j,2) * deg_to_rad
                    distance = earth_radius * acos(cos(lat1) * cos(lat2) * cos(lon1-lon2) + sin(lat1)*sin(lat2))
                end if
                if (distance .le. corr_length) then
                    ! j goes over all possible nearest neighbors of i, and this block is only executed for the neighbors that count
                    sigma_sum = sigma_sum + in_errors(j)
                    var_sum = var_sum + in_errors(j)**2
                    num_samples(i) = num_samples(i) + 1
                end if
            end do
            out_errors(i) = in_errors(i) * sqrt(sigma_sum**2/var_sum)

        case ('k', 'K')

            ! only need the number of nearest neighbors, their errors are not important
            do j = w_beg, w_end
                if (i == j) then
                    distance = 0.0
                else
                    lat2 = locations(j,1) * deg_to_rad
                    lon2 = locations(j,2) * deg_to_rad
                    distance = earth_radius * acos(cos(lat1) * cos(lat2) * cos(lon1-lon2) + sin(lat1)*sin(lat2))
                end if
                if (distance .le. corr_length) then
                    ! j goes over all possible nearest neighbors of i, and this block is only executed for the neighbors that count
                    num_samples(i) = num_samples(i) + 1
                end if
            end do
            random_error = max(in_errors(i)**2 - err_floor**2, 0.0)
            out_errors(i) = sqrt(random_error + err_floor**2 * num_samples(i))

        end select

    end do
    !$omp end parallel do

end subroutine errorInflation

end module sat_utils

module interpolate_fields

implicit none

public :: rebin, rebin_2D, rebin_3D, dailyAverage, regridDailyField, inv_rebin_arr

character(len=*), parameter :: mname = 'interpolate_fields'

contains

subroutine movingAverage(num_records, times, in_errors, locations, corr_time, corr_length, out_errors, num_samples)

    ! We need to inflate satellite measurement errors based on how many other observations are near a certain observation
    ! This will be used in lieu of using the full observation covariance matrix
    ! num_records -- total number of samples
    ! times       -- array of times, in ascending order
    ! in_errors   -- array of standard deviations
    ! locations   -- (lat, lon) for each observation
    ! corr_length -- floating point number, radius of circle inside which to average
    ! corr_time   -- floating point number, half the window width
    ! out_errors  -- array of inflated errors
    ! num_samples -- number of samples from the original time array used to inflate each error
    integer, intent(in)           :: num_records
    double precision, intent(in)  :: times(num_records), in_errors(num_records), locations(num_records, 2)
    double precision, intent(in)  :: corr_time, corr_length
    double precision, intent(out) :: out_errors(num_records)
    integer, intent(out)          :: num_samples(num_records)
    ! local variables
    logical                       :: input_ok
    integer                       :: i, j, w_beg, w_end
    double precision              :: distance, sigma_sum, var_sum, lat1, lat2, lon1, lon2
    double precision, parameter   :: earth_radius = 6371.0 ! km, to be consistent with what's in TM5
    double precision, parameter   :: pi = 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214 ! http://oeis.org/A000796/constant
    double precision, parameter   :: deg_to_rad = pi/180.0
    character(len=*), parameter   :: rname = mname//'/movingAverage'

    input_ok = .true.
    ! first check whether the time array is sorted or not
    do i = 1, num_records-1
        if (times(i) > times(i+1)) input_ok = .false.
    end do
    if ( .not. input_OK ) then
        write(0,*) 'Times are not strictly increasing'
        stop
    end if
    ! then check if all the errors are positive or not
    do i = 1, num_records
        if (in_errors(i) .lt. 0.0) input_ok = .false.
    end do
    if ( .not. input_OK ) then
        write(0,*) 'Some errors are negative'
        stop
    end if
    ! now check if the time_window is positive or not
    if (corr_length .le. 0.0) input_ok = .false.
    if (corr_time .le. 0.0) input_ok = .false.
    if ( .not. input_OK ) then
        write(0,*) 'Correlation length and time are not positive'
        stop
    end if

    !$omp parallel do schedule(guided) &
    !$omp private(i, w_beg, w_end, sigma_sum, var_sum, lat1, lon1, j, distance, lat2, lon2)
    do i = 1, num_records
        w_beg = 1
        do while (times(w_beg) .lt. times(i) - corr_time)
            w_beg = w_beg + 1
        end do
        w_end = num_records
        do while (times(w_end) .ge. times(i) + corr_time)
            w_end = w_end - 1
        end do

        ! select by spatial distance
        sigma_sum = 0.0
        var_sum = 0.0
        lat1 = locations(i,1) * deg_to_rad
        lon1 = locations(i,2) * deg_to_rad
        num_samples(i) = 0
        do j = w_beg, w_end
            if (i == j) then
                distance = 0.0
            else
                lat2 = locations(j,1) * deg_to_rad
                lon2 = locations(j,2) * deg_to_rad
                distance = earth_radius * acos(cos(lat1) * cos(lat2) * cos(lon1-lon2) + sin(lat1)*sin(lat2))
            end if
            if (distance .le. corr_length) then
                sigma_sum = sigma_sum + in_errors(j)
                var_sum = var_sum + in_errors(j)**2
                num_samples(i) = num_samples(i) + 1
            end if
        end do
        out_errors(i) = in_errors(i) * sqrt(sigma_sum**2/var_sum)
    end do
    !$omp end parallel do

end subroutine movingAverage

subroutine dailyAverage(tm, lm_in, lats, lons, isobars, in_mixing, out_mixing, out_ps)

    ! I'm moving the code for daily averaging to Fortran
    ! So far only for CarbonTracker output, since that's the largest
    integer, parameter            :: archive_levels = 26
    integer, intent(in)           :: tm, lm_in, lats, lons
    double precision, intent(in)  :: isobars(tm,lm_in+1,lats,lons), in_mixing(tm,lm_in,lats,lons)
    double precision, intent(out) :: out_mixing(archive_levels-1,lats,lons), out_ps(archive_levels,lats,lons)
    integer                       :: i,j,k
    double precision              :: surface_pressure(lats,lons), pressure_column(archive_levels)
    double precision              :: avg_isobars(lm_in+1,lats,lons), avg_mixing(lm_in,lats,lons)
    double precision, parameter   :: at_coeffs(archive_levels) = (/ 0.000000, 95.636963, 298.495789, 713.218079, 1680.640259, &
                                  3960.291504, 6018.019531, 8765.053711, 12077.446289, 15379.805664, 18045.183594, 19755.109375, &
                                  20429.863281, 20097.402344, 18864.750000, 16899.468750, 14411.124023, 11632.758789, &
                                  8802.356445, 6144.314941, 3850.913330, 2063.779785, 855.361755, 210.393890, 7.367743, &
                                  0.000000 /) ! from top to surface
    double precision, parameter   :: bt_coeffs(archive_levels) = (/ 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, &
                                  0.000000, 0.000076, 0.001815, 0.011143, 0.034121, 0.073534, 0.130023, 0.202476, 0.288323, &
                                  0.383892, 0.484772, 0.586168, 0.683269, 0.771597, 0.847375, 0.907884, 0.951822, 0.979663, &
                                  0.994019, 1.000000 /) ! from top to surface
    character(len=*), parameter   :: rname = mname//'/dailyAverage'
    ! we assume that the arrays passed in have sane indexing, i.e., they start from 1
    ! From python, we need to pass (isobars, in_mixing)

    ! Keep in mind that the python code sorts levels from the ground up
    ! whereas rebin() needs them from the top down
    ! do the daily average first
    avg_isobars = sum(isobars,1)/tm
    avg_mixing = sum(in_mixing,1)/tm
    if (lm_in+1 == archive_levels) then
        out_ps = avg_isobars
        out_mixing = avg_mixing
    else
        ! get the surface pressure
        surface_pressure = avg_isobars(1,:,:)
        do i=1,lats
            do j=1,lons
                pressure_column = at_coeffs + bt_coeffs * surface_pressure(i,j)
                out_ps(:,i,j) = pressure_column(archive_levels:1:-1)
                call rebin(lm_in, archive_levels-2, avg_isobars(lm_in+1:1:-1,i,j), avg_mixing(lm_in:1:-1,i,j), pressure_column(2:archive_levels-1), out_mixing(:,i,j))
            end do
        end do
        ! now invert out_mixing to be from surface to top
        out_mixing = out_mixing(archive_levels-1:1:-1,:,:)
    end if

end subroutine dailyAverage

subroutine columnAverage(tm, lm, nlevs, ndays, mixing_ratio, isobars, averaging_kernel, averaging_kernel_pres, prior_mixing, prior_pres, doy_array, total_column)

    integer, intent(in)           :: tm, lm, nlevs, ndays
    double precision, intent(in)  :: mixing_ratio(tm,nlevs), isobars(tm,nlevs+1)
    double precision, intent(in)  :: averaging_kernel(tm,lm), averaging_kernel_pres(lm+1)
    double precision, intent(in)  :: prior_mixing(ndays,lm), prior_pres(ndays,lm+1)
    integer, intent(in)           :: doy_array(tm)
    double precision, intent(out) :: total_column(tm)
    integer                       :: day, obs
    double precision              :: mixing_for_ak(lm), rebinned_prior(ndays,lm), ak(lm), i_ak(lm)
    double precision, allocatable :: faux_isobars(:)
    double precision              :: slope, intercept
    character(len=*), parameter   :: rname = mname//'/columnAverage'

    ! Say there are 67123 observations. Then for Lamont/CarbonTracker, mixing_ratio(67123,25), isobars(67123,26), averaging_kernel(67123,70),
    ! averaging_kernel_pres(71), prior_mixing(366,70), prior_pres(366,71), doy_array(67123) and total_column(67123). The pressure levels
    ! averaging_kernel_pres and prior_pres are 71 layers thick because they represent pressures at the boundaries. For the same reason,
    ! isobars is 26 layers thick.

    ! First, rebin the prior profiles onto the averaging kernel pressure levels
    allocate(faux_isobars(lm+1))
    do day=1,ndays
        ! align the ends of the prior pressure levels with the averaging kernel pressure levels
        slope = (averaging_kernel_pres(1)-averaging_kernel_pres(lm+1))/(prior_pres(day,1)-prior_pres(day,lm+1))
        intercept = averaging_kernel_pres(1) - slope * prior_pres(day,1)
        faux_isobars = intercept + slope * prior_pres(day,:)
        call rebin(lm+1, lm-1, prior_pres(day,lm+1:1:-1), prior_mixing(day,lm:1:-1), faux_isobars(lm:2:-1), rebinned_prior(day,:))
        rebinned_prior(day,:) = rebinned_prior(day,lm:1:-1)
    end do
    deallocate(faux_isobars)

    allocate(faux_isobars(nlevs+1))
    do obs=1,tm
        day = doy_array(obs)
        ! if the ends of the mixing ratio isobars are not the same as the averaging kernel isobars, make them equal
        slope = (averaging_kernel_pres(1)-averaging_kernel_pres(lm+1))/(isobars(obs,1)-isobars(obs,nlevs+1))
        intercept = averaging_kernel_pres(1) - slope * isobars(obs,1)
        faux_isobars = intercept + slope * isobars(obs,:)
        call rebin(nlevs+1, lm-1, faux_isobars(nlevs+1:1:-1), mixing_ratio(obs,nlevs:1:-1), averaging_kernel_pres(lm:2:-1), mixing_for_ak)
        ! invert mixing_for_ak so that 1 corresponds to surface
        mixing_for_ak = mixing_for_ak(lm:1:-1)
        ! now compute the total column mixing ratio
        ak = - averaging_kernel(obs,:) * diff(averaging_kernel_pres)
        i_ak = - diff(averaging_kernel_pres) - ak
        total_column(obs) = (dot_product(ak, mixing_for_ak) + dot_product(i_ak, rebinned_prior(day,:)))/(averaging_kernel_pres(1)-averaging_kernel_pres(lm+1))
    end do
    deallocate(faux_isobars)

end subroutine columnAverage

subroutine takeColumn(lm, nlevs, mixing_ratio, isobars, averaging_kernel, averaging_kernel_pres, prior_mixing, prior_pres, total_column)

    integer, intent(in)           :: lm, nlevs
    double precision, intent(in)  :: mixing_ratio(nlevs), isobars(nlevs+1) ! nlevs is the number of model levels, e.g., 25/34 for CarbonTracker, 25/60 for TM5
    double precision, intent(in)  :: averaging_kernel(lm), averaging_kernel_pres(lm+1) ! lm is the number of levels over which the averaging kernel is defined, e.g., 70 for TCCON
    double precision, intent(in)  :: prior_mixing(lm), prior_pres(lm+1)
    double precision, intent(out) :: total_column
    integer                       :: day, obs
    double precision              :: mixing_for_ak(lm), rebinned_prior(lm), ak(lm), i_ak(lm)
    double precision, allocatable :: faux_isobars(:)
    double precision              :: slope, intercept, p_tot, prior_column
    character(len=*), parameter   :: rname = mname//'/takeColumn'

    ! There is a single column observation. Then for Lamont/CarbonTracker, mixing_ratio(25), isobars(26), averaging_kernel(70),
    ! averaging_kernel_pres(71), prior_mixing(70), prior_pres(71) and total_column is a single real scalar. The pressure levels
    ! averaging_kernel_pres and prior_pres are 71 layers thick because they represent pressures at the boundaries. For the same reason,
    ! isobars is 26 layers thick.

    ! First, rebin the prior profiles onto the averaging kernel pressure levels
    allocate(faux_isobars(lm+1))
    ! align the ends of the prior pressure levels with the averaging kernel pressure levels
    slope = (averaging_kernel_pres(1)-averaging_kernel_pres(lm+1))/(prior_pres(1)-prior_pres(lm+1))
    intercept = averaging_kernel_pres(1) - slope * prior_pres(1)
    faux_isobars = intercept + slope * prior_pres
    call rebin(lm+1, lm-1, prior_pres(lm+1:1:-1), prior_mixing(lm:1:-1), faux_isobars(lm:2:-1), rebinned_prior(:))
    ! rebinned_prior goes from top to bottom now, so invert that to go from the ground to the top of the atmosphere
    rebinned_prior = rebinned_prior(lm:1:-1)
    deallocate(faux_isobars)

    allocate(faux_isobars(nlevs+1))
    ! if the ends of the mixing ratio isobars are not the same as the averaging kernel isobars, make them equal
    slope = (averaging_kernel_pres(1)-averaging_kernel_pres(lm+1))/(isobars(1)-isobars(nlevs+1))
    intercept = averaging_kernel_pres(1) - slope * isobars(1)
    faux_isobars = intercept + slope * isobars
    call rebin(nlevs+1, lm-1, faux_isobars(nlevs+1:1:-1), mixing_ratio(nlevs:1:-1), averaging_kernel_pres(lm:2:-1), mixing_for_ak)
    ! invert mixing_for_ak so that 1 corresponds to surface
    mixing_for_ak = mixing_for_ak(lm:1:-1)
    deallocate(faux_isobars)

    ! now compute the total column mixing ratio:
    !   p_tot = averaging_kernel_pres(1)-averaging_kernel_pres(lm+1)
    !   prior_column = dot_product(inv_diff(averaging_kernel_pres), rebinned_prior)/p_tot
    !   total_column = prior_column + dot_product(averaging_kernel, mixing_for_ak - rebinned_prior)

    ! Commented out old code
!    ak = - averaging_kernel * diff(averaging_kernel_pres)
!    i_ak = - diff(averaging_kernel_pres) - ak
!    total_column = (dot_product(ak, mixing_for_ak) + dot_product(i_ak, rebinned_prior))/(averaging_kernel_pres(1)-averaging_kernel_pres(lm+1))

    ! New code
    p_tot = averaging_kernel_pres(1)-averaging_kernel_pres(lm+1)
    prior_column = dot_product(inv_diff(averaging_kernel_pres), rebinned_prior)/p_tot
    total_column = prior_column + dot_product(averaging_kernel, mixing_for_ak - rebinned_prior)

end subroutine takeColumn

subroutine regridDailyField(tm, lm_in, lats, lons, isobars, in_mixing, out_mixing, out_ps)

    ! re-grids 34-level (or however many level) 3D field to 25 vertical levels
    integer, parameter            :: archive_levels = 26
    integer, intent(in)           :: tm, lm_in, lats, lons
    double precision, intent(in)  :: isobars(tm,lm_in+1,lats,lons), in_mixing(tm,lm_in,lats,lons)
    double precision, intent(out) :: out_mixing(tm,archive_levels-1,lats,lons), out_ps(tm,archive_levels,lats,lons)
    integer                       :: i,j,k
    double precision, parameter   :: at_coeffs(archive_levels) = (/ 0.000000, 95.636963, 298.495789, 713.218079, 1680.640259, &
                                  3960.291504, 6018.019531, 8765.053711, 12077.446289, 15379.805664, 18045.183594, 19755.109375, &
                                  20429.863281, 20097.402344, 18864.750000, 16899.468750, 14411.124023, 11632.758789, &
                                  8802.356445, 6144.314941, 3850.913330, 2063.779785, 855.361755, 210.393890, 7.367743, &
                                  0.000000 /) ! from top to surface
    double precision, parameter   :: bt_coeffs(archive_levels) = (/ 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, &
                                  0.000000, 0.000076, 0.001815, 0.011143, 0.034121, 0.073534, 0.130023, 0.202476, 0.288323, &
                                  0.383892, 0.484772, 0.586168, 0.683269, 0.771597, 0.847375, 0.907884, 0.951822, 0.979663, &
                                  0.994019, 1.000000 /) ! from top to surface
    double precision              :: surface_pressure(tm,lats,lons), pressure_column(archive_levels)
    character(len=*), parameter   :: rname = mname//'/regridDailyField'

    ! we assume that the arrays passed in have sane indexing, i.e., they start from 1
    ! From python, we need to pass (isobars, in_mixing)
    ! Keep in mind that the python code sorts levels from the ground up
    ! whereas rebin() needs them from the top down
    if (lm_in+1 == archive_levels) then
        out_ps = isobars
        out_mixing = in_mixing
    else
        ! get the surface pressure
        surface_pressure = isobars(:,1,:,:)
        do i=1,tm
            do j=1,lats
                do k=1,lons
                    pressure_column = at_coeffs + bt_coeffs * surface_pressure(i,j,k)
                    out_ps(i,:,j,k) = pressure_column(archive_levels:1:-1)
                    call rebin(lm_in, archive_levels-2, isobars(i,lm_in+1:1:-1,j,k), in_mixing(i,lm_in:1:-1,j,k), pressure_column(2:archive_levels-1), out_mixing(i,:,j,k))
                end do
            end do
        end do
        ! now invert out_mixing to be from surface to top
        out_mixing = out_mixing(:,archive_levels-1:1:-1,:,:)
    end if

end subroutine regridDailyField

subroutine rebin_3D(lm_in, lats, lons, lm_out, in_levels, in_mixing, out_levels, out_mixing)

    ! in_mixing is a 3D field, with the first axis refering to levels
    integer, intent(in)           :: lm_in, lats, lons, lm_out
    double precision, intent(in)  :: in_levels(lm_in,lats,lons), in_mixing(lm_in-1,lats,lons), out_levels(lm_out,lats,lons)
    double precision, intent(out) :: out_mixing(lm_out+1,lats,lons)
    integer                       :: i, j
    character(len=*), parameter   :: rname = mname//'/rebin_3D'
    ! we assume that the arrays passed in have sane indexing, i.e., they start from 1

    do i=1,lats
        do j=1,lons
            call rebin(lm_in, lm_out, in_levels(:,i,j), in_mixing(:,i,j), out_levels(:,i,j), out_mixing(:,i,j))
        end do
    end do

end subroutine rebin_3D

subroutine rebin_2D(lm_in, time_in, lm_out, in_levels, in_mixing, out_levels, out_mixing)

    ! in_mixing is a 2D field with the second axis refering to levels
    integer, intent(in)           :: lm_in, time_in, lm_out
    double precision, intent(in)  :: in_levels(time_in,lm_in), in_mixing(time_in,lm_in-1), out_levels(time_in,lm_out)
    double precision, intent(out) :: out_mixing(time_in,lm_out+1)
    integer                       :: i
    character(len=*), parameter   :: rname = mname//'/rebin_2D'
    ! we assume that the arrays passed in have sane indexing, i.e., they start from 1

    do i=1,time_in
        call rebin(lm_in, lm_out, in_levels(i,:), in_mixing(i,:), out_levels(i,:), out_mixing(i,:))
    end do

end subroutine rebin_2D

subroutine rebin(in_levels, out_levels, in_x_grid, in_y_vals, out_x_grid, out_y_vals)

    ! out_x_grid is the output grid except the end-points, since
    ! the end-points of out_x_grid are the same as those of in_x_grid
    ! in_x_grid and out_x_grid are in ascending order
    integer, intent(in)         :: in_levels, out_levels
    double precision, intent(in)  :: in_x_grid(in_levels), in_y_vals(in_levels-1), out_x_grid(out_levels)
    double precision, intent(out) :: out_y_vals(out_levels+1)
    logical                     :: input_OK
    integer                     :: i, breakpoints(out_levels)
    double precision              :: cumul(out_levels+2)
    double precision              :: in_x(in_levels), out_x(out_levels+2)
    character(len=*), parameter   :: rname = mname//'/rebin'

    ! first some error checking
    input_OK = .true.
    ! check whether x_grid's are sorted and unique
    do i=lbound(in_x_grid,1),ubound(in_x_grid,1)-1
        if (in_x_grid(i) >= in_x_grid(i+1)) then
            input_OK = .false.
            write(0,*) 'Input grid is not in ascending order'
        end if
    end do
    do i=lbound(out_x_grid,1),ubound(out_x_grid,1)-1
        if (out_x_grid(i) >= out_x_grid(i+1)) then
            input_OK = .false.
            write(0,*) 'Output grid is not in ascending order'
        end if
    end do
    ! check whether the output grid is contained inside the input grid
    if (out_x_grid(lbound(out_x_grid,1)) < in_x_grid(lbound(in_x_grid,1))) then
        input_OK = .false.
        write(0,*) 'Output grid starts earlier than input grid'
    end if
    if (out_x_grid(ubound(out_x_grid,1)) > in_x_grid(ubound(in_x_grid,1))) then
        input_OK = .false.
        write(0,*) 'Output grid ends later than input grid'
    end if
    if ( .not. input_OK ) then
        write(0,*) 'There is a problem with the input arrays'
        stop
    end if
    ! now that we know the input is OK, proceed with rebinning :-)

    ! copy over data
    in_x(:) = in_x_grid(:)
    out_x(2:size(out_x)-1) = out_x_grid(:)
    out_x(1) = in_x(1)
    out_x(size(out_x)) = in_x(size(in_x))

    breakpoints = searchSorted(in_x, out_x_grid)
    cumul(1) = 0.0
    do i=1,size(breakpoints)
        cumul(i+1) = sum(in_y_vals(:breakpoints(i)-1) * diff(in_x(:breakpoints(i))))
        cumul(i+1) = cumul(i+1) + in_y_vals(breakpoints(i)) * (out_x(i+1) - in_x(breakpoints(i)))
    end do
    cumul(size(cumul)) = sum(in_y_vals * diff(in_x))
    out_y_vals = diff(cumul)/diff(out_x)

end subroutine rebin

subroutine inv_rebin_arr(num_obs, in_levels, out_levels, in_x_grid, in_y_vals, out_x_grid, out_y_vals, distrib_matrix)

    integer, intent(in)             :: in_levels, out_levels, num_obs
    double precision, intent(in)    :: in_x_grid(num_obs, in_levels), in_y_vals(num_obs, in_levels-1), out_x_grid(num_obs, out_levels)
    double precision, intent(out)   :: out_y_vals(num_obs, out_levels+1), distrib_matrix(num_obs, out_levels+1, in_levels-1)

    logical                         :: input_OK
    integer                         :: i_obs, i, j
    double precision, allocatable   :: in_x(:,:), out_x(:,:), cumul(:)
    integer, allocatable            :: breakpoints(:), ext_breakpoints(:)
    double precision                :: mass_error, dummy(out_levels+1), lb, ub
    character(len=1024)             :: error_msg
    character(len=*), parameter     :: rname = mname//'/inv_rebin_arr'

    ! first some error checking
    input_OK = .true.
    ! check whether x_grid's are sorted in descending order and unique
    do i = lbound(in_x_grid, 2), ubound(in_x_grid, 2)-1
        if (any(in_x_grid(:, i) <= in_x_grid(:, i+1))) then
            input_OK = .false.
            error_msg = 'Input grid is not in descending order'
        end if
    end do
    do i = lbound(out_x_grid, 2), ubound(out_x_grid, 2)-1
        if (any(out_x_grid(:, i) <= out_x_grid(:, i+1))) then
            input_OK = .false.
            error_msg = 'Output grid is not in descending order'
        end if
    end do

    ! check whether the output grid is contained inside the input grid
    if (any(out_x_grid(:, lbound(out_x_grid,2)) > in_x_grid(:, lbound(in_x_grid,2)))) then
        input_OK = .false.
        error_msg = 'First value of output grid is higher than first value of input grid'
    end if
    if (any(out_x_grid(:, ubound(out_x_grid,2)) < in_x_grid(:, ubound(in_x_grid, 2)))) then
        input_OK = .false.
        error_msg = 'Last value of output grid is less than last value of input grid'
    end if
    if ( .not. input_OK ) then
        PRINT_ERROR(error_msg)
        stop
    end if
    ! now that we know the input is OK, proceed with rebinning :-)

    ! copy over data
    allocate(in_x(num_obs, in_levels), out_x(num_obs, out_levels+2))
    in_x = in_x_grid
    out_x(:, 2:out_levels+1) = out_x_grid(:, :)
    out_x(:, 1) = in_x(:, 1)
    out_x(:, out_levels+2) = in_x(:, in_levels)

    distrib_matrix = 0.0
    !$omp parallel private(i_obs, breakpoints, ext_breakpoints, cumul, i, j, ub, lb)
    allocate(breakpoints(out_levels), ext_breakpoints(out_levels+2), cumul(out_levels+2))
    !$omp do schedule(guided)
    do i_obs = 1, num_obs
        breakpoints = inv_searchSorted(in_x(i_obs, :), out_x_grid(i_obs, :))
        ext_breakpoints(2:out_levels+1) = breakpoints
        ext_breakpoints(1) = 1
        ext_breakpoints(out_levels+2) = in_levels-1

        cumul(1) = 0.0
        do i=1,size(breakpoints)
            cumul(i+1) = sum(in_y_vals(i_obs, :breakpoints(i)-1) * inv_diff(in_x(i_obs, :breakpoints(i))))
            cumul(i+1) = cumul(i+1) + in_y_vals(i_obs, breakpoints(i)) * (in_x(i_obs, breakpoints(i)) - out_x(i_obs, i+1))
        end do ! i
        cumul(size(cumul)) = sum(in_y_vals(i_obs,:) * inv_diff(in_x(i_obs,:)))
        out_y_vals(i_obs,:) = diff(cumul)/inv_diff(out_x(i_obs,:))

        ! now calculate the distribution matrix
        ! out_y_vals(i) has contributions from in_y_vals(ext_breakpoints(i)) to in_y_vals(ext_breakpoints(i+1))
        ! in_y_vals(j) lies between in_x(j) and in_x(j+1)
        ! hence, out_y_vals(i), which lines between out_x(i) and out_x(i+1), will need to examine levels in_x(ext_breakpoints(i)) and in_x(ext_breakpoints(i+1)+1)
        do i=1,out_levels+1
            do j=1,in_levels-1
                ! output level i is between out_x(i) and out_x(i+1)
                ! input level j is between in_x(j) and in_x(j+1)
                ! do something only if there is some overlap
                if ((in_x(i_obs,j+1) .ge. out_x(i_obs,i)) .or. (out_x(i_obs,i+1) .ge. in_x(i_obs,j))) then
                    continue
                else
                    ub = min(out_x(i_obs,i), in_x(i_obs,j))
                    lb = max(in_x(i_obs,j+1), out_x(i_obs,i+1))
                    distrib_matrix(i_obs,i,j) = (ub-lb)/(out_x(i_obs,i)-out_x(i_obs,i+1))
                end if
            end do
        end do
    end do ! i_obs
    !$omp end do
    deallocate(breakpoints, ext_breakpoints, cumul)
    !$omp end parallel

    deallocate(in_x, out_x)

end subroutine inv_rebin_arr

subroutine inv_rebin(in_levels, out_levels, in_x_grid, in_y_vals, out_x_grid, out_y_vals, distrib_matrix)

    ! out_x_grid is the output grid except the end-points, since
    ! the end-points of out_x_grid are the same as those of in_x_grid
    ! in_x_grid and out_x_grid are in descending order, and distrib_matrix is such that
    !     out_y_vals = (distrib_matrix) (in_y_vals)
    integer, intent(in)           :: in_levels, out_levels
    double precision, intent(in)  :: in_x_grid(in_levels), in_y_vals(in_levels-1), out_x_grid(out_levels)
    double precision, intent(out) :: out_y_vals(out_levels+1), distrib_matrix(out_levels+1,in_levels-1)
    logical                       :: input_OK
    integer                       :: i, breakpoints(out_levels), j, ext_breakpoints(out_levels+2)
    double precision              :: cumul(out_levels+2)
    double precision              :: in_x(in_levels), out_x(out_levels+2)
    double precision              :: mass_error, dummy(out_levels+1), lb, ub
    character(len=1024)           :: error_msg
    character(len=*), parameter   :: rname = mname//'/inv_rebin'

    ! first some error checking
    input_OK = .true.
    ! check whether x_grid's are sorted in descending order and unique
    do i=lbound(in_x_grid,1),ubound(in_x_grid,1)-1
        if (in_x_grid(i) <= in_x_grid(i+1)) then
            input_OK = .false.
            error_msg = 'Input grid is not in descending order'
        end if
    end do
    do i=lbound(out_x_grid,1),ubound(out_x_grid,1)-1
        if (out_x_grid(i) <= out_x_grid(i+1)) then
            input_OK = .false.
            error_msg = 'Output grid is not in descending order'
        end if
    end do
    ! check whether the output grid is contained inside the input grid
    if (out_x_grid(lbound(out_x_grid,1)) > in_x_grid(lbound(in_x_grid,1))) then
        input_OK = .false.
        error_msg = 'First value of output grid is higher than first value of input grid'
    end if
    if (out_x_grid(ubound(out_x_grid,1)) < in_x_grid(ubound(in_x_grid,1))) then
        input_OK = .false.
        error_msg = 'Last value of output grid is less than last value of input grid'
    end if
    if ( .not. input_OK ) then
        PRINT_ERROR(error_msg)
        stop
    end if
    ! now that we know the input is OK, proceed with rebinning :-)

    ! copy over data
    in_x(:) = in_x_grid(:)
    out_x(2:size(out_x)-1) = out_x_grid(:)
    out_x(1) = in_x(1)
    out_x(size(out_x)) = in_x(size(in_x))

    breakpoints = inv_searchSorted(in_x, out_x_grid)
    ext_breakpoints(2:out_levels+1) = breakpoints(:)
    ext_breakpoints(1) = 1
    ext_breakpoints(out_levels+2) = in_levels-1
    !print *, ext_breakpoints
    cumul(1) = 0.0
    do i=1,size(breakpoints)
        cumul(i+1) = sum(in_y_vals(:breakpoints(i)-1) * inv_diff(in_x(:breakpoints(i))))
        cumul(i+1) = cumul(i+1) + in_y_vals(breakpoints(i)) * (in_x(breakpoints(i)) - out_x(i+1))
    end do
    cumul(size(cumul)) = sum(in_y_vals * inv_diff(in_x))
    !write(*,'(a,f20.8)') 'Total input mass is ', sum(in_y_vals * inv_diff(in_x))
    out_y_vals = diff(cumul)/inv_diff(out_x)
    !write(*,'(a,f20.8)') 'Total output mass is ', sum(out_y_vals * inv_diff(out_x))

    ! now calculate the distribution matrix
    distrib_matrix = 0.0 ! set all out_levels+1 x in_levels-1 elements to 0
    ! out_y_vals(i) has contributions from in_y_vals(ext_breakpoints(i)) to in_y_vals(ext_breakpoints(i+1))
    ! in_y_vals(j) lies between in_x(j) and in_x(j+1)
    ! hence, out_y_vals(i), which lines between out_x(i) and out_x(i+1), will need to examine levels in_x(ext_breakpoints(i)) and in_x(ext_breakpoints(i+1)+1)
    do i=1,out_levels+1
        do j=1,in_levels-1
            ! output level i is between out_x(i) and out_x(i+1)
            ! input level j is between in_x(j) and in_x(j+1)
            ! do something only if there is some overlap
            if ((in_x(j+1) .ge. out_x(i)) .or. (out_x(i+1) .ge. in_x(j))) then
                continue
            else
                ub = min(out_x(i),in_x(j))
                lb = max(in_x(j+1),out_x(i+1))
                distrib_matrix(i,j) = (ub-lb)/(out_x(i)-out_x(i+1))
            end if
        end do
    end do

end subroutine inv_rebin

pure function searchSorted(in_array, values) result (index_array)

    ! searchSorted = return the indices in array in_array such that
    ! in_array(idx(i)) <= value(i) < in_array(idx(i)+1). If value(i) is outside
    ! in_array then either lbound(in_array)-1 or ubound(in_array) is returned.

    double precision, dimension(:), intent(in)  :: in_array, values
    integer, dimension(size(values))            :: index_array
    integer                                     :: ub, lb, i, j, idx
    character(len=*), parameter                 :: rname = mname//'/searchSorted'

    lb = lbound(in_array,1)
    ub = ubound(in_array,1)
    do j=lbound(values,1),ubound(values,1)
        idx = lb-1
        do i=lb,ub
            if (values(j) .lt. in_array(i)) exit
            idx = i
        end do
        index_array(j+1-lbound(values,1)) = idx
    end do

end function searchSorted

function inv_searchSorted(in_array, values) result (index_array)

    ! searchSorted = return the indices in array in_array such that
    ! in_array(idx(i)) >= value(i) > in_array(idx(i)+1). If value(i) is outside
    ! in_array then either lbound(in_array)-1 or ubound(in_array) is returned.

    double precision, dimension(:), intent(in) :: in_array, values
    integer, dimension(size(values))         :: index_array
    integer                                  :: ub, lb, i, j, idx
    character(len=*), parameter              :: rname = mname//'/inv_searchSorted'

    lb = lbound(in_array,1)
    ub = ubound(in_array,1)
    do j=lbound(values,1),ubound(values,1)
        idx = lb-1
        do i=lb,ub
            if (values(j) .ge. in_array(i)) exit
            idx = i
        end do
        index_array(j+1-lbound(values,1)) = idx
    end do

end function inv_searchSorted

pure function diff(in_array)

    double precision, intent(in), dimension(:)    :: in_array
    double precision, dimension(size(in_array)-1) :: diff
    integer                                     :: i, lb
    character(len=*), parameter                 :: rname = mname//'/diff'

    lb = lbound(in_array,1)
    do i=1,size(diff)
        diff(i) = in_array(i+lb)-in_array(i+lb-1)
    end do

end function diff

pure function inv_diff(in_array)

    double precision, intent(in), dimension(:)      :: in_array
    double precision, dimension(size(in_array)-1)   :: inv_diff
    integer                                         :: i, lb
    character(len=*), parameter                     :: rname = mname//'/inv_diff'

    lb = lbound(in_array,1)
    do i=1,size(inv_diff)
        inv_diff(i) = in_array(i+lb-1)-in_array(i+lb)
    end do

end function inv_diff

end module interpolate_fields
