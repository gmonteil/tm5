!
! 3d grid transforms
!
!
!  call Fill3D( lli, levi, ps, field, &
!               lliX, leviX, fieldX, &
!               combkey, status )
!
!    o lli,levi   : output field defintions  (in)
!      ps         : surface pressure in hPa  (in)
!      field      : output field             (out)
!
!    o lliX, leviX : input field definitons  (in)
!      fieldX      : input field             (in)
!
!    o combkey : 'mass-aver', 'sum', 'aver'   (in)
!
!    o note that psX is not required ...
!
!  call FillMass      ( m    , lli, levi, sp          , status )
!
!    ! Fill 3D mass array given surface pressure and grid cell aera's.
!
!  call FillMassChange( dm_dt, lli, levi, sp_t1, sp_t2, status )
!
!    ! Fill 3D mass change between two times given different
!    ! surface pressures;
!    ! note that ( a + b sp2 ) - ( a + b sp1 )
!    ! is not the same as  a + b (sp2 - sp1) ...
!
!
!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_NOTOK_STOP if (status/=0) then; TRACEBACK; stop; end if
!
!###############################################################################

module grid_3d

  use GO, only : gol, goPr, goErr

  implicit none

  ! --- in/out -----------------------------

  private

  public  ::  Fill3D
  public  ::  FillMass
  public  ::  FillMassChange
  public  ::  rebin_levels, regrid_2d, regrid_3d_mix, regrid_3d_PCTM

  ! --- const ---------------------------------

  character(len=*), parameter  ::  mname = 'grid_3d'

contains

  subroutine regrid_2d(ip_lats, ip_lons, ip_data, ip_per_area, op_lats, op_lons, op_data, op_per_area, status, print_diag)

    ! Arjo's regridding routines assume that the grids are multiples or submultiples of each other.
    ! This is not always the case. This routine can redistribute a 2D array from one grid to another,
    ! without even assuming that the grids are regular. The arguments are:
    !
    ! ip_lats, ip_lons (in) : Array of latitudes and longitudes on which the data are defined. These denote
    !                         boundaries and not midpoints. For example, for a global 1x1 grid, ip_lats will
    !                         be -90., -89., ... 90., and ip_lons will be -180., -179., ... 180.
    ! ip_data (in)          : 2D data, dimension (len(ip_lons)-1) x (len(ip_lats)-1), containing the data to be redistributed.
    ! ip_per_area (in)      : Logical variable, true if the input data are per area, false if they are per grid cell
    ! op_lats, op_lons (in) : Array of latitudes and longitudes defining the output grid. Again, these are boundaries.
    ! op_data (out)         : 2D array of regridded data.
    ! op_per_area (in)      : Logical variable, true => output data are per area, false => output data are per grid cell

    use binas,  only : R_e => ae ! radius of the Earth
    use binas,  only : deg2rad

    implicit none

    real, intent(in)        :: ip_lats(:), ip_lons(:), op_lats(:), op_lons(:)
    real, intent(in)        :: ip_data(:,:)
    real, intent(out)       :: op_data(:,:)
    logical, intent(in)     :: ip_per_area, op_per_area
    integer, intent(out)    :: status
    logical, intent(in), optional   :: print_diag

    character(len=*), parameter ::  rname = mname//'/regrid_2d'

    integer                 :: im_ip, im_op, jm_ip, jm_op
    real, allocatable       :: op_surf_area(:,:), ip_surf_area(:,:), ip_data_per_area(:,:)
    real, allocatable       :: dlon_ip(:), dlat_ip(:), dlon_op(:), dlat_op(:), op_data_stripe(:)
    integer                 :: i,j,k,l
    real                    :: ovlap_lat_up, ovlap_lat_dn, ovlap_lon_up, ovlap_lon_dn, ovlap_area
    logical                 :: diag
    real                    :: ip_total, op_total

    diag = .false.
    if (present(print_diag)) diag = print_diag

    jm_ip = size(ip_data, 2)
    im_ip = size(ip_data, 1)
    jm_op = size(op_data, 2)
    im_op = size(op_data, 1)

    ! check sizes
    if (size(ip_lats) /= jm_ip+1 .or. size(ip_lons) /= im_ip+1) then
        write(*,'("Error in ", a)') rname
        write(*,'("Size of input data = ", i3, " x ", i3)') im_ip, jm_ip
        write(*,'("Length of input boundaries = ", i3, ", ", i3)') size(ip_lons), size(ip_lats)
        status=1
        return
    end if

    if (size(op_lats) /= jm_op+1 .or. size(op_lons) /= im_op+1) then
        write(*,'("Error in ", a)') rname
        write(*,'("Size of output data = ", i3, " x ", i3)') im_op, jm_op
        write(*,'("Length of output boundaries = ", i3, ", ", i3)') size(op_lons), size(op_lats)
        status=1
        return
    end if

    ! In the routine below, it's easier to handle the input data if it is defined per area. If it is not,
    ! we will convert it to per area.
    allocate(ip_data_per_area(im_ip, jm_ip))

    if (diag .or. (.not. ip_per_area)) then
        allocate(dlon_ip(im_ip))
        allocate(dlat_ip(jm_ip))
        allocate(ip_surf_area(im_ip,jm_ip))
        dlon_ip = deg2rad * (ip_lons(2:im_ip+1) - ip_lons(1:im_ip))
        dlat_ip = sin(deg2rad*ip_lats(2:jm_ip+1)) - sin(deg2rad*ip_lats(1:jm_ip))
        do j = 1, jm_ip
            do i = 1, im_ip
                ip_surf_area(i,j) = dlon_ip(i) * dlat_ip(j)
            end do
        end do
        ip_surf_area = R_e * R_e * ip_surf_area
        deallocate(dlon_ip, dlat_ip)
    end if

    if (.not. ip_per_area) then
        ip_data_per_area = ip_data/ip_surf_area
    else
        ip_data_per_area = ip_data
    end if

    if (op_per_area) then
        allocate(dlon_op(im_op))
        allocate(dlat_op(jm_op))
        allocate(op_surf_area(im_op,jm_op))
        dlon_op = deg2rad * (op_lons(2:im_op+1) - op_lons(1:im_op))
        dlat_op = sin(deg2rad*op_lats(2:jm_op+1)) - sin(deg2rad*op_lats(1:jm_op))
        do j = 1, jm_op
            do i = 1, im_op
                op_surf_area(i,j) = dlon_op(i) * dlat_op(j)
            end do
        end do
        op_surf_area = R_e * R_e * op_surf_area
        deallocate(dlon_op, dlat_op)
    end if

    ! time for redistribution
    op_data = 0.0
    !$omp parallel default(shared) private(i, j, k, l, op_data_stripe) &
    !$omp private(ovlap_lat_up, ovlap_lat_dn, ovlap_lon_up, ovlap_lon_dn, ovlap_area)
    allocate(op_data_stripe(jm_op))
    !$omp do schedule(static)
    do k = 1, im_op ! lon limits are op_lons(k) and op_lons(k+1)
        op_data_stripe = 0.0
        do l = 1, jm_op ! lat limits are op_lats(l) and op_lats(l+1)
            ! the box with output data is between latitudes op_lats(l) and op_lats(l+1),
            ! and between longitudes op_lons(k) and op_lons(k+1)
            do i = 1, im_ip
                if ((op_lons(k) .ge. ip_lons(i+1)) .or. (op_lons(k+1) .le. ip_lons(i))) cycle
                do j = 1, jm_ip
                    if ((op_lats(l) .ge. ip_lats(j+1)) .or. (op_lats(l+1) .le. ip_lats(j))) cycle
                    ! the box with input data is between latitudes ip_lats(j) and ip_lats(j+1),
                    ! and between longitudes ip_lons(i) and ip_lons(i+1)
                    ovlap_lat_up = min(ip_lats(j+1),op_lats(l+1))
                    ovlap_lat_dn = max(ip_lats(j),op_lats(l))
                    ovlap_lon_up = min(ip_lons(i+1),op_lons(k+1))
                    ovlap_lon_dn = max(ip_lons(i),op_lons(k))
                    ovlap_area = R_e * R_e * deg2rad * (ovlap_lon_up-ovlap_lon_dn) * (sin(deg2rad*ovlap_lat_up) - sin(deg2rad*ovlap_lat_dn))
                    !op_data(k,l) = op_data(k,l) + ip_data_per_area(i,j) * ovlap_area
                    op_data_stripe(l) = op_data_stripe(l) + ip_data_per_area(i,j) * ovlap_area
                end do ! j
            end do ! i
        end do ! l
        op_data(k,:) = op_data_stripe
    end do ! k
    !$omp end do
    deallocate(op_data_stripe)
    !$omp end parallel

    if (op_per_area) then
        op_data = op_data/op_surf_area
    end if

    if (diag) then
        if (ip_per_area) then
            ip_total = sum(ip_data * ip_surf_area)
        else
            ip_total = sum(ip_data)
        end if
        if (op_per_area) then
            op_total = sum(op_data * op_surf_area)
        else
            op_total = sum(op_data)
        end if
        write(gol,'(a," : total before regridding = ", es20.12)') rname, ip_total ; call goPr
        write(gol,'(a," :  total after regridding = ", es20.12)') rname, op_total ; call goPr
        write(gol,'(a," :        fractional error = ", es20.12)') rname, (op_total-ip_total)/ip_total ; call goPr
    end if

    deallocate(ip_data_per_area)
    if (allocated(ip_surf_area)) deallocate(ip_surf_area)
    if (allocated(op_surf_area)) deallocate(op_surf_area)

    status = 0

  end subroutine regrid_2d

  subroutine rebin_levels(ip_levels, ip_data, ip_mix, op_levels, op_data, status)

    ! This routine rebins a column defined on one vertical to another vertical grid. The boolean 'ip_mix'
    ! controls whether the input is an array of mixing ratios, or whether they are absolute number of molecules.
    !
    ! ip_levels : If ip_mix is true, then this is an array of pressure boundaries defining the input grid.
    !             If ip_mix is false, then this is an array of altitudes. In either case, they start from the
    !             ground and finish at the top, i.e., they're descending in value.
    ! ip_data   : If ip_mix is true, then this is an array of mixing ratios or 3D densities. If ip_mix is false,
    !             then this is an array of 2D densities or absolute numbers.
    ! ip_mix    : Boolean which controls whether we are dealing with mixing ratios or absolute numbers. In either
    !             case, the output data have the same meaning as the input data, i.e., if the input data are
    !             mixing ratios, so are the output data, and so on.
    ! op_levels : Similar to ip_levels, an array of either pressure boundaries or altitudes, from the ground up.
    !             There is no requirement that the end-points match with the ip_levels, although if they do,
    !             total mass will be conserved.
    ! op_data   : Similar to ip_data, rebinned mixing ratios or numbers.

    implicit none

    real, intent(in)        :: ip_levels(:), ip_data(:), op_levels(:)
    logical, intent(in)     :: ip_mix
    real, intent(out)       :: op_data(:)
    integer, intent(out)    :: status

    character(len=*), parameter ::  rname = mname//'/rebin_levels'

    integer                 :: nlay_ip, nlay_op, i, j
    real                    :: ovlap_lev_up, ovlap_lev_dn

    ! check the sizes
    nlay_ip = size(ip_data)
    if (nlay_ip /= size(ip_levels)-1) then
        write(*,'("Error in ", a)') rname
        write(*,'("Size of input data = ", i3)') nlay_ip
        write(*,'("Length of input boundaries = ", i3)') size(ip_levels)
        status=1
        return
    end if

    nlay_op = size(op_data)
    if (nlay_op /= size(op_levels)-1) then
        write(*,'("Error in ", a)') rname
        write(*,'("Size of input data = ", i3)') nlay_op
        write(*,'("Length of input boundaries = ", i3)') size(op_levels)
        status=1
        return
    end if

    op_data = 0.0

    if (ip_mix) then
        ! Levels are going down in value as we go up in indices
        do i = 1, nlay_op
            ! We will assemble the total mass between op_levels(i) and op_levels(i+1)
            do j = 1, nlay_ip
                ! Calculate the contribution from the input layer between ip_levels(j) and ip_levels(j+1)
                if (ip_levels(j+1) .ge. op_levels(i) .or. ip_levels(j) .le. op_levels(i+1)) cycle
                ! What is the overlap between the i-th output level and the j-th input level?
                ovlap_lev_up = max(op_levels(i+1), ip_levels(j+1)) ! max chooses the higher pressure, i.e., lower altitude
                ovlap_lev_dn = min(op_levels(i), ip_levels(j)) ! min chooses the lower pressure, i.e., higher altitude
                op_data(i) = op_data(i) + ip_data(j) * (ovlap_lev_dn - ovlap_lev_up)
            end do
        end do
        ! Divide the output data by the level thickness to get back to mixing ratios
        do i = 1, nlay_op
            op_data(i) = op_data(i)/(op_levels(i)-op_levels(i+1))
        end do
    else
        ! Levels are going up in value as we go up in indices
        do i = 1, nlay_op
            ! We will assemble the total mass between op_levels(i) and op_levels(i+1)
            do j = 1, nlay_ip
                ! Calculate the contribution from the input layer between ip_levels(j) and ip_levels(j+1)
                if (ip_levels(j+1) .le. op_levels(i) .or. ip_levels(j) .ge. op_levels(i+1)) cycle
                ! What is the overlap between the i-th output level and the j-th input level?
                ovlap_lev_up = min(op_levels(i+1), ip_levels(j+1))
                ovlap_lev_dn = max(op_levels(i), ip_levels(j))
                op_data(i) = op_data(i) + ip_data(j) * (ovlap_lev_up-ovlap_lev_dn)/(ip_levels(j+1)-ip_levels(j))
            end do
        end do
    end if

    status = 0

  end subroutine

  subroutine regrid_3d_PCTM(lli_ip, lli_op, ilev_op, pmids_in, psurf_in, in_mix, out_mix, status, print_diag)
    ! Abhishek has given me some PCTM initial conditions. I do not know if those are defined on
    ! sigma-pressure coordinates, so I can't use regrid_3d_mix. Instead, he has given me arrays
    ! with the surface pressure as well as mid-level pressures.

    use grid_type_ll , only : TllGridInfo
    use grid_type_hyb, only : TLevelInfo
    use binas,         only : grav, amu, xmair

    type(TllGridInfo), intent(in)   :: lli_ip, lli_op
    type(TLevelInfo), intent(in)    :: ilev_op
    real, intent(in)                :: pmids_in(:,:,:) ! mid-level pressures, in Pascals
    real, intent(in)                :: psurf_in(:,:)   ! surface pressures, in Pascals
    real, intent(in)                :: in_mix(:,:,:)   ! input mixing ratio
    real, intent(out)               :: out_mix(:,:,:)  ! output mixing ratio
    integer, intent(out)            :: status
    logical, intent(in), optional   :: print_diag

    ! --- const ---------------------------------

    character(len=*), parameter     ::  rname = mname//'/regrid_3d_PCTM'

    ! --- local ---------------------------------

    real, allocatable               :: inv_mat(:,:), ps_out(:,:), pres_in(:,:,:), pres_thick_in(:,:,:), pmid_vec(:), pres_out(:,:,:)
    real, allocatable               :: regrid_pthick_in(:,:), intermediate_mix(:,:,:), pb_out(:), pb_in(:,:,:), pres_thick_out(:,:,:)
    real, allocatable               :: pctm_lons_shifted(:), pctm_sp_shifted(:,:), pctm_pmid_shifted(:,:,:), pctm_mix_shifted(:,:,:)
    integer                         :: nlev_in, nlat_in, nlon_in, nlay_in, nlev_out, nlay_out, nlat_out, nlon_out, i, j
    logical                         :: diag_msg
    real                            :: total_mass_before, total_mass_after, filler

    if (present(print_diag)) then
        diag_msg = print_diag
    else
        diag_msg = .false.
    end if

    nlon_in = size(in_mix, 1)
    nlat_in = size(in_mix, 2)
    nlay_in = size(in_mix, 3)
    nlev_in = nlay_in + 1

    ! check sizes
    if (lli_ip%nlon /= nlon_in .or. lli_ip%nlat /= nlat_in .or. size(pmids_in,3) /= nlay_in) then
        write(*,'("Error in ", a)') rname
        write(*,'("Size of input mixing ratio array = ", i3, " x ", i3, " x ", i3)') nlon_in, nlat_in, nlay_in
        write(*,'("Size according to input grid = ", i3, " x ", i3, " x ", i3)') lli_ip%nlon, lli_ip%nlat, size(pmids_in,3)
        status=1
        return
    end if

    nlay_out = size(out_mix, 3)
    nlat_out = size(out_mix, 2)
    nlon_out = size(out_mix, 1)
    nlev_out = nlay_out + 1

    ! check sizes
    if (lli_op%nlon /= nlon_out .or. lli_op%nlat /= nlat_out .or. ilev_op%nlev /= nlay_out) then
        write(*,'("Error in ", a)') rname
        write(*,'("Size of output mixing ratio array = ", i3, " x ", i3, " x ", i3)') nlon_out, nlat_out, nlay_out
        write(*,'("Size according to output grid = ", i3, " x ", i3, " x ", i3)') lli_op%nlon, lli_op%nlat, ilev_op%nlev
        status=1
        return
    end if

    ! also the input surface pressure
    if (nlon_in /= size(psurf_in, 1) .or. nlat_in /= size(psurf_in, 2)) then
        write(*,'("Error in ", a)') rname
        write(*,'("Size of input surface pressure array = ", i3, " x ", i3)') size(psurf_in, 1), size(psurf_in, 2)
        write(*,'("Size according to input grid = ", i3, " x ", i3)') lli_ip%nlon, lli_ip%nlat
        status=1
        return
    end if

    ! create the output surface pressure array
    allocate(ps_out(nlon_out, nlat_out))
    ! PCTM grid is strange, with longitude boundaries
    ! -180.625, -179.375, -178.125, ... 179.375
    ! This creates a problem, because if the output grid is -180 to 180, then the regridding
    ! routine doesn't know that the last cell should contain stuff from PCTM's first cell. To
    ! solve this, we 'extend' the PCTM grid to the right by one cell, and copy over data from
    ! the first cell to this additional cell.
    allocate(pctm_lons_shifted(nlon_in+2))
    pctm_lons_shifted(1:nlon_in+1) = lli_ip%blon_deg
    pctm_lons_shifted(nlon_in+2) = pctm_lons_shifted(nlon_in+1) + lli_ip%dlon_deg
    allocate(pctm_sp_shifted(nlon_in+1, nlat_in))
    pctm_sp_shifted(1:nlon_in,:) = psurf_in
    pctm_sp_shifted(nlon_in+1,:) = psurf_in(1,:)
    allocate(pctm_pmid_shifted(nlon_in+1, nlat_in, nlay_in))
    pctm_pmid_shifted(1:nlon_in,:,:) = pmids_in
    pctm_pmid_shifted(nlon_in+1,:,:) = pmids_in(1,:,:)
    allocate(pctm_mix_shifted(nlon_in+1, nlat_in, nlay_in))
    pctm_mix_shifted(1:nlon_in,:,:) = in_mix
    pctm_mix_shifted(nlon_in+1,:,:) = in_mix(1,:,:)

    call regrid_2d(lli_ip%blat_deg, pctm_lons_shifted, pctm_sp_shifted, .true., lli_op%blat_deg, lli_op%blon_deg, &
        ps_out, .true., status)
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    if (diag_msg) then

        total_mass_before = 0.0
        do j = 1, nlat_in
            total_mass_before = total_mass_before + sum(psurf_in(:,j) * lli_ip%area_m2(j))
        end do
        total_mass_before = total_mass_before/grav

        total_mass_after = 0.0
        do j = 1, nlat_out
            total_mass_after = total_mass_after + sum(ps_out(:,j) * lli_op%area_m2(j))
        end do
        total_mass_after = total_mass_after/grav

        write(gol,'("Total air mass before regridding = ", es20.12, " Kg")') total_mass_before ; call goPr
        write(gol,'(" Total air mass after regridding = ", es20.12, " Kg")') total_mass_after ; call goPr
        write(gol,'("   Airmass error from regridding = ", es20.12, " Kg")') total_mass_after - total_mass_before ; call goPr

    end if

    ! first we need to do the lateral redistribution, keeping the vertical grid the same
    ! create the 3D input pressure grid
    allocate(pres_in(nlon_in+1, nlat_in, nlev_in))
    allocate(pres_thick_in(nlon_in+1, nlat_in, nlay_in))

    ! create the inversion matrix to go from the mid-level pressures to the pressure boundaries
    allocate(inv_mat(nlay_in, nlay_in))
    inv_mat = 0.0
    do i = 1, nlay_in
        if (mod(i,2) == 1) then
            filler = 2.0
        else
            filler = -2.0
        end if

        do j = 1, i
            inv_mat(i,j) = filler
            filler = -filler
        end do
    end do

    allocate(pmid_vec(nlay_in))
    do i = 1, nlon_in+1
        do j = 1, nlat_in
            pmid_vec = pctm_pmid_shifted(i,j,:)
            pmid_vec(1) = pmid_vec(1) - 0.5*pctm_sp_shifted(i,j)
            pres_in(i,j,2:nlay_in) = matmul(inv_mat, pmid_vec)
            pres_in(i,j,1) = pctm_sp_shifted(i,j)
        end do
    end do
    deallocate(pmid_vec, inv_mat)

    do i = 1, nlay_in
        pres_thick_in(:,:,i) = pres_in(:,:,i) - pres_in(:,:,i+1)
    end do

    allocate(regrid_pthick_in(nlon_out, nlat_out))
    allocate(intermediate_mix(nlon_out, nlat_out, nlay_in))
    do i = 1, nlay_in

!        call regrid_2d(lli_ip%blat_deg, lli_ip%blon_deg, in_mix(:,:,i)*pres_thick_in(:,:,i), .true., &
!            lli_op%blat_deg, lli_op%blon_deg, intermediate_mix(:,:,i), .false., status)
        call regrid_2d(lli_ip%blat_deg, pctm_lons_shifted, pctm_mix_shifted(:,:,i)*pres_thick_in(:,:,i), .true., &
            lli_op%blat_deg, lli_op%blon_deg, intermediate_mix(:,:,i), .false., status)
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
!        call regrid_2d(lli_ip%blat_deg, lli_ip%blon_deg, pres_thick_in(:,:,i), .true., &
!            lli_op%blat_deg, lli_op%blon_deg, regrid_pthick_in, .false., status)
        call regrid_2d(lli_ip%blat_deg, pctm_lons_shifted, pres_thick_in(:,:,i), .true., &
            lli_op%blat_deg, lli_op%blon_deg, regrid_pthick_in, .false., status)
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
        intermediate_mix(:,:,i) = intermediate_mix(:,:,i)/regrid_pthick_in

    end do
    deallocate(regrid_pthick_in)
    ! now intermediate_mix is the mixing ratio on the input vertical grid but the output lateral grid, so redistribute vertically

    ! pb_out are the pressure boundaries on the lateral output grid and vertical output grid
    allocate(pb_out(nlev_out))
    ! pb_in is the pressure grid on the output lateral grid and the input vertical grid
    allocate(pb_in(nlon_out, nlat_out, nlev_in))
    do i = 1, nlev_in
        !call regrid_2d(lli_ip%blat_deg, lli_ip%blon_deg, pres_in(:,:,i), .true., lli_op%blat_deg, lli_op%blon_deg, pb_in(:,:,i), .true., status)
        call regrid_2d(lli_ip%blat_deg, pctm_lons_shifted, pres_in(:,:,i), .true., lli_op%blat_deg, &
            lli_op%blon_deg, pb_in(:,:,i), .true., status)
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    end do

    do i = 1, nlon_out
        do j = 1, nlat_out
            pb_out = ilev_op%a + ilev_op%b * ps_out(i,j)
            call rebin_levels(pb_in(i,j,:), intermediate_mix(i,j,:), .true., &
                pb_out, out_mix(i,j,:), status)
            if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
        end do
    end do

    ! Check whether the total mass has been conserved, sort of a diagnostic
    if (diag_msg) then
        total_mass_before = 0.0
        do i = 1, nlon_in
            do j = 1, nlat_in
                total_mass_before = total_mass_before + lli_ip%area_m2(j) * sum(in_mix(i,j,:) * pres_thick_in(i,j,:))
            end do
        end do
        total_mass_before = total_mass_before/(grav * amu * xmair)

        allocate(pres_out(nlon_out, nlat_out, nlev_out))
        allocate(pres_thick_out(nlon_out, nlat_out, nlay_out))
        do i = 1, nlev_out
            pres_out(:,:,i) = ilev_op%a(i-1) + ilev_op%b(i-1) * ps_out(:,:)
        end do
        do i = 1, nlay_out
            pres_thick_out(:,:,i) = pres_out(:,:,i) - pres_out(:,:,i+1)
        end do

        total_mass_after = 0.0
        do i = 1, nlon_out
            do j = 1, nlat_out
                total_mass_after = total_mass_after + lli_op%area_m2(j) * sum(out_mix(i,j,:) * pres_thick_out(i,j,:))
            end do
        end do
        total_mass_after = total_mass_after/(grav * amu * xmair)

        deallocate(pres_out, pres_thick_out)

        write(*,'("Total number of tracer molecules before regridding = ", es20.12)') total_mass_before
        write(*,'(" Total number of tracer molecules after regridding = ", es20.12)') total_mass_after
        write(*,'("                                  Regridding error = ", es20.12)') total_mass_after - total_mass_before
    end if

    deallocate(ps_out, pres_in, pres_thick_in, pb_out, pb_in, intermediate_mix)
    deallocate(pctm_lons_shifted, pctm_sp_shifted, pctm_pmid_shifted, pctm_mix_shifted)

    status = 0

  end subroutine regrid_3d_PCTM

  subroutine regrid_3d_mix(lli_ip, lli_op, ilev_ip, ilev_op, ps_in, in_mix, out_mix, status, print_diag)

    ! Given a 3D mixing ratio field defined on the lateral grid lli_ip and the vertical grid ilev_ip,
    ! redistribute it on the lateral grid lli_op and the vertical grid ilev_op. For now, this handles
    ! only mixing ratios. In the future, we'll add support for redistributing number fields.

    use grid_type_ll , only : TllGridInfo
    use grid_type_hyb, only : TLevelInfo
    use binas,         only : grav, amu, xmair

    type(TllGridInfo), intent(in)   :: lli_ip, lli_op
    type(TLevelInfo), intent(in)    :: ilev_ip, ilev_op
    real, intent(in)                :: in_mix(:,:,:), ps_in(:,:)
    real, intent(out)               :: out_mix(:,:,:)
    integer, intent(out)            :: status
    logical, intent(in), optional   :: print_diag

    ! --- const ---------------------------------

    character(len=*), parameter     ::  rname = mname//'/regrid_3d_mix'

    real, allocatable               :: ps_out(:,:), intermediate_mix(:,:,:), pres_in(:,:,:), pres_thick_in(:,:,:)
    real, allocatable               :: pres_out(:,:,:), pres_thick_out(:,:,:), pb_out(:), regrid_pthick_in(:,:), pb_in(:)
    integer                         :: i, j, nlat_in, nlon_in, nlev_in, nlon_out, nlat_out, nlev_out
    real                            :: total_mass_before, total_mass_after
    logical                         :: diag_msg

    if (present(print_diag)) then
        diag_msg = print_diag
    else
        diag_msg = .false.
    end if

    nlev_in = size(in_mix, 3)
    nlat_in = size(in_mix, 2)
    nlon_in = size(in_mix, 1)

    ! check sizes
    if (lli_ip%nlon /= nlon_in .or. lli_ip%nlat /= nlat_in .or. ilev_ip%nlev /= nlev_in) then
        write(*,'("Error in ", a)') rname
        write(*,'("Size of input mixing ratio array = ", i3, " x ", i3, " x ", i3)') nlon_in, nlat_in, nlev_in
        write(*,'("Size according to input grid = ", i3, " x ", i3, " x ", i3)') lli_ip%nlon, lli_ip%nlat, ilev_ip%nlev
        status=1
        return
    end if

    nlev_out = size(out_mix, 3)
    nlat_out = size(out_mix, 2)
    nlon_out = size(out_mix, 1)

    ! check sizes
    if (lli_op%nlon /= nlon_out .or. lli_op%nlat /= nlat_out .or. ilev_op%nlev /= nlev_out) then
        write(*,'("Error in ", a)') rname
        write(*,'("Size of output mixing ratio array = ", i3, " x ", i3, " x ", i3)') nlon_out, nlat_out, nlev_out
        write(*,'("Size according to output grid = ", i3, " x ", i3, " x ", i3)') lli_op%nlon, lli_op%nlat, ilev_op%nlev
        status=1
        return
    end if

    ! also the input surface pressure
    if (nlon_in /= size(ps_in, 1) .or. nlat_in /= size(ps_in, 2)) then
        write(*,'("Error in ", a)') rname
        write(*,'("Size of input surface pressure array = ", i3, " x ", i3)') size(ps_in, 1), size(ps_in, 2)
        write(*,'("Size according to input grid = ", i3, " x ", i3)') lli_ip%nlon, lli_ip%nlat
        status=1
        return
    end if

    ! create the output surface pressure array
    allocate(ps_out(nlon_out, nlat_out))
    call regrid_2d(lli_ip%blat_deg, lli_ip%blon_deg, ps_in, .true., lli_op%blat_deg, lli_op%blon_deg, ps_out, .true., status, print_diag=diag_msg)
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    if (diag_msg) then

        total_mass_before = 0.0
        do j = 1, nlat_in
            total_mass_before = total_mass_before + sum(ps_in(:,j) * lli_ip%area_m2(j))
        end do
        total_mass_before = total_mass_before/grav

        total_mass_after = 0.0
        do j = 1, nlat_out
            total_mass_after = total_mass_after + sum(ps_out(:,j) * lli_op%area_m2(j))
        end do
        total_mass_after = total_mass_after/grav

        write(gol,'(a, " : total air mass before regridding = ", es20.12, " Kg")') rname, total_mass_before ; call goPr
        write(gol,'(a, " :  total air mass after regridding = ", es20.12, " Kg")') rname, total_mass_after ; call goPr
        write(gol,'(a, " :         fractional airmass error = ", es20.12, " Kg")') rname, (total_mass_after - total_mass_before)/total_mass_before ; call goPr
    end if

    ! first we need to do the lateral redistribution, keeping the vertical grid the same
    ! create the 3D input pressure grid
    allocate(pres_in(nlon_in, nlat_in, nlev_in+1))
    allocate(pres_thick_in(nlon_in, nlat_in, nlev_in))
    do i = 1, nlev_in+1
        pres_in(:,:,i) = ilev_ip%a(i-1) + ilev_ip%b(i-1) * ps_in(:,:)
    end do
    do i = 1, nlev_in
        pres_thick_in(:,:,i) = pres_in(:,:,i) - pres_in(:,:,i+1)
    end do

    allocate(regrid_pthick_in(nlon_out, nlat_out))
    allocate(intermediate_mix(nlon_out, nlat_out, nlev_in))
    do i = 1, nlev_in
        call regrid_2d(lli_ip%blat_deg, lli_ip%blon_deg, in_mix(:,:,i)*pres_thick_in(:,:,i), .true., &
            lli_op%blat_deg, lli_op%blon_deg, intermediate_mix(:,:,i), .false., status)
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
        call regrid_2d(lli_ip%blat_deg, lli_ip%blon_deg, pres_thick_in(:,:,i), .true., &
            lli_op%blat_deg, lli_op%blon_deg, regrid_pthick_in, .false., status)
        if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
        intermediate_mix(:,:,i) = intermediate_mix(:,:,i)/regrid_pthick_in

    end do
    deallocate(regrid_pthick_in)
    ! now intermediate_mix is the mixing ratio on the input vertical grid but the output lateral grid, so redistribute vertically

    allocate(pb_out(nlev_out+1))
    allocate(pb_in(nlev_in+1))
    ! pb_in are the pressure boundaries on the lateral output grid but vertical input grid
    ! pb_out are the pressure boundaries on the lateral output grid and vertical output grid
    do i = 1, nlon_out
        do j = 1, nlat_out
            pb_in = ilev_ip%a + ilev_ip%b * ps_out(i,j)
            pb_out = ilev_op%a + ilev_op%b * ps_out(i,j)
            call rebin_levels(pb_in, intermediate_mix(i,j,:), .true., &
                pb_out, out_mix(i,j,:), status)
            if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
        end do
    end do

    ! Check whether the total mass has been conserved, sort of a diagnostic
    if (diag_msg) then
        total_mass_before = 0.0
        do i = 1, nlon_in
            do j = 1, nlat_in
                total_mass_before = total_mass_before + lli_ip%area_m2(j) * sum(in_mix(i,j,:) * pres_thick_in(i,j,:))
            end do
        end do
        total_mass_before = total_mass_before/(grav * amu * xmair)

        allocate(pres_out(nlon_out, nlat_out, nlev_out+1))
        allocate(pres_thick_out(nlon_out, nlat_out, nlev_out))
        do i = 1, nlev_out+1
            pres_out(:,:,i) = ilev_op%a(i-1) + ilev_op%b(i-1) * ps_out(:,:)
        end do
        do i = 1, nlev_out
            pres_thick_out(:,:,i) = pres_out(:,:,i) - pres_out(:,:,i+1)
        end do

        total_mass_after = 0.0
        do i = 1, nlon_out
            do j = 1, nlat_out
                total_mass_after = total_mass_after + lli_op%area_m2(j) * sum(out_mix(i,j,:) * pres_thick_out(i,j,:))
            end do
        end do
        total_mass_after = total_mass_after/(grav * amu * xmair)

        deallocate(pres_out, pres_thick_out)

        write(gol,'(a, " : total number of tracer molecules before regridding = ", es20.12)') rname, total_mass_before ; call goPr
        write(gol,'(a, " :  total number of tracer molecules after regridding = ", es20.12)') rname, total_mass_after ; call goPr
        write(gol,'(a, " :                   fractional tracer molecule error = ", es20.12)') rname, (total_mass_after - total_mass_before)/total_mass_before ; call goPr
    end if

    deallocate(ps_out, pres_in, pres_thick_in, pb_out, pb_in)

    status = 0

  end subroutine regrid_3d_mix

  ! =========================================


  subroutine Fill3D( lli, levi, nw, ps, field, &
                     lliX, leviX, fieldX, &
                     combkey, status )

    use grid_type_ll , only : TllGridInfo, FillGrid
    use grid_type_hyb, only : TLevelInfo, FillLevels

    ! --- in/out --------------------------------

    type(TllGridInfo), intent(in)        ::  lli
    type(TLevelInfo), intent(in)         ::  levi
    character(len=*), intent(in)         ::  nw
    real, intent(in)                     ::  ps(:,:)        ! Pa
    real, intent(out)                    ::  field(:,:,:)
    type(TllGridInfo), intent(in)        ::  lliX
    type(TLevelInfo), intent(in)         ::  leviX
    real, intent(in)                     ::  fieldX(:,:,:)
    character(len=*), intent(in)         ::  combkey
    integer, intent(out)                 ::  status

    ! --- const ---------------------------------

    character(len=*), parameter  ::  name = mname//'/Fill3D'

    ! --- local -------------------------

    real, allocatable     ::  field_ll(:,:,:)
    integer               ::  l

    ! --- begin -------------------------

    ! output horizontal grid, input levels
    allocate( field_ll(lli%nlon,lli%nlat,leviX%nlev) )

    select case ( combkey )

      !
      ! mass average
      !

      case ( 'mass-aver' )

        ! horizontal
        do l = 1, leviX%nlev
          call FillGrid( lli , 'n', field_ll(:,:,l), &
                         lliX, 'n', fieldX(:,:,l), 'area-aver', status )
          if (status<0) then
            write (*,'("ERROR - only part of target grid filled")')
            write (*,'("ERROR in ",a)') name; status=1; return
          end if
          if (status/=0) then; write (*,'("ERROR in ",a)') name; status=1; return; end if
        end do

        ! vertical
        call FillLevels( levi, nw, ps, field, leviX, field_ll, 'mass-aver', status )
        if (status/=0) then; write (*,'("ERROR in ",a)') name; status=1; return; end if

      !
      ! other (should be supported by FillGrid and FillLevels)
      !

      case default

        ! horizontal
        do l = 1, leviX%nlev
          call FillGrid( lli , 'n', field_ll(:,:,l), &
                         lliX, 'n', fieldX(:,:,l), combkey, status )
          if (status<0) then
            write (*,'("ERROR - only part of target grid filled")')
            write (*,'("ERROR in ",a)') name; status=1; return
          end if
          if (status/=0) then; write (*,'("ERROR in ",a)') name; status=1; return; end if
        end do

        ! vertical
        call FillLevels( levi, nw, ps, field, leviX, field_ll, combkey, status )
        if (status/=0) then; write (*,'("ERROR in ",a)') name; status=1; return; end if

    end select

    ! done
    deallocate( field_ll )

    ! ok
    status = 0

  end subroutine Fill3D



  ! **************************************************************


  !
  !   p = a + b sp
  !


  subroutine FillMass( m, lli, levi, sp, status )

    use Binas        , only : grav
    use grid_type_ll , only : TllGridInfo, AreaOper
    use grid_type_hyb, only : TLevelInfo

    ! --- begin ---------------------------------

    real, intent(out)                    ::  m(:,:,:)   ! kg
    type(TllGridInfo), intent(in)        ::  lli
    type(TLevelInfo), intent(in)         ::  levi
    real, intent(in)                     ::  sp(:,:)       ! Pa
    integer, intent(out)                 ::  status

    ! --- const ---------------------------------

    character(len=*), parameter  ::  rname = mname//'/FillMass'

    ! --- local -------------------------

    integer               ::  l

    ! --- begin -------------------------

    ! check shape of target grid:
    if ( (size(m,1) /= lli%nlon ) .or. (size(m,2) /= lli%nlat) .or. &
         (size(m,3) /= levi%nlev) ) then
      write (*,'("ERROR - target array does not match with grid definition:")')
      write (*,'("ERROR -   lli    : ",i3," x ",i3         )') lli%nlon, lli%nlat
      write (*,'("ERROR -   levi   : ",i3                  )') levi%nlev
      write (*,'("ERROR -   ll     : ",i3," x ",i3," x ",i3)') shape(m)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! Pa = kg g / A   ->  kg = A * Pa/g

    ! loop over levels
    do l = 1, levi%nlev
      m(:,:,l) = levi%da(l) + levi%db(l) * sp / grav    ! Pa/g = kg/m2
      call AreaOper( lli, m(:,:,l), '*', 'm2', status )          ! kg
      if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    end do

    ! ok
    status = 0

  end subroutine FillMass


  ! ***


  !
  ! p1 = a + b sp1
  ! p2 = a + b sp2
  !
  ! p2 - p1 = b (sp2 - sp1)
  !
  ! m = (p2 - p1)/g * A
  !


  subroutine FillMassChange( dm, lli, levi, sp1, sp2, status )

    use Binas        , only : grav
    use grid_type_ll , only : TllGridInfo, AreaOper
    use grid_type_hyb, only : TLevelInfo

    ! --- begin ---------------------------------

    real, intent(out)                    ::  dm(:,:,:)   ! kg
    type(TllGridInfo), intent(in)        ::  lli
    type(TLevelInfo), intent(in)         ::  levi
    real, intent(in)                     ::  sp1(:,:)       ! Pa
    real, intent(in)                     ::  sp2(:,:)       ! Pa
    integer, intent(out)                 ::  status

    ! --- const ---------------------------------

    character(len=*), parameter  ::  rname = mname//'/FillMassChange'

    ! --- local -------------------------

    integer               ::  l

    ! --- begin -------------------------

    ! check shape of target grid:
    if ( (size(dm,1) /= lli%nlon ) .or. (size(dm,2) /= lli%nlat) .or. &
         (size(dm,3) /= levi%nlev) ) then
      write (*,'("ERROR - target array does not match with grid definition:")')
      write (*,'("ERROR -   lli    : ",i3," x ",i3         )') lli%nlon, lli%nlat
      write (*,'("ERROR -   levi   : ",i3                  )') levi%nlev
      write (*,'("ERROR -   ll     : ",i3," x ",i3," x ",i3)') shape(dm)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! Pa = kg g / A   ->  kg = A * Pa/g

    ! loop over levels
    do l = 1, levi%nlev
      dm(:,:,l) = abs(levi%db(l)) * ( sp2 - sp1 ) / grav    ! Pa/g = kg/m2
      call AreaOper( lli, dm(:,:,l), '*', 'm2', status )          ! kg
      if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    end do


    ! ok
    status = 0

  end subroutine FillMassChange




end module grid_3d
