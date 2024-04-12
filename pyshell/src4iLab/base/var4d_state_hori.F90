!###############################################################################
!
! Storage for 4D-var state data for horizontal grids.
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

module Var4D_State_Hori

    use GO, only : gol, goPr, goErr

    use Var4D_Covariance_Types, only : Tcov_data
    use os_specs, only : MAX_FILENAME_LEN, MAX_RCKEY_LEN
    use chem_param, only : ntracet, names

    implicit none


    ! --- in/out -----------------------------------

    private

    public  ::  Var4D_State_Hori_Init
    public  ::  Var4D_State_Hori_Done

    public  ::  n_hor
    public  ::  vec2ll
    public  ::  outdir_state

    ! Background error correlation stuff
    integer                             ::  n_hor_cor               ! # different spatial corr matrices

    ! --- const ------------------------------------

    character(len=*), parameter   ::  mname = 'Var4D_State_Hori'

    ! maximum number of entries in horizontal correletations:
    integer, parameter  ::  max_hor_cor = 10

    ! --- types ------------------------------------

    ! Mapping of state vector space on lat-lon grid
    type Tvec2ll
        integer                    :: region    ! region
        integer                    :: i         ! longitude grid cell number
        integer                    :: j         ! latitude grid cell number
        real                       :: lon       ! longitude
        real                       :: lat       ! latitude
        logical                    :: enter_region  ! true if entering region
        logical                    :: leave_region  ! true if leaving region
        integer                    :: enter_region2 ! 1 if entering region
        integer                    :: leave_region2 ! 1 if leaving region
    end type Tvec2ll

    ! Store a region mask per category (optional), so that we can enforce complete correlation within a region
    ! and zero correlation between regions
    type Tregion_mask
        real(4), allocatable            :: reg_mask(:,:,:) ! nreg x nlat x nlon
        integer, allocatable            :: reg_index(:,:)  ! nlat x nlon
        real                            :: dlat, dlon      ! resolution of the mask
        logical                         :: use_it          ! should this category be clumped by region?
        character(len=MAX_FILENAME_LEN) :: mask_filename   ! which file contains the region masks?
    end type Tregion_mask

    ! --- var ------------------------------------

    ! number of grid points horizontally
    integer                :: n_hor

    ! mapping from state vec to latlon grid
    type(Tvec2ll), allocatable    :: vec2ll(:)     ! (n_hor)

    ! Storing optional region masks for clumping each category
    type(Tregion_mask), allocatable, target   :: reg_mask(:) ! dimension cat_all

    character(len=MAX_FILENAME_LEN)   :: outdir_state            ! outpur dir for state variables
    logical                           :: optimize_emission       ! whether to optimize emissions or not



contains


    ! ============================================================================


    subroutine Var4D_State_Hori_Init( status )

        use dims,                      only : nregions
        use dims,                      only : isr, ier, jsr, jer, parent, children
        use dims,                      only : xbeg, dx, xref, ybeg, dy, yref
        use dims,                      only : jbeg, jend, ibeg, iend, region_name
        use dims,                      only : im,jm,xend,yend,tref
        use global_data,               only : region_dat
        use file_netcdf
        use netcdf,                    only : nf90_put_att, NF90_NOERR, NF90_GLOBAL
        use GO                       , only : TrcFile, ReadRc
        use global_data,               only : rcF

        ! --- in/out ----------------------------------------------

        integer, intent(out)               :: status

        ! --- const -----------------------------------------------

        character(len=*), parameter        :: rname = mname//'/Var4D_State_Hori_Init'

        ! --- local -----------------------------------------------

        character(len=MAX_FILENAME_LEN)    :: outdir, fname, lli_name
        character(len=MAX_FILENAME_LEN)    :: correlation_dir, my_zoom_name
        integer, pointer                   :: zoomed(:,:)
        integer                            :: i, j
        integer                            :: region, region_old
        integer                            :: i_hor
        integer                            :: output_fid, group_id
        character(len=20)                  :: child_name

        integer                            :: cat_all, cat_region, i_cat, optim, i_catb, reg, nhor_from_file, itr
        character(len=MAX_RCKEY_LEN)       :: rc_key
        character(len=MAX_RCKEY_LEN)       :: cor_line
        character(len=20)                  :: cat_name, tcor, reg_name
        character(len=1)                   :: choice
        real                               :: corlen
        logical                            :: file_exist, calculate_Bh

        character(len=20),dimension(:), allocatable    :: catnam_all
        character(len=1), dimension(:), allocatable    :: choice_all
        real            , dimension(:), allocatable    :: corlen_all
        character(len=20),dimension(:), allocatable    :: tcorli_all
        integer         , dimension(:), allocatable    :: nregion_all
        integer         , dimension(:,:), allocatable  :: regions_all
        character(len=20),dimension(:), allocatable    :: regions_name
        character(len=10)                              :: reg_list_string

        ! --- begin -----------------------------------------------

        ! changes to allow for optimizing certain sources only in certain regions.
        ! first we have to determine the configuration of the sources:
        ! 1. a source that is optimized on all the grids, needs the old procedure
        ! 2. a source can be optimized ONLY in the parent, and not in the children
        ! 3. a source can be optimized in ONLY the child and not the parent.
        ! so, for each case, we need a different horizontal correlation matrix.
        ! To check:  do not allow for mixed correlation lengths (e.g. smaller in the
        ! zoom regions? Should be possible
        ! Also: we will identify a specific category by name (e.g. anthropogenic).
        ! The correspoding vec2ll (how to loop through regions) and matrices are
        ! stored in anthropogenic_region1_region2.nc

        call ReadRc( rcF, 'optimize.emission', optimize_emission, status )
        if (.not. optimize_emission) return

        ! first count the categories to be optimized:
        cat_all = 0
        do region = 1, nregions
            do itr = 1, ntracet
                rc_key = 'emission.'//trim(names(itr))//'.'//trim(region_name(region))//'.categories'
                call ReadRc( rcF, rc_key, cat_region, status )
                IF_NOTOK_RETURN(status=1)
                do i_cat = 1, cat_region
                    if (i_cat < 10) then
                        write(rc_key, '(a,i1)') 'emission.'//trim(names(itr))//'.'//trim(region_name(region))//'.category', i_cat
                    else
                        write(rc_key, '(a,i2)') 'emission.'//trim(names(itr))//'.'//trim(region_name(region))//'.category', i_cat
                    endif
                    call ReadRc( rcF, rc_key, cor_line, status )
                    IF_NOTOK_RETURN(status=1)
                    cor_line = cor_line(index(cor_line,';')+1:)
                    cor_line = cor_line(index(cor_line,';')+1:)
                    cor_line = cor_line(index(cor_line,';')+1:)
                    cor_line = cor_line(index(cor_line,';')+1:)
                    read(cor_line,*) optim
                    if (optim == 1) cat_all = cat_all + 1
                end do ! i_cat
            end do ! itr
        end do ! region

        if (cat_all > 0) then
            allocate(catnam_all(cat_all))
            allocate(choice_all(cat_all))
            allocate(corlen_all(cat_all))
            allocate(tcorli_all(cat_all))
            allocate(nregion_all(cat_all))
            nregion_all(:) = 0
            allocate(regions_all(cat_all,nregions))
            allocate(regions_name(cat_all))
        else
            write(gol, *) '========================================================================================='  ; call goPr
            write(gol, *)  'WARNING :: You are not optimising any category for any tracer' ; call goPr
            write(gol, *) '========================================================================================='  ; call goPr
            return
        end if
        ! second: collect info about emission categories
        ! a unique caterogie is defined by name, time-correlation. The horizontal
        ! correlation length might differ, but this has not yet been implemented!

        cat_all = 0    ! number of unique categories to be optimized:

        do region = 1, nregions
            do itr = 1, ntracet
                rc_key = 'emission.'//trim(names(itr))//'.'//trim(region_name(region))//'.categories'
                call ReadRc( rcF, rc_key, cat_region, status )
                IF_NOTOK_RETURN(status=1)
                cats: do i_cat = 1, cat_region
                    write(rc_key, '(a,i1)') 'emission.'//trim(names(itr))//'.'//trim(region_name(region))//'.category',i_cat
                    call ReadRc( rcF, rc_key, cor_line, status )
                    IF_NOTOK_RETURN(status=1)
                    reg_name = region_name(region)
                    cat_name = cor_line(:index(cor_line,';')-1)
                    cor_line = cor_line(index(cor_line,';')+1:)
                    cor_line = cor_line(index(cor_line,';')+1:)
                    read(cor_line,'(f7.1,1x,a1)') corlen,choice
                    cor_line = cor_line(index(cor_line,';')+1:)
                    tcor = cor_line(:index(cor_line,';')-1)
                    cor_line = cor_line(index(cor_line,';')+1:)
                    read(cor_line,*) optim
                    if (optim == 1) then
                        do i_catb = 1, cat_all    ! check if idential to earlier category:
                            if ((trim(cat_name) == trim(catnam_all(i_catb))) .and. &   !same name
                                    (trim(choice  ) == trim(choice_all(i_catb))) .and. &   !
                                    (abs(corlen - corlen_all(i_catb)) < 1e-5)    .and. &
                                    (trim(tcor)     == trim(tcorli_all(i_catb)))) then
                                if (reg_name /= regions_name(i_catb)) then
                                    nregion_all(i_catb) = nregion_all(i_catb) + 1
                                    regions_all(i_catb,nregion_all(i_catb)) = region
                                end if
                                cycle cats
                            end if
                        end do

                        cat_all = cat_all + 1
                        catnam_all(cat_all) = cat_name
                        choice_all(cat_all) = choice
                        corlen_all(cat_all) = corlen
                        tcorli_all(cat_all) = trim(tcor)
                        nregion_all(cat_all) = 1
                        regions_all(cat_all,1) = region
                        regions_name(cat_all) = reg_name
                    end if
                end do cats
            end do ! itr
        end do ! region

        write(gol, *) '========================================================================================='  ; call goPr
        write(gol, *) 'Categories marked for optimization:'  ; call goPr
        write(gol, *) '========================================================================================='  ; call goPr
        do i_cat = 1, cat_all
            write(reg_list_string, '(i1)') regions_all(i_cat,1)
            do i = 2, nregion_all(i_cat)
                write(reg_list_string, '(a,i2)') trim(reg_list_string), regions_all(i_cat,i)
            end do
            write(gol, '(a24, " nregion: ", i1, " regions: ", a, " corlen: ", f8.2, " choice: ", a1, " time: ", a)') &
                    catnam_all(i_cat), nregion_all(i_cat), trim(reg_list_string), corlen_all(i_cat), &
                    choice_all(i_cat), trim(tcorli_all(i_cat))
            call goPr
        end do
        write(gol, *) '========================================================================================='  ; call goPr


        call ReadRc( rcF, 'correlation.inputdir', correlation_dir, status )
        IF_NOTOK_RETURN(status=1)
        ! now loop over cat_all
        allocate(reg_mask(cat_all))
        reg_mask(:)%use_it = .false.

        do i_cat = 1, cat_all
            n_hor = 0
            do reg = 1, nregion_all(i_cat)
                region = regions_all(i_cat,reg)
                zoomed => region_dat(region)%zoomed
                do j = jsr(region), jer(region)
                    do i = isr(region), ier(region)
                        if( zoomed(i,j) /= region ) cycle
                        n_hor = n_hor + 1
                    enddo
                enddo
                nullify(zoomed)
            enddo


            allocate( vec2ll(n_hor) )
            i_hor = 0
            region_old = 0
            do reg = 1, nregion_all(i_cat)
                region = regions_all(i_cat,reg)
                zoomed => region_dat(region)%zoomed
                do j = jsr(region), jer(region)
                    do i = isr(region), ier(region)
                        if ( zoomed(i,j) /= region ) cycle
                        i_hor = i_hor + 1
                        vec2ll(i_hor)%region = region
                        vec2ll(i_hor)%i = i
                        vec2ll(i_hor)%lon = xbeg(region) + (i-0.5)*dx/xref(region)
                        vec2ll(i_hor)%j = j
                        vec2ll(i_hor)%lat = ybeg(region) + (j-0.5)*dy/yref(region)
                        vec2ll(i_hor)%enter_region = .false.
                        vec2ll(i_hor)%leave_region = .false.
                        vec2ll(i_hor)%enter_region2 = 0
                        vec2ll(i_hor)%leave_region2 = 0
                        if ( region /= region_old ) then
                            vec2ll(i_hor)%enter_region = .true.
                            vec2ll(i_hor)%enter_region2 = 1
                            if ( i_hor > 1 ) then
                                vec2ll(i_hor-1)%leave_region = .true.
                                vec2ll(i_hor-1)%leave_region2 = 1
                            end if
                        end if
                        region_old = region
                    end do
                end do
                nullify(zoomed)
            end do
            vec2ll(i_hor)%leave_region = .true.
            vec2ll(i_hor)%leave_region2 = 1

            ! Print dimensions
            write (gol,'("category, name, n_hor    ",i3, 1x, a20, i6)') i_cat,catnam_all(i_cat), n_hor; call goPr

            ! write to file: filename =
            ! Bh:my.zoom:region1-region2_corlen_choice.nc:

            call ReadRc( rcF, 'my.zoom', my_zoom_name, status )
            IF_NOTOK_RETURN(status=1)
            lli_name = region_name(regions_all(i_cat,1))
            do region = 2, nregion_all(i_cat)
                write(lli_name,'(a,"_",a)') trim(lli_name), trim(region_name(regions_all(i_cat,region)))
            end do

            ! write vec2ll to a file, which should be specific to this combination of regions
            call ReadRc( rcF, 'my.run.dir', outdir, status )
            write(fname,"(a,a,a,a)") trim(outdir),'/vec2ll_',trim(lli_name),'.nc' ! Why is the vec2ll file written for each category to be optimized?
            ! Write only if the file does not exist
            inquire(file=fname, exist=file_exist)
            if (.not. file_exist) then
                output_fid = nc_open(fname, 'c', status)
                IF_NOTOK_RETURN(status=1)

                call nc_create_dim(output_fid, 'n_hor', n_hor)
                call nc_dump_var(output_fid, 'vec2ll_region', (/'n_hor'/), vec2ll(:)%region, (/'long_name'/), (/'Per grid box, the region is given'/))
                call nc_dump_var(output_fid, 'vec2ll_lon',    (/'n_hor'/), vec2ll(:)%lon, (/'long_name'/),    (/'Per grid box, the longitude of the center is given'/))
                call nc_dump_var(output_fid, 'vec2ll_i',      (/'n_hor'/), vec2ll(:)%i, (/'long_name'/),      (/'Per grid box, the i-box number within that region is given'/))
                call nc_dump_var(output_fid, 'vec2ll_lat',    (/'n_hor'/), vec2ll(:)%lat, (/'long_name'/),    (/'Per grid box, the latitude of the center is given'/))
                call nc_dump_var(output_fid, 'vec2ll_j',      (/'n_hor'/), vec2ll(:)%j, (/'long_name'/),      (/'Per grid box, the j-box number within that region is given'/))
                call nc_dump_var(output_fid, 'vec2ll_enter_region', (/'n_hor'/), vec2ll(:)%enter_region2,  (/'long_name'/), (/'Per grid box, the logical value for enter region is given'/))
                call nc_dump_var(output_fid, 'vec2ll_leave_region', (/'n_hor'/), vec2ll(:)%leave_region2, (/'long_name'/), (/'Per grid box, the logical value for leave region is given'/))
                call nc_close(output_fid)
            end if
            deallocate(vec2ll)
        end do ! i_cat

        ! Write region_dat(region)%zoomed to file
        call ReadRc( rcF, 'my.run.dir', outdir, status )
        IF_NOTOK_RETURN(status=1)
        write(fname,'(a,a)') trim(outdir),'/Zoomed.nc4'
        inquire(file=fname, exist=file_exist)
        if (.not. file_exist) then
            output_fid = nc_open(fname, 'c', status)
            IF_NOTOK_RETURN(status=1)

            call nc_create_dim(output_fid, 'nregions', nregions)
            call nc_create_dim(output_fid, 'nregions_plus', nregions+1)
            ! Write region begin and end coordinates.

            call nc_dump_var(output_fid, 'xbeg', (/'nregions'/), xbeg(:), (/'long_name'/), &
                    (/'List with most western longitude for the regions'/))
            call nc_dump_var(output_fid, 'xend', (/'nregions'/), xend(:), (/'long_name'/),  &
                    (/'List with most eastern longitude for the regions'/))
            call nc_dump_var(output_fid, 'ybeg', (/'nregions'/), ybeg(:), (/'long_name'/),  &
                    (/'List with most southern latitude for the regions'/))
            call nc_dump_var(output_fid, 'yend', (/'nregions'/), yend(:), (/'long_name'/),  &
                    (/'List with most northern latitude for the regions'/))
            ! Write refinement factors
            print*,'size(xref), value(xref): ',size(xref), xref(:), xref(0), xref(1)
            call nc_dump_var(output_fid, 'xref', (/'nregions_plus'/), xref(:), (/'long_name'/),&
                    (/'List with refinement factors in longitude direction'/))
            call nc_dump_var(output_fid, 'yref', (/'nregions_plus'/), yref(:), (/'long_name'/),&
                    (/'List with refinement factors in latitude direction'/))
            call nc_dump_var(output_fid, 'tref', (/'nregions_plus'/), tref(:), (/'long_name'/),&
                    (/'List with refinement factors in time'/))
            ! Write number of boxes per region
            call nc_dump_var(output_fid, 'im', (/'nregions'/), im(:), (/'long_name'/),&
                    (/'List with number of boxes in longitude direction'/))
            call nc_dump_var(output_fid, 'jm', (/'nregions'/), jm(:), (/'long_name'/),&
                    (/'List with number of boxes in latitude direction'/))
            ! Write start and end box per region
            call nc_dump_var(output_fid, 'isr', (/'nregions'/), isr(:), (/'long_name'/),&
                    (/'List with start box in longitude direction'/))
            call nc_dump_var(output_fid, 'jsr', (/'nregions'/), jsr(:), (/'long_name'/),&
                    (/'List with start box in latitude direction'/))
            call nc_dump_var(output_fid, 'ier', (/'nregions'/), ier(:), (/'long_name'/),&
                    (/'List with end box in longitude direction'/))
            call nc_dump_var(output_fid, 'jer', (/'nregions'/), jer(:), (/'long_name'/),&
                    (/'List with end box in latitude direction'/))
            ! Write start and end box relative to parent
            call nc_dump_var(output_fid, 'ibeg', (/'nregions'/), ibeg(:), (/'long_name'/),&
                    (/'List with start box in longitude direction relative to parent'/))
            call nc_dump_var(output_fid, 'iend', (/'nregions'/), iend(:), (/'long_name'/),&
                    (/'List with end box in longitude direction relative to parent'/))
            call nc_dump_var(output_fid, 'jbeg', (/'nregions'/), jbeg(:), (/'long_name'/),&
                    (/'List with start box in latitude direction relative to parent'/))
            call nc_dump_var(output_fid, 'jend', (/'nregions'/), jend(:), (/'long_name'/),&
                    (/'List with end box in latitude direction relative to parent'/))
            ! Write parent of region: glb6x4 has no parent => write 0
            call nc_dump_var(output_fid, 'parents', (/'nregions'/), parent(:), (/'long_name'/), (/'List with parents'/))
            ! Write children
            do region = 1, nregions
                write(child_name,"('children_',i2.2)") region
                call nc_dump_var(output_fid, trim(child_name), (/'nregions'/), children(region,1:), (/'long_name'/), (/'List with children per region'/))
            end do

            status = nf90_put_att(output_fid, NF90_GLOBAL, 'dx', dx)
            status = nf90_put_att(output_fid, NF90_GLOBAL, 'dy', dy)
            status = nf90_put_att(output_fid, NF90_GLOBAL, 'zoom_def', trim(my_zoom_name))

            do region = 1, nregions
                print*,'region = ',region
                print*,'shape(region)', size(region_dat(region)%zoomed,dim=1)
                print*,'shape(region)', size(region_dat(region)%zoomed,dim=2)
                group_id = nc_create_group(output_fid, region_name(region))
                call nc_create_dim(group_id, 'lon', size(region_dat(region)%zoomed,dim=1))
                call nc_create_dim(group_id, 'lat', size(region_dat(region)%zoomed,dim=2))
                call nc_dump_var(group_id, 'zoomed', (/'lon   ', 'lat   '/), region_dat(region)%zoomed(:,:), (/'long_name'/), (/'Array zoomed for this region'/))
                print*,'shape(edge)', size(region_dat(region)%edge,dim=1)
                print*,'shape(edge)', size(region_dat(region)%edge,dim=2)
                call nc_create_dim(group_id, 'lon_edge', size(region_dat(region)%edge,dim=1))
                call nc_create_dim(group_id, 'lat_edge', size(region_dat(region)%edge,dim=2))
                call nc_dump_var(group_id, 'edge', (/'lon_edge   ', 'lat_edge   '/), region_dat(region)%edge(:,:), (/'long_name'/), (/'Array edge for this region'/))
            end do
            call nc_close(output_fid)
        end if ! Zoomed.nc4 does not exist


        deallocate(catnam_all)
        deallocate(choice_all)
        deallocate(corlen_all)
        deallocate(tcorli_all)
        deallocate(regions_all)
        deallocate(nregion_all)
        deallocate(regions_name)

        ! ok:
        status = 0

    end subroutine Var4D_State_Hori_Init


    ! ***


    subroutine Var4D_State_Hori_Done( status )


        ! --- in/out ----------------------------------------------

        integer, intent(out)               :: status

        ! --- const -----------------------------------------------

        character(len=*), parameter        :: rname = mname//'/Var4D_State_Hori_Done'

        ! --- local -----------------------------------------------

        integer         ::  i

        ! --- begin -----------------------------------------------

        ! mapping of latlon grid on vector
        !   done already in new frame_work: deallocate( vec2ll )

        if (allocated(reg_mask)) deallocate(reg_mask)

        ! ok:
        status = 0

    end subroutine Var4D_State_Hori_Done

    ! We might want to specify the correlation in terms of regions. E.g., we might want to divide the world into
    ! ecoregions, and enforce perfect correlation within each ecoregion and perfect decorrelation between regions.
    ! This will ensure that the adjustments to the flux within one region strictly follows the standard deviation
    ! distribution within that region.
    ! This needs an input file with the regions. In the file, there will be a variable called 'mask', with shape
    ! nregions x nlat x nlon.

    subroutine read_region_mask(i_cat, status)

        ! --- modules -----------------------------------------------
        use global_data,            only : rcf
        use Go,                     only : ReadRc
        use file_netcdf

        ! --- in/out -----------------------------------------------
        integer, intent(in)                 :: i_cat
        integer, intent(out)                :: status

        ! --- const -----------------------------------------------
        character(len=*), parameter         :: rname = mname//'/read_region_mask'

        ! --- local -----------------------------------------------
        integer                             :: nc_id, nlat, nlon, nreg, i, j
        real(4), allocatable                :: check_array(:,:)

        ! --- start -----------------------------------------------
        nc_id = nc_open(reg_mask(i_cat)%mask_filename, 'r', status)
        IF_NOTOK_RETURN(status=1)

        reg_mask(i_cat)%reg_mask = nc_read_var(nc_id, 'mask') ! Thanks to Fortran, will be nlon x nlat x nregions
        call nc_close(nc_id)

        nreg = size(reg_mask(i_cat)%reg_mask, 3)
        nlat = size(reg_mask(i_cat)%reg_mask, 2)
        nlon = size(reg_mask(i_cat)%reg_mask, 1)
        write(*,'(a, " has ", i3, " regions and is ", i4, " x ", i4)') trim(reg_mask(i_cat)%mask_filename), nreg, nlat, nlon

        reg_mask(i_cat)%dlat = 180.0/nlat
        reg_mask(i_cat)%dlon = 360.0/nlon

        ! the array 'reg_index' is an integer array which has a different integer for every region mask
        ! first, check if a box is covered by more than one region mask
        allocate(check_array(nlon,nlat))
        check_array = sum(reg_mask(i_cat)%reg_mask, 3)

        if (any(check_array .ge. 1.1)) then
            do i = 1, nlon
                do j = 1, nlat
                    if (check_array(i,j) .ge. 1.1) then
                        write(*,'("Location ", i4, ", ", i4, " is covered by more than one mask")') j, i
                        print *, reg_mask(i_cat)%reg_mask(i,j,:)
                    end if
                end do
            end do
        end if

        if (any(check_array .le. 0.9)) then
            do i = 1, nlon
                do j = 1, nlat
                    if (check_array(i,j) .le. 0.9) then
                        write(*,'("Location ", i4, ", ", i4, " is not covered by any mask")') j, i
                        print *, reg_mask(i_cat)%reg_mask(i,j,:)
                    end if
                end do
            end do
        end if

        deallocate(check_array)

        allocate(reg_mask(i_cat)%reg_index(nlon, nlat))
        reg_mask(i_cat)%reg_index = 0
        do i = 1, nreg
            reg_mask(i_cat)%reg_index = reg_mask(i_cat)%reg_index + i * int(reg_mask(i_cat)%reg_mask(:,:,i))
        end do

        status = 0

    end subroutine read_region_mask


end module Var4D_State_Hori
