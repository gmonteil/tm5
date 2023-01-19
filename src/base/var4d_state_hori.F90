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
!    if (.not. optimize_emission) return

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
       write(fname,'(a,"/Bh:",a,":",a,"_",i5.5,"_",a1".nc")') trim(correlation_dir), trim(my_zoom_name), &
           trim(lli_name), nint( corlen_all(i_cat) ), choice_all(i_cat)

       inquire( file=fname, exist=file_exist )
       ! We only need to calculate the horizontal correlation matrices if we want to do 4DVAR
       calculate_Bh = .false.
       if (optimize_emission) then
        if ( .not. file_exist) then
            calculate_Bh = .true.
        else
            ! check if file has the correct n_hor, and if not, delete file and recalculate
            output_fid = nc_open(fname, 'r', status)
            IF_NOTOK_RETURN(status=1)

            nhor_from_file = nc_get_dim(output_fid, 'n_hor')
            call nc_close(output_fid)
            if (nhor_from_file == n_hor) then
                calculate_Bh = .false.
            else
                calculate_Bh = .true.
            end if
        end if ! Bh_* file does not exist
       end if ! optimize emissions

       if (calculate_Bh) then
          ! Calculate and write:
          write(gol,*) 'Calculating and writing horizontal correlation file:' ; call goPr
          write(gol,*) trim(fname) ; call goPr
          if (choice_all(i_cat) == 'r' .or. choice_all(i_cat) == 'h') then
            ! Clumping by region needed
            reg_mask(i_cat)%use_it = .true.
            ! The key specifying the region file should be, e.g., fossil fuel.region.filename
            rc_key = trim(catnam_all(i_cat))//'.region.filename'
            call ReadRc(rcF, rc_key, reg_mask(i_cat)%mask_filename, status)
            IF_NOTOK_RETURN(status=1)
            call read_region_mask(i_cat, status)
            IF_NOTOK_RETURN(status=1)
          end if
          call calc_latlon_covariance( fname, corlen_all(i_cat), choice_all(i_cat), i_cat, status )
          ! to be done: also write vec2ll to this file!
          IF_NOTOK_RETURN(status=1)
       end if

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

  ! ***


  !===========================================================================================================


  !-----------------------------------------------------------------------
  ! Calculate spatial correlation matrix C on zoomed grid
  !
  ! Perform eigen decomposition, and return eigenvectors and eigenvalues
  !
  ! Input:
  !  hc              correlation matrix
  !                    o hc%n            dimension of C (number of grid points)
  !                    o hc%corlen       correlation length in km
  !                        If zero, no correlations are assumed (P = I, D = I)
  !                    o hc%choice       type of correlation
  !                                         'g'  gaussian
  !                                         'e'  exponential
  !
  ! Output:
  !  hc              correlation matrix
  !                    o hc%sqrt_lam     contains eigenvalues of correlation matrix
  !                    o hc%P            contains the eigenvectors
  !                                         such that:
  !                          C = P D P^-1 , with P^-1 = P' and D_ii = sqrt_lam_i^2
  ! slightly restructured by MK, June, 2011
  !-----------------------------------------------------------------------

  subroutine calc_latlon_covariance( fname, corlen, choice, icat, status )

    ! --- modules ---------------------------------------------

    use misctools,                 only : dist, check_dir!, eigen_symm
    use linalg_interface,          only : eigvals, matmul_fast
    use dims,                      only : im,jm,lm,dx,dy,xref,yref,tref
    use dims,                      only : xbeg,ybeg,xend,yend,nregions
    use file_netcdf
    use netcdf,                    only : nf90_put_att, NF90_NOERR, NF90_GLOBAL
    use global_data,               only : rcF

    ! --- in/out ----------------------------------------------

    integer, intent(out)               :: status
    integer, intent(in)                :: icat
    character(len=1), intent(in)       :: choice
    real, intent(in)                   :: corlen
    character(len=*),intent(in)        :: fname

    ! --- const -----------------------------------------------

    character(len=*), parameter        :: rname = mname//'/calc_latlon_covariance'

    ! --- local -----------------------------------------------

    integer                            :: i_hor1, i_hor2, j, i, l
    real                               :: lon1, lat1, lon2, lat2, dst, cor, eval
    real, dimension(:,:), allocatable  :: P,B,P_diag
    real, dimension(:), allocatable    :: lam,lam_sqrt
    integer                            :: n, mat_size
    integer                            :: iexp, nc_id
    logical                            :: comp_save, diagonalized
    integer, pointer                   :: reg_index(:,:)
    integer                            :: lat_idx_1, lat_idx_2, lon_idx_1, lon_idx_2
    integer, allocatable               :: region_num(:)

    ! --- begin -----------------------------------------------
    n = n_hor
    ! Pointers
    allocate(P(n_hor,n_hor), P_diag(n_hor,n_hor))
    allocate(lam(n_hor), lam_sqrt(n_hor))

    write(*,'(a,": corlen = ",f7.1," ; choice = ",a1)') rname, corlen, choice
    write(*,'(a,a)') "calculating covariance matrix (may take a while)", fname

    ! Function for correlations
    select case ( choice )
    case ( 'e' )
       iexp = 1
    case ( 'g' )
       iexp = 2
    case ( 'r' )
        iexp = 0 ! Who cares!
    case ( 'h' )
        iexp = 0 ! Again, who cares!
    case default
       write (*,'("ERROR - invalid choice for spatial correlation: ",a1)') choice
       status = 1
       IF_NOTOK_RETURN(status=1)
    end select

    if ( corlen < 0 ) then

       write (*,'("ERROR - correlation length should be >= 0")')
       TRACEBACK; return

    end if

    ! Need to handle the correlation choice 'r' differently
    diagonalized = .false. ! P has not been diagonalized yet
    select case (choice)
    case ('r')
        !
        ! In this case, if all the i_hors belonging to the same region appeared consecutively, then P would be block diagonal
        ! with as many blocks as number of regions, and each block would be a matrix of ones. An MxM matrix of ones has
        ! eigenvalue M with an eigenvector of all ones, and eigenvalue 0 repeated M-1 times. However, the ones are dispersed
        ! around instead of being clustered in block diagonals. Consider we have only 2 regions, and i_hor runs to 5, with
        ! regions [2 1 1 2 1]. The correlation matrix is, in that case,
        !
        ! 1 0 0 1 0
        ! 0 1 1 0 1
        ! 0 1 1 0 1                                                  (1)
        ! 1 0 0 1 0
        ! 0 1 1 0 1
        !
        ! Consider the array [2 3 5 1 4]. This is the redistribution array to rearrange the i_hors so that cells belonging to
        ! the same region are together, i.e., [2 1 1 2 1][2 3 5 1 4] = [1 1 1 2 2]. The inverse permutation is [4 1 2 5 3], i.e.,
        ! [1 1 1 2 2][4 1 2 5 3] = [2 1 1 2 1], the original array (inverse permutation from
        ! http://stackoverflow.com/questions/2483696/undo-or-reverse-argsort-python). If we did cluster the i_hors in the same
        ! region together, the correlation matrix would look like
        !
        ! 1 1 1 0 0
        ! 1 1 1 0 0
        ! 1 1 1 0 0                                                  (2)
        ! 0 0 0 1 1
        ! 0 0 0 1 1
        !
        ! Say [v1 v2 v3] is an eigenvector of the 3x2 ones-matrix with eigenvalue e. Then [v1 v2 v3 0 0] is an eigenvector of
        ! matrix (2) with the same eigenvalue e. If we now do the permutation [v1 v2 v3 0 0][4 1 2 5 3] = [0 v1 v2 0 v3], then
        ! the resultant vector is an eigenvector of matrix (1) with the same eigenvalue e. So to find the eigenvector
        ! decomposition in this case,
        !
        ! 1. Going from 1 to n_hor, construct an array of the region number each i_hor belongs to, e.g. [2 1 1 2 1].
        ! 2. Find the permutation which will put all cells in the same region together, e.g., [2 3 5 1 4]
        ! 3. Find the inverse permutation of the above permutation, e.g., [4 1 2 5 3].
        ! 4. Find the number of regions, and cells per region, e.g., region 1 has M1 cells, region 2 has M2 cells...
        ! 5. For each region i, diagonalize an Mi x Mi array of ones. If the eigenvector is [v_1 ... v_Mi], construct
        !    the full eigenvector of the (2)-like matrix by padding it with zeros on both sides. E.g., for region 2,
        !    the full eigenvector will be [0 0 ... (M1 times) v_1 ... v_M2 0 0 ...].
        ! 6. Permute the full eigenvector with the permutation found in step 3 to get the eigenvector of the (1)-like array.
        !
        reg_index => reg_mask(icat)%reg_index
        allocate(region_num(n))
        do i = 1, n
            lon1 = vec2ll(i)%lon
            lat1 = vec2ll(i)%lat
            lat_idx_1 = int((lat1 + 90.0)/reg_mask(icat)%dlat + 0.999999)
            lon_idx_1 = int((lon1 + 180.0)/reg_mask(icat)%dlon + 0.999999)
            region_num(i) = reg_index(lon_idx_1, lat_idx_1)
        end do
        call block_correlation_eig(region_num, P, lam)
        ! P now has the eigenvectors and lam the eigenvalues
        !deallocate(region_num)
        lam_sqrt = sqrt(lam)
        P_diag = 0.0
        do i=1,n_hor
            P_diag(i,i) = lam(i)
        end do
        diagonalized = .true.

    case ('h')
        reg_index => reg_mask(icat)%reg_index
        allocate(region_num(n))
        do i = 1, n
            lon1 = vec2ll(i)%lon
            lat1 = vec2ll(i)%lat
            lat_idx_1 = int((lat1 + 90.0)/reg_mask(icat)%dlat + 0.999999)
            lon_idx_1 = int((lon1 + 180.0)/reg_mask(icat)%dlon + 0.999999)
            region_num(i) = reg_index(lon_idx_1, lat_idx_1)
        end do
        call block_correlation_hybrid(region_num, P, lam, corlen, 1)
        ! P now has the eigenvectors and lam the eigenvalues
        !deallocate(region_num)
        lam_sqrt = sqrt(lam)
        P_diag = 0.0
        do i=1,n_hor
            P_diag(i,i) = lam(i)
        end do
        diagonalized = .true.

    !case ('h')
        !! This is a hybrid choice, a combination of 'r' and 'e'. The idea is that there are regions defined, within
        !! which the correlation falls off exponentially, whereas across regions the correlation is zero, which
        !! avoids spurious long-range correlations.
        !reg_index => reg_mask(icat)%reg_index
        !! loop first index over latlon grid
        !do i_hor1 = 1, n
            !lon1 = vec2ll(i_hor1)%lon
            !lat1 = vec2ll(i_hor1)%lat
            !lat_idx_1 = int((lat1 + 90.0)/reg_mask(icat)%dlat + 0.999999)
            !lon_idx_1 = int((lon1 + 180.0)/reg_mask(icat)%dlon + 0.999999)
            !! loop second index over latlon grid
            !do i_hor2 = i_hor1, n
                !lon2 = vec2ll(i_hor2)%lon
                !lat2 = vec2ll(i_hor2)%lat
                !lat_idx_2 = int((lat2 + 90.0)/reg_mask(icat)%dlat + 0.999999)
                !lon_idx_2 = int((lon2 + 180.0)/reg_mask(icat)%dlon + 0.999999)
                !if (reg_index(lon_idx_1, lat_idx_1) == reg_index(lon_idx_2, lat_idx_2)) then ! in the same region
                    !call dist( lon1, lat1, lon2, lat2, dst )
                    !cor = exp( -dst/corlen )
                    !P(i_hor1,i_hor2) = cor
                    !P(i_hor2,i_hor1) = cor
                !else ! in different regions, so no correlation
                    !P(i_hor1,i_hor2) = 0.0
                    !P(i_hor2,i_hor1) = 0.0
                !end if
            !end do ! i_hor2
        !end do ! i_hor1

    case ('e', 'g')
        ! put stuff here to construct P for exponential or gaussian decay
        if ( corlen == 0.0 ) then
            P = 0.
            P_diag = 0.0
            do i_hor1 = 1, n
                P(i_hor1,i_hor1) = 1.
                P_diag(i_hor1,i_hor1) = 1.
            end do
            lam = 1.
            lam_sqrt = 1. ! Sean Crowell pointed out this is needed here, otherwise lam_sqrt stays zero and 4DVAR does not work
            diagonalized = .true.
        else
            ! loop first index over latlon grid
            do i_hor1 = 1, n
                lon1 = vec2ll(i_hor1)%lon
                lat1 = vec2ll(i_hor1)%lat
                ! loop second index over latlon grid
                do i_hor2 = i_hor1, n
                    lon2 = vec2ll(i_hor2)%lon
                    lat2 = vec2ll(i_hor2)%lat
                    call dist( lon1, lat1, lon2, lat2, dst )
                    cor = exp( -(dst/corlen)**iexp )
                    P(i_hor1,i_hor2) = cor
                    P(i_hor2,i_hor1) = cor
                end do
            end do
        end if ! corlen == 0.0
    end select ! choice

    if (.not. diagonalized) then
        ! Eigen decomposition of symmetric matrix
        mat_size = size(P, 1)
        call eigvals(mat_size, P, lam, .true.)
        call make_positive_semidef(lam)

        lam_sqrt = sqrt( lam )   ! Return square root
        P_diag = 0.0
        do i=1,n_hor
            P_diag(i,i) = lam(i)
        end do
        diagonalized = .true.
    end if

    ! write to file for later use:
    ! turn on compression for the variables
    comp_save = nc_variables_deflate
    nc_variables_deflate = .true.
    write(*,'("Writing covariance file to disk: ", a)') trim(fname)

    ! If the folder with the Bh files does not exist, create it
    call check_dir(fname)

    nc_id = nc_open(fname, 'c', status)
    IF_NOTOK_RETURN(status=1)

    call nc_create_dim(nc_id, 'n_hor', n_hor)
    ! write all the attributes
    status = nf90_put_att(nc_id, NF90_GLOBAL, 'nregions', nregions)
    status = nf90_put_att(nc_id, NF90_GLOBAL, 'im', im(1:nregions))
    status = nf90_put_att(nc_id, NF90_GLOBAL, 'jm', jm(1:nregions))
    status = nf90_put_att(nc_id, NF90_GLOBAL, 'lm', lm(1:nregions))
    status = nf90_put_att(nc_id, NF90_GLOBAL, 'dx', dx)
    status = nf90_put_att(nc_id, NF90_GLOBAL, 'dy', dy)
    status = nf90_put_att(nc_id, NF90_GLOBAL, 'xref', xref(0:nregions))
    status = nf90_put_att(nc_id, NF90_GLOBAL, 'yref', yref(0:nregions))
    status = nf90_put_att(nc_id, NF90_GLOBAL, 'tref', tref(0:nregions))
    status = nf90_put_att(nc_id, NF90_GLOBAL, 'xbeg', xbeg(1:nregions))
    status = nf90_put_att(nc_id, NF90_GLOBAL, 'ybeg', ybeg(1:nregions))
    status = nf90_put_att(nc_id, NF90_GLOBAL, 'xend', xend(1:nregions))
    status = nf90_put_att(nc_id, NF90_GLOBAL, 'yend', yend(1:nregions))
    status = nf90_put_att(nc_id, NF90_GLOBAL, 'corlen', corlen)
    status = nf90_put_att(nc_id, NF90_GLOBAL, 'corchoice', choice)
    if (choice == 'r' .or. choice == 'h') then
        ! write out the name of the region file as an attribute
        status = nf90_put_att(nc_id, NF90_GLOBAL, 'region_mask_file', trim(reg_mask(icat)%mask_filename))
        ! Also write out region_num for debugging
        call nc_dump_var(nc_id, 'region_num', (/'n_hor'/), region_num(:), (/'long_name'/), (/'region number each cell belongs to'/))
        deallocate(region_num)
    end if

    call nc_dump_var(nc_id, 'sqrt_lam', (/ 'n_hor' /), lam_sqrt(:), (/'long_name'/), (/'square root of the eigenvalues'/))
    call nc_dump_var(nc_id, 'lam', (/ 'n_hor' /), lam(:), (/'long_name'/), (/'eigenvalues of the spatial correlation matrix'/))
    call nc_dump_var(nc_id, 'P', (/'n_hor', 'n_hor'/), P(:,:), (/'long_name'/), (/'eigenvectors'/))

    ! JFM: also write B itself (first recalculate), for use in postprocessing
    write(*,'(a,": recalculating B matrix from eigenvectors/eigenvalues")') rname
    ! P_diag is the diagonalized P, and P itself has the eigenvectors
    ! P (reconstituted) = P * P_diag * P.T
    allocate(B(n,n))
    call matmul_fast(P, P_diag, B) ! B <- P * P_diag
    ! since P_diag will not be used any more, use it to hold P.T
    P_diag = transpose(P)
    call matmul_fast(B, P_diag, P) ! P <- B * P.T

    call nc_dump_var(nc_id, 'B', (/'n_hor', 'n_hor'/), P(:,:), (/'long_name'/), (/'Horizontal correlation matrix'/))

    deallocate( B )

    call nc_close(nc_id)
    ! restore global compression setting
    nc_variables_deflate = comp_save

    deallocate ( lam, lam_sqrt )
    deallocate ( P, P_diag )

    status = 0

  end subroutine calc_latlon_covariance

  subroutine block_correlation_hybrid(idx, evecs, evals, corlen, iexp)
    ! modules
    use sorting,            only : argsort
    use orderpack,          only : multiplicity
    use linalg_interface,   only : eigvals
    use misctools,          only : dist

    ! in/out
    integer, intent(in)     :: idx(:)
    real, intent(out)       :: evecs(size(idx), size(idx))
    real, intent(out)       :: evals(size(idx))
    real, intent(in)        :: corlen
    integer, intent(in)     :: iexp

    ! local
    integer                 :: N, nregs, ireg, M, i_hor1, i_hor2, i, j, start_point, i_global, i_local
    real                    :: lat1, lat2, lon1, lon2, dst, cor
    integer, allocatable   :: sort_order(:), inv_sort_order(:), elem_per_reg(:), sorted_idx(:), uniq_reg_idx(:)
    real, allocatable      :: small_P(:,:), small_lam(:), padded_evec(:), Q(:,:)

    ! begin
    N = size(idx)
    allocate(sort_order(N))
    allocate(inv_sort_order(N))
    allocate(sorted_idx(N))
    allocate(padded_evec(N))

    ! idx is an array of size n_hor, where each element is an integer pointing to the region number of that grid cell. So
    ! if there are 5 regions in total, idx could look like [2, 1, 3, 3, 1, 4, 4, 5, 1, 3, 5, 1, 2]
    sort_order = argsort(idx) ! [2, 5, 9, 12, 1, 13, 3, 4, 10, 6, 7, 8, 11]
    ! so idx[sort_order] = [1, 1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 5, 5]
    inv_sort_order = argsort(sort_order) ! [5, 1, 7, 8, 2, 10, 11, 12, 3, 9, 13, 4, 6]
    sorted_idx = idx(sort_order)
    ! sorted_idx[inv_sort_order] = [2, 1, 3, 3, 1, 4, 4, 5, 1, 3, 5, 1, 2], i.e., the original array

    evecs = 0.0
    evals = 0.0

    ! how many distinct regions are there?
    call multiplicity(sorted_idx, uniq_reg_idx, elem_per_reg, .true.)
    nregs = size(elem_per_reg)

    start_point = 0
    i_global = 1
    do ireg = 1, nregs
        write(*, '("Region ", i2, " has ", i5, " of ", i5, " pixels")') uniq_reg_idx(ireg), elem_per_reg(ireg), N
        M = elem_per_reg(ireg)
        allocate(small_P(M,M))
        allocate(small_lam(M))

        ! now fill up small_P
        do i = 1, M
            i_hor1 = sort_order(start_point+i)
            lon1 = vec2ll(i_hor1)%lon
            lat1 = vec2ll(i_hor1)%lat
            do j = i, M
                i_hor2 = sort_order(start_point+j)
                lon2 = vec2ll(i_hor2)%lon
                lat2 = vec2ll(i_hor2)%lat
                call dist( lon1, lat1, lon2, lat2, dst )
                cor = exp( -(dst/corlen)**iexp )
                small_P(i,j) = cor
                small_P(j,i) = cor
            end do
        end do

        call eigvals(M, small_P, small_lam, .true.)
        call make_positive_semidef(small_lam)

        do i_local = 1, M
            padded_evec = 0.0
            padded_evec(start_point+1:start_point+M) = small_P(:,i_local)
            ! permute this eigenvector back to the full matrix
            evecs(:, i_global) = padded_evec(inv_sort_order)
            evals(i_global) = small_lam(i_local)
            i_global = i_global + 1
        end do
        start_point = start_point + M
        deallocate(small_P, small_lam)
    end do

    ! sort the eigenvalues and eigenvectors in ascending order of eigenvalues
    sort_order = argsort(evals)
    allocate(small_P(N,N))
    do i_global = 1, N
        small_P(:,i_global) = evecs(:, sort_order(i_global))
    end do
    evecs = small_P
    evals = evals(sort_order)
    deallocate(small_P)

    ! make sure the norm is 1 for all eigenvectors, just reuse the variable cor
    do i = 1, N
        cor = sqrt(sum(evecs(:,i) * evecs(:,i)))
        evecs(:,i) = evecs(:,i)/cor
    end do

    deallocate(sort_order, inv_sort_order, sorted_idx, padded_evec)
    deallocate(uniq_reg_idx, elem_per_reg)

  end subroutine block_correlation_hybrid

  subroutine block_correlation_eig(idx, evecs, evals)
    ! modules
    use sorting,            only : argsort
    use orderpack,          only : multiplicity
    use linalg_interface,   only : eigvals, orthonormalize

    ! in/out
    integer, intent(in)     :: idx(:)
    real, intent(out)       :: evecs(size(idx), size(idx))
    real, intent(out)       :: evals(size(idx))

    ! local
    integer                 :: N, nregs, ireg, M, i_global, i_local, start_point
    integer, allocatable    :: sort_order(:), inv_sort_order(:), elem_per_reg(:), sorted_idx(:), uniq_reg_idx(:)
    real, allocatable       :: small_P(:,:), small_lam(:), padded_evec(:), Q(:,:)

    ! begin
    N = size(idx)
    allocate(sort_order(N))
    allocate(inv_sort_order(N))
    allocate(sorted_idx(N))
    allocate(padded_evec(N))

    sort_order = argsort(idx)
    inv_sort_order = argsort(sort_order)
    sorted_idx = idx(sort_order)

    evecs = 0.0
    evals = 0.0

    ! how many distinct regions are there?
    call multiplicity(sorted_idx, uniq_reg_idx, elem_per_reg, .true.)
    nregs = size(elem_per_reg)

    start_point = 0
    i_global = 1
    do ireg = 1, nregs
        write(*, '("Region ", i2, " has ", i5, " pixels")') uniq_reg_idx(ireg), elem_per_reg(ireg)
        M = elem_per_reg(ireg)
        allocate(small_P(M,M))
        small_P = 1.0
        allocate(small_lam(M))
        !
        ! call eigvals(M, small_P, small_lam, .false.)
        !
        ! For a MxM matrix of ones, there is one eigenvalue M with eigenvector [1 1 ... 1], and (M-1) eigenvalues 0.
        ! The eigenvectors of the latter (M-1) eigenvalues can be constructed as [1 -1 0 0 ... 0],
        ! [1 0 -1 0 ... 0], [1 0 0 -1 ... 0] ... [1 0 0 0 ... -1] (idea to construct the M-1 eigenvectors this way from
        ! http://answers.yahoo.com/question/index?qid=20110228210510AAIfJ8a). After that, the first eigenvector needs to
        ! be normalized and the rest need to be orthonormalized.
        !
        small_lam(1) = M
        small_P(:,1) = 1.0/sqrt(float(M))
        allocate(Q(M,M-1))
        Q = 0.0
        Q(1,:) = 1.0
        do i_local = 2, M
            Q(i_local, i_local-1) = -1.0
        end do
        call orthonormalize(Q)
        small_lam(2:M) = 0.0
        small_P(:,2:M) = Q
        deallocate(Q)

        do i_local = 1, M
            padded_evec = 0.0
            padded_evec(start_point+1:start_point+M) = small_P(:,i_local)
            ! permute this eigenvector back to the full matrix
            evecs(:, i_global) = padded_evec(inv_sort_order)
            evals(i_global) = small_lam(i_local)
            i_global = i_global + 1
        end do
        start_point = start_point + M
        deallocate(small_P, small_lam)
    end do

    ! The eigenvalues of a bunch of ones is either M or 0, so reduce noise by setting the small ones to zero
    evals = merge(0.0, evals, abs(evals) < 0.01)

    ! sort the eigenvalues and eigenvectors in ascending order of eigenvalues
    sort_order = argsort(evals)
    allocate(small_P(N,N))
    do i_global = 1, N
        small_P(:,i_global) = evecs(:, sort_order(i_global))
    end do
    evecs = small_P
    evals = evals(sort_order)
    deallocate(small_P)

    deallocate(sort_order, inv_sort_order, sorted_idx, padded_evec)
    deallocate(uniq_reg_idx, elem_per_reg)

  end subroutine block_correlation_eig

    subroutine make_positive_semidef(lam)

        use GO,          only : ReadRc
        use global_data, only : rcF

        real, intent(inout) :: lam(:)
        integer              :: i, n, n_neg, status
        real                 :: min_eigval

        call ReadRc( rcF, 'var4d.horcor.min_eigval', min_eigval, status )

        if (min_eigval .gt. 1.0E-10) then
            min_eigval = min_eigval * minval((/ 1.0, maxval(lam) /))
        end if

        write(*,'("Maximum eigenvalue = ", e10.3, ", minimum eigenvalue = ", e10.3)') maxval(lam), minval(lam)

        n_neg = 0
        n = size(lam,1)
        do i=1,n
            if (lam(i) .lt. min_eigval) then
                lam(i) = min_eigval
                n_neg = n_neg + 1
            end if
        end do
        if (n_neg .gt. 0) write(*,'("Set ", i5, " eigenvalues to ", f15.11)') n_neg, min_eigval

    end subroutine make_positive_semidef

end module Var4D_State_Hori
