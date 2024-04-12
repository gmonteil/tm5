!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module redgridZoom

  use GO   , only : gol, goErr, goPr
  use binas, only : pi
  use dims, only : nregions, im,jm

  implicit none

  private

  !public :: uni2red_mfs
  public :: uni2red_mf, initredgrid, calc_pdiff
  public :: uni2red, red2uni, red2uni_em
#ifdef secmom
  public :: uni2red_2nd, red2uni_em_2nd
#endif
  public :: nred, nredmax, clustsize, grid_reduced
  public :: jred, imred, imredj

  ! --- const --------------------------------------

  character(len=*), parameter     ::  mname = 'redgridZoom'

  ! public parameters

  integer,parameter                     :: nredmax=30
  logical,parameter                     :: grid_reduced = .true.
  ! the number of reduced latitude circles per region
  integer,dimension(nregions)           :: nred
  ! number of joined cells per latitude circle and region
  integer,dimension(nredmax,nregions)   :: clustsize
  ! reduced dimension......
  integer,dimension(nredmax,nregions)   :: imredj
  ! reduced dimension over the full jm
  integer,dimension(180    ,nregions)   :: imred
  ! latitude numbers where the reduction applies...
  integer,dimension(nredmax,nregions)   :: jred

  ! private parameters

  real,parameter                        :: dl = 2.0*pi/im(1)
  real,parameter                        :: dp = pi/jm(1)
  real,parameter                        :: epsx = epsilon(0.0)
  !logical,parameter                     :: fromfile = .true.


contains

  !  Fortran code for Split-rg (reduced grid) scheme for advection on the
  !  globe.
  !
  !  VERSION:        1.1
  !  DATE:           November 12, 1998
  !
  !  Version 1.0     Dated October 19, 1998.
  !  Version 1.1     Corrected red2uni routine and added new uni2red_m
  !                  routine.
  !                  (Orography effects were erroneously not taken into
  !                  account in the reduced to uniform grid conversion)
  !
  !  This code solves the discrete advection equation for the tracer
  !  mass, with the wind field specified in air mass fluxes:
  !
  !  dtracm/dt = (
  !   mu chi (i-half) - mu chi (i+half)
  !   mv chi (j-half) - mv chi (j+half)
  !   mw chi (k-half) - mw chi (k+half) )
  !
  !  where
  !
  !  tracm = tracer mass
  !  mu, mv, mw = air mass fluxes
  !  chi = tracm / mass = tracer mixing ratio
  !  mass = air mass
  !
  !  The integer values of indices i, j, and k correspond to cell centers.
  !
  !  See:
  !
  !  Petersen, A.C., E.J. Spee, H. van Dop, and W. Hundsdorfer,
  !  An evaluation and intercomparison of four new advection schemes
  !  for use in global chemistry models,
  !  Journal of Geophysical Research, 103, 19253-19269, 1998.
  !
  !  Written by Edwin Spee, CWI, Amsterdam, The Netherlands
  !
  !
  !  Parts of the code are written by
  !  Joke Blom and Willem Hundsdorfer (CWI) and
  !
  !  The advection routine is called in the following way:
  !
  !     call advect_split(tracm, mu, mv, mw, mass, dt)
  !
  !  tracm = tracer mass            real(im,jm,lm,ntracet) [kg]
  !
  !  mu = longitudinal mass flux    real(im,jm,lm) [kg/s]
  !                                 (defined on eastward side of grid
  !                                 cells; for all grid cells)
  !
  !  mv = latitudinal mass flux     real(im,jm,lm) [kg/s]
  !                                 (defined on southward side of grid
  !                                 cells; only j=2,..,jm are used)
  !
  !  mw = vertical massflux         real(im,jm,lm) [kg/s]
  !                                 (defined on upper side of grid cell;
  !                                 only l=1,..,lm-1 are used since it is
  !                                 assumed that there is no flux across
  !                                 the top of the model)
  !
  !  The number of latitudinal elements actually used for mu, mv, and mw
  !  depends on the advection grid; for mv this number is the determined
  !  by the latitude row with the largest number of cells bordering a
  !  certain cell-bounding latitude circle.
  !
  !  mass = air mass                real(im,jm,lm) [kg]
  !
  !  dt   = timestep                real [s]
  !
  !  The physical units were chosen for compatibility with the TM3 model
  !  but can be changed for other global tracer models, provided that
  !  the calculation of the fluxes is done with an equivalent of tracer
  !  mixing ratios and not of tracer concentrations (make for other tracer
  !  models, if necessary, the appropriate changes in the routines
  !  split_lab, split_phi, and split_eta, where now tracer mixing ratios
  !  are calculated as tracm / mass; also the grid transformations as done
  !  in advect_split are different for tracer mass and concentration).
  !
  !  The longitude is counted from west to east, the latitude from south
  !  to north, and the height from surface to top.
  !
  !  The variables tracm, mu, mv, mw and mass should be declared in
  !  the calling program under any name; dt is a  parameter. The
  !  parameters im, jm, and lm denote the number of grid cells in
  !  longitudinal, latitudinal, and vertical direction, respectively.
  !  ntracet is the number of transported tracers. The letters 'u', 'v',
  !  and 'w' for the three dimensions are in the code often replaced by
  !  'l', 'p', and 'e', denoting the coordinates lambda, phi, and eta,
  !  respectively (see Petersen et al. 1998). The advection scheme works
  !  in principle with any grid and coordinate system of the main model.
  !  However, in the current implementation (made for the TM3 model) the
  !  variables tracm and mass are assumed to be defined in the main
  !  program on a uniform grid with polar cap, and the variables mu, mv,
  !  and mw on the advection grid (i.e., reduced grid with polar cap).
  !  Locally in the advection routines a transformation of tracm and mass,
  !  if necessary, is made to the advection grid. Please contact us for
  !  information on the changes that have to be made to work with other
  !  main model grids. Depending on the interests of users, we intend to
  !  generalize the code with respect to the main model grid in the next
  !  version.
  !
  !  The following parameters and variables, related to the specification
  !  of the advection grid are used by the routines of the advection
  !  scheme:
  !
  !     integer, parameter :: nredmax = 10
  !     real, parameter    :: dl = 2.0 * pi / im
  !     real, parameter    :: dp = pi / jm
  !     integer               nred
  !     integer               clustsize(0:nredmax),
  !                           jred(-nredmax:nredmax+1),imredj(jm)
  !     real                  maxtau
  !     integer               advsub
  !
  !   nredmax   = maximum number of grid reductions per hemisphere
  !               (in this implementation the minimum number of grid
  !               reductions per hemisphere equals one, since use is
  !               made of polar caps)
  !   dl        = longitudinal uniform grid cell size in lambda
  !               coordinates
  !   dp        = latitudinal grid cell size in phi coordinates
  !   nred      = number of grid reductions
  !   clustsize = array with the ratio im / imred, as a function of the
  !               absolute value of the latitude zone index, where imred
  !               is the actually used number of cells in the longitudinal
  !               direction for the given latitude zone
  !   jred      = latitude number where grid reduction starts (seen from
  !               SP), as a function of the reduction zone index
  !   imredj(j) = number of cells for latitude #j
  !   maxtau    = time step such that CFL <= 1 in all directions (given
  !               the directional splitting, see routine advect_split)
  !   advsub    = largest integer such that dtadv / advsub <= maxtau,
  !               with dtadv the advection time step in operator splitting
  !
  !  The advection grid variables for the reduced grid are automatically
  !  set by the following line:
  !
  !     call initredgrid(clustsize, jred, imredj, nred, nredmax, im, jm,
  !   .                  pi, .false.)
  !
  !  The variables maxtau and advsub should be set each time new mass
  !  fluxes (wind fields) are available in the model through
  !
  !     call split_maxtau(mu, mv, mw, mass, dtadv)
  !
  !  where dtadv is the advection time step in the main model.
  !
  !  Furthermore, some routines in the advection scheme need information
  !  on the machine precision. The following variable is used for this
  !  purpose:
  !
  !     real, parameter    :: epsx = epsilon(0.0)
  !
  !  The value of epsx is dependent on the computing platform that is used.
  !  epsilon is a Fortran 90 intrinsic function. For example, for our Cray
  !  C90 eps is about 1e-14.
  !
  !  The above-mentioned parameters and variables are assumed to be
  !  declared in the module global_data (Fortran 90) or in some common
  !  block (Fortran 77). The current implementation (Fortran 90) uses the
  !  statements
  !
  !     use global_data
  !     implicit none
  !
  !  If Fortran 77 is preferred one should define a common block or
  !  use an existing common block to include the above-given parameters
  !  and variables, e.g.
  !
  !     implicit none
  !     (declarations as given in the above)
  !     common/advgrid/nredmax,dl,dp,nred,clustsize,jred,maxtau,advsub
  !
  !  Also the parameters im, jm, lm, ntracet, pi, twopi (= 2 * pi), and
  !  eps are assumed to be declared in the module global_data.
  !
  !  --------------------------
  !  Example of a reduced grid
  !  --------------------------
  !
  !  NP=North Pole
  !  SP=South Pole
  !  EQ=Equator
  !  | grid boundary in longitudal direction
  !  o cell center
  !                                 j   imredj   clustsize   latitude zone
  !   NP|           o           |  13     1        12             3
  !     |   o   |   o   |   o   |  12     3         4             2
  !     | o | o | o | o | o | o |  11     6         2             1
  !     | o | o | o | o | o | o |  10     6         2             1
  !     |o|o|o|o|o|o|o|o|o|o|o|o|   9    12         1             0
  !     |o|o|o|o|o|o|o|o|o|o|o|o|   8    12         1             0
  !   EQ|o|o|o|o|o|o|o|o|o|o|o|o|   7    12         1             0
  !     |o|o|o|o|o|o|o|o|o|o|o|o|   6    12         1             0
  !     |o|o|o|o|o|o|o|o|o|o|o|o|   5    12         1             0
  !     | o | o | o | o | o | o |   4     6         2            -1
  !     | o | o | o | o | o | o |   3     6         2            -1
  !     |   o   |   o   |   o   |   2     3         4            -2
  !   SP|           o           |   1     1        12            -3
  !
  !  This grid would have
  !  im = 12, jm = 13, nred = 3
  !
  !  jred(-3) = 1 by definition
  !  jred(-2) = 2
  !  jred(-1) = 3
  !  jred( 0) = 5
  !  jred( 1) =10
  !  jred( 2) =12
  !  jred( 3) =13
  !  jred( 4) =14=jm+1 by definition
  !
  !  clustsize(3) = 12
  !  clustsize(2) =  4
  !  clustsize(1) =  2
  !  clustsize(0) =  1

  !CMK total new definition to allow zoom regios..
  !  nred(nregions)    :: number of reduced latitudes circles   < jm(region)
  !  nredmax           :: maximum number of reduced latitude circles...
  !  clustsize(nredmax.nregions)  :: number of joined cells for n'th circle
  !  jred(nredmax.nregions)       :: corresponding j value in 1...jm(region) array
  !  imredj(nredmax.nregions)     :: 'new' im value for the reduced grid...
  !CMK

!  subroutine uni2red_mfs(region,afl,bfl,cfl)
!
!    !converts am,bm, and cm to reduced grid.
!    !am is compressed
!    !bm is averaged (all mass-fluxes will be the same)
!    !cm is added    (since we just loop over the reduced grid)
!    use dims, only : im, jm, lm
!    use dims, only : okdebug
!    use dims, only : one
!#ifdef MPI
!    use mpi_const,only : lmloc
!#endif
!    implicit none
!
!    ! in/out
!    integer,intent(in) :: region
!    real,dimension(0:im(region)+1,0:jm(region)+1,0:lm(region)) :: afl
!    real,dimension(0:im(region)+1,0:jm(region)+1,0:lm(region)) :: bfl
!    real,dimension(0:im(region)+1,0:jm(region)+1,0:lm(region)) :: cfl
!
!    ! local
!    integer :: clust_size, i, j, l, i2, imax, jreg, is, ie
!    integer :: lmr
!
!    if (okdebug) then
!       print *,'uni2red_mfs: region =',region
!    end if
!
!#ifdef MPI
!    lmr=lmloc
!#else
!    lmr = lm(region)
!#endif
!
!    do l=1,lmr
!       do jreg = 1,nred(region)   !loop over the reduced latitudes....
!          j = jred(jreg,region)   !latitide
!          clust_size = clustsize(jreg,region)   !clustersize....
!          do i=1,imredj(jreg,region)
!             afl(i,j,l) = afl(i*clust_size,j,l)   !compress afl
!             is = (i-1)*clust_size + 1
!             ie = i*clust_size
!             cfl(i,j,l) = sum(cfl(is:ie,j,l))            !put the total flux
!             !        bfl(i,j+1,l) = sum(bfl(is:ie,j+1,l))/(ie-is+1)
!             ! calculate the average of norther flux!
!          end do
!          do i=imredj(jreg,region)+1,im(region)
!             afl(i,j,l) = -one
!             cfl(i,j,l) = -one
!          end do
!          do i=imredj(jreg,region),1,-1
!             is = (i-1)*clust_size + 1
!             ie = i*clust_size
!             !        bfl(is:ie,j+1,l)=bfl(i,j+1,l)
!             !  put the same value everywhere on northern edge..
!          end do
!          ! the next statement is only meaningful for region cyclic regions;
!          ! this subroutine should be applied only in this region
!          ! set eastern boundary condition...
!          afl(0,j,l)=afl(imredj(jreg,region),j,l)
!       end do
!    end do
!
!  end subroutine uni2red_mfs



  subroutine uni2red_mf(region)

    ! Converts mass fluxes mfl from
    ! uniform grid to reduced grid
    !
    ! Written by Edwin Spee, CWI, Amsterdam, The Netherlands
    !
    ! On input:
    ! mfl = mass fluxes on uniform grid
    !
    ! On output:
    ! mfl = mass fluxes on reduced grid
    !cmk changed.....
    use dims, only : lm
    use dims, only : okdebug
#ifdef MPI
    use mpi_const, only: lmloc
#endif
    use global_data, only: wind_dat

    implicit none

    ! in/out
    integer,intent(in) :: region

    ! local
    real,dimension(:,:,:),pointer :: mfl
    integer :: clust_size, i, j, l, i2, imax, jreg
    integer :: lmr

    mfl => wind_dat(region)%am_t

    lmr=lm(region)

    if ( okdebug ) then
       print *,'uni2red_mf: region =',region
    end if

    ! selection for mfl
    do l=1,lmr
       do jreg = 1,nred(region)   !loop over the reduced latitudes....
          j = jred(jreg,region)   !latitide
          clust_size = clustsize(jreg,region)   !clustersize....
          do i=1,imredj(jreg,region)
             mfl(i,j,l) = mfl(i*clust_size,j,l)   !compress mfl
          end do
          do i=imredj(jreg,region)+1,im(region)
             mfl(i,j,l) = -1.0
          end do
          ! the next statement is only meaningful for region cyclic regions;
          ! this subroutine should be applied only in this region
          ! set eastern boundary condition...
          mfl(0,j,l)=mfl(imredj(jreg,region),j,l)
       end do
    end do

    nullify(mfl)

  end subroutine uni2red_mf


  !-----------------------------------------------------------------------
  ! Below are some routines related to the initialization of the advection
  ! grid data structure (i.e. the reduced grid).
  !-----------------------------------------------------------------------

  subroutine initredgrid( region, imr, jmr, status )

    ! The subroutine initredgrid fills the arrays clustsize, jred, and
    ! imredj with appropriate values for the use of a reduced grid.
    !
    ! Written by Edwin Spee, CWI, Amsterdam, The Netherlands
    !
    ! This implementation is based on work of Joke Blom.
    !
    ! On input:
    ! clustsize,jred,imredj,nred undefined
    ! nredmax
    ! im,jm   number of grids in lon, lat direction
    ! pi      3.141....
    !
    ! !if fromfile=true: read grid definition from file, else generate one.
    !
    ! On output:
    ! array clustsize,jred,imredj contain grid definition
    ! integer nred is number of grid reductions

    use dims, only : okdebug
    use toolbox,     only : escape_tm
    implicit none

    ! in/out
    integer,intent(in) :: region
    integer,intent(in) :: imr
    integer,intent(in) :: jmr
    integer, intent(out)          ::  status

    ! --- const ------------------------------

    character(len=*), parameter ::  rname = mname//'/initredgrid'

    ! --- local ------------------------------

    integer   :: j, lrg

    ! --- begin -----------------------------

    call read_redgrid( region, status )
    IF_NOTOK_RETURN(status=1)

    imred(1:jm(region),region) = im(region)    ! the full jm array...
    do lrg=1,nred(region)
       imredj(lrg,region)=imr/clustsize(lrg,region)
       imred(jred(lrg,region),region)=imredj(lrg,region)
    end do

    if ( okdebug ) then
       do j=1,nred(region)
          print *,'initredgrid: imredj(',jred(j,region),')=',imredj(j,region)
       end do
    end if

    ! ok
    status = 0

  end subroutine initredgrid


  ! *


  subroutine read_redgrid( region, status )

    use GO         , only : ReadRc
    use global_data, only : rcF
    use dims   , only : okdebug
    use dims   , only : region_name, children
    use dims   , only : ybeg, yend, dy, yref

    ! --- in/out ----------------------------------

    integer, intent(in)           ::  region
    integer, intent(out)          ::  status

    ! --- const ------------------------------

    character(len=*), parameter ::  rname = mname//'/read_redgrid'

    ! --- local ----------------------------------------

    integer               ::  jmr
    integer               ::  i,j
    integer               ::  xam, nim
    integer               ::  nc, child,latps,latpn
    integer               ::  jband, ncomb
    integer, allocatable  ::  ncombs(:)
    integer               ::  nred_nh
    integer, allocatable  ::  ncombs_nh(:)
    integer               ::  nred_sh
    integer, allocatable  ::  ncombs_sh(:)

    ! --- begin ----------------------------------------

    ! number of lats
    jmr = jm(region)

    ! number of reduced lat bands:
    call ReadRc( rcF, 'region.'//trim(region_name(region))//'.redgrid.nh.n', nred_nh, status )
    IF_NOTOK_RETURN(status=1)
    ! number of reduced lat bands:
    call ReadRc( rcF, 'region.'//trim(region_name(region))//'.redgrid.sh.n', nred_sh, status )
    IF_NOTOK_RETURN(status=1)

    ! total number of reduced bands:
    nred(region) = nred_nh + nred_sh

    ! reduced grid ?
    if ( nred(region) > 0 ) then

      if ( okdebug ) then
        write (gol,'("read_redgrid:  Reading parameters from file rcfile.")'); call goPr
        write (gol,'("read_redgrid:  Region ",a)') trim(region_name(region)); call goPr
      end if
      if ( nred(region) > nredmax ) then
        write (gol,'("problem with reduced grid. nred > nredmax")'); call goErr
        TRACEBACK; status=1; return
      end if

      ! storage for combined cells for each lat; 1 by default
      allocate( ncombs(jmr) )
      ncombs = 1

      ! counter of reduced lat bands from south-pole to north-pole:
      i = 0

      ! southern hemisphere ?
      if ( nred_sh > 0 ) then
        ! temporary storage:
        allocate( ncombs_sh(jmr) )
        ! number of combined cells per lat band:
        call ReadRc( rcF, 'region.'//trim(region_name(region))//'.redgrid.sh.comb', ncombs_sh(1:nred_sh), status )
        IF_NOTOK_RETURN(status=1)
        ! loop over lat bands:
        do jband = 1, nred_sh
          ! increase counter:
          i = i + 1
          ! select:
          ncomb = ncombs_sh(jband)
          ! store band:
          jred(i,region) = jband
          ! clustsize is number of original cells combined into one reduced cell:
          clustsize(i,region) = ncomb
          ! idem for testing
          ncombs(jmr+1-jband) = ncomb
          ! now check if the latitude circle overlaps with
          ! southernmost latitude of gridcel parent ...
          latps = ybeg(region) + (jred(i,region)-1)*nint(dy)/yref(region)
          ! northermost latitude ...
          latpn = ybeg(region) + (jred(i,region))*nint(dy)/yref(region)
          do nc = 1, children(region,0)
            child = children(region,nc)
            if (latps < yend(child) .and. latpn > ybeg(child) ) then
              write (gol,'("overlap detected between child region and reduced grid parent")'); call goErr
              write (gol,'("  parent : ",i4)') region; call goErr
              write (gol,'("  child  : ",i4)') child; call goErr
              TRACEBACK; status=1; return
            end if
          end do
        end do  ! lat bands sh
        ! clear:
        deallocate( ncombs_sh )
      end if  ! sh

      ! northern hemisphere ?
      if ( nred_nh > 0 ) then
        ! temporary storage:
        allocate( ncombs_nh(jmr) )
        ! number of combined cells per lat band:
        call ReadRc( rcF, 'region.'//trim(region_name(region))//'.redgrid.nh.comb', ncombs_nh(1:nred_nh), status )
        IF_NOTOK_RETURN(status=1)
        ! loop over lat bands; read in order from pole to equator, thus reverse for original storage:
        do jband = nred_nh, 1, -1
          ! increase counter:
          i = i + 1
          ! select:
          ncomb = ncombs_nh(jband)
          ! store band:
          jred(i,region) = jm(region) + 1 - jband
          ! clustsize is number of original cells combined into one reduced cell:
          clustsize(i,region) = ncomb
          ! idem for testing
          ncombs(jband) = ncomb
          ! now check if the latitude circle overlaps with
          ! southermost latitude of gridcel parent ...
          latps = ybeg(region) + (jred(i,region)-1)*nint(dy)/yref(region)
          ! northermost latitude ...
          latpn = ybeg(region) + (jred(i,region))*nint(dy)/yref(region)
          do nc = 1, children(region,0)
            child = children(region,nc)
            if (latps < yend(child) .and. latpn > ybeg(child) ) then
              write (gol,'("overlap detected between child region and reduced grid parent")'); call goErr
              write (gol,'("  parent : ",i4)') region; call goErr
              write (gol,'("  child  : ",i4)') child; call goErr
              TRACEBACK; status=1; return
            end if
          end do
        end do  ! lat bands nh
        ! clear:
        deallocate( ncombs_nh )
      end if  ! nh

      ! testing ...
      do j = 1, jmr
        ! check if number of combined cells matches with grid:
        if ( modulo(im(region),ncombs(j)) /= 0 ) then
          write (gol,'("number of combined cells not ok:")'); call goErr
          write (gol,'("  region         : ",i2," ",a)') region, trim(region_name(region)); call goErr
          write (gol,'("  lat band       : ",i4)') j; call goErr
          write (gol,'("  combined cells : ",i4)') ncombs(j); call goErr
          write (gol,'("  im             : ",i4)') im(region); call goErr
          TRACEBACK; status=1; return
        end if
        ! check with previous ...
        if ( j > 1 ) then
          xam = max(ncombs(j-1),ncombs(j))
          nim = min(ncombs(j-1),ncombs(j))
          if ( modulo( xam, nim ) /= 0 ) then
            write (gol,'("number of combined cells does match with previous:")'); call goErr
            write (gol,'("  region         : ",i2," ",a)') region, trim(region_name(region)); call goErr
            write (gol,'("  lat band       : ",i4)') j; call goErr
            write (gol,'("  combined cells : ",i4)') ncombs(j); call goErr
            write (gol,'("  previous       : ",i4)') ncombs(j-1); call goErr
            TRACEBACK; status=1; return
          end if
        end if
      end do  ! lats

      !! info ...
      !write (gol,'("")'); call goPr
      !write (gol,'("region         : ",i2," ",a)') region, trim(region_name(region)); call goPr
      !write (gol,'("  band clustsize ncomb")'); call goPr
      !do j = 1, jmr
      !  write (gol,'(i6,i10,i6)') j, clustsize(j,region), ncombs(j); call goPr
      !end do

      ! clear
      deallocate( ncombs )

    else

      write (gol,'("WARNING - No reduced grid for region :",a)') trim(region_name(region)); call goPr
      nred(region) = 0   !default no reduced grid...

    end if

    ! ok
    status = 0

  end subroutine read_redgrid


  ! *


  subroutine calc_pdiff(region,p,pold,pdiffmax)
    !----------------------------------
    !
    !----------------------------------
    use dims, only : xref, yref, parent, xcyc, touch_sp, touch_np
    implicit none

    ! in/out
    integer,intent(in)                                         :: region
    real,dimension(-1:im(region)+2,-1:jm(region)+2),intent(in) :: pold,p
    real,intent(out)                                           :: pdiffmax

    ! local
    integer :: lrg,i,j,ratio,idx,iu,xref_,yref_
    real    :: work,work1

    pdiffmax = 0.0
    yref_ = yref(region)/yref(parent(region))
    xref_ = xref(region)/xref(parent(region))

    if ( xcyc(region) /= 1 .or. nred(region) == 0 ) then   ! no reduced grid:
       jloop: do j=1,jm(region)
          if (yref_>1) then    !there is a refinement in the y direction
             if(  (j < yref_ + 1) .and. &
                  (touch_sp(region) /= 1) ) cycle jloop
             if(  (j > jm(region) - yref_) .and. &
                  (touch_np(region) /= 1) ) cycle jloop
          end if
          iloop: do i=1,im(region)
             if ( xref_ > 1 .and. xcyc(region) /= 1 ) then
                !there is refinement in the x direction and region not cyclic.
                if( (i < xref_+1) .or. ( i>im(region)-xref_) ) cycle iloop
             end if
             pdiffmax = max(pdiffmax,abs(p(i,j)-pold(i,j)))
          end do iloop
       end do jloop
    else   ! gid is reduced: do something special:
       do lrg = 1,nred(region)   !first maximum over the reduced grid...!
          ratio = clustsize(lrg,region)
          j = jred(lrg,region)
          do i =  1,imredj(lrg,region)
             work = 0.0
             work1 = 0.0
             do iu = 1,ratio
                idx  = (i-1)*ratio+iu
                work = work + p(idx,j)
                work1 = work1 + pold(idx,j)
             end do
             pdiffmax=max(pdiffmax,abs(work-work1)/ratio)
          end do
       end do
       ! now the pressure difference in the non-reduced part...
       yref_ = yref(region)/yref(parent(region))
       j2loop: do j=1,jm(region)
          if(imred(j,region).ne.im(region)) cycle j2loop   !if reduced...skip
          if ( yref_ > 1 ) then
             ! there is a refinement in the y direction: skip the edges!
             if(  (j < yref_ + 1) .and. &
                  (touch_sp(region) /= 1) ) cycle j2loop
             if(  (j > jm(region) - yref_) .and. &
                  (touch_np(region) /= 1) ) cycle j2loop
          end if
          i2loop: do i=1,im(region)
             pdiffmax = max(pdiffmax,abs(p(i,j)-pold(i,j)))
          end do i2loop
       end do j2loop
    end if
  end subroutine calc_pdiff


  subroutine uni2red(region)
    !
    ! transforms data from uniform grid to reduced grid
    ! written by mike botchev, march-june 1999
    ! modified by Maarten Krol, dec 2002
    !
    use dims
    use global_data, only: mass_dat
    use MeteoData  , only : m_dat
#ifdef MPI
    use mpi_const,only : ntracetloc
#else
    use chem_param,only : ntracet
#endif
    implicit none

    ! input
    integer,intent(in)  :: region

    ! local
    real,dimension(:,:,:,:),pointer         :: rm, rxm,rym,rzm
    real,dimension(:,:,:),  pointer         :: m

    integer i,ie,is,j,l,lrg,redfact,n,lmr
    real summ,sumrm

    ! start
    m => m_dat(region)%data
    rm => mass_dat(region)%rm_t
    rxm => mass_dat(region)%rxm_t
    rym => mass_dat(region)%rym_t
    rzm => mass_dat(region)%rzm_t

    lmr=lm(region)

    do lrg=1,nred(region)
       redfact=clustsize(lrg,region)
       j = jred(lrg,region)
       do l=1,lmr
          do i =  1,imredj(lrg,region)
             ! the is:ie  array section will be reduced to i
             is = (i-1)*redfact + 1
             ie = i*redfact
             summ = sum(m(is:ie,j,l))
             m(i,j,l) = summ
             !cmkm_uni(is:ie,j,l) = m_uni(is:ie,j,l)/summ
             !   use as distribution function
             ! when transferring back from reduced--->uniform grid...
             ! these summations mean that mixing ratio and the
             ! the slopes are averaged out within the is:ie section
             ! with m(is:ie,...) taken as the weights
#ifdef MPI
             do n=1,ntracetloc
#else
             do n=1,ntracet
#endif
                sumrm = sum(rm(is:ie,j,l,n))
                rm(i,j,l,n)  = sumrm
                sumrm = sum(rxm(is:ie,j,l,n))
                rxm(i,j,l,n) = sumrm
                sumrm = sum(rym(is:ie,j,l,n))
                rym(i,j,l,n) = sumrm
                sumrm = sum(rzm(is:ie,j,l,n))
                rzm(i,j,l,n) = sumrm
             end do   !n
          end do  !i
          ! JFM: set remaining masses to zero
          !      for consistency with adjoint
          do i = imredj(lrg,region)+1, im(region)
#ifdef MPI
             do n=1,ntracetloc
#else
             do n=1,ntracet
#endif
                rm(i,j,l,n) = 0.
                rxm(i,j,l,n) = 0.
                rym(i,j,l,n) = 0.
                rzm(i,j,l,n) = 0.
             end do  !n
          end do  !i

          !put periodic boundary...
       end do  !l

    end do   !redgrid...

    nullify(m)
    nullify(rm)
    nullify(rxm)
    nullify(rym)
    nullify(rzm)

  end subroutine uni2red

#ifdef secmom
  subroutine uni2red_2nd(region)
    !
    ! transforms data from uniform grid to reduced grid
    ! written by mike botchev, march-june 1999
    ! modified by Maarten Krol, dec 2002
    !
    use dims
    use global_data, only: mass_dat
#ifdef MPI
    use mpi_const,only : ntracetloc
#else
    use chem_param,only : ntracet
#endif
    implicit none

    ! input
    integer,intent(in)  :: region

    ! local
    real,dimension(:,:,:,:),pointer         :: rm, rxm,rym,rzm, rxxm, rxym, rxzm, ryym, ryzm, rzzm
    real,dimension(:,:,:),  pointer         :: m

    integer i,ie,is,j,l,lrg,redfact,n,lmr
    real summ,sumrm

    ! start
    m => mass_dat(region)%m_t
    rm => mass_dat(region)%rm_t
    rxm => mass_dat(region)%rxm_t
    rym => mass_dat(region)%rym_t
    rzm => mass_dat(region)%rzm_t
     rxxm => mass_dat(region)%rxxm_t
     rxym => mass_dat(region)%rxym_t
     rxzm => mass_dat(region)%rxzm_t
     ryym => mass_dat(region)%ryym_t
     ryzm => mass_dat(region)%ryzm_t
     rzzm => mass_dat(region)%rzzm_t

    lmr=lm(region)

    do lrg=1,nred(region)
       redfact=clustsize(lrg,region)
       j = jred(lrg,region)
       do l=1,lmr
          do i =  1,imredj(lrg,region)
             ! the is:ie  array section will be reduced to i
             is = (i-1)*redfact + 1
             ie = i*redfact
             summ = sum(m(is:ie,j,l))
             m(i,j,l) = summ
             !cmkm_uni(is:ie,j,l) = m_uni(is:ie,j,l)/summ
             !   use as distribution function
             ! when transferring back from reduced--->uniform grid...
             ! these summations mean that mixing ratio and the
             ! the slopes are averaged out within the is:ie section
             ! with m(is:ie,...) taken as the weights
#ifdef MPI
             do n=1,ntracetloc
#else
             do n=1,ntracet
#endif
                sumrm = sum(rm(is:ie,j,l,n))
                rm(i,j,l,n)  = sumrm
                sumrm = sum(rxm(is:ie,j,l,n))
                rxm(i,j,l,n) = sumrm
                sumrm = sum(rym(is:ie,j,l,n))
                rym(i,j,l,n) = sumrm
                sumrm = sum(rzm(is:ie,j,l,n))
                rzm(i,j,l,n) = sumrm
                sumrm = sum(rxxm(is:ie,j,l,n))
                rxxm(i,j,l,n) = sumrm
                sumrm = sum(rxym(is:ie,j,l,n))
                rxym(i,j,l,n) = sumrm
                sumrm = sum(rxzm(is:ie,j,l,n))
                rxzm(i,j,l,n) = sumrm
                sumrm = sum(ryym(is:ie,j,l,n))
                ryym(i,j,l,n) = sumrm
                sumrm = sum(ryzm(is:ie,j,l,n))
                ryzm(i,j,l,n) = sumrm
                sumrm = sum(rzzm(is:ie,j,l,n))
                rzzm(i,j,l,n) = sumrm
             end do   !n
          end do  !i
          !put periodic boundary...
       end do  !l

    end do   !redgrid...

    nullify(m)
    nullify(rm)
    nullify(rxm)
    nullify(rym)
    nullify(rzm)
    nullify(rxxm)
    nullify(rxym)
    nullify(rxzm)
    nullify(ryym)
    nullify(ryzm)
    nullify(rzzm)

  end subroutine uni2red_2nd
#endif


  subroutine red2uni(region,m_uni)
    !
    ! transforms data from reduced grid back to uniform grid
    ! written by mike botchev, march-june 1999
    !
    use dims
    use global_data, only : mass_dat
    use MeteoData  , only : m_dat
#ifdef MPI
    use mpi_const,   only : ntracetloc
#else
    use chem_param,  only : ntracet
#endif
    implicit none

    ! input
    integer,intent(in) :: region
    real,dimension(-1:im(region)+2,-1:jm(region)+2,lm(region)), &
         intent(in) :: m_uni
    ! m_uni: same declaration as m!

    ! local
    real,dimension(:,:,:,:),pointer    :: rm,rxm,rym,rzm
    real,dimension(:,:,:)  ,pointer    :: m
    integer :: i, ie, ii, is, j, l, lrg, n, redfact, lmr
    real    :: hi, mass, mass_coord, rmm, slope, m_old
    character(len=5) :: distr_mode

    ! start
    lmr=lm(region)

    m => m_dat(region)%data
    rm => mass_dat(region)%rm_t
    rxm => mass_dat(region)%rxm_t
    rym => mass_dat(region)%rym_t
    rzm => mass_dat(region)%rzm_t

    distr_mode = 'unfrm' ! 'slope' or 'unfrm'

    do lrg=1,nred(region)
       redfact=clustsize(lrg,region)
       j = jred(lrg,region)
       do l=1,lmr
          do i =  imredj(lrg,region),1,-1
             ! the i cell will be distributed within the is:ie array section
             is = (i-1)*redfact + 1
             ie = i*redfact

             !m_uni is the mass-distribution in the non-reduced grid/divided by
             !the reduced_grid mass. This is used as distribution function!....
             mass=m(i,j,l); m(is:ie,j,l)= m_uni(is:ie,j,l)

             if (distr_mode=='unfrm') then
                ! mixing ratio and x-slope will be UNiFoRMly distributed
                ! within is:ie
#ifdef MPI
                do n=1,ntracetloc
#else
                do n=1,ntracet
#endif
                   rmm = rm(i,j,l,n)
                   rm(is:ie,j,l,n)= m_uni(is:ie,j,l)/mass* rmm
                   rmm = rxm(i,j,l,n)
                   rxm(is:ie,j,l,n)= m_uni(is:ie,j,l)/mass * rmm
                   ! rym and rzm are always distributed uniformly:
                   rmm = rym(i,j,l,n)
                   rym(is:ie,j,l,n)= m_uni(is:ie,j,l)/mass * rmm
                   rmm = rzm(i,j,l,n)
                   rzm(is:ie,j,l,n)= m_uni(is:ie,j,l)/mass * rmm
                end do
             end if
             !cmkelseif(distr_mode=='slope') then
          end do
          ! update cell(0,...) according to the periodic bc's:
#ifdef MPI
          rm(0,j,l,1:ntracetloc) = rm(im(region),j,l,1:ntracetloc)
#else
          rm(0,j,l,1:ntracet) = rm(im(region),j,l,1:ntracet)
#endif
          rxm(0,j,l,:) = rxm(im(region),j,l,:)
          rym(0,j,l,:) = rym(im(region),j,l,:)
          rzm(0,j,l,:) = rzm(im(region),j,l,:)
          m(0,j,l) = m(im(region),j,l)

       end do
    end do

    nullify(m)
    nullify(rm)
    nullify(rxm)
    nullify(rym)
    nullify(rzm)

  end subroutine red2uni



  subroutine red2uni_em(region)
    !
    ! transforms data from reduced grid back to uniform grid
    ! written by mike botchev, march-june 1999
    !
    use dims
    use global_data, only: mass_dat
    use MeteoData  , only : m_dat
#ifdef MPI
    use mpi_const, only : ntracetloc
#else
    use chem_param,only : ntracet
#endif

    implicit none

    ! input
    integer,intent(in) :: region

    ! local
    real,dimension(:,:,:,:),pointer    :: rm,rxm,rym,rzm
    real,dimension(:,:,:)  ,pointer    :: m
    integer i,ie,ii,is,j,l,lrg,n,redfact,lmr
    real hi,mass,mass_coord,rmm,slope,m_old
    character(len=5) distr_mode

    ! start

    lmr=lm(region)

    m => m_dat(region)%data
    rm => mass_dat(region)%rm_t
    rxm => mass_dat(region)%rxm_t
    rym => mass_dat(region)%rym_t
    rzm => mass_dat(region)%rzm_t

    distr_mode = 'unfrm' ! 'slope' or 'unfrm'

    do lrg=1,nred(region)
       redfact=clustsize(lrg,region)
       j = jred(lrg,region)
       do l=1,lmr
          do i =  imredj(lrg,region),1,-1
             ! the i cell will be distributed within the is:ie array section
             is = (i-1)*redfact + 1
             ie = i*redfact

             !m_uni is the mass-distribution in the non-reduced grid/divided by
             !the reduced_grid mass. This is used as distribution function!....
             mass=m(i,j,l); m(is:ie,j,l)= mass/(ie-is+1)

             if (distr_mode=='unfrm') then
                ! mixing ratio and x-slope will be UNiFoRMly distributed
                ! within is:ie
#ifdef MPI
                do n=1,ntracetloc
#else
                do n=1,ntracet
#endif
                   rmm = rm(i,j,l,n)
                   rm(is:ie,j,l,n)= m(is:ie,j,l)/mass* rmm
                   rmm = rxm(i,j,l,n)
                   rxm(is:ie,j,l,n)= m(is:ie,j,l)/mass * rmm
                   ! rym and rzm are always distributed uniformly:
                   rmm = rym(i,j,l,n)
                   rym(is:ie,j,l,n)= m(is:ie,j,l)/mass * rmm
                   rmm = rzm(i,j,l,n)
                   rzm(is:ie,j,l,n)= m(is:ie,j,l)/mass * rmm
                end do
             end if
             !cmkelseif(distr_mode=='slope') then
          end do
          ! update cell(0,...) according to the periodic bc's:
#ifdef MPI
          rm(0,j,l,1:ntracetloc) = rm(im(region),j,l,1:ntracetloc)
#else
          rm(0,j,l,1:ntracet) = rm(im(region),j,l,1:ntracet)
#endif
          rxm(0,j,l,:) = rxm(im(region),j,l,:)
          rym(0,j,l,:) = rym(im(region),j,l,:)
          rzm(0,j,l,:) = rzm(im(region),j,l,:)
          m(0,j,l) = m(im(region),j,l)

       end do
    end do

    nullify(m)
    nullify(rm)
    nullify(rxm)
    nullify(rym)
    nullify(rzm)

  end subroutine red2uni_em

#ifdef secmom
  subroutine red2uni_em_2nd(region)
    !
    ! transforms data from reduced grid back to uniform grid
    ! written by mike botchev, march-june 1999
    !
    use dims
    use global_data, only: mass_dat
#ifdef MPI
    use mpi_const, only : ntracetloc
#else
    use chem_param,only : ntracet
#endif

    implicit none

    ! input
    integer,intent(in) :: region

    ! local
    real,dimension(:,:,:,:),pointer    :: rm,rxm,rym,rzm, rxxm, rxym, rxzm, ryym, ryzm, rzzm
    real,dimension(:,:,:)  ,pointer    :: m
    integer i,ie,ii,is,j,l,lrg,n,redfact,lmr
    real hi,mass,mass_coord,rmm,slope,m_old
    character(len=5) distr_mode

    ! start

    lmr=lm(region)

    m => mass_dat(region)%m_t
    rm => mass_dat(region)%rm_t
    rxm => mass_dat(region)%rxm_t
    rym => mass_dat(region)%rym_t
    rzm => mass_dat(region)%rzm_t
     rxxm => mass_dat(region)%rxxm_t
     rxym => mass_dat(region)%rxym_t
     rxzm => mass_dat(region)%rxzm_t
     ryym => mass_dat(region)%ryym_t
     ryzm => mass_dat(region)%ryzm_t
     rzzm => mass_dat(region)%rzzm_t

    distr_mode = 'unfrm' ! 'slope' or 'unfrm'

    do lrg=1,nred(region)
       redfact=clustsize(lrg,region)
       j = jred(lrg,region)
       do l=1,lmr
          do i =  imredj(lrg,region),1,-1
             ! the i cell will be distributed within the is:ie array section
             is = (i-1)*redfact + 1
             ie = i*redfact

             !m_uni is the mass-distribution in the non-reduced grid/divided by
             !the reduced_grid mass. This is used as distribution function!....
             mass=m(i,j,l); m(is:ie,j,l)= mass/(ie-is+1)

             if (distr_mode=='unfrm') then
                ! mixing ratio and x-slope will be UNiFoRMly distributed
                ! within is:ie
#ifdef MPI
                do n=1,ntracetloc
#else
                do n=1,ntracet
#endif
                   rmm = rm(i,j,l,n)
                   rm(is:ie,j,l,n)= m(is:ie,j,l)/mass* rmm
                   rmm = rxm(i,j,l,n)
                   rxm(is:ie,j,l,n)= m(is:ie,j,l)/mass * rmm
                   ! rym and rzm are always distributed uniformly:
                   rmm = rym(i,j,l,n)
                   rym(is:ie,j,l,n)= m(is:ie,j,l)/mass * rmm
                   rmm = rzm(i,j,l,n)
                   rzm(is:ie,j,l,n)= m(is:ie,j,l)/mass * rmm
                   rmm = rxxm(i,j,l,n)
                   rxxm(is:ie,j,l,n)= m(is:ie,j,l)/mass * rmm
                   rmm = rxym(i,j,l,n)
                   rxym(is:ie,j,l,n)= m(is:ie,j,l)/mass * rmm
                   rmm = rxzm(i,j,l,n)
                   rxzm(is:ie,j,l,n)= m(is:ie,j,l)/mass * rmm
                   rmm = ryym(i,j,l,n)
                   ryym(is:ie,j,l,n)= m(is:ie,j,l)/mass * rmm
                   rmm = ryzm(i,j,l,n)
                   ryzm(is:ie,j,l,n)= m(is:ie,j,l)/mass * rmm
                   rmm = rzzm(i,j,l,n)
                   rzzm(is:ie,j,l,n)= m(is:ie,j,l)/mass * rmm
                end do
             end if
             !cmkelseif(distr_mode=='slope') then
          end do
          ! update cell(0,...) according to the periodic bc's:
#ifdef MPI
          rm(0,j,l,1:ntracetloc) = rm(im(region),j,l,1:ntracetloc)
#else
          rm(0,j,l,1:ntracet) = rm(im(region),j,l,1:ntracet)
#endif
          rxm(0,j,l,:) = rxm(im(region),j,l,:)
          rym(0,j,l,:) = rym(im(region),j,l,:)
          rzm(0,j,l,:) = rzm(im(region),j,l,:)
          m(0,j,l) = m(im(region),j,l)

       end do
    end do

    nullify(m)
    nullify(rm)
    nullify(rxm)
    nullify(rym)
    nullify(rzm)
     nullify(rxxm)
     nullify(rxym)
     nullify(rxzm)
     nullify(ryym)
     nullify(ryzm)
     nullify(rzzm)

  end subroutine red2uni_em_2nd
#endif

end module redgridZoom
