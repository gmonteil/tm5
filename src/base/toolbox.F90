!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module toolbox

  ! --- modules ------------------------------

  use GO,         only : gol, goErr, goPr

  implicit none

  private

  public :: zfarr
  public :: ltropo, lvlpress, print_totalmass
  public :: coarsen_emission
  public :: escape_tm
  public :: distribute_emis2D
  public :: distribute1x1
  public :: distribute1x1b
  public :: iter_2d_s, iter_2d_m, setup_iterator_2d

  interface coarsen_emission
     module procedure coarsen_emission_1d
     module procedure coarsen_emission_2d
     module procedure coarsen_emission_3d
     module procedure coarsen_emission_d23
  end interface

  type Tij_iterator
    integer                 :: n_iter
    integer, allocatable   :: iter(:,:)
  end type Tij_iterator

  type(Tij_iterator), allocatable   :: iter_2d_s(:), iter_2d_m(:)

  character(len=*), parameter   :: mname='toolbox'

contains

  subroutine setup_iterator_2d
    !-----------------------------------------------------------------------
    !
    ! Sometimes we need to iterate over i=isr(region),ier(region), j=jsr(region),jer(region). If we want to
    ! parallelize such a double loop, it is useful to have a single iterator. This subroutine creates such
    ! single iterators, two for each region, (a) from isr to ier, (b) from 1 to im.
    !
    !-----------------------------------------------------------------------
    use dims,   only : im, jm, nregions
    use dims,   only : isr, ier, jsr, jer

    implicit none

    integer         :: region, i, j, i_iter, n_iter

    allocate(iter_2d_s(nregions))
    allocate(iter_2d_m(nregions))

    do region = 1, nregions

        ! first the isr --> ier pairs
        n_iter = (ier(region) - isr(region) + 1) * (jer(region) - jsr(region) + 1)
        allocate(iter_2d_s(region)%iter(n_iter,2))
        iter_2d_s(region)%n_iter = n_iter

        i_iter = 1
        do i = isr(region), ier(region)
            do j = jsr(region), jer(region)
                iter_2d_s(region)%iter(i_iter, 1) = i
                iter_2d_s(region)%iter(i_iter, 2) = j
                i_iter = i_iter + 1
            end do
        end do
        if (i_iter-1 /= n_iter) then
            write(0,'("ERROR :: setting up iter_2d_s for region ", i1)') region
            write(0,'("ERROR :: conflicting values ", i6, " and ", i6, " for total number of iterations")') i_iter-1, n_iter
            call escape_tm('Run!')
        end if

        ! now the 1 --> im pairs
        n_iter = im(region) * jm(region)
        allocate(iter_2d_m(region)%iter(n_iter,2))
        iter_2d_m(region)%n_iter = n_iter

        i_iter = 1
        do i = 1, im(region)
            do j = 1, jm(region)
                iter_2d_m(region)%iter(i_iter, 1) = i
                iter_2d_m(region)%iter(i_iter, 2) = j
                i_iter = i_iter + 1
            end do
        end do
        if (i_iter-1 /= n_iter) then
            write(0,'("ERROR :: setting up iter_2d_m for region ", i1)') region
            write(0,'("ERROR :: conflicting values ", i6, " and ", i6, " for total number of iterations")') i_iter-1, n_iter
            call escape_tm('Run!')
        end if

    end do ! region

  end subroutine setup_iterator_2d

  integer function ltropo(region,tx,gphx,lmx)
    !-----------------------------------------------------------------------
    !
    !****   tropo determine tropopause level
    !
    !   programmed by        fd 01-11-1996
    !       changed to function      fd 12-07-2002
    !
    !   purpose
    !   -------
    !
    !   interface
    !   ---------
    !   i=ltropo(region,tx,gphx,lmx)
    !
    !   method
    !   ------
    !   determine tropopause level
    !       from temperature gradient
    !
    !   externals
    !   ---------
    !       function lvlpress
    !
    !   reference
    !   ---------
    !   WMO /Bram Bregman
    !----------------------------------------------------------------------
    implicit none

    ! in/out
    integer,intent(in)               :: region
    integer,intent(in)               :: lmx
    real,dimension(lmx),intent(in)   :: tx
    real,dimension(lmx+1),intent(in) :: gphx

    ! local
    integer :: ltmin,ltmax,l
    real    :: dz,dt

    !   ltropo is highest tropospheric level
    !   above is defined as stratosphere

    ! min tropopause level 450 hPa (ca. 8 km)
    ltmin=lvlpress(region,45000.,98400.)
    ! max tropopause level 70 hPa (ca. 20 km)
    ltmax=lvlpress(region,7000.,98400.)
    ltropo=ltmin
    do l=ltmin,ltmax
       dz=(gphx(l)-gphx(l-2))/2.
       dt=(tx(l)-tx(l-1))/dz
       if ( dt < -0.002) ltropo=l !wmo 85 criterium for stratosphere
       ! increase upper tropospheric level
    end do !l

  end function ltropo



  integer function lvlpress(region,press,press0)
    !-----------------------------------------------------------------------
    !
    !****   lvlpress determines the index of the level with a pressure
    !       greater than press i.e. in height below press
    !
    !       programmed by           Peter van Velthoven, 6-november-2000
    !
    !       purpose
    !       -------
    !       determines level independent of vertical resolution lm
    !
    !       interface
    !       ---------
    !       i = lvlpress(press,press0)
    !
    !       method
    !       ------
    !       use hybrid level coefficients to determine lvlpress
    !
    !       input
    !       ---------
    !       press  : pressure in Pa
    !       press0 : surface pressure in Pa
    !
    !----------------------------------------------------------------------
    use dims, only: at,bt,lm
    implicit none

    ! in/out
    integer, intent(in) :: region
    real, intent(in)    :: press
    real, intent(in)    :: press0

    ! local
    real             :: zpress0, zpress
    integer          :: l,lvl

    if ( press0 == 0.0 ) then
       zpress0 = 98400.
    else
       zpress0 = press0
    endif
    lvl = 1
    ! l increases from bottom (l=1) to top (l=lm)
    do l = 1, lm(region)
       zpress =  (at(l)+at(l+1) + (bt(l)+bt(l+1))*zpress0)*0.5
       if ( press < zpress ) then
          lvl = l
       endif
    end do
    lvlpress = lvl

  end function lvlpress



!  subroutine calc_region_coordinates(region,xlat,xlon)
!    !
!    !
!    !
!    use dims, only: im, jm, dx, dy, xref, yref, xbeg, ybeg
!    implicit none
!
!    ! in/out
!    integer, intent(in)                       :: region
!    real, dimension(im(region)+1),intent(out) :: xlon
!    real, dimension(jm(region)+1),intent(out) :: xlat
!
!    ! local
!    integer :: i,j,l,n
!    real    :: dx1,dy1   ! gridsize in degrees longitude and latitude
!
!    dx1=dx/xref(region) !longitude
!    dy1=dy/yref(region) !latitude
!
!
!    do j=1,jm(region)+1
!       ! southern border,last value is North Border
!       xlat(j)=ybeg(region)+dy1*(j-1)
!    enddo
!
!    do i=1,im(region)+1
!       ! western border, last value is East Border
!       xlon(i)=xbeg(region)+dx1*(i-1)
!    enddo
!
!  end subroutine calc_region_coordinates



  real function zfarr(rx1,er,ztrec)
    !------------------------------------------------------------------
    !
    !****  ZFARR calculation of Arrhenius expression for rate constants
    !
    !------------------------------------------------------------------
    !
    implicit none
    real,intent(in) :: rx1,er,ztrec
    !
    zfarr=rx1*exp(er*ztrec)
    !
  end function zfarr



  subroutine print_totalmass(region,tracer)

    use dims,           only : im, jm, lm, region_name
    use global_data,    only : mass_dat
    use chem_param,     only : emis_unit, emis_unit_name, names

    implicit none

    ! in/out
    integer,intent(in)            :: region, tracer

    ! local
    integer       :: i
    real          :: mass
    character(len=*), parameter :: rname = mname//'/print_totalmass'

    mass = emis_unit(tracer) * sum(mass_dat(region)%rm_t(1:im(region),1:jm(region),1:lm(region),tracer))
    write(gol, '(a, ": Total mass of tracer ", a, " in region ", a, " = ", es11.4, " ", a)') rname, trim(names(tracer)), trim(region_name(region)), mass, trim(emis_unit_name(tracer))
    call goPr

  end subroutine print_totalmass


  subroutine coarsen_emission_1d(name_field,jm_emis,field_in,field,avg)
    !
    ! Purpose: Transform the 1D field available on e.g. 1 degree resolution
    !          to the desired zoom geometry
    ! name_field : name, only for printing reasons
    ! jm_emis    : the dimension of the input field
    ! field_in   : the input field
    ! field      : type d2_data  (defined in chem_param)
    ! avg        : flags the mode: avg = 1 means that a surface area
    !              weighted area is calulated (e.g. land fraction)
    !              avg = 0 means that the sum of the underlying field
    !              is taken (e.g. emissions in kg/month).
    !
    use dims,       only: nregions, dy, yref, ybeg, jm, dxy11, idate
    use chem_param, only: d2_data
    implicit none

    ! in/out
    character(len=*),intent(in)               :: name_field
    integer,intent(in)                        :: avg  ! 0=no 1=yes
    integer,intent(in)                        :: jm_emis
    real,dimension(jm_emis),intent(in)        :: field_in
    type(d2_data),dimension(nregions),target  :: field
    ! contains per region field d2(jm)

    ! local
    real,dimension(:),pointer    :: coarse_field
    integer :: region
    real    :: scale,w,wtot
    integer :: ystart
    integer :: yres
    integer :: j,j_index,jj
    integer :: jm_start

    ! start

    do region=1,nregions      ! main loop over regions

       yres=dy/yref(region)

       if ( yres < 1 ) then
          call escape_tm('coarsen_emission_1d:'//&
               ' 1 degree minimum resolution in emissions')
       end if
       if ( jm_emis /= 180 ) then
          call escape_tm( &
               'coarsen_emission_1d: input resolution should be 1 degree')
       end if

       !WP! find index of southmost gridpoint based on 1x1 degree
       jm_start=ybeg(region)+91

       if ( jm_start <= -1 ) then
          call escape_tm('coarsen_emission_1d:'//&
               ' requested emission fields outside valid region')
       end if

       coarse_field => field(region)%d2
       do j=1,jm(region)
          !cycle through 1x1 grid with steps for this region
          j_index=jm_start+(j-1)*yres
          coarse_field(j) = 0.0
          wtot = 0.0
          do jj=0,yres-1
             if(avg==1) then
                w = dxy11(j_index+jj)
                wtot = wtot + w
                coarse_field(j)=coarse_field(j) + field_in(j_index+jj)*w
             else
                coarse_field(j)=coarse_field(j) + field_in(j_index+jj)
             end if
          end do
          if ( avg == 1 )  coarse_field(j) = coarse_field(j)/wtot
       end do

       if ( avg == 0 ) then
          write(*,'(a,i1,x,a12,a,i2,a,1pe14.3)') &
               ' coarsen_emission_1d: region ',region,name_field, &
               ' Field coarsened month :',idate(2),'total:',&
               sum(coarse_field)
       else
          write(*,'(a,i1,x,a12,a,i2,a,1pe14.3)') &
               ' coarsen_emission_1d: region:',region,name_field, &
               ' Field averaged month :',idate(2)
       end if

       nullify(coarse_field)

    end do ! loop over regions....

  end subroutine coarsen_emission_1d



  !CMK
  subroutine coarsen_emission_2d(name_field,im_emis,jm_emis,field_in,field,avg)
    !
    ! Purpose: Transform the 2D emissions available on
    !          e.g. 1x1 degree resolution
    !          to the desired zoom geometry
    ! name_field : name, only for printing reasons
    ! im_emis    : the x-dimension of the input field
    ! jm_emis    : the y-dimension of the input field
    ! field_in   : the input field
    ! field      : type emis_data  (defined in chem_param)
    ! avg        : flags the mode: avg = 1 means that a surface area
    !              weighted area is calulated (e.g. land fraction)
    !              avg = 0 means that the sum of the underlying field
    !              is taken (e.g. emissions in kg/month).
    !
    use dims,       only: nregions, dx, dy, xref, yref, xbeg, ybeg
    use dims,       only: im, jm, dxy11, idate
    use chem_param, only: emis_data
    implicit none

    ! in/out
    character(len=*),intent(in)                :: name_field
    integer,intent(in)                         :: avg  !0=no 1=yes
    integer,intent(in)                         :: im_emis,jm_emis
    real,dimension(im_emis,jm_emis),intent(in) :: field_in
    type(emis_data),dimension(nregions),target :: field
    !                  contains per region field surf(im,jm)

    ! local
    real,dimension(:,:),pointer    :: coarse_field
    integer :: region
    real    :: scale,w,wtot
    integer :: xstart,ystart
    integer :: xres,yres
    integer :: i,j,i_index,j_index,ii,jj
    integer :: im_start,jm_start

    ! start

    do region=1,nregions      ! main loop over regions

       xres=dx/xref(region)
       yres=dy/yref(region)

       if ( xres < 1 .or. yres < 1 ) then
          call escape_tm( &
               'coarsen_emission_2d: 1 degree minimum resolution in emissions')
       end if
       if ( im_emis /= 360 .or. jm_emis /= 180 ) then
          call escape_tm('coarsen_emission_2d: input resolution should be 1x1')
       end if

       !WP! find index of westmost gridpoint based on 1x1 degree
       im_start=xbeg(region)+181
       !WP! find index of southmost gridpoint based on 1x1 degree
       jm_start=ybeg(region)+91

       if( im_start <= -1 .or. jm_start <= -1 ) then
          call escape_tm( 'coarsen_emission_2d: '// &
               'requested emission fields outside valid region')
       end if

       coarse_field => field(region)%surf

       do j=1,jm(region)
          do i=1,im(region)
             !cycle through 1x1 grid with steps for this region
             i_index=im_start+(i-1)*xres
             !cycle through 1x1 grid with steps for this region
             j_index=jm_start+(j-1)*yres
             coarse_field(i,j) = 0.0
             wtot = 0.0
             do ii=0,xres-1
                do jj=0,yres-1
                   if(avg==1) then
                      w = dxy11(j_index+jj)
                      wtot = wtot + w
                      coarse_field(i,j) = coarse_field(i,j) + &
                           field_in(i_index+ii,j_index+jj)*w
                   else
                      coarse_field(i,j) = coarse_field(i,j) + &
                           field_in(i_index+ii,j_index+jj)
                   end if
                end do
             end do
             if ( avg == 1 )  coarse_field(i,j) = coarse_field(i,j)/wtot
          end do
       end do

       !if ( avg == 0 ) then
       !   write(*,'(a,i1,x,a12,a,i2,a,1pe14.3)') &
       !        ' coarsen_emission_2d: region ',region,name_field, &
       !        ' Field coarsened month ',idate(2),' total ',&
       !        sum(coarse_field)
       !else
       !   write(*,'(a,i1,x,a12,a,i2,a,1pe14.3)') &
       !        ' coarsen_emission_2d: region ',region,name_field, &
       !        ' Field averaged month ',idate(2)
       !end if

       nullify(coarse_field)

    end do ! loop over regions....

  end subroutine coarsen_emission_2d



  !CMK
  subroutine coarsen_emission_3d( name_field, im_emis, jm_emis, lm_emis, &
                                  field_in, field, avg )
    !
    ! Purpose: Transform the 3D emissions available on
    !          e.g. 1x1 degree resolution
    !          to the desired zoom geometry
    ! name_field : name, only for printing reasons
    ! im_emis    : the x-dimension of the input field
    ! jm_emis    : the y-dimension of the input field
    ! lm_emis    : the z-dimension of the input field
    ! field_in   : the input field
    ! field      : type d3_data  (defined in chem_param)
    ! avg        : flags the mode: avg = 1 means that a surface area
    !              weighted area is calulated (e.g. land fraction)
    !              avg = 0 means that the sum of the underlying field
    !              is taken (e.g. emissions in kg/month).
    !
    use dims,       only: nregions, dx, dy, xref, yref, xbeg, ybeg
    use dims,       only: im, jm, lm, dxy11, idate
    use chem_param, only: d3_data
    implicit none

    ! in/out
    character(len=*),intent(in)           :: name_field
    integer,intent(in)                    :: avg  !0=no 1=yes
    integer,intent(in)                    :: im_emis,jm_emis,lm_emis
    real,dimension(im_emis,jm_emis,lm_emis),intent(in) :: field_in
    type(d3_data),dimension(nregions),target           :: field
    !             contains per region 3d field d3(im,jm,lm)

    ! local
    real,dimension(:,:,:),pointer    :: coarse_field
    integer :: region
    real    :: scale,w,wtot
    integer :: xstart,ystart
    integer :: xres,yres
    integer :: i,j,i_index,j_index,ii,jj,l
    integer :: im_start,jm_start

    ! start

    do region=1,nregions      ! main loop over regions

       xres=dx/xref(region)
       yres=dy/yref(region)

       if ( xres < 1 .or. yres < 1 ) then
          call escape_tm('coarsen_emission_3d: '//&
               '1 degree minimum resolution in emissions')
       end if
       if ( im_emis /= 360 .or. jm_emis /= 180 ) then
          call escape_tm('coarsen_emission_3d: input resolution should be 1x1')
       end if

       !WP! find index of westmost gridpoint based on 1x1 degree
       im_start=xbeg(region)+181
       !WP! find index of southmost gridpoint based on 1x1 degree
       jm_start=ybeg(region)+91

       if ( im_start <= -1 .or. jm_start <= -1 ) then
          call escape_tm( 'coarsen_emission_3d: '//&
               'requested emission fields outside valid region')
       end if

       coarse_field => field(region)%d3
       do l=1,lm(region)
          do j=1,jm(region)
             do i=1,im(region)
                !cycle through 1x1 grid with steps for this region
                i_index=im_start+(i-1)*xres
                !cycle through 1x1 grid with steps for this region
                j_index=jm_start+(j-1)*yres
                coarse_field(i,j,l) = 0.0
                wtot = 0.0
                do ii=0,xres-1
                   do jj=0,yres-1
                      if(avg==1) then
                         w = dxy11(j_index+jj)
                         wtot = wtot + w
                         coarse_field(i,j,l)=coarse_field(i,j,l) + &
                              field_in(i_index+ii,j_index+jj,l)*w
                      else
                         coarse_field(i,j,l)=coarse_field(i,j,l) + &
                              field_in(i_index+ii,j_index+jj,l)
                      end if
                   end do
                end do
                if ( avg == 1 )  coarse_field(i,j,l) = coarse_field(i,j,l)/wtot
             end do
          end do
       end do  !levels

       !if (avg == 0) then
       !
       !   write(*,'(a,i1,x,a12,a,i2,a,1pe14.3)') &
       !        ' coarsen_emission_3d: region ',region,name_field, &
       !        ' Field coarsened month ',idate(2),'  total ',&
       !        sum(coarse_field)
       !else
       !   write(*,'(a,i1,x,a12,a,i2,a,1pe14.3)') &
       !        ' coarsen_emission_3d: region ',region,name_field, &
       !        ' Field averaged month ',idate(2)
       !end if

       nullify(coarse_field)

    end do ! loop over regions....

  end subroutine coarsen_emission_3d



  subroutine coarsen_emission_d23( name_field, jm_emis, lm_emis, &
                                   field_in, field, avg )
    !
    ! Purpose: Transform the 2D emissions available on lat-pressure resolution
    !          to the desired zoom geometry
    ! name_field : name, only for printing reasons
    ! jm_emis    : the y-dimension of the input field
    ! lm_emis    : the z-dimension of the input field
    ! field_in   : the input field
    ! field      : type d23_data (defined in chem_param)
    ! avg        : flags the mode: avg = 1 means that a surface area
    !              weighted area is calulated (e.g. land fraction)
    !              avg = 0 means that the sum of the underlying field
    !              is taken (e.g. emissions in kg/month).
    !
    use dims,       only: nregions, dy, yref, ybeg, jm, lm, dxy11, idate
    use chem_param, only: d23_data
    implicit none

    ! in/out
    character(len=*),intent(in)                        :: name_field
    integer,intent(in)                                 :: avg  !0=no 1=yes
    integer,intent(in)                                 :: jm_emis,lm_emis
    real,dimension(jm_emis,lm_emis),intent(in)         :: field_in
    type(d23_data),dimension(nregions),target          :: field
    !                contains per region field d23(jm,lm)

    ! local
    real,dimension(:,:),pointer    :: coarse_field
    integer :: region
    real    :: scale,w,wtot
    integer :: ystart
    integer :: yres
    integer :: j,j_index,jj,l
    integer :: jm_start

    ! start

    do region=1,nregions      ! main loop over regions

       yres=dy/yref(region)

       if ( yres < 1 ) then
          call escape_tm('coarsen_emission_d23:'//&
               ' 1 degree minimum resolution in emissions')
       end if
       if ( jm_emis /= 180 ) then
          call escape_tm( 'coarsen_emission_d23:'//&
               'input resolution should be 1 degree')
       end if

       !WP! find index of southmost gridpoint based on 1x1 degree
       jm_start=ybeg(region)+91

       if(jm_start<=(-1)) then
          call escape_tm( 'coarsen_emission_d23:'//&
               'requested emission fields outside valid region')
       end if

       coarse_field => field(region)%d23
       do l=1,lm(region)
          do j=1,jm(region)
             !cycle through 1x1 grid with steps for this region
             j_index=jm_start+(j-1)*yres
             coarse_field(j,l) = 0.0
             wtot = 0.0
             do jj=0,yres-1
                if(avg==1) then
                   w = dxy11(j_index+jj)
                   wtot = wtot + w
                   coarse_field(j,l)=coarse_field(j,l) + &
                        field_in(j_index+jj,l)*w
                else
                   coarse_field(j,l)=coarse_field(j,l) + &
                        field_in(j_index+jj,l)
                end if
             end do
             if ( avg == 1 )  coarse_field(j,l) = coarse_field(j,l)/wtot
          end do  !j
       end do  !l

       if ( avg == 0 ) then

          write(*,'(a,i1,x,a12,a,i2,a,1pe14.3)') &
               ' coarsen_emission_d23: region ',region,name_field// &
               ' Field coarsened month ',idate(2),'  total ',&
               sum(coarse_field)
       else
          write(*,'(a,i1,x,a12,a,i2,a,1pe14.3)') &
               ' coarsen_emission_d23: region ',region,name_field// &
               ' Field averaged month ',idate(2)
       end if

       nullify(coarse_field)

    end do ! loop over regions....

  end subroutine coarsen_emission_d23



  subroutine escape_tm(msg)
    !--------------------------------------------------------------
    !
    ! abnormal termination of a model run
    !
    ! msg       character string to be printed on unit kmain
    !
    !-------------------------------------------------------------
    use dims,        only : okdebug, kmain, itau
    use datetime,    only : tstamp
#ifdef MPI
    use mpi_const,   only : myid, root
    use mpi_comm,    only : abortmpi
#endif
    implicit none
    !
    character(*),intent(in) :: msg
    !
#ifdef MPI
    if ( myid == root ) then
#endif
       write(kmain,'(/2(1x,39("--")/))')
       call tstamp(kmain,itau,'ERROR - ESCAPE_TM')
       write(kmain,'(1x,a)') msg
       write(kmain,'(/2(1x,39("--")/))')
#ifdef MPI
    endif
#endif

#ifdef MPI
    call abortmpi
#else
    stop
#endif

  end subroutine escape_tm

  subroutine distribute_emis2D(emis2D, emis3D, hlow, hhigh, xfrac)
  ! subroutine to distribute the emissions given as a TM5 2D field
  ! into a TM5 3D emission field. Hlow and Hhigh are the bounds of
  ! the 2D emission fields. They give the height RELATIVE to oro
  ! From that the distribution is calculated
  ! employing the geopotential height (relative to oro)
  use dims,        only: lm, nregions, im, jm
  use MeteoData  , only : gph_dat
  use chem_param,  only: d3_data, emis_data
  implicit none
  !______________________IO______________________________________________________
  type(emis_data),intent(in),dimension(nregions),target   :: emis2D    ! 2D emission field (kg/gridbox/month)
  type(d3_data),intent(inout),dimension(nregions),target  :: emis3D    ! 3D emission field (kg/gridbox/month)
  real,intent(in)                       :: hlow      ! lowest level (m)
  real,intent(in)                       :: hhigh     ! highest level (m)
  real,intent(in), optional             :: xfrac     ! fraction of emissions to put
  !______________________local___________________________________________________
  real, dimension(:,:,:),pointer        :: gph       ! geopotential height (m)
  real, dimension(:,:,:),pointer        :: e3d       ! 3D emissions
  real, dimension(:,:),pointer          :: e2d       ! 2D emissions
  integer                               :: region, i,j,l
  real, dimension(lm(1)+1)              :: height
  real, dimension(lm(1))                :: dz
  real                                  :: dze
  real                                  :: totw, f, tot2d, tot3db, tot3da, fraction
  integer, dimension(3)                 :: ubound_e3d
  integer                               :: lmmax
  real                                  :: hhighb
  !______________________start___________________________________________________

  if (present(xfrac)) then
     fraction = xfrac
  else
     fraction = 1.0
  endif
  if ( hhigh <= 0.0) call escape_tm(' Routine distribute_emis2D hhigh <= 0')
  if ( hlow   < 0.0) call escape_tm(' Routine distribute_emis2D hlow < 0')
  dze = hhigh - hlow
  if ( dze <= 0.0) call escape_tm(' Routine distribute_emis2D hhigh-hlow <= 0')
  do region = 1, nregions
  nullify(gph)
  nullify(e2d)
  nullify(e3d)
  gph => gph_dat(region)%data
  e2d => emis2d(region)%surf
  e3d => emis3d(region)%d3
  ubound_e3d = ubound(e3d)
  lmmax = ubound_e3d(3)
  tot2d = sum(e2d*fraction)
  tot3db = sum(e3d)
  do j=1,jm(region)
     do i=1,im(region)
        height(1) = 0.0
        do l=1, lm(region)
           dz(l) = gph(i,j,l+1) - gph(i,j,l)
           height(l+1) = height(l) + dz(l)
        enddo
        totw = 0.0
        if(hhigh > height(lmmax+1) ) then
           print *, 'hhigh, height(lmmax+1)', hhigh, height(lmmax+1)
           call escape_tm('distribute Emis2D: try to put emissions higher than array allows')
        endif
        zz: do l=1, lmmax
            if (hhigh > height(l+1)) then
               if ( hlow < height(l) ) then
                  f =  dz(l)/dze
                  totw = totw + f
                  e3d(i,j,l) = e3d(i,j,l) + f*fraction*e2d(i,j)
               else if( hlow < height(l+1)) then
                  f = (height(l+1)-hlow)/dze
                  totw = totw + f
                  e3d(i,j,l) = e3d(i,j,l) + f*fraction*e2d(i,j)
               endif
            else
               if ( hlow < height(l)) then
                  f = (hhigh - height(l))/dze
                  totw = totw + f
                  e3d(i,j,l) = e3d(i,j,l) + f*fraction*e2d(i,j)
               else
                  totw = totw + 1.0
                  e3d(i,j,l) = e3d(i,j,l) + fraction*e2d(i,j)
               endif
               exit zz
            endif
        enddo zz
        if ( abs(totw-1.0) > 1e-14 ) &
             call escape_tm(' Routine distribute_emis2D : sum weight /= 1')
     enddo   !i
  enddo   !j
  tot3da = sum(e3d)
  if (abs((tot3da-tot3db)-tot2d) > 1e-8*tot2d ) &
      call escape_tm(' Routine distribute_emis2D : emissions have not been distributed mass-conserving')
  nullify(gph)
  nullify(e2d)
  nullify(e3d)
  enddo ! regions
  end subroutine distribute_emis2D

  subroutine distribute1x1(emi1x1, hlow, hhigh, emis3d, xfrac)
  !
  ! subroutine to distribute 1x1 emissions liniairly between
  ! hlow and hhigh. The vertical level is determined by
  ! the orography which is read from the surface file...
  ! A simple scale height vertical structure is assumed.
  !
  use dims,          only: at, bt, nlon360, nlat180, lm, itau
  use Binas,         only: grav
  use Dims         , only : iglbsfc
  use MeteoData    , only : oro_dat
  implicit none
  !____________________________________IO_____________________________________________________________________
  real, dimension(nlon360,nlat180), intent(in)            :: emi1x1  ! (kg/1x1 gridbox) 2D field of emissions
  real, dimension(nlon360,nlat180), intent(in)            :: hlow    ! (m) lower bound of emission
  real, dimension(nlon360,nlat180), intent(in)            :: hhigh   ! (m) higher bound of emission
  real, dimension(nlon360,nlat180,lm(1)), intent(inout)   :: emis3d  ! (kg/box) distributed in height
  real,intent(in), optional                               :: xfrac   ! fraction of emissions to put
  !____________________________________locals_________________________________________________________________
  real, parameter                                :: scalh = 8000.0
  real                                           :: p0, pt
  real,dimension(lm(1))                          :: height, dz
  integer                                        :: i,j,l
  real                                           :: hh,hl,dze,totw,f,e3da,e3db, fract
  !____________________________________start_________________________________________________________________
  if (present(xfrac)) then
     fract = xfrac
  else
     fract = 1.0
  endif
  e3db = sum(emis3d)
  do j=1,nlat180
    do i=1,nlon360
       if (fract*emi1x1(i,j) > 1e-14) then
          height(1) = oro_dat(iglbsfc)%data(i,j,1) / grav
          p0 = 1e5*exp(-height(1)/scalh)
          do l=1,lm(1)-1   ! bug reported by FD: alog(0) crashes!
             pt = at(l+1) + bt(l+1)*p0
             height(l+1) = height(1)-scalh*alog(pt/p0)
             dz(l) = height(l+1)-height(l)
          enddo
          hl = max(height(1),hlow(i,j))
          hh = max(hhigh(i,j),height(1))
          dze = hh-hl
          if(dze < 0.0) then
             call escape_tm('In distribute1x1: dze <= 0')
          else if ( dze == 0.0) then   ! this somehow happens!
             hh = height(1)+1.0
             hl = height(1)
          endif
          totw = 0.0
          zz: do l=1, lm(1)
              if (hh > height(l+1)) then
                 if ( hl < height(l) ) then
                    f =  dz(l)/dze
                    totw = totw + f
                    emis3d(i,j,l) = emis3d(i,j,l) + f*fract*emi1x1(i,j)
                 else if( hl < height(l+1)) then
                    f = (height(l+1)-hl)/dze
                    totw = totw + f
                    emis3d(i,j,l) = emis3d(i,j,l) + f*fract*emi1x1(i,j)
                 endif
              else
                 if ( hl < height(l)) then
                    f = (hh - height(l))/dze
                    totw = totw + f
                    emis3d(i,j,l) = emis3d(i,j,l) + f*fract*emi1x1(i,j)
                 else
                    totw = totw + 1.0
                    emis3d(i,j,l) = emis3d(i,j,l) + fract*emi1x1(i,j)
                 endif
                 exit zz
              endif
          enddo zz
          if ( abs(totw-1.0) > 1e-14 ) &
               call escape_tm(' Routine distribute1x1 : sum weight /= 1')
       endif
    enddo
  enddo
  e3da = sum(emis3d)
  if (abs(e3da-e3db-sum(fract*emi1x1)) > e3da*1e-8 ) &
                call escape_tm(' Routine distribute1x1 : distributed amount differs!')
  end subroutine distribute1x1

  subroutine distribute1x1b(emi1x1, hlow, hhigh, emis3d, xfrac)
  !
  ! subroutine to distribute 1x1 emissions liniairly between
  ! hlow and hhigh. The vertical level is determined by
  ! the orography which is read from the surface file...
  ! A simple scale height vertical structure is assumed.
  ! same as distribute1x1 but hlow, hhigh now scalar
  ! ALSO: the height is now defined relative to the orography!!!
  !
  use dims,          only: at, bt, nlon360, nlat180, lm, itau
  use Binas,         only: grav
  use Dims         , only : iglbsfc
  use MeteoData    , only : oro_dat
  implicit none
  !____________________________________IO_____________________________________________________________________
  real, dimension(nlon360,nlat180), intent(in)            :: emi1x1  ! (kg/1x1 gridbox) 2D field of emissions
  real,  intent(in)                                       :: hlow    ! (m) lower bound of emission
  real,  intent(in)                                       :: hhigh   ! (m) higher bound of emission
  real, dimension(nlon360,nlat180,lm(1)), intent(inout)   :: emis3d  ! (kg/box) distributed in height
  real,intent(in), optional                               :: xfrac   ! fraction of emissions to put
  !____________________________________locals_________________________________________________________________
  real, parameter                                :: scalh = 8000.0
  real                                           :: p0, pt
  real,dimension(lm(1))                          :: height, dz
  integer                                        :: i,j,l
  real                                           :: hh,hl,dze,totw,f,e3da,e3db, fract, hlow_oro, hhigh_oro
  !____________________________________start_________________________________________________________________
  if (present(xfrac)) then
     fract = xfrac
  else
     fract = 1.0
  endif
  e3db = sum(emis3d)
  do j=1,nlat180
    do i=1,nlon360
       if (fract*emi1x1(i,j) > 1e-14) then
          height(1) = oro_dat(iglbsfc)%data(i,j,1) / grav
          hlow_oro = hlow + height(1)
          hhigh_oro = hhigh + height(1)
          p0 = 1e5*exp(-height(1)/scalh)
          do l=1,lm(1)-1   ! bug reported by FD: alog(0) crashes!
             pt = at(l+1) + bt(l+1)*p0
             height(l+1) = height(1)-scalh*alog(pt/p0)
             dz(l) = height(l+1)-height(l)
          enddo
          hl = max(height(1),hlow_oro)
          hh = max(hhigh_oro,height(1))
          dze = hh-hl
          if(dze < 0.0) then
             call escape_tm('In distribute1x1b: dze <= 0')
          else if ( dze == 0.0) then   ! this somehow happens!
             hh = height(1)+1.0
             hl = height(1)
          endif
          totw = 0.0
          zz: do l=1, lm(1)
              if (hh > height(l+1)) then
                 if ( hl < height(l) ) then
                    f =  dz(l)/dze
                    totw = totw + f
                    emis3d(i,j,l) = emis3d(i,j,l) + f*fract*emi1x1(i,j)
                 else if( hl < height(l+1)) then
                    f = (height(l+1)-hl)/dze
                    totw = totw + f
                    emis3d(i,j,l) = emis3d(i,j,l) + f*fract*emi1x1(i,j)
                 endif
              else
                 if ( hl < height(l)) then
                    f = (hh - height(l))/dze
                    totw = totw + f
                    emis3d(i,j,l) = emis3d(i,j,l) + f*fract*emi1x1(i,j)
                 else
                    totw = totw + 1.0
                    emis3d(i,j,l) = emis3d(i,j,l) + fract*emi1x1(i,j)
                 endif
                 exit zz
              endif
          enddo zz
          if ( abs(totw-1.0) > 1e-14 ) &
               call escape_tm(' Routine distribute1x1 : sum weight /= 1')
       endif
    enddo
  enddo
  e3da = sum(emis3d)
  if (abs(e3da-e3db-sum(fract*emi1x1)) > e3da*1e-8 ) &
                call escape_tm(' Routine distribute1x1 : distributed amount differs!')
  end subroutine distribute1x1b

end module toolbox
