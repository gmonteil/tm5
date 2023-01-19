!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module geometry

  !*****************************
  !* geometry related routines *
  !*****************************
  !
  implicit none

  private

  public :: geomtryh, geomtryv, calc_dxy, calc_dxy0505

contains
  
     
  subroutine geomtryh(region)
    !--------------------------------------------------------------------
    ! construct arrays of horizontal model geometry
    !
    !                             mh, 27-jun-1989
    !                             mh, 26-sep-1992
    !                             v 9.1
    ! Changed vertical definition to hybrid levels according to v 9.1 in new
    ! subroutine geomtryv
    !
    !                                       aj, 23-may-1995
    !                                       mk, 5-nov-1999 zoom version
    !--------------------------------------------------------------------

    use binas,       only: ae, pi
    use dims,        only: dx, gtor, xref, dy, yref, im, jm, ybeg, areag
    use global_data, only: region_dat

    implicit none

    ! input/output
    integer,intent(in) :: region

    ! local
    real,dimension(:),pointer :: dxyp
    integer :: j
    real    :: dxx,dyy,lat,area

    ! start

    dxyp => region_dat(region)%dxyp

    !-------------------------------------------
    ! Standard horizontal geometry parent region (always global)
    !-------------------------------------------

    dxx = dx*gtor/xref(region)
    dyy = dy*gtor/yref(region)
    lat = ybeg(region)*gtor
    area =0.0
    do j=1,jm(region)
       dxyp(j) = dxx * (sin(lat+dyy)-sin(lat))*ae**2
       lat = lat+dyy
       area = area + dxyp(j)*im(region)
    end do

    print *,'geomtryh: global area ',4*pi*ae**2
    print *,'geomtryh: area region ',region,area

    areag(region) = area

    nullify(dxyp)

  end subroutine geomtryh



  subroutine calc_dxy(dxy,nlat)
    !
    ! compute area of the model grid cells
    !

    use binas,  only : ae, pi
    use dims,   only : gtor, nlon360

    implicit none

    integer,intent(in)   :: nlat
    real,dimension(nlat),intent(out) :: dxy

    real                 :: dxx,dyy,lat,areag2
    integer              :: j
    dxx = 1.0*gtor
    dyy = 1.0*gtor
    lat = -90.0*gtor
    areag2 =0.0
    do j=1,nlat
       dxy(j) = dxx * (sin(lat+dyy)-sin(lat))*ae**2
       lat = lat+dyy
       areag2 = areag2 + dxy(j)*nlon360
    enddo
    print *,'calc_dxy: global area                   ',4*pi*ae**2
    print *,'calc_dxy: area calculated for 1x1 grid..',areag2

  end subroutine calc_dxy

  subroutine calc_dxy0505(dxy,nlat)
    !
    ! compute area of the model grid cells
    !

    use binas,  only : ae, pi
    use dims,   only : gtor, nlon360

    implicit none

    integer,intent(in)   :: nlat
    real,dimension(nlat),intent(out) :: dxy

    real                 :: dxx,dyy,lat,areag2
    integer              :: j
    dxx = 0.5*gtor
    dyy = 0.5*gtor
    lat = -90.0*gtor
    areag2 =0.0
    do j=1,nlat
       dxy(j) = dxx * (sin(lat+dyy)-sin(lat))*ae**2
       lat = lat+dyy
       areag2 = areag2 + dxy(j)*nlon360*2
    enddo
    print *,'calc_dxy: global area                       ',4*pi*ae**2
    print *,'calc_dxy: area calculated for 0.5x0.5 grid..',areag2

  end subroutine calc_dxy0505




  subroutine geomtryv(region)
    !--------------------------------------------------------------------
    !  define the vertical geometry of the tm model grid v9.1knmi
    !                                   aj, 30-8-1995
    !--------------------------------------------------------------------

    use binas,       only: grav
    use dims,        only: im, jm, lm, at, bt
    use global_data, only: mass_dat, region_dat
    use MeteoData  , only : sp_dat, phlb_dat, m_dat
    implicit none

    ! input/output
    integer,intent(in) :: region

    ! local
    real,dimension(:,:,:),pointer    :: m
    real,dimension(:,:,:),pointer    :: phlb
    real,dimension(:,:,:),pointer    :: ps
    real,dimension(:),pointer        :: dxyp
    integer                          :: i,j,l

    ! start

    !c hybrid coordinate coefficients specifying 19 TM vertical boundaries
    !c (from bottom to top), extracted from the 19 hybrid level coefficients
    !c of ECHAM
    ! coefficient a is given in [Pa]
    !
    !mk    in the new system, the array m is updated by the advection routine.
    !mk    the surface pressure is compared to the surface pressure 
    !mk    stored on disk in routine rwind.    

    phlb => phlb_dat(region)%data
    m    =>    m_dat(region)%data
    ps   =>   sp_dat(region)%data
    dxyp => region_dat(region)%dxyp

    do l=1,lm(region)+1
       do j=1,jm(region)
          do  i=1,im(region)
             phlb(i,j,l) = at(l)+bt(l)*ps(i,j,1)
          end do
       end do
    end do

    !
    !cmk ----
    !     compute m (kg), the mass of air in each box.  (at the poles, m
    !     is the air mass of a full cylindrical grid box. This same mass is
    !     placed in every cell for j=1 or j=jm)
    !----

    do l=1,lm(region)
       do j=1,jm(region)
          do  i=1,im(region)
             m(i,j,l)=(phlb(i,j,l)-phlb(i,j,l+1))*dxyp(j)/grav
          end do
       end do
    end do

  end subroutine geomtryv


end module geometry



