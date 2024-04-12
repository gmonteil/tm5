!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module TMPhys_Convec

  implicit none


  ! --- in/out ----------------------------------

  private

  public  ::  ConvCloudDim


  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'module TMPhys_Convec'


contains


  ! ==================================================

  !
  ! updo : level order
  !  'u'   :  upwards   :  1=surface, .., n=top
  !  'd'   :  downwards :  1=top, ..., n=surface
  !

  subroutine ConvCloudDim( detu, entd, updo, &
                           iclbas, icltop, icllfs, &
                           status )

    ! --- in/out ------------------------

    ! 3-D fields on TM grid
    real, intent(in)              ::  detu(:)
    real, intent(in)              ::  entd(:)

    character(len=1), intent(in)  ::  updo

    ! cloud base, top, level of free sinking
    integer, intent(out)          ::  iclbas
    integer, intent(out)          ::  icltop
    integer, intent(out)          ::  icllfs

    integer, intent(out)          ::  status

    ! --- const ----------------------------------------

    character(len=*), parameter  ::  rname = mname//', ConvCloudDim'

    ! --- local ----------------------------------------

    integer         ::  l, lm
    integer         ::  lbot, ltop, ldir
    integer         ::  nocloud

    ! --- begin -----------------------------------------

    ! number of levels
    lm = size(detu)

    select case ( updo )
      case ( 'u', 'U' )
        lbot = 1
        ltop = lm
        ldir = +1
      case ( 'd', 'D' )
        lbot = lm
        ltop = 1
        ldir = -1
    end select

    ! check ...
    if ( size(entd) /= lm ) then
      write (*,'("ERROR - input arrays should have save size:")')
      write (*,'("ERROR -   size(detu) : ",i3)') size(detu)
      write (*,'("ERROR -   size(entd) : ",i3)') size(entd)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! determine cloud top level
    ! (cloud top level is the highest TM model level where detrainment
    ! is greater than 0)
    ! no cloud top present by default:
    icltop = 0
    do l = ltop, lbot, -ldir
      if ( detu(l) > 0.0 ) then
        icltop = l
        exit
      end if
    end do


    ! determine cloud base level
    ! (cloud base level is the lowest TM model level where detrainment
    ! is greater than 0)
    ! no cloud base present by default:
    iclbas = 0
    do l = lbot, ltop, ldir
      if ( detu(l) > 0.0 ) then
        iclbas = l
        exit
      end if
    end do

    ! determine level of free sinking (start of cumulus downdraft)
    ! (level of free sinking is the highest TM model level where
    ! entrainment (downdraft) is greater than 0)
    ! no cumulus downdraft present by default
    icllfs = 0
    do l = ltop, lbot, -ldir
      if ( entd(l) > 0.0 ) then
        icllfs = l
        exit
      end if
    end do

    ! ok
    status = 0

  end subroutine ConvCloudDim



end module TMPhys_Convec


