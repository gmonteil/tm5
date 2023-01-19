module Phys_Convec_Clouds

  implicit none
  
  
  ! --- in/out ----------------------------------
  
  private
  
  public  ::  ConvCloudDim
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'module Phys_Convec_Clouds'


contains


  ! ==============================================================
  ! ===
  ! === convective clouds
  ! ===
  ! ==============================================================
  
  !
  ! updo : level order
  !  'u'   :  upwards   :  1=surface, .., n=top
  !  'd'   :  downwards :  1=top, ..., n=surface
  !
  
  subroutine ConvCloudDim( updo, lm, detu, entd, &
                           iclbas, ictop, icllfs, &
                           status )

    ! --- in/out ------------------------

    character(len=1), intent(in)   ::  updo
    integer, intent(in)            ::  lm

    real, intent(in)              ::  detu(lm)
    real, intent(in)              ::  entd(lm)
    
    ! cloud base, top, level of free sinking
    integer, intent(out)          ::  iclbas
    integer, intent(out)          ::  ictop
    integer, intent(out)          ::  icllfs
 
    integer, intent(out)          ::  status
    
    ! --- const ----------------------------------------
    
    character(len=*), parameter  ::  rname = mname//', ConvCloudDim'

    ! --- local ----------------------------------------

    integer         ::  l
    integer         ::  bot, top, one

    ! --- begin -----------------------------------------
    
    select case ( updo )
      case ( 'u', 'U' )
        bot  = 1
        top  = lm
        one = +1
      case ( 'd', 'D' )
        bot  = lm
        top  = 1
        one = -1
      case default
        write (*,'("ERROR - updo should be `u` or `d` ...")')
        write (*,'("ERROR in ",a)') rname; status=1; return
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
    ictop = 0
    do l = top, bot, -one
      if ( detu(l) > 0.0 ) then
        ictop = l
        exit
      end if
    end do


    ! determine cloud base level
    ! (cloud base level is the lowest TM model level where detrainment
    ! is greater than 0)
    ! no cloud base present by default:
    iclbas = 0
    do l = bot, top, one
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
    do l = top, bot, -one
      if ( entd(l) > 0.0 ) then
        icllfs = l
        exit
      end if
    end do
    
    ! ok
    status = 0

  end subroutine ConvCloudDim
  
  


end module Phys_Convec_Clouds


