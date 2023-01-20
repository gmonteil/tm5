!###############################################################################
!
! Main program from TM5 4D-var.
!
!### macro's ###################################################################
!
#include "tm5.inc"
!
!###############################################################################

program Main

  use GO      , only : goExit
  use TM5var4D, only : T_TM5var4D
  use TM5var4D, only : TM5var4D_Init, TM5var4D_Run, TM5var4D_Done
  
  ! --- local ------------------------------------
  
  type(T_TM5var4D)   ::  tm5var4d_data
  integer            ::  status
    
  ! --- const ------------------------------------

  character(len=*), parameter   ::  rname = 'Main'
  
  ! --- begin ------------------------------------
  
  ! initialise:
  call TM5var4D_Init( tm5var4d_data, status )
  if (status/=0) then
    write (*,'("ERROR - non-zero return status from TM5var4D_Init : ",i6)') status
    write (*,'("ERROR - in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__
    status=1; call goExit(1)
  end if
  
  ! main routine:
  call TM5var4D_Run( tm5var4d_data, status )
  if (status/=0) then
    write (*,'("ERROR - non-zero return status from TM5var4D_Run : ", i6)') status
    write (*,'("ERROR - in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__
    status=1; call goExit(1)
  end if
  
  ! done with everything:
  call TM5var4D_Done( tm5var4d_data, status )
  if (status/=0) then
    write (*,'("ERROR - non-zero return status from  TM5var4D_Done : ", i6)') status
    write (*,'("ERROR - in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__
    status=1; call goExit(1)
  end if

end program Main

