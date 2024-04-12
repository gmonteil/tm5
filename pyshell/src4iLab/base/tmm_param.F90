!###############################################################################
!
! NAME
!   tmm_param  -  parameter specfic keys
!
! USAGE
!
!   use tmm_param
!
!   ! Character keys for vertical combination,
!   ! to be used as input for 'CopyAndProcess' routine of 'Grid' module.
!   ! Specifies how to combine multiple levels:
!   !
!   !  'bottom'     :  use bottom value (most close to the ground)
!   !  'top'        :  use top value (most close to the model top)
!   !  'sum'        :  sum values
!   !  'aver'       :  average of all levels
!   !  'mass-aver'  : mass weighted average
!   !
!   call CombineKey( combkey, 'pu'|'T'|... )
!
!
!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tmm.inc"
!
!###############################################################################

module tmm_param

  implicit none
  
  ! --- in/out ------------------------------
  
  private
  
  public   ::  GetCombineKeys
  
  ! --- const ---------------------------------------
  
  character(len=*), parameter  ::  mname = 'module tmm_param'
  
  
  
contains


  ! ==============================================================


  subroutine GetCombineKeys( hcomb, vcomb, paramkey, status )
  
    ! --- in/out ---------------------------------
    
    character(len=*), intent(out)    ::  hcomb
    character(len=*), intent(out)    ::  vcomb
    character(len=*), intent(in)     ::  paramkey
    integer, intent(out)             ::  status
    
    ! --- const --------------------------------------
    
    character(len=*), parameter  ::  name = mname//', GetCombineKeys'
    
    ! --- begin ---------------------------------
    
    select case ( paramkey )
    
      ! ~~~ tm fields ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      case ( 'sp', 'spm', 'sp_sfc' )
        hcomb = 'area-aver'
        vcomb = 'none'

      case ( 'oro' , 'lsm' , 'sr', 'srols', 'sr_ols', 'sr_mer', &
             'cvl', 'cvh', &
             'tv01', 'tv02', 'tv03', 'tv04', 'tv05', &
             'tv06', 'tv07', 'tv08', 'tv09', 'tv10', &
             'tv11', 'tv12', 'tv13', 'tv14', 'tv15', &
             'tv16', 'tv17', 'tv18', 'tv19', 'tv20', &
             'swvl1', &
             'al', 'albedo', &
             'lsrh', 'ci'  , '10fg', 'g10m', 'u10m', 'v10m', 'sd', &
             'lsp' , 'cp'  , 'sf'  , 'sshf', 'slhf', 'blh' , &
             't2m' , 'd2m' , 'ssr' , 'sstr', 'ewss', 'nsss', &
             'src' , 'raero', 'ustar', &
             'sst' , 'skt' )
        hcomb = 'area-aver'
        vcomb = 'none'

      case ( 'pu', 'pv', 'mfu', 'mfv' )
        hcomb = 'sum'
        vcomb = 'sum'

      case ( 'pw', 'mfw' )
        hcomb = 'sum'
        vcomb = 'bottom'

      case ( 'T', 'Tv', 'Q', 'PVo', 'PV', 'theta' )
        hcomb = 'mass-aver'
        vcomb = 'mass-aver'

      case ( 'CLWC', 'CIWC', 'CC', &
             'clwc', 'ciwc', 'cc' )
        hcomb = 'mass-aver'
        vcomb = 'mass-aver'

      case ( 'eu', 'du', 'ed', 'dd' )
        hcomb = 'area-aver'
        vcomb = 'sum'

      ! ecmwf convective fields:
      !  o mass fluxes in kg/m2, upper half level
      !  o detrainments in kg/m2/m
      case ( 'UDMF', 'DDMF' )
        hcomb = 'area-aver'
        vcomb = 'top'
      case ( 'UDDR', 'DDDR' )
        hcomb = 'area-aver'
        vcomb = 'height-aver'

      case ( 'u', 'v' )
        hcomb = 'aver'
        vcomb = 'aver'

      case ( 'Kz' )
        hcomb = 'none'
        vcomb = 'bottom'
        
      case ( 'cco', 'CCO' )
        hcomb = 'area-aver'
        vcomb = 'bottom'
        
      case ( 'ccu', 'CCU' )
        hcomb = 'area-aver'
        vcomb = 'top'
        
      ! ~~~ sh/gg fields ~~~~~~~~~~~~~~~~~~~~~
      
      case ( 'lnsp' )
        hcomb = 'exp,aver'
        vcomb = 'none'
      
      ! ~~~ dummy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      case default
        
        !write (*,'("ERROR - do not know how to combine levels for parameter `",a,"`")') paramkey
        !write (*,'("ERROR in ",a)') name; status=1; return
        
        hcomb = 'unknown'
        vcomb = 'unknown'
        
    end select
    
    ! ok
    status = 0

  end subroutine GetCombineKeys

  
end module tmm_param
