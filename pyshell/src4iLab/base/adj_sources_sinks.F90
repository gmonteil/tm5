!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module adj_sources_sinks

  use GO, only : gol, goPr, goErr

  implicit none

  private

  ! --- const ------------------------------------

  character(len=*), parameter   ::  mname = 'adj_sources_sinks'

  !public routines
  public :: adj_trace0
  public :: adj_trace1
  public :: adj_trace_after_read
  public :: adj_source1
  public :: adj_source2
  public :: adj_trace_end

contains

  !=====================================================================================================
  !=====================================================================================================

  subroutine adj_trace0

    use toolbox,                  only : escape_tm

!    use adj_chemistry,            only : adj_chemistry_init
!#endif
!#ifndef without_dry_deposition
!    use dry_deposition,           only : dry_deposition_init
!#endif
    use dims,                     only : newsrun
    implicit none

    !__IO___________________________________________________________________

    !__CONST________________________________________________________________

    character(len=*), parameter       :: rname = mname//'/adj_trace0'

    !__LOCAL_VARIABLES______________________________________________________

    integer                           :: status

    !__START_SUBROUTINE______________________________________________________

    !=======================
    ! surface fields
    !=======================

!    PB: now from initexit
!    ! fields for surface processes (1x1 --> coarsened)
!    if ( newsrun ) call declare_surface_fields

    !=======================
    ! adjoint emissions
    ! JFM: this is now called from initexit: start_TM5
    !  because it depends on run_mode and iteration whether this has to be done.
    !=======================
!    if ( newsrun ) call adj_emission_init



!    !=======================
!    ! chemistry (OH sink)
!    !=======================
!    if (newsrun) then
!        call adj_chemistry_init( status )
!        IF_NOTOK_RETURN(status=1)
!    end if
!#endif
!
!
!#ifndef without_dry_deposition
!    call dry_deposition_init( status )
!    if ( status /= 0 ) call escape_tm('adj_trace0: initializing dry_deposition failed')
!#endif
!    <<<

  end subroutine adj_trace0


  !=====================================================================================================
  !=====================================================================================================

  subroutine adj_trace1

    !
    ! Initialize adjoint tracer masses
    !

    use dims,        only : im, jm, lm
    use dims,        only : adv_scheme
    use dims,        only : nregions
    use global_data, only : mass_dat
    use chem_param,  only : ntracet
    use ParTools  ,  only : ntracetloc

    ! testing ...
    use dims       , only : istart

    implicit none

    !__IO___________________________________________________________________

    !__LOCAL_VARIABLES______________________________________________________

    real,dimension(:,:,:,:),pointer   :: rm, rxm, rym, rzm
    integer                           :: region
    integer                           :: n

    !__START_SUBROUTINE______________________________________________________


    do region=1,nregions

       rm => mass_dat(region)%rm_t
       rxm => mass_dat(region)%rxm_t
       rym => mass_dat(region)%rym_t
       rzm => mass_dat(region)%rzm_t

       do n=1,ntracetloc

          rm(:,:,:,n) = 0.0

          if ( adv_scheme == 'slope' ) then
             rxm(:,:,:,n) = 0.0
             rym(:,:,:,n) = 0.0
             rzm(:,:,:,n) = 0.0
          end if
       end do

       ! adhoc for testing changes in adj model ...
       if (istart==88) then
         print *, 'WARNING - init adjoint rm to 1.0 for region ', region
         rm = 1.0
       end if

       !?       allocate(chem_dat(region)%rm_k(im(region),jm(region), &
       !?            lm(region),ntracet+1:ntracet+ntrace_chem))

       nullify(rm)
       nullify(rxm)
       nullify(rym)
       nullify(rzm)

    end do


  end subroutine adj_trace1

  !=====================================================================================================
  !=====================================================================================================

  subroutine adj_trace_after_read

    implicit none

    !__IO___________________________________________________________________

    !__LOCAL_VARIABLES______________________________________________________

    !__START_SUBROUTINE______________________________________________________


  end subroutine adj_trace_after_read

  !=====================================================================================================
  !=====================================================================================================

  subroutine adj_source1( region, tr, status )

    use GO               , only : TDate
    use Emission,          only : Emission_Adj_Apply
#ifndef without_dry_deposition
    use dry_deposition,    only : Dry_Deposition_Apply
#endif


    !__IO___________________________________________________________________

    integer, intent(in)        :: region
    type(TDate), intent(in)    :: tr(2)
    integer, intent(out)       :: status

    !__CONST________________________________________________________________

    character(len=*), parameter   ::  rname = mname//'/adj_source1'

    !__LOCAL_VARIABLES______________________________________________________

    !__START_SUBROUTINE______________________________________________________

#ifndef without_dry_deposition
    call Dry_Deposition_Apply( region, tr, status )
    IF_NOTOK_RETURN(status=1)
#endif

    call Emission_Adj_Apply( region, tr, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine adj_source1

  !=====================================================================================================
  !=====================================================================================================

  subroutine adj_source2(region)

    implicit none

    !__IO___________________________________________________________________

    integer, intent(in)        :: region

    !__LOCAL_VARIABLES______________________________________________________

    !__START_SUBROUTINE______________________________________________________


  end subroutine adj_source2

  !=====================================================================================================
  !=====================================================================================================

  subroutine adj_trace_end(msg, file_name)

!    use adj_emission,    only : adj_emission_done


!#endif
!#ifndef without_dry_deposition
!    use dry_deposition,  only : Dry_Deposition_Done
!#endif
!
!
!    implicit none

    !__IO___________________________________________________________________

    character(len=*),intent(in) :: file_name
    character(len=*),intent(in) :: msg

    !__LOCAL_VARIABLES______________________________________________________

    integer          :: status

    !__CONST________________________________________________________________

    character(len=*), parameter   ::  rname = mname//'/adj_trace_end'

    !__START_SUBROUTINE______________________________________________________

!    call adj_emission_done
!#ifndef without_chemistry
!    call adj_chemistry_done
!#endif
!#ifndef without_dry_deposition
!    call  Dry_Deposition_Done( status )
!    IF_NOTOK_RETURN(status=1)
!#endif


  end subroutine adj_trace_end

  !=====================================================================================================
  !=====================================================================================================

end module adj_sources_sinks
