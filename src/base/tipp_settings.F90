!----------------------------------------------------------
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!----------------------------------------------------------

module TIPP_Settings

  use GO, only : gol, goPr, goErr
  use os_specs, only : MAX_FILENAME_LEN

  implicit none


  ! --- in/out -----------------------------

  public


  ! --- const ------------------------------

  character(len=*), parameter  ::  mname = 'TIPP_Settings'


  ! --- var -------------------------------

  character(len=MAX_FILENAME_LEN)      ::  tipp_ecmwf_dir
  character(len=MAX_FILENAME_LEN)      ::  em_gfed_dir
  character(len=MAX_FILENAME_LEN)      ::  em_gfed_8day_dir
  logical                 ::  em_gfed_hprof
  character(len=MAX_FILENAME_LEN)      ::  em_giss_dir
  character(len=MAX_FILENAME_LEN)      ::  em_geia_dir
  character(len=MAX_FILENAME_LEN)      ::  em_kaplan_dir
  character(len=MAX_FILENAME_LEN)      ::  em_edgar_32_dir, em_edgar_32ft_dir, em_edgar_40_dir
  character(len=MAX_FILENAME_LEN)      ::  em_matthews_dir
  character(len=MAX_FILENAME_LEN)      ::  em_tmm_dir


contains


  ! ====================================================================


  subroutine TIPP_Settings_Init( rcf, status )

    use GO, only : TRcFile, ReadRc

    ! --- in/out -------------------------------------

    type(TRcFile), intent(in)       ::  rcf
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  :: rname = mname//'/TIPP_Settings_Init'

    ! --- local --------------------------------------

    ! --- begin --------------------------------------

    ! TIPP settings:
    call ReadRc( rcf, 'tipp.ecmwf.dir', tipp_ecmwf_dir, status )
    IF_NOTOK_RETURN(status=1)

    ! GFED settings:
    call ReadRc( rcf, 'em.gfed.dir', em_gfed_dir, status )
    IF_NOTOK_RETURN(status=1)
    call ReadRc( rcf, 'em.gfed_8day.dir', em_gfed_8day_dir, status )
    IF_NOTOK_RETURN(status=1)
    call ReadRc( rcf, 'em.gfed.hprof', em_gfed_hprof, status )
    IF_NOTOK_RETURN(status=1)

    ! GISS settings:
    call ReadRc( rcf, 'em.giss.dir', em_giss_dir, status )
    IF_NOTOK_RETURN(status=1)

    ! GEIA settings:
    call ReadRc( rcf, 'em.geia.dir', em_geia_dir, status )
    IF_NOTOK_RETURN(status=1)

    ! Kaplan model settings:
    call ReadRc( rcf, 'em.kaplan.dir', em_kaplan_dir, status )
    IF_NOTOK_RETURN(status=1)

    ! EDGAR settings:
    call ReadRc( rcf, 'em.edgar_32.dir', em_edgar_32_dir, status )
    IF_NOTOK_RETURN(status=1)
    call ReadRc( rcf, 'em.edgar_32ft.dir', em_edgar_32ft_dir, status )
    IF_NOTOK_RETURN(status=1)
    call ReadRc( rcf, 'em.edgar_40.dir', em_edgar_40_dir, status )
    IF_NOTOK_RETURN(status=1)

    ! Matthews rice seasonality
    call ReadRc( rcf, 'em.matthews.dir', em_matthews_dir, status )
    IF_NOTOK_RETURN(status=1)

    ! TM meteo archive:
    call ReadRc( rcf, 'em.tmm.dir', em_tmm_dir, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine TIPP_Settings_Init


  ! ***


  subroutine TIPP_Settings_Done( status )

    use GO, only : TRcFile

    ! --- in/out -------------------------------------

    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  :: rname = mname//'/TIPP_Settings_Done'

    ! --- local --------------------------------------

    ! --- begin --------------------------------------

    ! ok
    status = 0

  end subroutine TIPP_Settings_Done


end module TIPP_Settings
