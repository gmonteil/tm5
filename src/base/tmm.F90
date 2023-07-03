!###############################################################################
!
!ProTeX: 1.14-AJS
!
!BOI
!
! !TITLE:        TMM - TM Meteo
! !AUTHORS:      Arjo Segers
! !AFFILIATION:  KNMI
! !DATE:         21/04/2004
!
! !INTRODUCTION: Introduction
!
!   Module to access TM meteo.
!
!   The main structure provides access to a list of opened meteo files.
!   If a new meteo field is required, the subroutines search in the
!   list wether the field is available.
!   If not, a new file is opened and added to the list.
!   Optionally, a shell script is invoked to search for a file
!   and to store it locally if necessary.
!
!
! !INTRODUCTION: Usage
!
! \bv
!
!   ! --- modules -----------------------------------------
!
!   use GO, only : TDate, NewDate
!
!   use grid, only : TllGridInfo, Init, Done
!   use grid, only : TLevelInfo, Init, Done
!
!   use TMM, only : TTmMeteo
!   use TMM, only : Init, Done
!   use TMM, only : ReadField, ReadUVSP
!
!   ! --- local -----------------------------------------
!
!   type(TllGridInfo)           ::  lli
!   type(TLevelInfo)            ::  levi_ec, levi
!
!   type(TTmMeteo)              ::  tmmd
!
!   type(TDate)                 ::  tday, t1, t2
!
!   real                        ::  psurf(120,90)
!   real                        ::  temper(120,90,25)
!   real                        ::  pu(0:120,90,25)
!
!   ! --- begin -----------------------------------------
!
!   ! define horizontal grid
!   call Init( lli, -178.5, 3.0, 120, -89.0, 2.0, 90, status )
!
!   ! define vertical hybride levels
!   call Init( levi_ec, 'ec60', status )            ! ecmwf levels
!   call Init( levi, levi_ec, (/60,..,0/), status ) ! tm half level selection
!
!   ! setup TM meteo access:
!   call Init( tmmd, 'tm5.rc', status )
!
!   ! define time range for field:
!   tday = NewDate( year=2003, month=01, day=01 )
!   t1   = NewDate( year=2002, month=12, day=31, hour=21 )
!   t2   = NewDate( year=2003, month=01, day=01, hour=03 )
!
!   ! Read meteo arrays; specify grid, parameter, time, etc.
!   !
!   ! Type of grid is defined by nuv key:
!   !  'n' = normal grid (cell centers)
!   !  'u' = u-grid (east/west boundaries)
!   !  'v' = v-grid (north/south boundaries)
!   !
!   ! Requested grid (lli/levi) might be different from the grid in the file.
!   ! Horizontal:
!   !   o if file contains same resolutions as defined by lli,
!   !     a part of the data in the file is selected;
!   !   o if file contains higher resolution fields,
!   !     the ll array is filled with combined values from the file;
!   !     depending on the parameter, the result is summed/avaraged/etc.
!   ! Vertical:
!   !   o if a file contains a supperset of the levels in levi,
!   !     some levels are combined;
!   !     depending on the parameter, the result is summed/avaraged/etc;
!   !   o the file might contain fields with reversed level order.
!   !
!   ! 3D fields require surface pressure for the level definition.
!   ! It should be valid for [t1,t2] !
!   ! Best is to use the spm from ReadUVSP.
!
!   call ReadUVSP ( tmmd, 'tmpp:od-fc-ml60-glb3x2', tday, t1, t2, lli, levi, sp1, spm, sp2, pu, pv, status )
!
!   call ReadField( tmmd, 'tmpp:od-fc-ml60-glb3x2', 'T' , tday, t1, t2, lli, 'n', levi, spm, temper, status )
!   call ReadField( tmmd, 'tmpp:od-fc-ml60-glb3x2', 'pu', tday, t1, t2, lli, 'u', levi, spm, pu    , status )
!
!   !
!   ! TMPP surface fields can be called using the 1x1 as well as the 3x2 key.
!   !
!   call ReadField( tmmd, 'tmpp:od-fc-sfc-glb1x1' , 'ci', tday, t1, t2, lli, 'n', ci, status )
!   call ReadField( tmmd, 'tmpp:od-fc-ml60-glb3x2', 'ci', tday, t1, t2, lli, 'n', ci, status )
!   !
!   ! Similar for spm surface pressures:
!   !
!   call ReadField( tmmd, 'tmpp:od-fc-ml1-glb3x2' , 'spm', tday, t1, t1, lli, 'n', spm, status )
!   call ReadField( tmmd, 'tmpp:od-fc-ml60-glb3x2', 'spm', tday, t1, t1, lli, 'n', spm, status )
!
!
!   ! *** output
!
!   call WriteField( tmmd, 'tmpp:od-fc-ml60-glb3x2', 'T', 'K', t1, t2, &
!                          lli, 'n', levi, spm, temper, status )
!
!   ! *** finish
!
!   call Done( tmmd, status )
!   call Done( levi, status    )
!   call Done( levi_ec, status )
!   call Done( lli )
!
!   ! -- end
!
! \ev
!
!
! !INTRODUCTION: Rcfile
!
! \bv
!
!   !
!   ! Meteo files are linked to or unpacked in a buffer directory.
!   !  o Set the clean flag (T|F) such that files that have not been accessed
!   !    for a long time are removed if a maximum buffer usage is exceeded.
!   !  o specify a maximum size in Mb
!   !
!   tmm.dir        : ${RUNDIR}/tmm-buf
!   tmm.dir.clean  : T
!   tmm.dir.size   : 500
!
!   !
!   ! TMM requires keys on how to form meteo for a certain region.
!   ! A key should be defined for each region, names are in 'dims_grid.F90'
!   ! For example:
!   !
!   !   tmpp:od-fc-ml60-glb3x2
!   !     Read global 3x2, 60 level files produced by TMPP.
!   !     Optionally, the meteo is combined over levels or grid cells.
!   !     The files are expected to be present in the buffer directory
!   !     specified below after 'tmm.buf.dir' .
!   !     To have the appropriate files installed at the begin of a run,
!   !     use the 'tmm.setup.*' stuff below.
!   !
!   !   tmppS:od-fc-ml60-glb3x2
!   !     Idem, but also calls a script to search for an appropriate file
!   !     from within the fortran program.
!   !     The system call to this script turned out to be rather slow.
!   !     This source type should be avoided therefore, but might be
!   !     very usefull in case of limitted disk space.
!   !
!   !   prism:
!   !     Receive meteo from the prism coupler.
!   !
!   tmm.sourcekey.glb6x4  : tmpp:od-fc-ml60-glb3x2
!   tmm.sourcekey.eur3x2  : tmpp:od-fc-ml60-glb3x2
!   tmm.sourcekey.eur1x1  : tmpp:od-fc-tropo25-eur1x1
!
!   !
!   ! Meteo files could be setup before the actual program is started.
!   ! Fill the following settings:
!   !  o set the apply flag apply this feature or not (T|F)
!   !  o specify a list of meteo files to be installed (spm,uvsp, etc)
!   !  o specify a list of meteo sources (od-fc-ml60-glb3x2 etc)
!   !  o specify wether message are printend or not (T|F)
!   !
!   tmm.setup.apply    :  T
!   tmm.setup.files    :  spm uvsp t q cld sub surf
!   tmm.setup.sources  :  od-fc-ml60-glb3x2
!   tmm.setup.verbose  :  T
!
!   !
!   ! Archive(s) to be searched for monthly tar files.
!   ! If more than one is specified (space seperated list),
!   ! multiple directories are examined.
!   !  o disk archives
!   !  o tape archives ecfs/mos ('massive-storage-system')
!   !
!   tmm.search.disk   : ${DATADIR}/meteo
!   tmm.search.mss    : /nlh/TM/meteo
!
! \ev
!
!
! !INTRODUCTION: Source and scripts
!
! \bv
!
!   tmm.f90                : Main routines and collecting data structure.
!                            Provides access to a list of open meteo files.
!
!     tmm_mf.f90           : Search, open, close a meteo file.
!                            Calls shell script 'bin/tmm_getmeteo'.
!                            Calls specific routines for hdf/etc files.
!
!       tmm_mf_hdf.f90     : Read fields from hdf files.
!
!       tmm_param.f90      : Parameter specific stuff:
!                            what to do with temperature fields,
!                            what to do with mass fluxes etc.
!
! \ev
!
!EOI
!

!
! The spm/spmid problem ...
!
!   Glossary:
!     spm      : surface pressures for 00, 06, ...;
!                read from 'spm_' files, computed from ecmwf lnsp
!     spmid    : average of surface pressures at begin/end interval
!
!   To have the same algorithm as in TMPP,
!   the following procedure should be used:
!
!    1) use ReadUVSP to get sp1, spm, sp2
!       spm is now in fact a spmid field (average of sp1 and sp2);
!       in future it should be the real spm:
!
!         call ReadUVSP( tmm, archivekey, tday, t1, t2, &
!                             lli, levi, sp1, spm, sp2, pu, pv, status )
!
!    2) use this spm in calls to ReadField:
!
!         call ReadField( tmm, archivekey, paramkey, tday, t1, t2, &
!                              lli, nuv, levi, spm, ll, status )
!
!   Examples of current implementation:
!
!     1) Temperature 6x4x25 from 3x2x60
!          read spm 3x2 from file
!          horizontal mass average to 6x4x60 using spm 3x2 from
!          vertical combination to 6x4x25 using provided sp 3x2
!
!     2) Temperature 6x4x25 from N80x60
!          read spm N80
!          horizontal mass average to 6x4x60 using spm N80 and provided sp 3x2
!            --> should be spm 3x2 from spm N80
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tmm.inc"
!
!###############################################################################

module TMM

  use GO      , only : gol, goPr, goErr, goBug
  use GO      , only : TDate
  use Grid    , only : TllGridInfo, TggGridInfo, TshGridInfo, TLevelInfo
  use tmm_mf  , only : TMeteoFile
  use tmm_info, only : TMeteoInfo, SetHistory, AddHistory
  use os_specs, only : MAX_FILENAME_LEN, MAX_RCKEY_LEN

  implicit none

  ! --- in/out -------------------------------------

  private

  public  ::  TTmMeteo
  public  ::  Init, Done

  public  ::  ReadField
  public  ::  Read_SP, Read_MFUV, Read_MFW, Read_TQ
  public  ::  Read_Convec, Read_CloudCovers
  public  ::  Read_SR_OLS
!  public  ::  ReadEqvLat

  public  ::  WriteField

  public  ::  TMeteoInfo, SetHistory, AddHistory

  ! --- const ---------------------------------------

  character(len=*), parameter  ::  mname = 'tmm'

  ! maximum number of opened meteo files
  integer, parameter           ::  maxmf = 200


  ! --- types ---------------------------------------

  type TTmMeteo
    ! rc file with paths etc
    character(len=MAX_FILENAME_LEN)            ::  rcfilename
    ! paths
    character(len=MAX_FILENAME_LEN)            ::  input_dir, output_dir, input_1x1_dir
    ! list of meteo files
    integer                       ::  nmf
    type(TMeteoFile)              ::  mf(maxmf)
    !
    ! buffers for latest read raw field
    !
    logical                    ::  buf_filled
    !
    character(len=10)          ::  buf_archivekey
    character(len=10)          ::  buf_paramkey
    type(TDate)                ::  buf_tday, buf_t1, buf_t2
    character(len=1)           ::  buf_nuv
    character(len=1)           ::  buf_nw
    type(TMeteoInfo)           ::  buf_tmi
    !
    character(len=2)           ::  buf_gridtype
    !
    type(TllGridInfo)          ::  buf_lli
    real, pointer              ::  buf_ll(:,:,:)
    real, pointer              ::  buf_sp_ll(:,:)
    !
    type(TggGridInfo)          ::  buf_ggi
    real, pointer              ::  buf_gg(:,:)
    real, pointer              ::  buf_sp_gg(:)
    !
    type(TshGridInfo)          ::  buf_shi
    complex, pointer           ::  buf_sh(:,:)
    complex, pointer           ::  buf_lnsp_sh(:)
    !
    type(TLevelInfo)           ::  buf_levi
    !
    ! storage for horizontal interpol
    real, pointer              ::  llX(:,:,:)
  end type TTmMeteo


  ! --- interfaces -------------------------------------

  interface Init
    module procedure tmm_Init
  end interface

  interface Done
    module procedure tmm_Done
  end interface

  interface SelectMF
    module procedure tmm_SelectMF
  end interface

  interface ReadField
    module procedure tmm_ReadField_2d
    module procedure tmm_ReadField_3d
  end interface

!  interface ReadUVSP
!    module procedure tmm_ReadUVSP
!  end interface

  interface Read_SP
    module procedure tmm_Read_SP
  end interface

  interface Read_MFUV
    module procedure tmm_Read_MFUV
  end interface

  interface Read_MFW
    module procedure tmm_Read_MFW
  end interface

  interface Read_TQ
    module procedure tmm_Read_TQ
  end interface

  interface Read_Convec
    module procedure tmm_Read_Convec
  end interface

  interface Read_CloudCovers
    module procedure tmm_Read_CloudCovers
  end interface

  interface Read_SR_OLS
    module procedure tmm_Read_SR_OLS
  end interface

!  interface ReadEqvLat
!    module procedure tmm_ReadEqvLat
!  end interface

  interface WriteField
    module procedure tmm_WriteField_2d
    module procedure tmm_WriteField_3d
  end interface


  ! --- var --------------------------------------

  ! timer id's:
  integer           ::  itim_fillbuffer
  integer           ::  itim_readfield_2d
  integer           ::  itim_readfield_3d
  integer           ::  itim_transform_2d
  integer           ::  itim_transform_3dh
  integer           ::  itim_transform_3dv


contains


  ! ===================================================================


  subroutine tmm_Init( tmm, rcF, status )

    use PArray      , only : pa_Init
    use GO          , only : goErr, gol, goPr
    use GO          , only : NewDate
    use GO          , only : TrcFile, ReadRc
    use GO          , only : GO_Timer_Def
#ifdef with_tmm_tm5
    use tmm_mf_tm5_nc , only : TMM_MF_TM5_NC_Init
#endif
#ifdef with_prism
    use tmm_mf_prism, only : mfPrism_Init
#endif

    ! --- in/out ----------------------------------

    type(TTmMeteo), intent(out)         ::  tmm
    type(TrcFile), intent(in)           ::  rcF
    integer, intent(out)                ::  status

    ! --- const ------------------------------------

    character(len=*), parameter ::  rname = mname//'/tmm_Init'

    ! --- local -----------------------------------

    ! --- begin -----------------------------------

    ! store rc file name
    tmm%rcfilename = trim(rcf%fname)

    ! read paths
    call ReadRc( rcF, 'tmm.dir', tmm%input_dir, status )
    IF_NOTOK_RETURN(status=1)
    call ReadRc( rcF, 'my.meteo.source.dir', tmm%input_1x1_dir, status)
    IF_NOTOK_RETURN(status=1)
    ! read output path: set to a dummy value if key not found in rcfiles,
    ! since probably no meteo output requested anyway in this case ...
    call ReadRc( rcF, 'tmm.output.dir', tmm%output_dir, status, &
                         default='/no/tmm/output/dir/specified/' )

    ! no files open yet
    tmm%nmf = 0

    ! buffer empty
    tmm%buf_filled     = .false.
    tmm%buf_archivekey = 'none'
    tmm%buf_paramkey   = 'none'
    tmm%buf_tday       = NewDate(time5=(/0001,01,01,00,00/))
    tmm%buf_t1         = NewDate(time5=(/0001,01,01,00,00/))
    tmm%buf_t2         = NewDate(time5=(/9999,12,31,00,00/))
    tmm%buf_nuv        = 'n'
    tmm%buf_nw         = 'n'
    tmm%buf_gridtype   = 'xx'
    call pa_Init( tmm%buf_ll )
    call pa_Init( tmm%buf_sp_ll )
    call pa_Init( tmm%buf_gg )
    call pa_Init( tmm%buf_sp_gg )
    call pa_Init( tmm%buf_sh )
    call pa_Init( tmm%buf_lnsp_sh )

    ! init temp grid on destination grid but raw levels
    call pa_Init( tmm%llX )

#ifdef with_prism
    ! init prism stuff:
    call mfPrism_Init( status )
    IF_NOTOK_RETURN(status=1)
#endif

#ifdef with_tmm_tm5
    ! init input of TM5/NetCDF files:
    call TMM_MF_TM5_NC_Init( rcf, status )
    IF_NOTOK_RETURN(status=1)
#endif

    ! init timers:
    call GO_Timer_Def( itim_fillbuffer   , 'tmm fill buffer', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_readfield_2d , 'tmm readfield 2D', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_readfield_3d , 'tmm readfield 3D', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_transform_2d , 'tmm transform 2D', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_transform_3dh, 'tmm transform 3D hor' , status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_transform_3dv, 'tmm transform 3D vert', status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine tmm_Init


  ! ***


  subroutine tmm_Done( tmm, status )

    use PArray      , only : pa_Done
    use GO          , only : goErr, gol, goPr
    use tmm_mf      , only : Opened, Done
#ifdef with_tmm_tm5
    use tmm_mf_tm5_nc , only : TMM_MF_TM5_NC_Done
#endif
#ifdef with_prism
    use tmm_mf_prism, only : mfPrism_Done
#endif

    ! --- in/out ----------------------------------

    type(TTmMeteo), intent(inout)          ::  tmm
    integer, intent(out)                   ::  status

    ! --- const ------------------------------------

    character(len=*), parameter ::  rname = mname//'/tmm_Done'

    ! --- local ----------------------------------------

    integer          ::  imf

    ! --- begin -------------------------------------------

#ifdef with_tmm_tm5
    ! init input of TM5/NetCDF files:
    call TMM_MF_TM5_NC_Done( status )
    IF_NOTOK_RETURN(status=1)
#endif

#ifdef with_prism
    ! done with prism stuff:
    call mfPrism_Done( status )
    IF_NOTOK_RETURN(status=1)
#endif

    ! loop over all available meteo files
    do imf = 1, tmm%nmf

      ! in use ?
      if ( .not. Opened( tmm%mf(imf) ) ) cycle

      ! close ...
      call Done( tmm%mf(imf), status )
      IF_NOTOK_RETURN(status=1)

    end do

    ! clear buffer
    call ClearBuffer( tmm, status )
    IF_NOTOK_RETURN(status=1)
    call pa_Done( tmm%buf_ll )
    call pa_Done( tmm%buf_sp_ll )
    call pa_Done( tmm%buf_gg )
    call pa_Done( tmm%buf_sp_gg )
    call pa_Done( tmm%buf_sh )
    call pa_Done( tmm%buf_lnsp_sh )

    ! clear temp grid
    call pa_Done( tmm%llX )

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! ok
    status = 0

  end subroutine tmm_Done


  ! ***


  ! Select meteo file for this param and time,
  ! search for a new file if necessary

  subroutine tmm_SelectMF( tmm, io, archivekey, paramkey, tday, t1, t2, imf, status )

    use GO    , only : goErr, gol, goPr
    use GO    , only : TDate, wrtgol, Pretty
    use tmm_mf, only : Init, Done, Opened, CheckTime, CheckParam
    use tmm_mf, only : SetupInput, SetupOutput
    use dims  , only : okdebug_tmm

    ! --- in/out ----------------------------------------

    type(TTmMeteo), intent(inout)           ::  tmm
    character(len=1), intent(in)            ::  io
    character(len=*), intent(in)            ::  archivekey
    character(len=*), intent(in)            ::  paramkey
    type(TDate), intent(in)                 ::  tday, t1, t2
    integer, intent(out)                    ::  imf
    integer, intent(out)                    ::  status

    ! --- const ------------------------------------

    character(len=*), parameter ::  rname = mname//'/tmm_SelectMF'

    ! --- local ----------------------------------------

    integer                  ::  i
    logical                  :: file_exists

    ! --- begin ----------------------------------------

    if (okdebug_tmm) then
        write(gol,'(a, ": paramkey = ", a)') rname, trim(paramkey); call gobug
        write(gol,'(a, ": archivekey = ", a)') rname, trim(archivekey); call gobug
        write(gol,'(a, ": tday = ", a)') rname, trim(Pretty(tday)); call gobug
        write(gol,'(a, ": t1 - t2 = ", a, " - ", a)') rname, trim(Pretty(t1)), trim(Pretty(t2)); call gobug
        !call wrtgol( ' DDD selectmf tday       : ', tday ); call gobug
        !call wrtgol( ' DDD selectmf t1 - t2    : ', t1, ' - ', t2 ); call gobug
    end if

    ! not found yet ...
    imf = -1

    ! loop over currently used meteo files
    do i = 1, tmm%nmf
      ! in use ?
      if ( .not. Opened( tmm%mf(i) ) ) cycle
      ! test on requested grid and param:
      call CheckParam( tmm%mf(i), io, archivekey, paramkey, status )
      !write (gol,'("DDD selectmf chk param  : ", i2, "  ", a, "  ", i2)') i, trim(tmm%mf(i)%filename), status; call gobug
      if ( status == 0 ) then
        ! param included: leave
        imf = i
        exit
      else if ( status < 0 ) then
        ! param not included, try next
        cycle
      else
        ! error ...
        TRACEBACK; status=1; return
      end if
    end do

    ! time ok ? close if data in file is too old:
    if ( imf > 0 ) then
      ! time ok ?
      call CheckTime( tmm%mf(imf), t1, t2, status )
      !write (gol,*) 'DDD selectmf chk time   : ', status; call gobug
      if ( status == 0 ) then
        ! mf includes [t1,t2]; return with this imf
        status = 0; return
      else if ( status < 0  ) then
        ! mf does not include [t1,t2]; close file
        !write (gol,*) 'DDD selectmf close      : ', imf, trim(tmm%mf(imf)%filename); call gobug
        call Done( tmm%mf(imf), status )
        IF_NOTOK_RETURN(status=1)
        ! not in use anymore ...
        imf = -1
      else
        TRACEBACK; status=1; return
      end if
    end if
    !write (gol,*) 'DDD selectmf imf        : ', imf; call gobug

    ! open new meteo file ?
    if ( imf < 0 ) then
      ! search first available empty mf
      do i = 1, tmm%nmf
        if ( .not. Opened(tmm%mf(i)) ) then
          imf = i
          exit
        end if
      end do
      ! next number ?
      if ( imf < 0 ) then
        tmm%nmf = tmm%nmf + 1
        if ( tmm%nmf > maxmf ) then
          write (gol,'("Tried to init meteo file beyond maximum number: ",i6)') maxmf; call goErr
          write (gol,'("Initialized files:")'); call goErr
          do i = 1, maxmf
            write (gol,'("  ",a)') trim(tmm%mf(i)%filename); call goErr
          end do
          write (gol,'("Please increase parameter `maxmf` in ",a)') mname; call goErr
          TRACEBACK; status=1; return
        end if
        imf = tmm%nmf
      end if
      ! start new mf ...
      call Init( tmm%mf(imf), io, status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! input or output ?
    select case ( io )

      case ( 'i' )  ! input

        ! open file, store description, etc
        call SetupInput( tmm%mf(imf), archivekey, paramkey, tday, t1, t2, tmm%rcfilename, tmm%input_dir, status )
        IF_NOTOK_RETURN(status=1)
        ! If the file does not exist, go back to the global 1x1 archive
        inquire( file=trim(tmm%mf(imf)%filename), exist=file_exists )
        if (.not. file_exists) then
            write(gol,'(a,": ",a," does not exist")') rname, trim(tmm%mf(imf)%filename)
            call goPr
            call SetupInput( tmm%mf(imf), archivekey, paramkey, tday, t1, t2, tmm%rcfilename, tmm%input_1x1_dir, status )
            IF_NOTOK_RETURN(status=1)
            write(gol,'(a,": ",a," will be read instead")') rname, trim(tmm%mf(imf)%filename)
            call goPr
        end if

      case ( 'o' )  ! output

        ! open file, store description, etc
        call SetupOutput( tmm%mf(imf), archivekey, paramkey, tday, t1, t2, &
                                       tmm%rcfilename, tmm%output_dir, status )
        IF_NOTOK_RETURN(status=1)

      case default
        write (gol,'("unsupported io `",a,"`")') io; call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine tmm_SelectMF


  ! ==================================================================
  ! ===
  ! ===  buffer
  ! ===
  ! ==================================================================


  subroutine FillBuffer( tmm, archivekey, paramkey, unit, tday, t1, t2, nuv, nw, status )

    use GO    , only : goErr, gol, goPr
    use GO    , only : TDate, operator(==), Pretty
    use GO    , only : GO_Timer_Start, GO_Timer_End
    use tmm_mf, only : ReadRecord
    use dims  , only : okdebug_tmm

    ! --- in/out --------------------------------

    type(TTmMeteo), intent(inout)           ::  tmm
    character(len=*), intent(in)            ::  archivekey
    character(len=*), intent(in)            ::  paramkey
    character(len=*), intent(in)            ::  unit
    type(TDate), intent(in)                 ::  tday, t1, t2
    character(len=1), intent(in)            ::  nuv
    character(len=1), intent(in)            ::  nw
    integer, intent(out)                    ::  status

    ! --- const ------------------------------------

    character(len=*), parameter ::  rname = mname//'/FillBuffer'

    ! --- local -------------------------------

    integer          ::  imf

    ! --- begin -------------------------------

    ! start timing:
    call GO_Timer_Start( itim_fillbuffer, status )
    IF_NOTOK_RETURN(status=1)

    if ( (archivekey == tmm%buf_archivekey) .and. (paramkey == tmm%buf_paramkey) .and. &
         (nuv == tmm%buf_nuv) .and. (nw == tmm%buf_nw) .and. &
         (tday == tmm%buf_tday) .and. (t1 == tmm%buf_t1) .and. (t2 == tmm%buf_t2) ) then
      ! requested field already in buffer ...
      ! ok
      status = 0; return
    end if

    ! select (eventually retrieve first) the meteo file with this param:
    call SelectMF( tmm, 'i', archivekey, paramkey, tday, t1, t2, imf, status )
    IF_NOTOK_RETURN(status=1)

    ! clear buffer
    call ClearBuffer( tmm, status )
    IF_NOTOK_RETURN(status=1)

    ! fill keys:
    tmm%buf_archivekey = archivekey
    tmm%buf_paramkey   = paramkey
    tmm%buf_t1         = t1
    tmm%buf_t2         = t2
    tmm%buf_tday       = tday
    tmm%buf_nuv        = nuv
    tmm%buf_nw         = nw

    if (okdebug_tmm) then
        write(gol,'(a, ": ", a, " between ", a, " and ", a, " to be read from ", a)') &
            rname, trim(paramkey), trim(Pretty(t1)), trim(Pretty(t2)), trim(tmm%mf(imf)%filename)
        call goPr
    end if
    ! read field:
    call ReadRecord( tmm%mf(imf), paramkey, unit, tday, t1, t2, nuv, nw, &
                     tmm%buf_gridtype, tmm%buf_levi, &
                     tmm%buf_lli, tmm%buf_ll, tmm%buf_sp_ll, &
                     tmm%buf_ggi, tmm%buf_gg, tmm%buf_sp_gg, &
                     tmm%buf_shi, tmm%buf_sh, tmm%buf_lnsp_sh, &
                     tmm%buf_tmi, status )
    IF_NOTOK_RETURN(status=1)

    ! some data present ..
    tmm%buf_filled   = .true.

    ! end timing:
    call GO_Timer_End( itim_fillbuffer, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine FillBuffer


  ! ***


  subroutine ClearBuffer( tmm, status )

    use GO      , only : goErr, gol, goPr
    use PArray  , only : pa_Done
    use GO      , only : goErr
    use Grid    , only : Done
    use tmm_info, only : Done

    ! --- in/out --------------------------------

    type(TTmMeteo), intent(inout)           ::  tmm
    integer, intent(out)                    ::  status

    ! --- const ------------------------------------

    character(len=*), parameter ::  rname = mname//'/ClearBuffer'

    ! --- begin -------------------------------

    if ( tmm%buf_filled ) then

      call Done( tmm%buf_tmi, status )
      IF_NOTOK_RETURN(status=1)

      if ( tmm%buf_gridtype == 'll' ) then
        call Done( tmm%buf_lli, status )
        IF_NOTOK_RETURN(status=1)
        call pa_Done( tmm%buf_ll    )
        call pa_Done( tmm%buf_sp_ll )
      end if

      if ( tmm%buf_gridtype == 'gg' ) then
        call Done( tmm%buf_ggi, status )
        IF_NOTOK_RETURN(status=1)
        call pa_Done( tmm%buf_gg    )
        call pa_Done( tmm%buf_sp_gg )
      end if

      if ( tmm%buf_gridtype == 'sh' ) then
        call Done( tmm%buf_shi )
        call pa_Done( tmm%buf_sh    )
        call pa_Done( tmm%buf_lnsp_sh )
      end if

      call Done( tmm%buf_levi, status )
      IF_NOTOK_RETURN(status=1)

    end if

    ! ok
    status = 0

  end subroutine ClearBuffer


  ! ***


  subroutine CheckBuffer( tmm, status, gridtype, np, shT, nlev )

    use GO, only : goErr, gol, goPr

    ! --- in/out --------------------------------

    type(TTmMeteo), intent(inout)           ::  tmm
    integer, intent(out)                    ::  status
    character(len=*), intent(in), optional  ::  gridtype
    integer, intent(in), optional           ::  np
    integer, intent(in), optional           ::  shT
    integer, intent(in), optional           ::  nlev

    ! --- const ------------------------------------

    character(len=*), parameter ::  rname = mname//'/CheckBuffer'

    ! --- begin -------------------------------

    if ( .not. tmm%buf_filled ) then
      write (gol,'("buffer not filled ...")'); call goErr
      TRACEBACK; status=1; return
    end if

    if ( present(gridtype) ) then
      if ( tmm%buf_gridtype /= gridtype ) then
        write (gol,'("unexpected gridtype in buffer:")'); call goErr
        write (gol,'("  buffer     : ",a)') tmm%buf_gridtype; call goErr
        write (gol,'("  expected   : ",a)') gridtype; call goErr
        TRACEBACK; status=1; return
      end if
    end if

    if ( present(nlev) ) then
      if ( tmm%buf_levi%nlev /= nlev ) then
        write (gol,'("unexpected number of levels in buffer:")'); call goErr
        write (gol,'("  buffer     : ",i4)') tmm%buf_levi%nlev ; call goErr
        write (gol,'("  expected   : ",i4)') nlev ; call goErr
        TRACEBACK; status=1; return
      end if
    end if

    if ( present(np) ) then
      if ( tmm%buf_ggi%np /= np ) then
        write (gol,'("unexpected grid size in buffer:")'); call goErr
        write (gol,'("  buffer ggi%np : ",i6)') tmm%buf_ggi%np; call goErr
        write (gol,'("  expected      : ",i6)') np; call goErr
        TRACEBACK; status=1; return
      end if
    end if

    if ( present(shT) ) then
      if ( tmm%buf_shi%T /= shT ) then
        write (gol,'("unexpected grid size in buffer:")'); call goErr
        write (gol,'("  buffer shi%shT : ",i6)') tmm%buf_shi%T; call goErr
        write (gol,'("  expected       : ",i6)') shT; call goErr
        TRACEBACK; status=1; return
      end if
    end if

    ! ok
    status = 0

  end subroutine CheckBuffer


  ! ***


  subroutine Transform2D( tmm, lli, nuv, tmi, status )

    use PArray       , only : pa_SetShape
    use GO           , only : goErr, gol, goPr
    use GO           , only : GO_Timer_Start, GO_Timer_End
    use Grid         , only : TllGridInfo, FillGrid
    use Grid         , only : Tgg2llFracs, Init, Done
    use Grid         , only : Interpol
    use Grid         , only : FracSum
    use tmm_param    , only : GetCombineKeys
    use tmm_info     , only : TMeteoInfo, Init, AddHistory

    ! --- in/out --------------------------------

    type(TTmMeteo), intent(inout)           ::  tmm
    type(TllGridInfo), intent(in)           ::  lli
    character(len=1), intent(in)            ::  nuv
    type(TMeteoInfo), intent(out)           ::  tmi
    integer, intent(out)                    ::  status

    ! --- const ------------------------------------

    character(len=*), parameter ::  rname = mname//'/Transform2D'

    ! --- local -------------------------------

    character(len=10)       ::  hcomb, vcomb
    type(Tgg2llFracs)       ::  gg2ll

    ! --- begin ----------------------------------

    ! start timing:
    call GO_Timer_Start( itim_transform_2d, status )
    IF_NOTOK_RETURN(status=1)

    ! copy info:
    tmi = tmm%buf_tmi

    ! set shape of target grid:
    select case ( nuv )
      case ( 'n' )
        call pa_SetShape( tmm%llX, lli%nlon  , lli%nlat  , 1 )
      case ( 'u' )
        call pa_SetShape( tmm%llX, lli%nlon+1, lli%nlat  , 1 )
      case ( 'v' )
        call pa_SetShape( tmm%llX, lli%nlon  , lli%nlat+1, 1 )
      case default
        write (gol,'("unsupported nuv `",a,"`")') nuv
        TRACEBACK; status=1; return
    end select

    ! fill with zero's for safety:
    tmm%llX = 0.0

    ! define how to combine horizontal cells and vertical levels:
    call GetCombineKeys( hcomb, vcomb, tmm%buf_paramkey, status )
    IF_NOTOK_RETURN(status=1)

    ! transform raw field to ll :
    select case ( tmm%buf_gridtype )

      ! ~~~ lat/lon ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'll' )

        ! copy or combine horizontal cells from buffer into lli/llX :
        call FillGrid( lli, nuv, tmm%llX(:,:,1), &
                       tmm%buf_lli, tmm%buf_nuv, tmm%buf_ll(:,:,1), &
                       hcomb, status )
        IF_NOTOK_RETURN(status=1)

      ! ~~~ gg ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'gg' )

        ! determine fraction
        call Init( gg2ll, tmm%buf_ggi, lli, status )
        IF_NOTOK_RETURN(status=1)

        ! deceide how to interpolate given hcomb key:
        select case ( hcomb )

          case ( 'area-aver' )

            ! take fractions of overlapping cells
            call FracSum( gg2ll, tmm%buf_gg(:,1), tmm%llX(:,:,1), status, 'area-aver' )
            IF_NOTOK_RETURN(status=1)

          case default

            write (gol,'("unsupported horizonal combination key for gg :",a)') hcomb; call goErr
            write (gol,'("TIP: hcomb not set for this variable in module tmm_param ?")'); call goErr
            TRACEBACK; status=1; return

        end select

        ! done
        call Done( gg2ll )

      ! ~~~ sh ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'sh' )

        ! interpolate field in shi/sh into lli/ll :
        call Interpol( tmm%buf_shi, tmm%buf_sh(:,1), lli, tmm%llX(:,:,1), status )
        IF_NOTOK_RETURN(status=1)

      case default

        write (gol,'("unsupported input grid type `",a,"`")') tmm%buf_gridtype; call goErr
        TRACEBACK; status=1; return

    end select

    ! fill history:
    call Init( tmi, tmm%buf_tmi, status )
    call AddHistory( tmi, 'oper==hcomb,'//trim(hcomb), status )

    ! end timing:
    call GO_Timer_End( itim_transform_2d, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Transform2D


  ! ***


  subroutine Transform3Dh( tmm, lli, nuv, levi, nw, sp, tmi, status )

    use PArray       , only : pa_SetShape
    use GO           , only : goErr, gol, goPr
    use GO           , only : GO_Timer_Start, GO_Timer_End
    use Grid         , only : TLevelInfo, FillLevels
    use Grid         , only : TllGridInfo, FillGrid, AreaOper
    use Grid         , only : Tgg2llFracs, Init, Done
    use Grid         , only : IntArea
    use Grid         , only : FracSum
    use tmm_param    , only : GetCombineKeys
    use tmm_info     , only : TMeteoInfo, Init, AddHistory

    ! --- in/out --------------------------------

    type(TTmMeteo), intent(inout)           ::  tmm
    type(TllGridInfo), intent(in)           ::  lli
    character(len=1), intent(in)            ::  nuv
    type(TLevelInfo), intent(in)            ::  levi
    character(len=1), intent(in)            ::  nw
    real, intent(out)                       ::  sp(:,:)
    type(TMeteoInfo), intent(out)           ::  tmi
    integer, intent(out)                    ::  status

    ! --- const ------------------------------------

    character(len=*), parameter ::  rname = mname//'/Transform3Dh'

    ! --- local -------------------------------

    integer                 ::  nlon, nlat, nlev
    character(len=10)       ::  hcomb, vcomb
    integer                 ::  l
    type(Tgg2llFracs)       ::  gg2ll

    real, allocatable       ::  dp_ll(:,:)
    real, allocatable       ::  dp_llX(:,:)
    real, allocatable       ::  dp_gg(:)

    ! --- begin ----------------------------------

    ! start timing:
    call GO_Timer_Start( itim_transform_3dh, status )
    IF_NOTOK_RETURN(status=1)

    ! set shape of target grid:
    select case ( nuv )
      case ( 'n' )
        nlon = lli%nlon
        nlat = lli%nlat
      case ( 'u' )
        nlon = lli%nlon+1
        nlat = lli%nlat
      case ( 'v' )
        nlon = lli%nlon
        nlat = lli%nlat+1
      case default
        write (gol,'("unsupported nuv `",a,"`")') nuv; call goErr
        TRACEBACK; status=1; return
    end select
    select case ( nw )
      case ( 'n' )
        nlev = tmm%buf_levi%nlev
      case ( 'w' )
        nlev = tmm%buf_levi%nlev+1
      case default
        write (gol,'("unsupported nw `",a,"`")') nw; call goErr
        TRACEBACK; status=1; return
    end select
    call pa_SetShape( tmm%llX, nlon, nlat, nlev )
    tmm%llX = 0.0

    ! define how to combine horizontal cells and vertical levels:
    call GetCombineKeys( hcomb, vcomb, tmm%buf_paramkey, status )
    IF_NOTOK_RETURN(status=1)

    ! transform raw field to ll :
    select case ( tmm%buf_gridtype )

      ! ~~~ lat/lon ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'll' )

        !! check size of buffer array
        !if ( size(tmm%buf_ll,3) /= nlev ) then
        !  write (gol,'("buffer array does not match with level definition:")')
        !  write (gol,'("  nw             : ",a )') nw
        !  write (gol,'("  buf_levi nlev  : ",i3)') tmm%buf_levi%nlev
        !  write (gol,'("  buf_ll levels  : ",i3)') size(tmm%buf_ll,3)
        !  TRACEBACK; status=1; return
        !end if
        ! convective fluxes have only lowest layers ...
        if ( size(tmm%buf_ll,3) > nlev ) then
          write (gol,'("buffer array has more layers than implied by nw and level definition:")'); call goErr
          write (gol,'("  nw             : ",a )') nw; call goErr
          write (gol,'("  buf_levi nlev  : ",i3)') tmm%buf_levi%nlev; call goErr
          write (gol,'("  buf_ll levels  : ",i3)') size(tmm%buf_ll,3); call goErr
          TRACEBACK; status=1; return
        end if

        ! surface pressure: aver horizontal cells from buffer into lli/llX :
        call FillGrid( lli, nuv, sp, &
                       tmm%buf_lli, tmm%buf_nuv, tmm%buf_sp_ll, &
                       'area-aver', status )
        IF_NOTOK_RETURN(status=1)

        ! deceide how to interpolate given hcomb key;
        ! for most keys, let 'FillGrid' determine what to do ...
        select case ( hcomb )

          case ( 'mass-aver' )

            ! check ...
            if ( nw /= 'n' ) then
              write (gol,'("nw should be `n` for mass aver ...")'); call goErr
              TRACEBACK; status=1; return
            end if

            !
            ! mass weighted fractions:
            !
            !    sum_i  f_i dp_i/g dA_i
            !    ----------------------
            !      sum_i  dp_i/g dA_i
            !

            ! temporary pressure field:
            allocate( dp_ll(tmm%buf_lli%nlon,tmm%buf_lli%nlat) )

            ! loop over layers
            !do l = 1, nlev
            do l = 1, size(tmm%buf_ll,3)
              ! pressure gradients in this layer:
              dp_ll  = tmm%buf_levi%da(l) + tmm%buf_levi%db(l) * tmm%buf_sp_ll
              call AreaOper( tmm%buf_lli, dp_ll, '*', 'm2', status )
              IF_NOTOK_RETURN(status=1)
              ! copy or combine horizontal cells from buffer into lli/llX;
              ! combining cells is weighted by 'mass' dp_ll :
              call FillGrid( lli, nuv, tmm%llX(:,:,l), &
                             tmm%buf_lli, tmm%buf_nuv, tmm%buf_ll(:,:,l), &
                             'weight', status, dp_ll )
              IF_NOTOK_RETURN(status=1)
            end do

            ! clear:
            deallocate( dp_ll )

          case default

            ! loop over layers in ll array:
            !do l = 1, nlev
            do l = 1, size(tmm%buf_ll,3)
              ! copy or combine horizontal cells from buffer into lli/llX;
              ! pass hcomb to FillGrid to define how cells are combined:
              call FillGrid( lli, nuv, tmm%llX(:,:,l), &
                             tmm%buf_lli, tmm%buf_nuv, tmm%buf_ll(:,:,l), &
                             hcomb, status )
              IF_NOTOK_RETURN(status=1)
            end do

        end select

      ! ~~~ gg ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'gg' )

        ! check size of buffer array
        if ( size(tmm%buf_gg,2) /= nlev ) then
          write (gol,'("buffer array does not match with level definition:")'); call goErr
          write (gol,'("  nw             : ",a )') nw; call goErr
          write (gol,'("  buf_levi nlev  : ",i3)') tmm%buf_levi%nlev; call goErr
          write (gol,'("  buf_gg levels  : ",i3)') size(tmm%buf_gg,2); call goErr
          TRACEBACK; status=1; return
        end if

        ! determine fraction
        call Init( gg2ll, tmm%buf_ggi, lli, status )
        IF_NOTOK_RETURN(status=1)

        ! surface pressure: take fractions of overlapping cells
        call FracSum( gg2ll, tmm%buf_sp_gg, sp, status, 'area-aver' )
        IF_NOTOK_RETURN(status=1)

        ! deceide how to interpolate given hcomb key:
        select case ( hcomb )

          case ( 'area-aver' )

            !
            ! area weighted fractions:
            !
            !    sum_i  f_i dA_i
            !    ---------------
            !      sum_i  dA_i
            !

            ! loop over layers
            do l = 1, nlev
              ! take fractions of overlapping cells
              call FracSum( gg2ll, tmm%buf_gg(:,l), tmm%llX(:,:,l), status, 'area-aver' )
              IF_NOTOK_RETURN(status=1)
            end do

          case ( 'mass-aver' )

            ! check ...
            if ( nw /= 'n' ) then
              write (gol,'("nw should be `n` for mass aver ...")'); call goErr
              TRACEBACK; status=1; return
            end if

            !
            ! mass weighted fractions:
            !
            !    sum_i  f_i dp_i/g dA_i
            !    ----------------------
            !      sum_i  dp_i/g dA_i
            !

            ! temporary pressure field:
            allocate( dp_gg(size(tmm%buf_sp_gg)) )
            allocate( dp_llX(size(sp,1),size(sp,2)) )

            ! loop over layers
            do l = 1, nlev
              ! pressure gradients in this layer:
              dp_gg  = tmm%buf_levi%da(l) + tmm%buf_levi%db(l) * tmm%buf_sp_gg
              dp_llX = tmm%buf_levi%da(l) + tmm%buf_levi%db(l) * sp
              ! take fractions of overlapping cells:
              !   sum_i  f_i  dp_i  dA_i
              call FracSum( gg2ll, tmm%buf_gg(:,l)*dp_gg, tmm%llX(:,:,l), status, 'area-sum' )
              IF_NOTOK_RETURN(status=1)
              ! weight by total dp dA
              tmm%llX(:,:,l) = tmm%llX(:,:,l) / dp_llX
              call AreaOper( lli, tmm%llX(:,:,l), '/', 'm2', status )
              IF_NOTOK_RETURN(status=1)
            end do

            ! clear:
            deallocate( dp_gg )
            deallocate( dp_llX )

          case default

            write (gol,'("unsupported horizonal combination key for gg :",a)') hcomb; call goErr
            TRACEBACK; status=1; return

        end select

        ! done
        call Done( gg2ll )

      ! ~~~ sh ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'sh' )

        ! check size of buffer array
        if ( size(tmm%buf_sh,2) /= nlev ) then
          write (gol,'("buffer array does not match with level definition:")'); call goErr
          write (gol,'("  nw             : ",a )') nw; call goErr
          write (gol,'("  buf_levi nlev  : ",i3)') tmm%buf_levi%nlev; call goErr
          write (gol,'("  buf_sh levels  : ",i3)') size(tmm%buf_sh,2); call goErr
          TRACEBACK; status=1; return
        end if

        ! surface pressure: interpolate lnsp to target horizontal resolution:
        call IntArea( 'exp,aver', tmm%buf_shi, tmm%buf_lnsp_sh, lli, sp, status )
        IF_NOTOK_RETURN(status=1)

        ! deceide how to interpolate given hcomb key:
        select case ( hcomb )

          !case ( 'area-aver' )
          !
          !  ! area average over sepectral fields:
          !  call IntArea( 'aver', tmm%buf_shi, tmm%buf_sh, lli, tmm%llX, status )
          !  IF_NOTOK_RETURN(status=1)

          case ( 'mass-aver' )

            ! check ...
            if ( nw /= 'n' ) then
              write (gol,'("nw should be `n` for mass aver ...")'); call goErr
              TRACEBACK; status=1; return
            end if

            ! mass average over sepectral fields:
            call IntArea( '[F*(da+db*exp(H))*cos]/[()*cos]', &
                          tmm%buf_shi, nlev, tmm%buf_sh, &
                          tmm%buf_lnsp_sh, tmm%buf_levi%da, tmm%buf_levi%db, &
                          lli, tmm%llX, status )
            IF_NOTOK_RETURN(status=1)

          case default

            write (gol,'("unsupported horizonal combination key for sh :",a)') hcomb; call goErr
            TRACEBACK; status=1; return

        end select

      ! ~~~ ?? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case default

        write (gol,'("unsupported input grid type `",a,"`")') tmm%buf_gridtype; call goErr
        TRACEBACK; status=1; return

    end select

    ! fill history:
    call Init( tmi, tmm%buf_tmi, status )
    call AddHistory( tmi, 'oper==hcomb,'//trim(hcomb), status )

    ! end timing:
    call GO_Timer_End( itim_transform_3dh, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Transform3Dh


  ! ***


  subroutine Transform3Dv( tmm, levi, nw, sp, ll, tmi, status )

    use PArray       , only : pa_SetShape
    use GO           , only : goErr, gol, goPr
    use GO           , only : GO_Timer_Start, GO_Timer_End
    use Grid         , only : TLevelInfo, FillLevels
    use Grid         , only : TllGridInfo, FillGrid, AreaOper
    use Grid         , only : Tgg2llFracs, Init, Done
    use Grid         , only : IntArea
    use Grid         , only : FracSum
    use tmm_param    , only : GetCombineKeys
    use tmm_info     , only : TMeteoInfo, AddHistory

    ! --- in/out --------------------------------

    type(TTmMeteo), intent(inout)           ::  tmm
    type(TLevelInfo), intent(in)            ::  levi
    character(len=1), intent(in)            ::  nw
    real, intent(in)                        ::  sp(:,:)
    real, intent(out)                       ::  ll(:,:,:)
    type(TMeteoInfo), intent(inout)         ::  tmi
    integer, intent(out)                    ::  status

    ! --- const ------------------------------------

    character(len=*), parameter ::  rname = mname//'/Transform3Dv'

    ! --- local -------------------------------

    character(len=10)       ::  hcomb, vcomb

    ! --- begin ----------------------------------

    ! start timing:
    call GO_Timer_Start( itim_transform_3dv, status )
    IF_NOTOK_RETURN(status=1)

    ! define how to combine horizontal cells and vertical levels:
    call GetCombineKeys( hcomb, vcomb, tmm%buf_paramkey, status )
    IF_NOTOK_RETURN(status=1)

    ! collect or copy levels, eventually reversed, from leviX/llX into levi/ll :
    !write (gol,'(a,": vertical ...")') rname
    call FillLevels( levi, nw, sp, ll, tmm%buf_levi, tmm%llX, vcomb, status )
    IF_NOTOK_RETURN(status=1)

    ! add history:
    call AddHistory( tmi, 'oper==vcomb,'//trim(vcomb), status )

    ! end timing:
    call GO_Timer_End( itim_transform_3dv, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Transform3Dv


  ! ==================================================================
  ! ===
  ! ===  read ll field
  ! ===
  ! ==================================================================


  subroutine tmm_ReadField_2d( tmm, archivekey, paramkey, unit, tday, t1, t2, &
                                    lli, nuv, ll, tmi, status )

    use PArray  , only : pa_Init, pa_Done
    use GO      , only : goErr, gol, goPr, Pretty
    use GO      , only : TDate, wrtgol
    use GO      , only : GO_Timer_Start, GO_Timer_End
    use Grid    , only : TllGridInfo
    use tmm_info, only : TMeteoInfo
    use dims,     only : okdebug_tmm

    ! --- in/out --------------------------------

    type(TTmMeteo), intent(inout)           ::  tmm
    character(len=*), intent(in)            ::  archivekey
    character(len=*), intent(in)            ::  paramkey
    character(len=*), intent(in)            ::  unit
    type(TDate), intent(in)                 ::  tday, t1, t2
    type(TllGridInfo), intent(in)           ::  lli
    character(len=1), intent(in)            ::  nuv
    real, intent(out)                       ::  ll(:,:)
    type(TMeteoInfo), intent(out)           ::  tmi
    integer, intent(out)                    ::  status

    ! --- const ------------------------------------

    character(len=*), parameter ::  rname = mname//'/tmm_ReadField_2d'

    ! --- local -------------------------------

    character(len=10)       ::  hcomb, vcomb

    ! --- begin ----------------------------------

    ! start timing:
    call GO_Timer_Start( itim_readfield_2d, status )
    IF_NOTOK_RETURN(status=1)

    if (okdebug_tmm) then
        write(gol, '(a, " : ", a, " from ", a, " to ", a)') rname, trim(paramkey), Pretty(t1), Pretty(t2)
        call goPr
    end if

    !call wrtgol( 'tmm read  : '//trim(paramkey)//' (', tday, ') ', t1, '  -  ', t2 ); call goPr

    ! check shape of target grid:
    if ( ((nuv == 'n') .and. ((size(ll,1) /= lli%nlon  ) .or. (size(ll,2) /= lli%nlat  ))) .or. &
         ((nuv == 'u') .and. ((size(ll,1) /= lli%nlon+1) .or. (size(ll,2) /= lli%nlat  ))) .or. &
         ((nuv == 'v') .and. ((size(ll,1) /= lli%nlon  ) .or. (size(ll,2) /= lli%nlat+1))) ) then
      write (gol,'("target array does not match with grid definition:")'); call goErr
      write (gol,'("  param  : ",a          )') paramkey; call goErr
      write (gol,'("  lli    : ",i3," x ",i3)') lli%nlon, lli%nlat; call goErr
      write (gol,'("  nuv    : ",a          )') nuv; call goErr
      write (gol,'("  ll     : ",i3," x ",i3)') shape(ll); call goErr
      TRACEBACK; status=1; return
    end if

    !
    ! convert grid
    !

    ! standard  values?
    if ( trim(archivekey) == 'standard' ) then

      ! dummy value
      ll = 0.0

    else

      ! read raw field in buffer
      call FillBuffer( tmm, archivekey, paramkey, unit, tday, t1, t2, nuv, 'n', status )
      IF_NOTOK_RETURN(status=1)

      ! horizontal interpolation:
      call Transform2D( tmm, lli, nuv, tmi, status )
      IF_NOTOK_RETURN(status=1)

      ! expecting single level ...
      if ( size(tmm%llX,3) /= 1 ) then
        write (gol,'("expecting single level:")'); call goErr
        write (gol,'("  paramkey  : ",a)') paramkey; call goErr
        write (gol,'("  levels    : ",a)') size(tmm%llX,3); call goErr
        TRACEBACK; status=1; return
      end if

      ! extract first level
      ll = tmm%llX(:,:,1)

    end if

    !
    ! unit conversions, truncations, etc
    !

    select case ( paramkey )
      case ( 'oro', 'cp', 'lsp' )
        ! set minium; some negative values if made from spectral field:
        ll = max( 0.0, ll )
    end select

    !
    ! done
    !

    ! end timing:
    call GO_Timer_End( itim_readfield_2d, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine tmm_ReadField_2d


  ! ***


  subroutine tmm_ReadField_3d( tmm, archivekey, paramkey, unit, tday, t1, t2, &
                                    lli, nuv, levi, nw, sp, ll, tmi, status )

    use GO      , only : goErr, gol, goPr, Pretty
    use GO      , only : TDate, wrtgol
    use GO      , only : GO_Timer_Start, GO_Timer_End
    use Binas   , only : p_global
    use Grid    , only : TllGridInfo, TLevelInfo
    use tmm_info, only : TMeteoInfo
    use dims,     only : okdebug_tmm

    ! --- in/out --------------------------------

    type(TTmMeteo), intent(inout)           ::  tmm
    character(len=*), intent(in)            ::  archivekey
    character(len=*), intent(in)            ::  paramkey, unit
    type(TDate), intent(in)                 ::  tday, t1, t2
    type(TllGridInfo), intent(in)           ::  lli
    character(len=1), intent(in)            ::  nuv
    type(TLevelInfo), intent(in)            ::  levi
    character(len=1), intent(in)            ::  nw
    real, intent(out)                       ::  sp(:,:)
    real, intent(out)                       ::  ll(:,:,:)
    type(TMeteoInfo), intent(out)           ::  tmi
    integer, intent(out)                    ::  status

    ! --- const ------------------------------------

    character(len=*), parameter ::  rname = mname//'/tmm_ReadField_3d'

    ! --- begin ----------------------------------

    ! start timing:
    call GO_Timer_Start( itim_readfield_3d, status )
    IF_NOTOK_RETURN(status=1)

    if (okdebug_tmm) then
        write(gol, '(a, " : ", a, " from ", a, " to ", a)') rname, trim(paramkey), Pretty(t1), Pretty(t2)
        call goPr
    end if

!    call wrtgol( 'tmm read  : '//trim(paramkey)//' ', t1, '  -  ', t2 ); call goPr

    ! check shape of target grid:
    if ( ((nuv == 'n') .and. ((size(ll,1) /= lli%nlon  ) .or. (size(ll,2) /= lli%nlat  ))) .or. &
         ((nuv == 'u') .and. ((size(ll,1) /= lli%nlon+1) .or. (size(ll,2) /= lli%nlat  ))) .or. &
         ((nuv == 'v') .and. ((size(ll,1) /= lli%nlon  ) .or. (size(ll,2) /= lli%nlat+1))) .or. &
         ((nuv == 'n') .and. ((size(sp,1) /= lli%nlon  ) .or. (size(sp,2) /= lli%nlat  ))) .or. &
         ((nuv == 'u') .and. ((size(sp,1) /= lli%nlon+1) .or. (size(sp,2) /= lli%nlat  ))) .or. &
         ((nuv == 'v') .and. ((size(sp,1) /= lli%nlon  ) .or. (size(sp,2) /= lli%nlat+1))) .or. &
         ((nw  == 'n') .and. (size(ll,3) /= levi%nlev  )) .or. &
         ((nw  == 'w') .and. (size(ll,3) /= levi%nlev+1)) ) then
      write (gol,'("target arrays do not match with grid definition:")'); call goErr
      write (gol,'("  param  : ",a          )') paramkey   ; call goErr
      write (gol,'("  lli    : ",i3," x ",i3         )') lli%nlon, lli%nlat; call goErr
      write (gol,'("  nuv    : ",a                   )') nuv; call goErr
      write (gol,'("  levi   : ",i3                  )') levi%nlev; call goErr
      write (gol,'("  nw     : ",a                   )') nw; call goErr
      write (gol,'("  sp     : ",i3," x ",i3         )') shape(sp); call goErr
      write (gol,'("  ll     : ",i3," x ",i3," x ",i3)') shape(ll); call goErr
      TRACEBACK; status=1; return
    end if

    !
    ! read, convert
    !

    ! standard  values?
    if ( trim(archivekey) == 'standard' ) then

      ! dummy values
      sp = p_global
      ll = 0.0

    else

      ! read raw field in buffer
      call FillBuffer( tmm, archivekey, paramkey, unit, tday, t1, t2, nuv, nw, status )
      IF_NOTOK_RETURN(status=1)

      ! 3d horizontal conversion:
      call Transform3Dh( tmm, lli, nuv, levi, nw, sp, tmi, status )
      IF_NOTOK_RETURN(status=1)

      ! 3d vertical conversion:
      call Transform3Dv( tmm, levi, nw, sp, ll, tmi, status )
      IF_NOTOK_RETURN(status=1)

    end if

    !
    ! unit conversion, extreme values
    !

    select case ( paramkey )
      case ( 'Q', 'CLWC', 'CIWC' )
        ! set minimum, stored values could be slightly negative
        ll = max( 0.0, ll )
        ll = max( 0.0, ll )
    end select

    !
    ! done
    !

    ! start timing:
    call GO_Timer_End( itim_readfield_3d, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine tmm_ReadField_3d



  ! ==================================================================
  ! ===
  ! ===  VO/D/LNSP -> pu/pv/sp
  ! ===
  ! ==================================================================



  subroutine tmm_Read_SP( tmm, archivekey, paramname, paramunit, tday, t1, t2, &
                                lli, sp, tmi, status )

    use binas   , only : p_global
    use GO      , only : goErr, gol, goPr
    use GO      , only : TDate, Get, wrtgol
    use GO      , only : goSplitLine
    use Grid    , only : TllGridInfo
    use Grid    , only : TshGridInfo, TshGrid, Init, Done, Set
    use Grid    , only : IntArea
    use Grid    , only : Tgg2llFracs, FracSum
    use tmm_info, only : TMeteoInfo, Init, AddHistory

    ! --- in/out --------------------------------

    type(TTmMeteo), intent(inout)     ::  tmm
    character(len=*), intent(in)      ::  archivekey, paramname, paramunit
    type(TDate), intent(in)           ::  tday, t1, t2
    type(TllGridInfo), intent(in)     ::  lli
    real, intent(out)                 ::  sp(:,:)       ! Pa
    type(TMeteoInfo), intent(out)     ::  tmi
    integer, intent(out)              ::  status

    ! --- const ------------------------------------

    character(len=*), parameter ::  rname = mname//'/tmm_Read_SP'

    ! --- local -------------------------------

    character(len=10)         ::  sourcetype
    character(len=MAX_RCKEY_LEN)        ::  sourcename
    integer                   ::  hour
    type(Tgg2llFracs)         ::  gg2ll

    ! --- begin -------------------------------

    ! split source key in type and name:
    call goSplitLine( archivekey, sourcetype, ':', sourcename, status )
    IF_NOTOK_RETURN(status=1)

    ! input TMPP fields or raw prism fields ?
    select case ( sourcetype )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! standard
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'standard' )

        !write (gol,'("tmm read  : sp standard global mean pressure")'); call goPr

        ! fill field with global mean pressure:
        sp = p_global

        ! set history info
        call Init( tmi, 'sp', 'Pa', status )
        call AddHistory( tmi, 'standard', status )


      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! read directly from tmpp hdf file
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'tmpp' )

        ! select parameter given hour:
        call Get( t1, hour=hour )
        if ( modulo(hour,6) == 3 ) then

          ! hour = 21/03/..  : uvsp files
          call ReadField( tmm, archivekey, 'sp', 'Pa', tday, t1, t1, &
                               lli, 'n', sp, tmi, status )
          IF_NOTOK_RETURN(status=1)

        else

          ! hour = 00/06/..  : spm files
          call ReadField( tmm, archivekey, 'spm', 'Pa', tday, t1, t1, &
                               lli, 'n', sp, tmi, status )
          IF_NOTOK_RETURN(status=1)

        end if

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! read directly from tm5 hdf file
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'tm5-hdf', 'tm5-nc' )

        ! sp always storred in 'sp' files ...
        call ReadField( tmm, archivekey, paramname, paramunit, tday, t1, t1, &
                             lli, 'n', sp, tmi, status )
        IF_NOTOK_RETURN(status=1)

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! convert spectral lnsp
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'ecmwf-tmpp', 'ncep-cdc', 'ncep-gfs', 'msc-data', 'prism' )

        !call wrtgol( 'tmm read  : sp ', t1, '  -  ', t2 ); call goPr

        ! read LNSP
        call FillBuffer( tmm, archivekey, 'LNSP', '1', tday, t1, t1, 'n', 'n', status )
        IF_NOTOK_RETURN(status=1)

        ! should be spectral ...
        if ( tmm%buf_gridtype /= 'sh' ) then
          write (gol,'("expecting sh field, not ",a)') tmm%buf_gridtype
          TRACEBACK; status=1; return
        end if

        ! aera average exp(lnsp)
        call IntArea( 'exp,aver', tmm%buf_shi, tmm%buf_sh(:,1), lli, sp, status )
        IF_NOTOK_RETURN(status=1)

        ! set history info
        call Init( tmi, 'sp', 'Pa', status, (/tmm%buf_tmi/) )
        call AddHistory( tmi, 'oper==exp,aver', status )


      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! convert gaussian sp or spectral lnsp
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'ecmwf-tm5' )

        select case ( paramname )

          case ( 'sp', 'SP' )

            !call wrtgol( 'tmm read  : sp ', t1, '  -  ', t2 ); call goPr

            ! read LNSP
            call FillBuffer( tmm, archivekey, 'LNSP', '1', tday, t1, t1, 'n', 'n', status )
            IF_NOTOK_RETURN(status=1)

            ! should be spectral ...
            if ( tmm%buf_gridtype /= 'sh' ) then
              write (gol,'("expecting sh field, not ",a)') tmm%buf_gridtype
              TRACEBACK; status=1; return
            end if

            ! aera average exp(lnsp)
            call IntArea( 'exp,aver', tmm%buf_shi, tmm%buf_sh(:,1), lli, sp, status )
            IF_NOTOK_RETURN(status=1)

            ! set history info
            call Init( tmi, 'sp', 'Pa', status, (/tmm%buf_tmi/) )
            call AddHistory( tmi, 'oper==exp,aver', status )

          case ( 'sps', 'SPS' )

            !call wrtgol( 'tmm read  : sps ', t1, '  -  ', t2 ); call goPr

            ! read gg field
            call FillBuffer( tmm, archivekey, 'sp', 'Pa', tday, t1, t1, 'n', 'n', status )
            IF_NOTOK_RETURN(status=1)

            ! should be gg ...
            if ( tmm%buf_gridtype /= 'gg' ) then
              write (gol,'("expecting gg field, not ",a)') tmm%buf_gridtype
              TRACEBACK; status=1; return
            end if

            ! determine fraction
            call Init( gg2ll, tmm%buf_ggi, lli, status )
            IF_NOTOK_RETURN(status=1)

            ! take fractions of overlapping cells
            call FracSum( gg2ll, tmm%buf_gg(:,1), sp, status, 'area-aver' )
            IF_NOTOK_RETURN(status=1)

            ! done
            call Done( gg2ll )

            ! set history info
            call Init( tmi, 'sp', 'Pa', status, (/tmm%buf_tmi/) )
            call AddHistory( tmi, 'oper==area-aver', status )

          case default

            write (gol,'("unsupported param `",a,"` for source type `",a,"`")') &
                             trim(paramname), trim(sourcetype); call goErr
            TRACEBACK; status=1; return

        end select

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! error ...
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case default

        write (gol,'("unsupported source type `",a,"`")') trim(sourcetype); call goErr
        TRACEBACK; status=1; return

    end select

    ! ok
    status = 0

  end subroutine tmm_Read_SP


  ! *****************************************************************


  subroutine tmm_Read_MFUV( tmm, archivekey, tday, t1, t2, &
                                lli, levi, &
                                spu, mfu, mfu_tmi, &
                                spv, mfv, mfv_tmi, &
                                status )

    use Binas   , only : R => ae
    use Binas   , only : g => grav
    use Binas   , only : p_global
    use GO      , only : goErr, gol, goPr
    use GO      , only : TDate, wrtgol
    use GO      , only : goSplitLine
    use Grid    , only : TllGridInfo, TLevelInfo
    use Grid    , only : TshGridInfo, TshGrid, Init, Done, Set
    use Grid    , only : IntLat, IntLon, vod2uv
    use Grid    , only : FillLevels
    use tmm_info, only : TMeteoInfo, Init, Done, AddHistory

    ! --- in/out --------------------------------

    type(TTmMeteo), intent(inout)     ::  tmm
    character(len=*), intent(in)      ::  archivekey
    type(TDate), intent(in)           ::  tday, t1, t2
    type(TllGridInfo), intent(in)     ::  lli
    type(TLevelInfo), intent(in)      ::  levi
    real, intent(out)                 ::  spu(:,:) , spv(:,:)    ! Pa
    real, intent(out)                 ::  mfu(:,:,:), mfv(:,:,:)   ! kg/s
    type(TMeteoInfo), intent(out)     ::  mfu_tmi, mfv_tmi
    integer, intent(out)              ::  status

    ! --- const ------------------------------------

    character(len=*), parameter ::  rname = mname//'/tmm_Read_MFUV'

    ! --- local -------------------------------

    character(len=10)        ::  sourcetype
    character(len=MAX_RCKEY_LEN)       ::  sourcename

    type(TshGridInfo)         ::  shi
    integer                   ::  shT
    integer                   ::  nlev

    ! spectral fields
    type(TshGrid)             ::  LNSP_sh

    complex , pointer         ::  D_sh(:,:), VO_sh(:,:)
    complex , pointer         ::  U_sh(:,:), V_sh(:,:)
    complex , pointer         ::  Help_LNSP_sh(:,:)

    ! mfu/mfv on source levels
    real, allocatable         ::  mfuX(:,:,:)
    real, allocatable         ::  mfvX(:,:,:)

    ! loops etc
    integer                   ::  l

    ! extra info
    type(TMeteoInfo)          ::  LNSP_tmi
    type(TMeteoInfo)          ::  vo_tmi, div_tmi, u_tmi, v_tmi

    ! temporary arrays
    real, allocatable         ::  tmp_sp(:,:,:)

    ! --- begin -------------------------------

    ! split source key in type and name:
    call goSplitLine( archivekey, sourcetype, ':', sourcename, status )
    IF_NOTOK_RETURN(status=1)

    ! input TMPP fields or raw prism fields ?
    select case ( sourcetype )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! standaard values
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'standard' )

        ! fill field with global mean pressure:
        spu = p_global
        spv = p_global

        ! dummy values:
        mfu = 0.0
        mfv = 0.0

        ! set history info:
        call Init( mfu_tmi, 'mfu', 'kg/s', status )
        call AddHistory( mfu_tmi, 'standard', status )
        call Init( mfu_tmi, 'mfv', 'kg/s', status )
        call AddHistory( mfv_tmi, 'standard', status )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! read directly from hdf/netcdf file
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'tmpp', 'tm5-hdf','tm5-nc' )

        ! flux through and surface pressure on west and east boundaries:
        call ReadField( tmm, archivekey, 'mfu', 'kg/s', tday, t1, t2, &
                               lli, 'u', levi, 'n', spu, mfu, mfu_tmi, status )
        IF_NOTOK_RETURN(status=1)

        ! flux through and surface pressure on south and north boundaries:
        call ReadField( tmm, archivekey, 'mfv', 'kg/s', tday, t1, t2, &
                               lli, 'v', levi, 'n', spv, mfv, mfv_tmi, status )
        IF_NOTOK_RETURN(status=1)

#if defined(with_tmm_ecmwf) || defined(with_tmm_ncep) || defined(with_msc) || defined(with_prism)
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! convert spectral fields
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'ecmwf-tmpp', 'ecmwf-tm5', 'ncep-cdc', 'ncep-gfs', 'msc-data', 'prism' )

        !call wrtgol( 'tmm read  : mfuv ', t1, '  -  ', t2 ); call goPr

        !
        ! read spectral fields
        !

        ! read LNSP
        call FillBuffer( tmm, archivekey, 'LNSP', 'Pa', tday, t1, t1, 'n', 'n', status )
        IF_NOTOK_RETURN(status=1)

          ! check ...
          call CheckBuffer( tmm, status, gridtype='sh' )
          IF_NOTOK_RETURN(status=1)

          ! extract grid size
          call Init( shi, tmm%buf_shi, status )
          IF_NOTOK_RETURN(status=1)
          shT  = tmm%buf_shi%T

          ! copy 1d spectral lnsp from buffer:
          call Init( LNSP_sh )
          call Set( LNSP_sh, tmm%buf_shi%T, tmm%buf_sh(:,1) )
          LNSP_tmi = tmm%buf_tmi

        ! read VO
        call FillBuffer( tmm, archivekey, 'VO', '1/s', tday, t1, t2, 'n', 'n', status )
        IF_NOTOK_RETURN(status=1)

          ! check ...
          call CheckBuffer( tmm, status, gridtype='sh', shT=shi%T )
          IF_NOTOK_RETURN(status=1)
          ! extract grid size
          nlev = tmm%buf_levi%nlev

          ! copy 3d spectral field from buffer:
          allocate( VO_sh(shi%np,nlev) )
          VO_sh = tmm%buf_sh
          vo_tmi = tmm%buf_tmi

        ! read D
        call FillBuffer( tmm, archivekey, 'D', '1/s', tday, t1, t2, 'n', 'n', status )
        IF_NOTOK_RETURN(status=1)

          ! check ...
          call CheckBuffer( tmm, status, gridtype='sh', shT=shi%T, nlev=nlev )
          IF_NOTOK_RETURN(status=1)

          ! copy 3d spectral field from buffer:
          allocate( D_sh(shi%np,nlev) )
          D_sh = tmm%buf_sh
          div_tmi = tmm%buf_tmi

        !
        ! compute U/V from VO/D
        !
        allocate( U_sh(shi%np,nlev) )
        allocate( V_sh(shi%np,nlev) )
        allocate( Help_LNSP_sh(shi%np,1) )

        !$OMP PARALLEL &
        !$OMP       default ( none ) &
        !$OMP       shared  ( nlev, shi, VO_sh, D_sh, U_sh, V_sh ) &
        !$OMP       private ( l )
        do l = 1, nlev
           call vod2uv( shi, VO_sh(:,l), D_sh(:,l), shi, U_sh(:,l), V_sh(:,l) )
        enddo
        !$OMP END PARALLEL
        ! history ...
        call Init( u_tmi, 'U', 'm/s', status, (/vo_tmi,div_tmi/) )
        call Init( v_tmi, 'V', 'm/s', status, (/vo_tmi,div_tmi/) )

        !
        !  mfu  =  R/g int  U (da+db*exp(LNSP)) / cos(lat)  dlat
        !

        allocate( mfuX(0:lli%nlon,lli%nlat,nlev) )
        allocate( tmp_sp(0:lli%nlon,lli%nlat,1) )

        ! integral over boundary:
        call IntLat(  '(da+exp*db)/cos', shi, nlev, U_sh, LNSP_sh%c, &
                     tmm%buf_levi%da, tmm%buf_levi%db, lli, mfuX, status )
        IF_NOTOK_RETURN(status=1)
        mfuX = mfuX * R/g

        ! average surface pressure on boundary:
        Help_LNSP_sh(:,1)=LNSP_sh%c(:)
        call IntLat( 'exp(H),aver',shi,1, Help_LNSP_sh, LNSP_sh%c, (/0.0/), (/0.0/), lli, tmp_sp, status )
        IF_NOTOK_RETURN(status=1)
        spu = tmp_sp(:,:,1)

        ! collect levels:
        call FillLevels( levi, 'n', spu, mfu, tmm%buf_levi, mfuX, 'sum', status )
        IF_NOTOK_RETURN(status=1)

        ! clear:
        deallocate( mfuX )
        deallocate( tmp_sp )

        ! info ...
        call Init( mfu_tmi, 'mfu', 'kg/s', status, (/u_tmi,LNSP_tmi/) )
        call AddHistory( mfu_tmi, 'oper==intlat', status )
        call AddHistory( mfu_tmi, 'oper==collectlevels', status )

        !
        !  mfv  =  R/g int  V (da+db*exp(LNSP)) dlon
        !

        allocate( mfvX(lli%nlon,0:lli%nlat,nlev) )
        allocate( tmp_sp(lli%nlon,0:lli%nlat,1) )

        ! integral over boundary:
        call IntLon( '(da+exp*db)',shi,nlev, V_sh, LNSP_sh%c, &
                     tmm%buf_levi%da, tmm%buf_levi%db, lli, mfvX, status )
        IF_NOTOK_RETURN(status=1)
        mfvX = mfvX * R/g

        ! average surface pressure on boundary:
        Help_LNSP_sh(:,1)=LNSP_sh%c(:)
        call IntLon( 'exp(H),aver', shi,1,Help_LNSP_sh, LNSP_sh%c, (/0.0/), (/0.0/), lli, tmp_sp, status )
        IF_NOTOK_RETURN(status=1)
        spv = tmp_sp(:,:,1)

        ! collect levels:
        call FillLevels( levi, 'n', spv, mfv, tmm%buf_levi, mfvX, 'sum', status )
        IF_NOTOK_RETURN(status=1)

        ! clear:
        deallocate( Help_LNSP_sh)
        deallocate( mfvX )
        deallocate( tmp_sp )

        ! info ...
        call Init( mfv_tmi, 'mfv', 'kg/s', status, (/v_tmi,LNSP_tmi/) )
        call AddHistory( mfv_tmi, 'oper==intlon', status )
        call AddHistory( mfv_tmi, 'oper==collectlevels', status )

        !
        ! done
        !

        call Done( LNSP_sh )
        deallocate( VO_sh )
        deallocate(  D_sh )
        deallocate(  U_sh )
        deallocate(  V_sh )

        call Done( lnsp_tmi, status )
        call Done(   vo_tmi, status )
        call Done(  div_tmi, status )
        call Done(    u_tmi, status )
        call Done(    v_tmi, status )
#endif

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! error
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case default

        write (gol,'("unsupported source type `",a,"`")') trim(sourcetype); call goErr
        TRACEBACK; status=1; return

    end select

    ! ok
    status = 0

  end subroutine tmm_Read_MFUV


  ! *****************************************************************


  subroutine tmm_Read_MFW( tmm, archivekey, tday, t1, t2, lli, levi, &
                                sp, mfw, tsp, tmi, status )

    use Binas   , only : R => ae
    use Binas   , only : g => grav
    use Binas   , only : p_global
    use GO      , only : goErr, gol, goPr
    use GO      , only : TDate, wrtgol
    use GO      , only : goSplitLine
    use Grid    , only : TllGridInfo, TLevelInfo
    use Grid    , only : TshGridInfo, TshGrid, Init, Done, Set
    use Grid    , only : vod2uv, Nabla, IntArea
    use Grid    , only : FillLevels
    use Grid    , only : AreaOper
    use tmm_info, only : TMeteoInfo, Init, Done, AddHistory

    ! --- in/out --------------------------------

    type(TTmMeteo), intent(inout)     ::  tmm
    character(len=*), intent(in)      ::  archivekey
    type(TDate), intent(in)           ::  tday, t1, t2
    type(TllGridInfo), intent(in)     ::  lli
    type(TLevelInfo), intent(in)      ::  levi
    real, intent(out)                 ::  sp(:,:)      ! Pa
    real, intent(out)                 ::  mfw(:,:,:)   ! kg/s
    real, intent(out)                 ::  tsp(:,:)     ! tendency of surface pressure Pa/s
    type(TMeteoInfo), intent(out)     ::  tmi
    integer, intent(out)              ::  status

    ! --- const ------------------------------------

    character(len=*), parameter ::  rname = mname//'/tmm_Read_MFW'

    ! --- local -------------------------------

    character(len=10)        ::  sourcetype
    character(len=MAX_RCKEY_LEN)       ::  sourcename

    type(TshGridInfo)         ::  shi
    integer                   ::  shT
    integer                   ::  nlev

    ! spectral fields
    type(TshGrid)             ::  LNSP_sh
    complex , pointer         ::  D_sh(:,:), VO_sh(:,:)
    complex , pointer         ::  U_sh(:,:), V_sh(:,:)

    ! nabla.lnps
    type(TshGrid)             ::  NabLNSP_sh(2)

    ! integrated Omega arrays:
    real, pointer             ::  IIOmega (:,:,:)
    real, pointer             ::  IIOmega2(:,:,:)

    ! mfw on source levels
    real, pointer             ::  mfwX(:,:,:)

    ! loops etc
    integer                   ::  l

    ! extra info
    type(TMeteoInfo)          ::  LNSP_tmi
    type(TMeteoInfo)          ::  vo_tmi, div_tmi, u_tmi, v_tmi

    ! --- begin -------------------------------

    ! split source key in type and name:
    call goSplitLine( archivekey, sourcetype, ':', sourcename, status )
    IF_NOTOK_RETURN(status=1)

    ! input TMPP fields or raw prism fields ?
    select case ( sourcetype )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! standaard values
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'standard' )

        ! fill field with global mean pressure:
        sp = p_global

        ! dummy values:
        mfw = 0.0

        ! set history info
        call Init( tmi, 'mfw', 'kg/s', status )
        call AddHistory( tmi, 'standard', status )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! read directly from hdf file
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'tmpp', 'tm5-hdf', 'tm5-nc' )

        ! downward mass flux through bottom:
        call ReadField( tmm, archivekey, 'mfw', 'kg/s', tday, t1, t2, &
                               lli, 'n', levi, 'w', sp, mfw, tmi, status )
        IF_NOTOK_RETURN(status=1)

        ! dummy ...
        tsp = 0.0   ! Pa/s

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! convert spectral fields
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'ecmwf-tmpp', 'ecmwf-tm5', 'ncep-cdc', 'ncep-gfs', 'msc-data', 'prism' )

        !call wrtgol( 'tmm read  : mfw ', t1, '  -  ', t2 ); call goPr

        !
        ! read spectral fields
        !

        ! read LNSP
        call FillBuffer( tmm, archivekey, 'LNSP', '1', tday, t1, t1, 'n', 'n', status )
        IF_NOTOK_RETURN(status=1)

          ! check ...
          call CheckBuffer( tmm, status, gridtype='sh' )
          IF_NOTOK_RETURN(status=1)

          ! extract grid size
          call Init( shi, tmm%buf_shi, status )
          IF_NOTOK_RETURN(status=1)
          shT  = tmm%buf_shi%T

          ! copy 1d spectral lnsp from buffer:
          call Init( LNSP_sh )
          call Set( LNSP_sh, tmm%buf_shi%T, tmm%buf_sh(:,1) )
          LNSP_tmi = tmm%buf_tmi

        ! read VO
        call FillBuffer( tmm, archivekey, 'VO', '1/s', tday, t1, t2, 'n', 'n', status )
        IF_NOTOK_RETURN(status=1)

          ! check ...
          call CheckBuffer( tmm, status, gridtype='sh', shT=shi%T )
          IF_NOTOK_RETURN(status=1)

          ! extract grid size
          nlev = tmm%buf_levi%nlev

          ! copy 3d spectral field from buffer:
          allocate( VO_sh(shi%np,nlev) )
          VO_sh = tmm%buf_sh
          vo_tmi = tmm%buf_tmi

        ! read D
        call FillBuffer( tmm, archivekey, 'D', '1/s', tday, t1, t2, 'n', 'n', status )
        IF_NOTOK_RETURN(status=1)

          ! check ...
          call CheckBuffer( tmm, status, gridtype='sh', shT=shi%T, nlev=nlev )
          IF_NOTOK_RETURN(status=1)

          ! copy 3d spectral field from buffer:
          allocate( D_sh(shi%np,nlev) )
          D_sh = tmm%buf_sh
          div_tmi = tmm%buf_tmi

        !
        ! compute U/V from VO/D
        !

        allocate( U_sh(shi%np,nlev) )
        allocate( V_sh(shi%np,nlev) )

        print *, '        vod2uv ...'
        !$OMP PARALLEL &
        !$OMP       default ( none ) &
        !$OMP       shared  ( nlev, shi, VO_sh, D_sh, U_sh, V_sh ) &
        !$OMP       private ( l )
        do l = 1, nlev
          call vod2uv( shi, VO_sh(:,l), D_sh(:,l), shi, U_sh(:,l), V_sh(:,l) )
        end do
        !$OMP END PARALLEL

        ! history ...
        call Init( u_tmi, 'U', 'm/s', status, (/vo_tmi,div_tmi/) )
        call Init( v_tmi, 'V', 'm/s', status, (/vo_tmi,div_tmi/) )


        !
        ! allocate 3d fields
        !

        allocate( IIOmega (lli%nlon,lli%nlat,nlev) )
        allocate( IIOmega2(lli%nlon,lli%nlat,nlev) )

        allocate( mfwX(lli%nlon,lli%nlat,nlev+1) )


        !
        !   int int   D  (da+db*exp(LNSP)) cos(lat)  dA
        !

        IIOmega = 0.0
        call IntArea( 'F*(da+db*exp(H))*cos', shi, nlev, D_sh, LNSP_sh%c, LNSP_sh%c, &
                      tmm%buf_levi%da, tmm%buf_levi%db, lli, IIOmega, status )
        IF_NOTOK_RETURN(status=1)

        !
        !   int int   U  dLNSP1  exp(LNSP) db / cos(lat)  dA
        !   int int   V  dLNSP2  exp(LNSP) db / cos(lat)  dA
        !

        ! allocate
        call Init( NabLNSP_sh(1) )
        call Init( NabLNSP_sh(2) )

        ! compute nabla.lnsp :
        call Nabla( LNSP_sh, NabLNSP_sh )

        ! integral over cell area; add contributions
        IIOmega2 = 0.0
        call IntArea( 'F*G*(db*exp(H))/cos', shi, nlev, U_sh, NabLNSP_sh(1)%c, LNSP_sh%c, &
                      tmm%buf_levi%da, tmm%buf_levi%db, lli, IIOmega2, status )
        IF_NOTOK_RETURN(status=1)

        call IntArea( 'F*G*(db*exp(H))/cos', shi, nlev, V_sh, NabLNSP_sh(2)%c, LNSP_sh%c, &
                      tmm%buf_levi%da, tmm%buf_levi%db, lli, IIOmega2, status )
        IF_NOTOK_RETURN(status=1)

        ! deallocate
        call Done( NabLNSP_sh(1) )
        call Done( NabLNSP_sh(2) )

        !
        ! column integrated Omega
        !

        ! parent levels downwards or upwards ?
        if ( tmm%buf_levi%updo == 'd' ) then

          ! loop from top to bottom:
          do l = 1, nlev
            ! replace with contribution of current level:
            IIOmega(:,:,l) = (R**2)*IIOmega(:,:,l)/g  +  R*IIOmega2(:,:,l)/g
            ! add contribution of previous levels:
            if ( l > 1 ) then
              IIOmega(:,:,l) = IIOmega(:,:,l) + IIOmega(:,:,l-1)
            end if
          end do

        else

          ! loop from top to bottom:
          do l = nlev, 1, -1
            ! replace with contribution of current level:
            IIOmega(:,:,l) = (R**2)*IIOmega(:,:,l)/g  +  R*IIOmega2(:,:,l)/g
            ! add contribution of levels above:
            if ( l < nlev ) then
              IIOmega(:,:,l) = IIOmega(:,:,l) + IIOmega(:,:,l+1)
            end if
          end do

        end if

        !
        ! tendency of surface pressure
        !

        !  dps        1              dp
        !  ---  =  - int  nabla ( v ---- ) deta  =  - IIOmega(:,:,bot)
        !  dt       eta=0           deta

        ! parent levels downwards or upwards ?
        if ( tmm%buf_levi%updo == 'd' ) then
          tsp = -1.0 * IIOmega(:,:,nlev)   ! kg/s
        else
          tsp = -1.0 * IIOmega(:,:,1)   ! kg/s
        end if

        ! convert to Pa/s :
        call AreaOper( lli, tsp, '/', 'm2', status )   ! kg/m2/s
        IF_NOTOK_RETURN(status=1)
        tsp = tsp * g   ! Pa/s

        !
        ! compute vertical flux:
        !

        ! parent levels downwards or upwards ?
        if ( tmm%buf_levi%updo == 'd' ) then

          ! top hlev:
          mfwX(:,:,1) = 0.0    ! kg/s
          ! loop from top to bottom layer:
          do l = 1, nlev
            ! lay l bot hlev    surflay bot hlev    lay l bot hlev          lay l bot hlev
            ! 2 .. nlev+1       nlev                0 .. nlev-1             1 .. nlev
            mfwX(:,:,l+1)   =   IIOmega(:,:,nlev) * tmm%buf_levi%b(l)   -   IIOmega(:,:,l)
          end do

        else

          ! top hlev:
          mfwX(:,:,nlev+1) = 0.0    ! kg/s
          ! loop from top to bottom layer:
          do l = nlev, 1, -1
            ! lay l bot hlev    surflay bot hlev   lay l bot hlev            lay l bot hlev
            ! 1 .. nlev         nlev               0 .. nlev-1               1 .. nlev
            mfwX(:,:,l)    =    IIOmega(:,:,1)  *  tmm%buf_levi%b(l-1)   -   IIOmega(:,:,l)
          end do

        end if

        ! check: fluxh through bottom should be zero ...
        if ( maxval(abs(mfwX(:,:,1))) > 1.0 ) then
          write (gol,'("vert.flux through bottom half level should be zero ...")'); call goErr
          write (gol,'("  max value  : ",es12.4)') maxval(abs(mfwX(:,:,1))); call goErr
          TRACEBACK; status=1; return
        end if

        ! create history
        call Init( tmi, 'mfw', 'kg/s', status, (/lnsp_tmi,div_tmi,u_tmi,v_tmi/) )


        !
        ! collect levels
        !

        ! average surface pressure
        call IntArea( 'exp,aver', tmm%buf_shi, LNSP_sh%c, lli, sp, status )
        IF_NOTOK_RETURN(status=1)

        ! combine levels etc
        call FillLevels( levi, 'w', sp, mfw, tmm%buf_levi, mfwX, 'bottom', status )
        IF_NOTOK_RETURN(status=1)

        call AddHistory( tmi, 'oper==collectlevels', status )

        !
        ! upward flux
        !

        ! flux should be upwards (in direction of increasing level):
        mfw = - mfw

        call AddHistory( tmi, 'oper==upwards', status )


        !
        ! done
        !

        call Done( LNSP_sh )

        deallocate( VO_sh )
        deallocate(  D_sh )
        deallocate(  U_sh )
        deallocate(  V_sh )

        deallocate( IIOmega  )
        deallocate( IIOmega2 )
        deallocate( mfwX      )

        call Done( lnsp_tmi, status )
        call Done(   vo_tmi, status )
        call Done(  div_tmi, status )
        call Done(    u_tmi, status )
        call Done(    v_tmi, status )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! error
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case default

        write (gol,'("unsupported source type `",a,"`")') trim(sourcetype); call goErr
        TRACEBACK; status=1; return

    end select

    ! ok
    status = 0

  end subroutine tmm_Read_MFW


  ! *****************************************************************


  subroutine tmm_Read_TQ( tmm, archivekey_T, archivekey_Q, tday, t1, t2, lli, levi, &
                                sp, T, T_tmi, Q, Q_tmi, status )

    use GO      , only : goErr, gol, goPr
    use GO      , only : TDate, wrtgol
    use GO      , only : goSplitLine
    use Binas   , only : p_global
    use Phys    , only : RealTemperature
    use Grid    , only : TllGridInfo, TLevelInfo
    use tmm_info, only : TMeteoInfo, Init, Done, AddHistory

    ! --- in/out --------------------------------

    type(TTmMeteo), intent(inout)     ::  tmm
    character(len=*), intent(in)      ::  archivekey_T
    character(len=*), intent(in)      ::  archivekey_Q
    type(TDate), intent(in)           ::  tday, t1, t2
    type(TllGridInfo), intent(in)     ::  lli
    type(TLevelInfo), intent(in)      ::  levi
    real, intent(out)                 ::  sp(:,:)      ! Pa
    real, intent(out)                 ::  T(:,:,:)     ! K
    type(TMeteoInfo), intent(out)     ::  T_tmi
    real, intent(out)                 ::  Q(:,:,:)     ! kg/kg
    type(TMeteoInfo), intent(out)     ::  Q_tmi
    integer, intent(out)              ::  status

    ! --- const ------------------------------------

    character(len=*), parameter ::  rname = mname//'/tmm_Read_TQ'

    ! --- local -------------------------------

    character(len=10)        ::  sourcetype
    character(len=MAX_RCKEY_LEN)       ::  sourcename

    ! --- begin -------------------------------

    ! split source key in type and name:
    call goSplitLine( archivekey_T, sourcetype, ':', sourcename, status )
    IF_NOTOK_RETURN(status=1)

    ! input TMPP fields or raw prism fields ?
    select case ( sourcetype )

      case ( 'standard' )

        ! fill field with global mean pressure:
        sp = p_global
        ! dummy values:
        Q = 0.0
        T = 0.0

        ! set history info
        call Init( Q_tmi, 'Q', 'kg/kg', status )
        call AddHistory( Q_tmi, 'standard', status )
        call Init( T_tmi, 'T', 'K', status )
        call AddHistory( T_tmi, 'standard', status )

      case ( 'tmpp', 'tm5-hdf', 'tm5-nc', 'ecmwf-tmpp', 'ecmwf-tm5', 'msc-data' )

        ! humidity:
        call ReadField( tmm, archivekey_Q, 'Q', 'kg/kg', tday, t1, t2, &
                               lli, 'n', levi, 'n', sp, Q, Q_tmi, status )
        IF_NOTOK_RETURN(status=1)

        ! temperature:
        call ReadField( tmm, archivekey_T, 'T', 'K', tday, t1, t2, &
                               lli, 'n', levi, 'n', sp, T, T_tmi, status )
        IF_NOTOK_RETURN(status=1)

      case ( 'ncep-cdc', 'ncep-gfs' )

        ! humidity:
        call ReadField( tmm, archivekey_Q, 'Q', 'kg/kg', tday, t1, t2, &
                               lli, 'n', levi, 'n', sp, Q, Q_tmi, status )
        IF_NOTOK_RETURN(status=1)

        ! virtual temperature:
        call ReadField( tmm, archivekey_T, 'Tv', 'K', tday, t1, t2, &
                               lli, 'n', levi, 'n', sp, T, T_tmi, status )
        IF_NOTOK_RETURN(status=1)

        ! convert:
        T = RealTemperature( T, Q )

        ! info:
        call AddHistory( T_tmi, 'convert from virtual temperature', status )

      case default

        write (gol,'("unsupported temper source type `",a,"`")') trim(sourcetype); call goErr
        TRACEBACK; status=1; return

    end select

    !
    ! done
    !

    ! ok
    status = 0

  end subroutine tmm_Read_TQ


  ! *****************************************************************


  subroutine tmm_Read_CloudCovers( tmm, archivekey, tday, t1, t2, lli, levi, &
                        sp, cc, cc_tmi, cco, cco_tmi, ccu, ccu_tmi, status )

    use GO       , only : goErr, gol, goPr
    use GO       , only : TDate, wrtgol, gol, goPr
    use GO       , only : goSplitLine
    use Binas    , only : p_global
    use Grid     , only : TllGridInfo, TLevelInfo, FillLevels
    use Phys     , only : cf_overhead
    use tmm_info , only : TMeteoInfo, Init, AddHistory
    use tmm_param, only : GetCombineKeys

    ! --- in/out --------------------------------

    type(TTmMeteo), intent(inout)     ::  tmm
    character(len=*), intent(in)      ::  archivekey
    type(TDate), intent(in)           ::  tday, t1, t2
    type(TllGridInfo), intent(in)     ::  lli
    type(TLevelInfo), intent(in)      ::  levi
    real, intent(out)                 ::  sp(:,:)      ! Pa
    real, intent(out)                 ::  cc(:,:,:), cco(:,:,:), ccu(:,:,:)   ! 0-1
    type(TMeteoInfo), intent(out)     ::  cc_tmi, cco_tmi, ccu_tmi
    integer, intent(out)              ::  status

    ! --- const ------------------------------------

    character(len=*), parameter ::  rname = mname//'/tmm_Read_CloudCovers'

    ! --- local -------------------------------

    character(len=10)        ::  sourcetype
    character(len=MAX_RCKEY_LEN)       ::  sourcename

    integer                  ::  i, j, l, lme
    real, allocatable        ::   cc_col(:)
    real, allocatable        ::  cco_col(:), ccoX(:,:,:)
    real, allocatable        ::  ccu_col(:), ccuX(:,:,:)
    character(len=10)        ::  hcomb, vcomb

    ! --- begin -------------------------------

    ! check ...
    if ( any(shape(cco)/=shape(cc)) .or. any(shape(ccu)/=shape(cc)) ) then
      write (gol,'("output arrays should have same shape:")'); call goErr
      write (gol,'("   cc : ",3i4)') shape( cc); call goErr
      write (gol,'("  cco : ",3i4)') shape(cco); call goErr
      write (gol,'("  ccu : ",3i4)') shape(ccu); call goErr
      TRACEBACK; status=1; return
    end if

    ! split source key in type and name:
    call goSplitLine( archivekey, sourcetype, ':', sourcename, status )
    IF_NOTOK_RETURN(status=1)

    ! input TMPP fields or raw prism fields ?
    select case ( sourcetype )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! standard values
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'standard' )

        ! fill field with global mean pressure:
        sp = p_global
        ! no clouds ...
        cc  = 0.0
        cco = 0.0
        ccu = 0.0

        ! set history info
        call Init( cc_tmi, 'CC', '0-1', status )
        call AddHistory( cc_tmi, 'standard', status )
        call Init( cco_tmi, 'CCO', '0-1', status )
        call AddHistory( cc_tmi, 'standard', status )
        call Init( ccu_tmi, 'CCU', '0-1', status )
        call AddHistory( cc_tmi, 'standard', status )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! read directly from hdf file
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'tmpp', 'tm5-hdf', 'tm5-nc', 'prism' )

        ! cloud cover
        call ReadField( tmm, archivekey, 'CC', '1', tday, t1, t2, &
                               lli, 'n', levi, 'n', sp, cc,  cc_tmi, status )
        IF_NOTOK_RETURN(status=1)

        ! cloud cover overhead
        call ReadField( tmm, archivekey, 'CCO', '1', tday, t1, t2, &
                               lli, 'n', levi, 'n', sp, cco, cco_tmi, status )
        IF_NOTOK_RETURN(status=1)

        ! cloud cover underfeet
        call ReadField( tmm, archivekey, 'CCU', '1', tday, t1, t2, &
                               lli, 'n', levi, 'n', sp, ccu, ccu_tmi, status )
        IF_NOTOK_RETURN(status=1)

        ! set extrema; stored values could be slightly outside [0,1]
         cc = max( 0.0, min(  cc, 1.0 ) )
        cco = max( 0.0, min( cco, 1.0 ) )
        ccu = max( 0.0, min( ccu, 1.0 ) )


      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! convert from raw gg fields
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'ecmwf-tmpp', 'ecmwf-tm5' )

        !call wrtgol( 'tmm read  : cc   ', t1, '  -  ', t2 ); call goPr

        !
        ! read cc: gg, all levels
        !

        call FillBuffer( tmm, archivekey, 'CC', '1', tday, t1, t2, 'n', 'n', status )
        IF_NOTOK_RETURN(status=1)

        ! 3d horizontal conversion; result is stored in tmm%llX
        call Transform3Dh( tmm, lli, 'n', levi, 'n', sp, cc_tmi, status )
        IF_NOTOK_RETURN(status=1)

        !
        ! compute cloudcover overhead using all levels
        !

        ! unreduced number of levels:
        lme = tmm%buf_levi%nlev

        ! storage:
        allocate(  cc_col(lme) )
        allocate( cco_col(lme) )
        allocate( ccoX(lli%nlon,lli%nlat,lme) )
        allocate( ccu_col(lme) )
        allocate( ccuX(lli%nlon,lli%nlat,lme) )

        ! loop over grid cells
        do j = 1, lli%nlat
          do i = 1, lli%nlon

            ! overhead cloud cover:
            cc_col = tmm%llX(i,j,:)
            call cf_overhead ( lme, cc_col, cco_col )
            ccoX(i,j,:) = cco_col

            ! underfeet cloud cover; first reverse layers:
            do l = 1, lme
              cc_col(l) = tmm%llX(i,j,lme+1-l)
            end do
            call cf_overhead( lme, cc_col, ccu_col )
            do l = 1, lme
              ccuX(i,j,l) = ccu_col(lme+1-l)
            end do

          end do
        end do

        ! clear
        deallocate(  cc_col )
        deallocate( cco_col )
        deallocate( ccu_col )

        ! info on this operation:
        call AddHistory( cco_tmi, 'oper==cf_overhead', status )

        !
        ! 3d vertical conversions
        !

        ! convert from tmm%buf_llX to cc
        call Transform3Dv( tmm, levi, 'n', sp, cc, cc_tmi, status )
        IF_NOTOK_RETURN(status=1)

        ! store ccoX in buffer, convert from tmm%llX to cco
        call GetCombineKeys( hcomb, vcomb, 'cco', status )
        IF_NOTOK_RETURN(status=1)
        call FillLevels( levi, 'n', sp, cco, tmm%buf_levi, ccoX, vcomb, status )
        IF_NOTOK_RETURN(status=1)
        call AddHistory( cco_tmi, 'oper==vcomb,'//trim(vcomb), status )

        ! store ccuX in buffer, convert from tmm%llX to ccu
        call GetCombineKeys( hcomb, vcomb, 'ccu', status )
        IF_NOTOK_RETURN(status=1)
        call FillLevels( levi, 'n', sp, ccu, tmm%buf_levi, ccuX, vcomb, status )
        IF_NOTOK_RETURN(status=1)
        call AddHistory( ccu_tmi, 'oper==vcomb,'//trim(vcomb), status )

        ! clear
        deallocate( ccoX )

        ! set extrema; stored values could be slightly outside [0,1]
         cc = max( 0.0, min(  cc, 1.0 ) )
        cco = max( 0.0, min( cco, 1.0 ) )
        ccu = max( 0.0, min( ccu, 1.0 ) )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! error
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case default

        write (gol,'("unsupported source type `",a,"`")') trim(sourcetype); call goErr
        TRACEBACK; status=1; return

    end select

    ! ok
    status = 0

  end subroutine tmm_Read_CloudCovers


  ! *****************************************************************


  subroutine tmm_Read_Convec( tmm, archivekey, tday, t1, t2, lli, levi, &
                                omega, omega_tmi, &
                                gph, gph_tmi, &
                                sp, &
                                entu, entu_tmi, entd, entd_tmi, &
                                detu, detu_tmi, detd, detd_tmi, &
                                status )

    use GO      , only : goErr, gol, goPr
    use GO      , only : TDate, wrtgol
    use GO      , only : operator(-), operator(+), rTotal
    use GO      , only : goSplitLine, goVarValue
    use Binas   , only : p_global
    use Grid    , only : TllGridInfo
    use Grid    , only : TLevelInfo
    use tmm_info, only : TMeteoInfo, Init, AddHistory

    ! --- in/out --------------------------------

    type(TTmMeteo), intent(inout)     ::  tmm
    character(len=*), intent(in)      ::  archivekey
    type(TDate), intent(in)           ::  tday, t1, t2
    type(TllGridInfo), intent(in)     ::  lli
    type(TLevelInfo), intent(in)      ::  levi
    real, intent(in)                  ::  omega(:,:,:)   ! Pa/s
    type(TMeteoInfo), intent(in)      ::  omega_tmi
    real, intent(in)                  ::  gph(:,:,:)   ! m
    type(TMeteoInfo), intent(in)      ::  gph_tmi
    real, intent(out)                 ::  sp(:,:)      ! Pa
    real, intent(out)                 ::  entu(:,:,:), entd(:,:,:)   ! kg/m2/s
    type(TMeteoInfo), intent(out)     ::  entu_tmi, entd_tmi
    real, intent(out)                 ::  detu(:,:,:), detd(:,:,:)   ! kg/m2/s
    type(TMeteoInfo), intent(out)     ::  detu_tmi, detd_tmi
    integer, intent(out)              ::  status

    ! --- const ------------------------------------

    character(len=*), parameter ::  rname = mname//'/tmm_Read_Convec'

    ! --- local -------------------------------

    character(len=10)        ::  sourcetype
    character(len=MAX_RCKEY_LEN)       ::  sourcename

    integer                  ::  lout
    real, allocatable        ::  ll(:,:,:)

    character(len=8)         ::  method

    ! --- begin -------------------------------

    !call wrtgol( 'tmm read  convec ', t1, '  -  ', t2 ); call goPr

    ! number of levels in output arrays:
    lout = size(entu,3)

    ! check ...
    if ( (size(entd,3)/=lout) .or. &
         (size(detu,3)/=lout) .or. (size(detd,3)/=lout) ) then
      write (gol,'("output arrays should have same number of levels:")'); call goErr
      write (gol,'("  entu : ",i4)') size(entu,3); call goErr
      write (gol,'("  entd : ",i4)') size(entd,3); call goErr
      write (gol,'("  detu : ",i4)') size(detu,3); call goErr
      write (gol,'("  detd : ",i4)') size(detd,3); call goErr
      TRACEBACK; status=1; return
    end if

    ! split source key in type and name:
    call goSplitLine( archivekey, sourcetype, ':', sourcename, status )
    IF_NOTOK_RETURN(status=1)

    ! input TMPP fields or raw prism fields ?
    select case ( sourcetype )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! standard values
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'standard')

        ! fill field with global mean pressure:
        sp = p_global
        ! no convection ...
        entu  = 0.0
        entd  = 0.0
        detu  = 0.0
        detd  = 0.0

        ! set history info
        call Init( entu_tmi, 'eu', 'kg/m2/s', status )
        call AddHistory( entu_tmi, 'standard', status )
        call Init( entd_tmi, 'ed', 'kg/m2/s', status )
        call AddHistory( entd_tmi, 'standard', status )
        call Init( detu_tmi, 'du', 'kg/m2/s', status )
        call AddHistory( detu_tmi, 'standard', status )
        call Init( entd_tmi, 'dd', 'kg/m2/s', status )
        call AddHistory( detd_tmi, 'standard', status )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! read directly from hdf file
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'tmpp', 'tm5-hdf', 'tm5-nc' )

        ! full level output:
        allocate( ll(lli%nlon,lli%nlat,levi%nlev) )

        ! entrement updraught
        call ReadField( tmm, archivekey, 'eu', 'kg/m2/s', tday, t1, t2, &
                               lli, 'n', levi, 'n', sp, ll, entu_tmi, status )
        IF_NOTOK_RETURN(status=1)
        !
        entu = ll(:,:,1:lout)

        ! entrement downdraught
        call ReadField( tmm, archivekey, 'ed', 'kg/m2/s', tday, t1, t2, &
                               lli, 'n', levi, 'n', sp, ll, entd_tmi, status )
        IF_NOTOK_RETURN(status=1)
        !
        entd = ll(:,:,1:lout)

        ! detrement updraught
        call ReadField( tmm, archivekey, 'du', 'kg/m2/s', tday, t1, t2, &
                               lli, 'n', levi, 'n', sp, ll, detu_tmi, status )
        IF_NOTOK_RETURN(status=1)
        !
        detu = ll(:,:,1:lout)

        ! detrement downdraught
        call ReadField( tmm, archivekey, 'dd', 'kg/m2/s', tday, t1, t2, &
                               lli, 'n', levi, 'n', sp, ll, detd_tmi, status )
        IF_NOTOK_RETURN(status=1)
        !
        detd = ll(:,:,1:lout)

        ! clear
        deallocate( ll )


      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! ecmwf convective stuff
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'ecmwf-tm5','prism' )

        method = 'raw'
#ifdef with_prism
        method = 'ec-ll'
#endif
        call goVarValue( sourcename, ';', 'method', '=', method, status )
        IF_ERROR_RETURN(status=1)

        select case ( method )

#ifdef with_tmm_convec_raw

          case ( 'raw' )

            call tmm_Read_Convec_raw( tmm, archivekey, tday, t1, t2, lli, levi, &
                                    omega, omega_tmi, &
                                    sp, &
                                    entu, entu_tmi, entd, entd_tmi, &
                                    detu, detu_tmi, detd, detd_tmi, &
                                    status )
            IF_NOTOK_RETURN(status=1)

#endif

#ifdef with_tmm_convec_ec_gg

          case ( 'ec-gg' )

            stop 'check implementation of tmm_Read_Convec_EC_gg'

            !! read ec flux/detr, convert to tm entr/detr, average to tm ll
            !call tmm_Read_Convec_EC_gg( tmm, archivekey, tday, t1, t2, lli, levi, &
            !                        omega, omega_tmi, &
            !                        sp, &
            !                        entu, entu_tmi, entd, entd_tmi, &
            !                        detu, detu_tmi, detd, detd_tmi, &
            !                        status )
            !IF_NOTOK_RETURN(status=1)

#endif

#ifdef with_tmm_convec_ec

         case ( 'ec-ll' )

            ! read ec flux/detr, aver to tm ll, convert to tm entr/detr
            ! note: gph instead of omega
            call tmm_Read_Convec_EC( tmm, archivekey, tday, t1, t2, lli, levi, &
                                    gph, gph_tmi, &
                                    sp, &
                                    entu, entu_tmi, entd, entd_tmi, &
                                    detu, detu_tmi, detd, detd_tmi, &
                                    status )
            IF_NOTOK_RETURN(status=1)

#endif

         case default

           write (gol,'("unsupported convec method : ",a)') trim(method); call goErr
           TRACEBACK; status=1; return

        end select


#ifdef with_tmm_convec_raw

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! convert from raw fields (sh,gg)
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'ecmwf-tmpp', 'ncep-cdc', 'ncep-gfs', 'msc-data' )

        call tmm_Read_Convec_raw( tmm, archivekey, tday, t1, t2, lli, levi, &
                                omega, omega_tmi, &
                                sp, &
                                entu, entu_tmi, entd, entd_tmi, &
                                detu, detu_tmi, detd, detd_tmi, &
                                status )
        IF_NOTOK_RETURN(status=1)

#endif


      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! error
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case default

        write (gol,'("unsupported source type `",a,"`")') trim(sourcetype); call goErr
        TRACEBACK; status=1; return

    end select

    ! ok
    status = 0

  end subroutine tmm_Read_Convec


  ! *****************************************************************


#ifdef with_tmm_convec_raw

  subroutine tmm_Read_Convec_raw( tmm, archivekey, tday, t1, t2, lli, levi, &
                                omega, omega_tmi, &
                                sp, &
                                entu, entu_tmi, entd, entd_tmi, &
                                detu, detu_tmi, detd, detd_tmi, &
                                status )

    use Binas   , only : grav
    use GO      , only : goErr, gol, goPr
    use GO      , only : TDate, wrtgol, Get, IncrDate
    use GO      , only : operator(-), operator(+), rTotal, iTotal
    use GO      , only : goSplitLine
    use Phys    , only : mid2bound_uv, mid2bound_w, mid2bound_t
    use Phys    , only : mid2bound_q, mid2bound_z, mid2bound_p
    use Phys    , only : subscal, phlev_ec_rev, geopot
    use Phys    , only : subscal_2d
    use Phys    , only : RealTemperature
    use Phys    , only : GeoPotentialHeightB
    use Phys    , only : ECconv_to_TMconv
    use Grid    , only : TllGridInfo
    use Grid    , only : TLevelInfo, HPressure, FPressure, FillLevels
    use Grid    , only : TShGrid, VoD2UV
    use Grid    , only : TshGridInfo
    use Grid    , only : TggGridInfo, InterpolMask, Divergence, Tgg2llFracs, FracSum, AreaOper
    use Grid    , only : Init, Done, Set, Interpol
    use Grid    , only : NewInterpol
    use tmm_info, only : TMeteoInfo, Init, AddHistory

    ! --- in/out --------------------------------

    type(TTmMeteo), intent(inout)     ::  tmm
    character(len=*), intent(in)      ::  archivekey
    type(TDate), intent(in)           ::  tday, t1, t2
    type(TllGridInfo), intent(in)     ::  lli
    type(TLevelInfo), intent(in)      ::  levi
    real, intent(in)                  ::  omega(:,:,:)   ! Pa/s
    type(TMeteoInfo), intent(in)      ::  omega_tmi
    real, intent(out)                 ::  sp(:,:)      ! Pa
    real, intent(out)                 ::  entu(:,:,:), entd(:,:,:)   ! kg/m2/s
    type(TMeteoInfo), intent(out)     ::  entu_tmi, entd_tmi
    real, intent(out)                 ::  detu(:,:,:), detd(:,:,:)   ! kg/m2/s
    type(TMeteoInfo), intent(out)     ::  detu_tmi, detd_tmi
    integer, intent(out)              ::  status

    ! --- const ------------------------------------

    character(len=*), parameter ::  rname = mname//'/tmm_Read_Convec_raw'

    ! --- local -------------------------------

    character(len=10)        ::  sourcetype
    character(len=MAX_RCKEY_LEN)       ::  sourcename

    integer                  ::  lout
    real, pointer            ::  ll(:,:,:)

    ! timing
    integer                  ::  hour
    real                     ::  dhour
    integer                  ::  dth
    type(TDate)              ::  t1s, t2s
    type(TDate)              ::  t1m, t2m
    logical                  ::  skip_second

    ! loops etc
    integer                  ::  l, nlev
    integer                  ::  i, i1, i2, j

    ! spectral fields
    type(TshGridInfo)        ::  shi
    complex, pointer         ::  LNSP_sh(:)
    complex, pointer         ::  D_sh(:,:), VO_sh(:,:)
    complex, pointer         ::  U_sh(:), V_sh(:)

    ! gaussian grids
    integer                  ::  ggN
    logical                  ::  reduced
    type(TggGridInfo)        ::  ggi
    real, pointer            ::  gg(:)
    real, pointer            ::  slhf_gg(:)
    real, pointer            ::   oro_gg(:)
    real, pointer            ::    sp_gg(:)
    real, pointer            ::     m_gg(:)
    real, pointer            ::  u_gg(:,:)
    real, pointer            ::  v_gg(:,:)
    real, pointer            ::  w_gg(:,:)
    real, pointer            ::  t_gg(:,:)
    real, pointer            ::  q_gg(:,:)
    real, pointer            ::  p_gg(:,:)
    real, pointer            ::  z_gg(:,:)
    logical, pointer         ::  ggmask(:)
    real, pointer            ::  qam_gg(:,:), qac_gg(:,:)
    real, pointer            ::  entu_gg(:,:)
    real, pointer            ::  detu_gg(:,:)
    real, pointer            ::  entd_gg(:,:)
    real, pointer            ::  detd_gg(:,:)
    type(Tgg2llFracs)        ::  gg2ll
    real, pointer            ::  llX(:,:,:)

    ! subgrid stuff
    real                     ::  dt_sec
    real                     ::  evap
    real, pointer            ::  p_hlev(:)

    ! extra info
    type(TMeteoInfo)         ::  slhf_tmi, sp_tmi, oro_tmi
    type(TMeteoInfo)         ::  vo_tmi, div_tmi, u_tmi, v_tmi
    type(TMeteoInfo)         ::  w_tmi, t_tmi, q_tmi

    ! reversed levels
    type(TLevelInfo)         ::  leviX
    real, pointer            ::  aX(:), bX(:)
    real, pointer            ::  tmp_gg(:,:)

    ! e4 convection fields
    real, pointer            ::  udmf_gg(:,:)
    real, pointer            ::  ddmf_gg(:,:)
    real, pointer            ::  uddr_gg(:,:)
    real, pointer            ::  dddr_gg(:,:)
    type(TMeteoInfo)         ::  udmf_tmi, ddmf_tmi, uddr_tmi, dddr_tmi
    real, pointer            ::  ph_ec(:)
    real, pointer            ::  zh_ec(:)

    ! --- begin -------------------------------

    ! number of levels in output arrays:
    lout = size(entu,3)

    ! check ...
    if ( (size(entd,3)/=lout) .or. &
         (size(detu,3)/=lout) .or. (size(detd,3)/=lout) ) then
      write (gol,'("output arrays should have same number of levels:")'); call goErr
      write (gol,'("  entu : ",i4)') size(entu,3); call goErr
      write (gol,'("  entd : ",i4)') size(entd,3); call goErr
      write (gol,'("  detu : ",i4)') size(detu,3); call goErr
      write (gol,'("  detd : ",i4)') size(detd,3); call goErr
      TRACEBACK; status=1; return
    end if

    ! split source key in type and name:
    call goSplitLine( archivekey, sourcetype, ':', sourcename, status )
    IF_NOTOK_RETURN(status=1)

    ! length of time interval:
    dth  = iTotal( t2 - t1, 'hour' )
    dt_sec = dth * 3600.0

    ! 3 hourly or 6 hourly ?
    select case ( dth )

      !
      ! ~~ 6 hourly
      !

      case ( 6 )

        ! only around hours 00/06/12/ etc yet ...
        call Get( t1, hour=hour )
        if ( modulo(hour,6) /= 3 ) then
          write (gol,'("6 hourly convection only for intervals [21,03], [03,09], ...")'); call goErr
          TRACEBACK; status=1; return
        end if

        ! times for model level fields is mid of requested interval:
        t1m = t1 + IncrDate(sec=nint(dt_sec/2))
        t2m = t1m

        ! use 6 hour interval around requested (instant!) time to read slhf;
        ! for 3 hourly request, this means double reading
        ! (no problem since slhf is time average anyway)

        !
        ! slhf (gg)
        !

        ! first interval: [t1-dt/2,t1]
        t1s = t1 - IncrDate(sec=int(dt_sec/2))
        t2s = t1

        !call wrtgol( 'tmm read  : slhf ', t1s, '  -  ', t2s ); call goPr

        ! read first slhf in buffer (W/m2 time aver):
        call FillBuffer( tmm, archivekey, 'slhf', 'W/m2', tday, t1s, t2s, 'n', 'n', status )
        IF_NOTOK_RETURN(status=1)

          ! check ...
          call CheckBuffer( tmm, status, gridtype='gg' )
          IF_NOTOK_RETURN(status=1)

          ! extract grid size
          ggN     = tmm%buf_ggi%N
          reduced = tmm%buf_ggi%reduced

          ! setup gg defintion from buffer:
          call Init( ggi, ggN, reduced, status )
          IF_NOTOK_RETURN(status=1)

          ! allocate storage:
          allocate(  oro_gg(ggi%np) )
          allocate( slhf_gg(ggi%np) )
          allocate(   sp_gg(ggi%np) )

          ! store first half of slhf (accumulated over 3 hr):
          slhf_gg  = tmm%buf_gg(:,1) * dt_sec/2 ! W/m2 s
          slhf_tmi = tmm%buf_tmi

        ! second interval: [t1,t1+dt/2]
        t1s = t1
        t2s = t1 + IncrDate(sec=int(dt_sec/2))

        !call wrtgol( 'tmm read  : slhf ', t1s, '  -  ', t2s ); call goPr

        ! read first slhf in buffer (W/m2 time aver):
        call FillBuffer( tmm, archivekey, 'slhf', 'W/m2', tday, t1s, t2s, 'n', 'n', status )
        IF_NOTOK_RETURN(status=1)

          ! check ...
          call CheckBuffer( tmm, status, gridtype='gg', np=ggi%np )
          IF_NOTOK_RETURN(status=1)

          ! add second half of slhf (accumulated over 3 hr):
          slhf_gg = slhf_gg + tmm%buf_gg(:,1) * dt_sec/2   ! W/m2 s

          ! copy info
          call AddHistory( slhf_tmi, tmm%buf_tmi, status )

      !
      ! ~~ 3 hourly
      !

      case ( 0 )

        ! convective stuff for instant time;
        ! assume this is the 3 hourly meteo ...

        ! only for times around 00, 03, ...
        call Get( t1, hour=hour )
        if ( modulo(hour,3) /= 0 ) then
          write (gol,'("3 hourly convection only for intervals [00,03], [03,06], ...")'); call goErr
          TRACEBACK; status=1; return
        end if

        ! times for model level fields is begin of requested interval
        ! (just a choice)
        t1m = t1
        t2m = t1m

        ! use 6 hour interval around requested (instant!) time to read slhf;
        ! for 3 hourly request, this means double reading
        ! (no problem since slhf is time average anyway)

        !
        ! slhf (gg), accumulated over [-3,3] hours
        !

        ! slhf will be valid for [-3,3]
        dt_sec = 6*3600.0

        ! forecast for 72 hour or later ? then slhf for [-6,6]
        if ( rTotal(t1-tday,'hour') >= 12.0+72.0 ) dt_sec = 12*3600.0

        ! forecast for fcday 10 is only for [00:00,12:00] where 12:00 is 12+240,
        ! thus no second interval; instead copy from first
        skip_second = rTotal(t1-tday,'hour') >= 12.0+240.0

           ! >>> adhoc: always skip second interval,
           ! because of problems with setting tday;
           ! for requested time, this could be previous day
           ! while slhf is taken from two different days
           ! (in fact routines should be called with two tdays rather than one)
           !skip_second = .true.
           !write (gol,'("WARNING - skipped second interval for slhf")'); call goPr

        ! first interval: [t1-dt,t1]
        t1s = t1 - IncrDate(sec=int(dt_sec/2))
        t2s = t1

        !call wrtgol( 'tmm read  : slhf ', t1s, '  -  ', t2s ); call goPr

        ! read first slhf in buffer (W/m2 time aver):
        call FillBuffer( tmm, archivekey, 'slhf', 'W/m2', tday, t1s, t2s, 'n', 'n', status )
        IF_NOTOK_RETURN(status=1)

          ! check ...
          call CheckBuffer( tmm, status, gridtype='gg' )
          IF_NOTOK_RETURN(status=1)

          ! extract grid size
          ggN     = tmm%buf_ggi%N
          reduced = tmm%buf_ggi%reduced

          ! setup gg defintion from buffer:
          call Init( ggi, ggN, reduced, status )
          IF_NOTOK_RETURN(status=1)

          ! allocate storage:
          allocate(  oro_gg(ggi%np) )
          allocate( slhf_gg(ggi%np) )
          allocate(   sp_gg(ggi%np) )

          ! store first half of slhf (accumulated over 3 hr):
          slhf_gg  = tmm%buf_gg(:,1) * dt_sec/2 ! W/m2 s
          slhf_tmi = tmm%buf_tmi

        ! second interval: [t1,t1+dt]
        t1s = t1
        t2s = t1 + IncrDate(sec=int(dt_sec/2))

        ! skip seond interval ?
        if ( skip_second ) then

          call wrtgol( 'tmm copy  : slhf ', t1s, '  -  ', t2s ); call goPr

          ! use for the second field the first:
          slhf_gg = slhf_gg * 2

        else

          !call wrtgol( 'tmm read  : slhf ', t1s, '  -  ', t2s ); call goPr

          ! read second slhf in buffer (W/m2 time aver):
          call FillBuffer( tmm, archivekey, 'slhf', 'W/m2', tday, t1s, t2s, 'n', 'n', status )
          IF_NOTOK_RETURN(status=1)

            ! check ...
            call CheckBuffer( tmm, status, gridtype='gg', np=ggi%np )
            IF_NOTOK_RETURN(status=1)

            ! add second half of slhf (accumulated over 3 hr):
            slhf_gg = slhf_gg + tmm%buf_gg(:,1) * dt_sec/2   ! W/m2 s

            ! copy info
            call AddHistory( slhf_tmi, tmm%buf_tmi, status )

      end if

      !
      ! ~~ error
      !

      case default

        write (gol,'("unsupported lenght of time interval : ",i4)') dth; call goErr
        TRACEBACK; status=1; return

    end select


    !
    ! orography (gg)
    !

    !write (gol,'("tmm read  : oro")'); call goPr

    ! read orography in buffer:
    call FillBuffer( tmm, archivekey, 'oro', 'm m/s2', tday, t1, t1, 'n', 'n', status )
    IF_NOTOK_RETURN(status=1)

      ! interpol from sh ?
      select case ( tmm%buf_gridtype )
        case ( 'gg' )
          ! copy from buffer:
          oro_gg  = tmm%buf_gg(:,1)  ! m*[g]
        case ( 'sh' )
          ! interpol from sh to gg :
          call Interpol( tmm%buf_shi, tmm%buf_sh(:,1), ggi, oro_gg, status )
          IF_NOTOK_RETURN(status=1)
        case default
          write (gol,'("unsupported grid type `",a,"` for raw oro")') tmm%buf_gridtype
          TRACEBACK; status=1; return
      end select

      ! store:
      oro_tmi = tmm%buf_tmi

    !
    ! read raw 3d fields
    !

    ! ~~~

    !call wrtgol( 'tmm read  : u,v  ', t1m, '  -  ', t2m ); call goPr

    ! read VO
    call FillBuffer( tmm, archivekey, 'VO', '1/s', tday, t1m, t2m, 'n', 'n', status )
    IF_NOTOK_RETURN(status=1)

      ! check ...
      call CheckBuffer( tmm, status, gridtype='sh' )
      IF_NOTOK_RETURN(status=1)

      ! extract other grid sizes
      call Init( shi, tmm%buf_shi, status )
      IF_NOTOK_RETURN(status=1)
      nlev = tmm%buf_levi%nlev

      ! allocate 3d storage:
      allocate(    u_gg(ggi%np,0:nlev) )
      allocate(    v_gg(ggi%np,0:nlev) )
      allocate(    w_gg(ggi%np,0:nlev) )
      allocate(    t_gg(ggi%np,0:nlev) )
      allocate(    q_gg(ggi%np,0:nlev) )
      allocate(    p_gg(ggi%np,0:nlev) )
      allocate(    z_gg(ggi%np,0:nlev) )

      ! copy 3d spectral field from buffer:
      allocate( VO_sh(shi%np,nlev) )
      VO_sh = tmm%buf_sh

      ! copy info
      vo_tmi = tmm%buf_tmi

    ! read D
    call FillBuffer( tmm, archivekey, 'D', '1/s', tday, t1m, t2m, 'n', 'n', status )
    IF_NOTOK_RETURN(status=1)

      ! check ...
      call CheckBuffer( tmm, status, gridtype='sh', shT=shi%T, nlev=nlev )
      IF_NOTOK_RETURN(status=1)

      ! copy 3d spectral field from buffer:
      allocate( D_sh(shi%np,nlev) )
      D_sh = tmm%buf_sh

      ! copy info
      div_tmi = tmm%buf_tmi

    ! convert from sh VO/D to gg U/V
    ! setup storage:
    allocate( U_sh(shi%np) )
    allocate( V_sh(shi%np) )

    ! loop over levels:
    !xOMP PARALLEL &
    !xOMP     default ( none ) &
    !xOMP     shared  ( nlev, shi, VO_sh, D_sh ) &
    !xOMP     shared  ( ggi, u_gg, v_gg ) &
    !xOMP     private ( l, U_sh, V_sh, status )
    !xOMP   DO
    do l = 1, nlev
      ! convert to U=u*cos(lat) and V=v*cos(lat) :
      call vod2uv( shi, VO_sh(:,l), D_sh(:,l), shi, U_sh, V_sh )
      ! convert to Gaussian grid:
      call NewInterpol( shi, U_sh, ggi, u_gg(:,l), status )
      !IF_NOTOK_RETURN(status=1)
      call NewInterpol( shi, V_sh, ggi, v_gg(:,l), status )
      !IF_NOTOK_RETURN(status=1)
    end do
    !xOMP   END DO
    !xOMP END PARALLEL

    ! clear
    call Done( shi )
    deallocate( U_sh )
    deallocate( V_sh )

    ! history ...
    call Init( u_tmi, 'u', 'm/s', status, (/vo_tmi,div_tmi/) )
    call Init( v_tmi, 'u', 'm/s', status, (/vo_tmi,div_tmi/) )
    call AddHistory( u_tmi, 'VoD to UV;;interpol to gg', status )
    call AddHistory( v_tmi, 'VoD to UV;;interpol to gg', status )

    ! clear
    deallocate( VO_sh )
    deallocate(  D_sh )

    ! remove cos(lat) factor:
    do j = 1, ggi%nlat
      i1 = ggi%i1(j)
      i2 = ggi%im(j)
      u_gg(i1:i2,1:nlev) = u_gg(i1:i2,1:nlev) / cos(ggi%lat(j))
      v_gg(i1:i2,1:nlev) = v_gg(i1:i2,1:nlev) / cos(ggi%lat(j))
    end do
    call AddHistory( u_tmi, 'divide by cos(lat)', status )
    call AddHistory( v_tmi, 'divide by cos(lat)', status )

    ! ~~~

    select case ( sourcetype )

      case ( 'ecmwf-tmpp', 'msc-data' )

        !call wrtgol( 'tmm read  : Q    ', t1m, '  -  ', t2m ); call goPr

        ! read humidity
        call FillBuffer( tmm, archivekey, 'Q', 'kg/kg', tday, t1m, t2m, 'n', 'n', status )
        IF_NOTOK_RETURN(status=1)

          ! check ...
          call CheckBuffer( tmm, status, gridtype='gg', np=ggi%np, nlev=nlev )
          IF_NOTOK_RETURN(status=1)

          ! copy 3d gg field from buffer; copy info:
          Q_gg(:,1:nlev) = tmm%buf_gg     ! kg h2o / kg air

          ! info ...
          call Init( Q_tmi, tmm%buf_tmi, status )

          ! copy surface pressure:
          sp_gg  = tmm%buf_sp_gg      ! Pa
          call Init( sp_tmi, 'sp', 'Pa', status, (/tmm%buf_tmi/) )

        !call wrtgol( 'tmm read  : T   ', t1m, '  -  ', t2m ); call goPr

        ! read T
        call FillBuffer( tmm, archivekey, 'T', 'K', tday, t1m, t2m, 'n', 'n', status )
        IF_NOTOK_RETURN(status=1)

          ! check ...
          call CheckBuffer( tmm, status, gridtype='sh', shT=shi%T, nlev=nlev )
          IF_NOTOK_RETURN(status=1)

          do l = 1, nlev
            call Interpol( tmm%buf_shi, tmm%buf_sh(:,l), ggi, T_gg(:,l), status )
            IF_NOTOK_RETURN(status=1)
          end do

          ! info ...
          call Init( T_tmi, tmm%buf_tmi, status )
          call AddHistory( T_tmi, 'interpol to gg', status )


      case ( 'ecmwf-tm5' )

        !call wrtgol( 'tmm read  : Q    ', t1m, '  -  ', t2m ); call goPr

        ! read humidity
        call FillBuffer( tmm, archivekey, 'Q', 'kg/kg', tday, t1m, t2m, 'n', 'n', status )
        IF_NOTOK_RETURN(status=1)

          ! check ...
          call CheckBuffer( tmm, status, gridtype='gg', np=ggi%np, nlev=nlev )
          IF_NOTOK_RETURN(status=1)

          ! copy 3d gg field from buffer::
          Q_gg(:,1:nlev) = tmm%buf_gg     ! kg h2o / kg air

          ! info ...
          call Init( Q_tmi, tmm%buf_tmi, status )

          ! copy surface pressure:
          sp_gg  = tmm%buf_sp_gg      ! Pa
          call Init( sp_tmi, 'sp', 'Pa', status, (/tmm%buf_tmi/) )

        !call wrtgol( 'tmm read  : T   ', t1m, '  -  ', t2m ); call goPr

        ! read T
        call FillBuffer( tmm, archivekey, 'T', 'K', tday, t1m, t2m, 'n', 'n', status )
        IF_NOTOK_RETURN(status=1)

          ! check ...
          call CheckBuffer( tmm, status, gridtype='gg', np=ggi%np, nlev=nlev )
          IF_NOTOK_RETURN(status=1)

          ! copy 3d gg field from buffer::
          T_gg(:,1:nlev) = tmm%buf_gg     ! K

          ! info ...
          call Init( T_tmi, tmm%buf_tmi, status )


      case ( 'ncep-cdc', 'ncep-gfs' )

        !call wrtgol( 'tmm read  : Q    ', t1m, '  -  ', t2m ); call goPr

        ! read humidity
        call FillBuffer( tmm, archivekey, 'Q', 'kg/kg', tday, t1m, t2m, 'n', 'n', status )
        IF_NOTOK_RETURN(status=1)

          ! check ...
          call CheckBuffer( tmm, status, gridtype='sh', np=ggi%np, nlev=nlev )
          IF_NOTOK_RETURN(status=1)

          ! check ...
          call CheckBuffer( tmm, status, gridtype='sh', shT=shi%T, nlev=nlev )
          IF_NOTOK_RETURN(status=1)

          ! convert from sh to gg
          do l = 1, nlev
            call Interpol( tmm%buf_shi, tmm%buf_sh(:,l), ggi, Q_gg(:,l), status )   ! kg h2o / kg air
            IF_NOTOK_RETURN(status=1)
          end do

          ! prevent negatives ...
          Q_gg = max( 0.0, Q_gg )

          ! info ...
          call Init( Q_tmi, tmm%buf_tmi, status )
          call AddHistory( Q_tmi, 'interpol to gg', status )
          call AddHistory( Q_tmi, 'truncate', status )

          ! interpolate surface pressure:
          call Interpol( tmm%buf_shi, tmm%buf_lnsp_sh, ggi, sp_gg, status )   ! ln(Pa)
          IF_NOTOK_RETURN(status=1)
          sp_gg = exp( sp_gg )        ! Pa
          ! info ...
          call Init( sp_tmi, 'sp', 'Pa', status, (/tmm%buf_tmi/) )

        !call wrtgol( 'tmm read  : T   ', t1m, '  -  ', t2m ); call goPr

        ! read virtual temperature
        call FillBuffer( tmm, archivekey, 'Tv', 'K', tday, t1m, t2m, 'n', 'n', status )
        IF_NOTOK_RETURN(status=1)

          ! check ...
          call CheckBuffer( tmm, status, gridtype='sh', shT=shi%T, nlev=nlev )
          IF_NOTOK_RETURN(status=1)

          ! convert from sh to gg
          do l = 1, nlev
            call Interpol( tmm%buf_shi, tmm%buf_sh(:,l), ggi, T_gg(:,l), status )
            IF_NOTOK_RETURN(status=1)
          end do

          ! convert from virtual to normal temperature:
          T_gg = RealTemperature( T_gg, Q_gg )

          ! info ...
          call Init( T_tmi, tmm%buf_tmi, status )
          call AddHistory( T_tmi, 'interpol to gg', status )


      case default

        write (gol,'("unsupported source type `",a,"` for raw surface fields")') trim(sourcetype); call goErr
        TRACEBACK; status=1; return

    end select


    ! ~~~

    write (gol,'("    interpol : W")'); call goPr

    ! omega on gg and TM levels:
    allocate( tmp_gg(ggi%np,levi%nlev+1) )

    ! convert from ll to gg:
    do l = 1, levi%nlev+1
      call Interpol( lli, omega(:,:,l), ggi, tmp_gg(:,l), status )
      IF_NOTOK_RETURN(status=1)
    end do

    ! convert to parent levels (buffer) from TM levels:
    call FillLevels( tmm%buf_levi, 'w', sp_gg, w_gg, &
                     levi, tmp_gg, 'bottom', status )
    IF_NOTOK_RETURN(status=1)

    ! info ...
    call Init( w_tmi, omega_tmi, status )
    call AddHistory( w_tmi, 'interpol to gg and raw levels', status )

    ! clear
    deallocate( tmp_gg )


    !
    ! ensure that layers are in ecmwf order (top->down)
    !

    select case ( sourcetype )

      case ( 'ecmwf-tmpp', 'ecmwf-tm5', 'msc-data' )

        ! copy level info from buffer:
        call Init( leviX, tmm%buf_levi%nlev, tmm%buf_levi%a, tmm%buf_levi%b, status )
        IF_NOTOK_RETURN(status=1)

      case ( 'ncep-cdc', 'ncep-gfs' )

        ! revert level info from buffer:
        allocate( aX(0:tmm%buf_levi%nlev) )
        allocate( bX(0:tmm%buf_levi%nlev) )
        do l = 0, tmm%buf_levi%nlev
          aX(l) = tmm%buf_levi%a(tmm%buf_levi%nlev-l)
          bX(l) = tmm%buf_levi%b(tmm%buf_levi%nlev-l)
        end do
        call Init( leviX, tmm%buf_levi%nlev, aX, bX, status )
        IF_NOTOK_RETURN(status=1)
        deallocate( aX )
        deallocate( bX )

        ! revert fields
        allocate( tmp_gg(ggi%np,0:nlev) )
        tmp_gg = u_gg
        do l = 1, nlev
          u_gg(:,l) = tmp_gg(:,nlev+1-l)
        end do
        tmp_gg = v_gg
        do l = 1, nlev
          v_gg(:,l) = tmp_gg(:,nlev+1-l)
        end do
        tmp_gg = w_gg
        do l = 1, nlev
          w_gg(:,l) = tmp_gg(:,nlev+1-l)
        end do
        tmp_gg = t_gg
        do l = 1, nlev
          t_gg(:,l) = tmp_gg(:,nlev+1-l)
        end do
        tmp_gg = q_gg
        do l = 1, nlev
          q_gg(:,l) = tmp_gg(:,nlev+1-l)
        end do
        deallocate( tmp_gg )

    end select


    !
    ! updraughts/downdraughts
    !

    write (gol,'("    tmm updr/downdr")'); call goPr

    !
    ! WARNING: the TMPP2/tmpp_conv-gg on bsgi59 is probably not bugfree:
    !  o oro not read
    !  o pw = w/g [kg/m2/s], but tmpp_conv_tiedtke expects Pa/s ?
    !  o factor to convert from lshf to evaporation
    !      here  :  1 m/s = 2.45e9 W/m2 (at 20 degrees C)
    !      binas :  lvap = 2.5 e6 J/kg
    !      binas :  Lc   = 2.26e6 J/kg  (at 0 deg C)
    !
    ! There are however good reasons to step over to the conv-gg code
    ! by Dirk Olivie
    !  o no unnecessary interpolations to half levels
    !  o removed stuff that was based on the ec 19 levels
    !

    ! extra fields:
    allocate( ggmask(ggi%np) )
    allocate(  qam_gg(ggi%np,0:nlev) )
    allocate(  qac_gg(ggi%np,0:nlev) )
    allocate( entu_gg(ggi%np,nlev) )
    allocate( detu_gg(ggi%np,nlev) )
    allocate( entd_gg(ggi%np,nlev) )
    allocate( detd_gg(ggi%np,nlev) )

    ! Set mask for averaging over ll;
    ! one extra row to compute derivatives.
    ! This routine changes ggi: flag for each row to be processed or not
    call InterpolMask( ggi, ggmask, lli, 1 )

    ! calculate geopotential (z)
    allocate( p_hlev(0:nlev) )
    do i = 1, ggi%np
      ! skip if not required for ll grid:
      if ( .not. ggmask(i) ) cycle
      ! compute pressure at half levels (surf -> top):
      call phlev_ec_rev( nlev, leviX%a, leviX%b, sp_gg(i), p_hlev )
      ! compute z for single column:
      call GeoPot( nlev, oro_gg(i), T_gg(i,1:nlev), Q_gg(i,1:nlev), &
                                                p_hlev, z_gg(i,1:nlev) )
    end do
    deallocate( p_hlev )

    ! interpolate variables from the center of parent-model layers to the
    ! boundaries of parent-model layers and save result in same memory location ..
    call mid2bound_uv( nlev, ggi%np, u_gg, v_gg, sp_gg, ggmask, leviX%a, leviX%b )
    call mid2bound_t ( nlev, ggi%np, t_gg, sp_gg, ggmask, leviX%a, leviX%b )
    call mid2bound_q ( nlev, ggi%np, q_gg, sp_gg, ggmask, leviX%a, leviX%b, t_gg )
    call mid2bound_z ( nlev, ggi%np, z_gg, sp_gg, ggmask, leviX%a, leviX%b, oro_gg )

    ! already on half levels, since interpolated from omega ...
    !!call mid2bound_w ( nlev, ggi%np, w_gg, sp_gg, ggmask, leviX%a, leviX%b )

    ! NOTE: p is not filled on input, but filled with pressures on output
    call mid2bound_p ( nlev, ggi%np, p_gg, sp_gg, ggmask, leviX%a, leviX%b )

    ! divergence fields
    do l = 0, nlev
      call Divergence( ggi, q_gg(:,l)*u_gg(:,l), q_gg(:,l)*v_gg(:,l), qac_gg(:,l) )
      call Divergence( ggi,           u_gg(:,l),           v_gg(:,l), qam_gg(:,l) )
    end do

    ! Convert from SLHF (W/m2*interval) to EVAP (m/s)
    ! 1 m/s = 2.45e9 W/m2 (at 20 degrees C).
    ! Don't forget to change sign from latent heatflux to evaporation!!!
    !    evap = - slhf_gg(i) / dt_sec / 2.45e9    ! m/s
    ! (apply in argument)

    ! work routine:
    call subscal_2d( ggi%np, nlev, leviX%a, leviX%b, &
                    z_gg, p_gg, w_gg, t_gg, &
                    q_gg, qac_gg, qam_gg, -1.0e3*slhf_gg/dt_sec/2.45e9, &
                    entu_gg, detu_gg, entd_gg, detd_gg )


    ! clear
    deallocate( ggmask )
    deallocate( slhf_gg )
    deallocate(  oro_gg )
    deallocate(    u_gg )
    deallocate(    v_gg )
    deallocate(    w_gg )
    deallocate(    t_gg )
    deallocate(    q_gg )
    deallocate(  qam_gg )
    deallocate(  qac_gg )
    deallocate(    p_gg )
    deallocate(    z_gg )

    ! history
    call Init( entu_tmi, 'entu', 'kg/m2/s', status, &
                 (/sp_tmi,slhf_tmi,oro_tmi,T_tmi,Q_tmi,u_tmi,v_tmi,w_tmi/) )
    call Init( entd_tmi, 'entd', 'kg/m2/s', status, &
                 (/sp_tmi,slhf_tmi,oro_tmi,T_tmi,Q_tmi,u_tmi,v_tmi,w_tmi/) )
    call Init( detu_tmi, 'detu', 'kg/m2/s', status, &
                 (/sp_tmi,slhf_tmi,oro_tmi,T_tmi,Q_tmi,u_tmi,v_tmi,w_tmi/) )
    call Init( detd_tmi, 'detd', 'kg/m2/s', status, &
                 (/sp_tmi,slhf_tmi,oro_tmi,T_tmi,Q_tmi,u_tmi,v_tmi,w_tmi/) )

    !
    ! convert from gg to ll
    !

    ! determine fraction
    call Init( gg2ll, ggi, lli, status )
    IF_NOTOK_RETURN(status=1)

    ! take fractions of overlapping cells
    call FracSum( gg2ll, sp_gg, sp, status, 'area-aver' )
    IF_NOTOK_RETURN(status=1)
    ! clear
    deallocate( sp_gg )

    ! full level output:
    allocate( llX(lli%nlon,lli%nlat,nlev     ) )
    allocate( ll (lli%nlon,lli%nlat,levi%nlev) )

    ! take fractions of overlapping cells
    do l = 1, nlev
      call FracSum( gg2ll, entu_gg(:,l), llX(:,:,l), status, 'area-aver' )
      IF_NOTOK_RETURN(status=1)
    end do
    ! integrated variables  might become slightly negative ...
    llX = max( 0.0, llX )
    ! combine levels from llX to ll:
    call FillLevels( levi, 'n', sp, ll, leviX, llX, 'sum', status )
    IF_NOTOK_RETURN(status=1)
    ! truncate to number of output levels:
    entu = ll(:,:,1:lout)
    ! info ..
    call AddHistory( entu_tmi, 'gg to ll, area aver;;sum levels', status )

    ! take fractions of overlapping cells
    do l = 1, nlev
      call FracSum( gg2ll, entd_gg(:,l), llX(:,:,l), status, 'area-aver' )
      IF_NOTOK_RETURN(status=1)
    end do
    ! integrated variables  might become slightly negative ...
    llX = max( 0.0, llX )
    ! combine levels from llX to ll:
    call FillLevels( levi, 'n', sp, ll, leviX, llX, 'sum', status )
    IF_NOTOK_RETURN(status=1)
    ! truncate to number of output levels:
    entd = ll(:,:,1:lout)
    ! info ..
    call AddHistory( entd_tmi, 'gg to ll, area aver;;sum levels', status )

    ! take fractions of overlapping cells
    do l = 1, nlev
      call FracSum( gg2ll, detu_gg(:,l), llX(:,:,l), status, 'area-aver' )
      IF_NOTOK_RETURN(status=1)
    end do
    ! integrated variables  might become slightly negative ...
    llX = max( 0.0, llX )
    ! combine levels from llX to ll:
    call FillLevels( levi, 'n', sp, ll, leviX, llX, 'sum', status )
    IF_NOTOK_RETURN(status=1)
    ! truncate to number of output levels:
    detu = ll(:,:,1:lout)
    ! info ..
    call AddHistory( detu_tmi, 'gg to ll, area aver;;sum levels', status )

    ! take fractions of overlapping cells
    do l = 1, nlev
      call FracSum( gg2ll, detd_gg(:,l), llX(:,:,l), status, 'area-aver' )
      IF_NOTOK_RETURN(status=1)
    end do
    ! integrated variables  might become slightly negative ...
    llX = max( 0.0, llX )
    ! combine levels from llX to ll:
    call FillLevels( levi, 'n', sp, ll, leviX, llX, 'sum', status )
    IF_NOTOK_RETURN(status=1)
    ! truncate to number of output levels:
    detd = ll(:,:,1:lout)
    ! info ..
    call AddHistory( detd_tmi, 'gg to ll, area aver;;sum levels', status )

    ! clear
    call Done( gg2ll )
    deallocate( entu_gg )
    deallocate( entd_gg )
    deallocate( detu_gg )
    deallocate( detd_gg )
    deallocate( llX )
    deallocate( ll  )
    call Done( leviX, status )
    IF_NOTOK_RETURN(status=1)


    !
    ! done
    !

    call Done( ggi, status )
    IF_NOTOK_RETURN(status=1)


    ! ok
    status = 0

  end subroutine tmm_Read_Convec_raw

#endif


  ! *****************************************************************


#ifdef with_tmm_convec_ec_gg

  subroutine tmm_Read_Convec_EC_gg( tmm, archivekey, tday, t1, t2, lli, levi, &
                                omega, omega_tmi, &
                                sp, &
                                entu, entu_tmi, entd, entd_tmi, &
                                detu, detu_tmi, detd, detd_tmi, &
                                status )

    use Binas   , only : grav
    use GO      , only : goErr, gol, goPr
    use GO      , only : TDate, wrtgol, Get, IncrDate
    use GO      , only : operator(-), operator(+), rTotal
    use GO      , only : goSplitLine
    use Phys    , only : mid2bound_uv, mid2bound_w, mid2bound_t
    use Phys    , only : mid2bound_q, mid2bound_z, mid2bound_p
    use Phys    , only : subscal, phlev_ec_rev, geopot
    use Phys    , only : RealTemperature
    use Phys    , only : GeoPotentialHeightB
    use Phys    , only : ECconv_to_TMconv
    use Grid    , only : TllGridInfo
    use Grid    , only : TLevelInfo, HPressure, FPressure, FillLevels
    use Grid    , only : TShGrid, VoD2UV
    use Grid    , only : TggGridInfo, InterpolMask, Divergence, Tgg2llFracs, FracSum, AreaOper
    use Grid    , only : Init, Done, Set, Interpol
    use tmm_info, only : TMeteoInfo, Init, AddHistory

    ! --- in/out --------------------------------

    type(TTmMeteo), intent(inout)     ::  tmm
    character(len=*), intent(in)      ::  archivekey
    type(TDate), intent(in)           ::  tday, t1, t2
    type(TllGridInfo), intent(in)     ::  lli
    type(TLevelInfo), intent(in)      ::  levi
    real, intent(in)                  ::  omega(:,:,:)   ! Pa/s
    type(TMeteoInfo), intent(in)      ::  omega_tmi
    real, intent(out)                 ::  sp(:,:)      ! Pa
    real, intent(out)                 ::  entu(:,:,:), entd(:,:,:)   ! kg/m2/s
    type(TMeteoInfo), intent(out)     ::  entu_tmi, entd_tmi
    real, intent(out)                 ::  detu(:,:,:), detd(:,:,:)   ! kg/m2/s
    type(TMeteoInfo), intent(out)     ::  detu_tmi, detd_tmi
    integer, intent(out)              ::  status

    ! --- const ------------------------------------

    character(len=*), parameter ::  rname = mname//'/tmm_Read_Convec_EC_gg'

    ! --- local -------------------------------

    character(len=10)        ::  sourcetype
    character(len=MAX_RCKEY_LEN)       ::  sourcename

    integer                  ::  lout
    real, allocatable        ::  ll(:,:,:)

    ! timing
    integer                  ::  hour
    real                     ::  dhour
    type(TDate)              ::  t1s, t2s

    ! loops etc
    integer                  ::  l, nlev
    integer                  ::  i, i1, i2, j

    ! gaussian grids
    integer                  ::  ggN
    logical                  ::  reduced
    type(TggGridInfo)        ::  ggi
    real, allocatable        ::  gg(:)
    real, allocatable        ::   oro_gg(:)
    real, allocatable        ::    sp_gg(:)
    real, allocatable        ::     m_gg(:)
    real, allocatable        ::  t_gg(:,:)
    real, allocatable        ::  q_gg(:,:)
    type(Tgg2llFracs)        ::  gg2ll
    real, allocatable        ::  llX(:,:,:)

    real, allocatable        ::  p_hlev(:)

    ! e4 convection fields
    real, allocatable        ::  udmf_gg(:,:)
    real, allocatable        ::  ddmf_gg(:,:)
    real, allocatable        ::  uddr_gg(:,:)
    real, allocatable        ::  dddr_gg(:,:)
    type(TMeteoInfo)         ::  udmf_tmi, ddmf_tmi, uddr_tmi, dddr_tmi
    real, allocatable        ::  ph_ec(:)
    real, allocatable        ::  zh_ec(:)

    ! --- begin -------------------------------

    ! number of levels in output arrays:
    lout = size(entu,3)

    ! check ...
    if ( (size(entd,3)/=lout) .or. &
         (size(detu,3)/=lout) .or. (size(detd,3)/=lout) ) then
      write (gol,'("output arrays should have same number of levels:")'); call goErr
      write (gol,'("  entu : ",i4)') size(entu,3); call goErr
      write (gol,'("  entd : ",i4)') size(entd,3); call goErr
      write (gol,'("  detu : ",i4)') size(detu,3); call goErr
      write (gol,'("  detd : ",i4)') size(detd,3); call goErr
      TRACEBACK; status=1; return
    end if

    ! split source key in type and name:
    call goSplitLine( archivekey, sourcetype, ':', sourcename, status )
    IF_NOTOK_RETURN(status=1)

    !
    ! Original parameters archived in MARS:
    !
    !     Parameter                      Surfaces           Code    Units   Units
    !                                                         1)      2)      3)
    !
    !     updraught   mass flux          half levels  UDMF  101     kg/m2   kg/m2/s
    !
    !     downdraught mass flux          half levels  DDMF  102     kg/m2   kg/m2/s
    !
    !     updraught   detrainment rate   full levels  UDDR  103     kg/m3   kg/m3/s
    !
    !     downdraught detrainment rate   full levels  DDDR  104     kg/m3   kg/m3/s
    !
    !      1) GRIB code table 2 version 128
    !      2) original units, accumulated
    !      3) time averaged after reading
    !

    ! only hours 00/03/06/etc yet ...
    call Get( t1, hour=hour )
    dhour = rTotal( t2 - t1, 'hour' )
    if ( (modulo(hour,3) /= 0) .or. (dhour /= 3.0) ) then
      write (gol,'("convection only for 3hr intervals [0,3] etc.")'); call goErr
      call wrtgol( '  requested : ', t1, '  -  ', t2 ); call goErr
      write (gol,'("  dhour     : ",f8.4)') dhour; call goErr
      TRACEBACK; status=1; return
    end if

    !
    ! read gg fields
    !

    !call wrtgol( 'tmm read  : UDMF ', t1, '  -  ', t2 ); call goPr

    ! read updraught mass flux
    call FillBuffer( tmm, archivekey, 'UDMF', tday, t1, t2, 'n', 'n', status )
    IF_NOTOK_RETURN(status=1)

      ! check ...
      call CheckBuffer( tmm, status, gridtype='gg' )
      IF_NOTOK_RETURN(status=1)

      ! extract grid sizes
      ggN     = tmm%buf_ggi%N
      reduced = tmm%buf_ggi%reduced

      ! copy level info from buffer:
      call Init( leviX, tmm%buf_levi%nlev, tmm%buf_levi%a, tmm%buf_levi%b, status )
      IF_NOTOK_RETURN(status=1)

      ! setup gg defintion from buffer:
      call Init( ggi, ggN, reduced, status )
      IF_NOTOK_RETURN(status=1)

      ! allocate storage:
      allocate( udmf_gg(ggi%np,0:leviX%nlev) )
      allocate( ddmf_gg(ggi%np,0:leviX%nlev) )
      allocate( uddr_gg(ggi%np,1:leviX%nlev) )
      allocate( dddr_gg(ggi%np,1:leviX%nlev) )
      allocate(   sp_gg(ggi%np) )
      allocate(  oro_gg(ggi%np) )
      allocate(    T_gg(ggi%np,1:leviX%nlev) )
      allocate(    Q_gg(ggi%np,1:leviX%nlev) )

      ! copy 3d field, levels top down, surface implicit zero, copy info:
      udmf_gg(:,0:nlev-1) = tmm%buf_gg   !  kg/m2/s
      udmf_gg(:,    nlev) = 0.0
      call AddHistory( udmf_tmi, tmm%buf_tmi, status )

    !call wrtgol( 'tmm read  : DDMF ', t1, '  -  ', t2 ); call goPr

    ! downdraught mass flux
    call FillBuffer( tmm, archivekey, 'DDMF', tday, t1, t2, 'n', 'n', status )
    IF_NOTOK_RETURN(status=1)

      ! check ...
      call CheckBuffer( tmm, status, gridtype='gg', np=ggi%np, nlev=leviX%nlev )
      IF_NOTOK_RETURN(status=1)

      ! copy 3d field, levels top down, surface implicit zero, copy info:
      ddmf_gg(:,0:nlev-1) = tmm%buf_gg   !  kg/m2/s
      ddmf_gg(:,    nlev) = 0.0
      call AddHistory( ddmf_tmi, tmm%buf_tmi, status )

    !call wrtgol( 'tmm read  : UDDR ', t1, '  -  ', t2 ); call goPr

    ! updraught detrainment rate
    call FillBuffer( tmm, archivekey, 'UDDR', tday, t1, t2, 'n', 'n', status )
    IF_NOTOK_RETURN(status=1)

      ! check ...
      call CheckBuffer( tmm, status, gridtype='gg', np=ggi%np, nlev=leviX%nlev )
      IF_NOTOK_RETURN(status=1)

      ! copy 3d field, levels top down, copy info:
      uddr_gg = tmm%buf_gg   !  kg/m3/s
      call AddHistory( uddr_tmi, tmm%buf_tmi, status )

    !call wrtgol( 'tmm read  : DDDR ', t1, '  -  ', t2 ); call goPr

    ! downdraught detrainment rate
    call FillBuffer( tmm, archivekey, 'DDDR', tday, t1, t2, 'n', 'n', status )
    IF_NOTOK_RETURN(status=1)

      ! check ...
      call CheckBuffer( tmm, status, gridtype='gg', np=ggi%np, nlev=leviX%nlev )
      IF_NOTOK_RETURN(status=1)

      ! copy 3d field, levels top down, copy info:
      dddr_gg = tmm%buf_gg   !  kg/m3/s
      call AddHistory( dddr_tmi, tmm%buf_tmi, status )

    ! temperature at begin of interval
    !call wrtgol( 'tmm read  : T    ', t1, '  -  ', t1 ); call goPr
    call FillBuffer( tmm, archivekey, 'T', tday, t1, t1, 'n', 'n', status )
    IF_NOTOK_RETURN(status=1)

      ! check ...
      call CheckBuffer( tmm, status, gridtype='gg', np=ggi%np, nlev=leviX%nlev )
      IF_NOTOK_RETURN(status=1)

      ! copy 3d field:
      T_gg = tmm%buf_gg   ! K

    ! temperature at end of interval
    !call wrtgol( 'tmm read  : T    ', t2, '  -  ', t2 ); call goPr
    call FillBuffer( tmm, archivekey, 'T', tday, t2, t2, 'n', 'n', status )
    IF_NOTOK_RETURN(status=1)

      ! check ...
      call CheckBuffer( tmm, status, gridtype='gg', np=ggi%np, nlev=leviX%nlev )
      IF_NOTOK_RETURN(status=1)

      ! add, average:
      T_gg = ( T_gg + tmm%buf_gg )/2.0   ! K

    ! humidity at begin of interval
    !call wrtgol( 'tmm read  : Q    ', t1, '  -  ', t1 ); call goPr
    call FillBuffer( tmm, archivekey, 'Q', tday, t1, t1, 'n', 'n', status )
    IF_NOTOK_RETURN(status=1)

      ! check ...
      call CheckBuffer( tmm, status, gridtype='gg', np=ggi%np, nlev=leviX%nlev )
      IF_NOTOK_RETURN(status=1)

      ! copy 3d field:
      Q_gg = tmm%buf_gg   ! kg/kg

    ! humidity at end of interval
    !call wrtgol( 'tmm read  : Q    ', t2, '  -  ', t2 ); call goPr
    call FillBuffer( tmm, archivekey, 'Q', tday, t2, t2, 'n', 'n', status )
    IF_NOTOK_RETURN(status=1)

      ! check ...
      call CheckBuffer( tmm, status, gridtype='gg', np=ggi%np, nlev=leviX%nlev )
      IF_NOTOK_RETURN(status=1)

      ! add, average:
      Q_gg = ( Q_gg + tmm%buf_gg )/2.0   ! kg/kg

    ! surface pressure at begin of interval
    !call wrtgol( 'tmm read  : sp   ', t1, '  -  ', t1 ); call goPr
    call FillBuffer( tmm, archivekey, 'sp', tday, t1, t1, 'n', 'n', status )
    IF_NOTOK_RETURN(status=1)

      ! check ...
      call CheckBuffer( tmm, status, gridtype='gg', np=ggi%np )
      IF_NOTOK_RETURN(status=1)

      ! copy 2d field:
      sp_gg = tmm%buf_gg(:,1)   ! Pa

    ! surface pressure at end of interval
    !call wrtgol( 'tmm read  : sp   ', t2, '  -  ', t2 ); call goPr
    call FillBuffer( tmm, archivekey, 'sp', tday, t2, t2, 'n', 'n', status )
    IF_NOTOK_RETURN(status=1)

      ! check ...
      call CheckBuffer( tmm, status, gridtype='gg', np=ggi%np )
      IF_NOTOK_RETURN(status=1)

      ! add, average:
      sp_gg = ( sp_gg + tmm%buf_gg(:,1) )/2.0   ! Pa

    ! read orography in buffer:
    !write (gol,'("    tmm read  : oro")'); call goPr
    call FillBuffer( tmm, archivekey, 'oro', tday, t1, t1, 'n', 'n', status )
    IF_NOTOK_RETURN(status=1)

      ! check ...
      call CheckBuffer( tmm, status, gridtype='gg', np=ggi%np )
      IF_NOTOK_RETURN(status=1)

      ! copy from buffer:
      oro_gg  = tmm%buf_gg(:,1)  ! m*[g]

    !
    ! convert
    !

    allocate( ph_ec(0:leviX%nlev) )
    allocate( zh_ec(0:leviX%nlev) )
    allocate( entu_gg(ggi%np,leviX%nlev) )
    allocate( detu_gg(ggi%np,leviX%nlev) )
    allocate( entd_gg(ggi%np,leviX%nlev) )
    allocate( detd_gg(ggi%np,leviX%nlev) )

    ! loop over gg cells:
    do i = 1, ggi%np

      ! half level pressure:
      call HPressure( leviX, sp_gg(i), ph_ec, status )
      IF_NOTOK_RETURN(status=1)

      ! gph at half levels:
      call GeoPotentialHeightB( nlev, ph_ec, T_gg(i,:), Q_gg(i,:), oro_gg(i)/grav, zh_ec )

      ! convert from fluxes to rates:
      call ECconv_to_TMconv( leviX%nlev, zh_ec, &
                             udmf_gg(i,:), uddr_gg(i,:), ddmf_gg(i,:), dddr_gg(i,:), &
                             entu_gg(i,:), detu_gg(i,:), entd_gg(i,:), detd_gg(i,:), &
                             status )
      IF_NOTOK_RETURN(status=1)

    end do

    deallocate(   ph_ec )
    deallocate(   zh_ec )

    ! clear
    deallocate( udmf_gg )
    deallocate( ddmf_gg )
    deallocate( uddr_gg )
    deallocate( dddr_gg )
    deallocate(  oro_gg )
    deallocate(    T_gg )
    deallocate(    Q_gg )

    ! history
    call Init( entu_tmi, 'entu', 'kg/m2/s', status, (/udmf_tmi,ddmf_tmi,uddr_tmi,dddr_tmi/) )
    call Init( entd_tmi, 'entd', 'kg/m2/s', status, (/udmf_tmi,ddmf_tmi,uddr_tmi,dddr_tmi/) )
    call Init( detu_tmi, 'detu', 'kg/m2/s', status, (/udmf_tmi,ddmf_tmi,uddr_tmi,dddr_tmi/) )
    call Init( detd_tmi, 'detd', 'kg/m2/s', status, (/udmf_tmi,ddmf_tmi,uddr_tmi,dddr_tmi/) )

    !
    ! convert from gg to ll
    !

    ! determine fraction
    call Init( gg2ll, ggi, lli, status )
    IF_NOTOK_RETURN(status=1)

    ! take fractions of overlapping cells
    call FracSum( gg2ll, sp_gg, sp, status, 'area-aver' )
    IF_NOTOK_RETURN(status=1)

    ! full level output:
    allocate( llX(lli%nlon,lli%nlat,nlev     ) )
    allocate( ll (lli%nlon,lli%nlat,levi%nlev) )

    ! take fractions of overlapping cells
    do l = 1, nlev
      call FracSum( gg2ll, entu_gg(:,l), llX(:,:,l), status, 'area-aver' )
      IF_NOTOK_RETURN(status=1)
    end do
    ! integrated variables  might become slightly negative ...
    llX = max( 0.0, llX )
    ! combine levels from llX to ll:
    call FillLevels( levi, 'n', sp, ll, leviX, llX, 'sum', status )
    IF_NOTOK_RETURN(status=1)
    ! truncate to number of output levels:
    entu = ll(:,:,1:lout)
    ! info ..
    call AddHistory( entu_tmi, 'gg to ll, area aver;;sum levels', status )

    ! take fractions of overlapping cells
    do l = 1, nlev
      call FracSum( gg2ll, entd_gg(:,l), llX(:,:,l), status, 'area-aver' )
      IF_NOTOK_RETURN(status=1)
    end do
    ! integrated variables  might become slightly negative ...
    llX = max( 0.0, llX )
    ! combine levels from llX to ll:
    call FillLevels( levi, 'n', sp, ll, leviX, llX, 'sum', status )
    IF_NOTOK_RETURN(status=1)
    ! truncate to number of output levels:
    entd = ll(:,:,1:lout)
    ! info ..
    call AddHistory( entd_tmi, 'gg to ll, area aver;;sum levels', status )

    ! take fractions of overlapping cells
    do l = 1, nlev
      call FracSum( gg2ll, detu_gg(:,l), llX(:,:,l), status, 'area-aver' )
      IF_NOTOK_RETURN(status=1)
    end do
    ! integrated variables  might become slightly negative ...
    llX = max( 0.0, llX )
    ! combine levels from llX to ll:
    call FillLevels( levi, 'n', sp, ll, leviX, llX, 'sum', status )
    IF_NOTOK_RETURN(status=1)
    ! truncate to number of output levels:
    detu = ll(:,:,1:lout)
    ! info ..
    call AddHistory( detu_tmi, 'gg to ll, area aver;;sum levels', status )

    ! take fractions of overlapping cells
    do l = 1, nlev
      call FracSum( gg2ll, detd_gg(:,l), llX(:,:,l), status, 'area-aver' )
      IF_NOTOK_RETURN(status=1)
    end do
    ! integrated variables  might become slightly negative ...
    llX = max( 0.0, llX )
    ! combine levels from llX to ll:
    call FillLevels( levi, 'n', sp, ll, leviX, llX, 'sum', status )
    IF_NOTOK_RETURN(status=1)
    ! truncate to number of output levels:
    detd = ll(:,:,1:lout)
    ! info ..
    call AddHistory( detd_tmi, 'gg to ll, area aver;;sum levels', status )

    ! clear
    call Done( gg2ll )
    deallocate( entu_gg )
    deallocate( entd_gg )
    deallocate( detu_gg )
    deallocate( detd_gg )
    deallocate( llX )
    deallocate( ll  )

    ! clear
    deallocate( sp_gg )

    !
    ! done
    !

    call Done( ggi, status )
    IF_NOTOK_RETURN(status=1)

    call Done( leviX, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine tmm_Read_Convec_EC_gg

#endif


  ! *****************************************************************


#ifdef with_tmm_convec_ec

  subroutine tmm_Read_Convec_EC( tmm, archivekey, tday, t1, t2, lli, levi, &
                                gph, gph_tmi, &
                                sp, &
                                entu, entu_tmi, entd, entd_tmi, &
                                detu, detu_tmi, detd, detd_tmi, &
                                status )

    use GO      , only : goErr, gol, goPr
    use GO      , only : TDate, wrtgol
    use Phys    , only : convec_mfldet_to_entdet
    use Grid    , only : TllGridInfo
    use Grid    , only : TLevelInfo
    use tmm_info, only : TMeteoInfo, Init, AddHistory

    ! --- in/out --------------------------------

    type(TTmMeteo), intent(inout)     ::  tmm
    character(len=*), intent(in)      ::  archivekey
    type(TDate), intent(in)           ::  tday, t1, t2
    type(TllGridInfo), intent(in)     ::  lli
    type(TLevelInfo), intent(in)      ::  levi
    real, intent(in)                  ::  gph(:,:,:)   ! m   (half levels)
    type(TMeteoInfo), intent(in)      ::  gph_tmi
    real, intent(out)                 ::  sp(:,:)      ! Pa
    real, intent(out)                 ::  entu(:,:,:), entd(:,:,:)   ! kg/m2/s (full levels)
    type(TMeteoInfo), intent(out)     ::  entu_tmi, entd_tmi
    real, intent(out)                 ::  detu(:,:,:), detd(:,:,:)   ! kg/m2/s (full levels)
    type(TMeteoInfo), intent(out)     ::  detu_tmi, detd_tmi
    integer, intent(out)              ::  status

    ! --- const ------------------------------------

    character(len=*), parameter ::  rname = mname//'/tmm_Read_Convec_EC'

    ! --- local -------------------------------

    integer                  ::  lout
    real, allocatable        ::  ll(:,:,:)
    real, allocatable        ::  udmf(:,:,:)
    real, allocatable        ::  ddmf(:,:,:)
    real, allocatable        ::  uddr(:)
    real, allocatable        ::  dddr(:)

    integer                  ::  i, j, k

    ! --- begin -------------------------------

    ! number of levels in output arrays:
    lout = size(entu,3)

    ! check ...
    if ( (size(entd,3)/=lout) .or. &
         (size(detu,3)/=lout) .or. (size(detd,3)/=lout) ) then
      write (gol,'("output arrays should have same number of levels:")'); call goErr
      write (gol,'("  entu : ",i4)') size(entu,3); call goErr
      write (gol,'("  entd : ",i4)') size(entd,3); call goErr
      write (gol,'("  detu : ",i4)') size(detu,3); call goErr
      write (gol,'("  detd : ",i4)') size(detd,3); call goErr
      TRACEBACK; status=1; return
    end if

    ! full level output:
    allocate(   ll(lli%nlon,lli%nlat,levi%nlev  ) )
    ! generalized for lout=lmax_conv < levi%nlev
    allocate( udmf(lli%nlon,lli%nlat,lout+1) )
    allocate( ddmf(lli%nlon,lli%nlat,lout+1) )
    ! detrainment rates:
    allocate( uddr(lout) )
    allocate( dddr(lout) )

    ! updraught mass flux, half levels !
    call ReadField( tmm, archivekey, 'UDMF', tday, t1, t2, &
                           lli, 'n', levi, 'n', sp, ll, entu_tmi, status )
    IF_NOTOK_RETURN(status=1)
    ! store ...
    udmf(:,:,2:lout+1) = ll(:,:,1:lout)
    udmf(:,:,1) = 0.

    ! updraught detraintment rate
    call ReadField( tmm, archivekey, 'UDDR', tday, t1, t2, &
                           lli, 'n', levi, 'n', sp, ll, detu_tmi, status )
    IF_NOTOK_RETURN(status=1)
    ! store ...
    detu = ll(:,:,1:lout)

    ! downdraught mass flux, half levels !
    call ReadField( tmm, archivekey, 'DDMF', tday, t1, t2, &
                           lli, 'n', levi, 'n', sp, ll, entd_tmi, status )
    IF_NOTOK_RETURN(status=1)
    ! store ...
    ddmf(:,:,2:lout+1) = ll(:,:,1:lout)
    ddmf(:,:,1) = 0.

    ! downdraught detraintment rate
    call ReadField( tmm, archivekey, 'DDDR', tday, t1, t2, &
                           lli, 'n', levi, 'n', sp, ll, detd_tmi, status )
    IF_NOTOK_RETURN(status=1)
    ! store ...
    detd = ll(:,:,1:lout)

    ! convert from flux/detr to entr/detr
    do j = 1, lli%nlat
      do i = 1, lli%nlon
        ! copy detrainment rates:
        uddr = detu(i,j,:)
        dddr = detd(i,j,:)
        ! convert:
        call convec_mfldet_to_entdet( 'u', lout, gph(i,j,1:lout+1), &
                udmf(i,j,1:lout+1), uddr       , ddmf(i,j,1:lout+1), dddr       , &
                entu(i,j,:       ), detu(i,j,:), entd(i,j,:       ), detd(i,j,:), status )
        IF_NOTOK_RETURN(status=1)
      end do
    end do

    ! info ...
    call AddHistory( entu_tmi, 'converted mflux/detr to entr/detr', status )
    call AddHistory( detu_tmi, 'converted mflux/detr to entr/detr', status )
    call AddHistory( entd_tmi, 'converted mflux/detr to entr/detr', status )
    call AddHistory( detd_tmi, 'converted mflux/detr to entr/detr', status )

    ! clear
    deallocate( ll )
    deallocate( udmf )
    deallocate( ddmf )
    deallocate( uddr )
    deallocate( dddr )

    ! ok
    status = 0

  end subroutine tmm_Read_Convec_EC

#endif


  ! ==================================================================
  ! ===
  ! ===  Olsson surface roughness
  ! ===
  ! ==================================================================


  subroutine tmm_Read_SR_OLS( tmm, archivekey, tday, t1, t2, &
                                lli, ll, tmi, status )

    use GO      , only : goErr, gol, goPr
    use GO      , only : TDate, Get, wrtgol
    use GO      , only : goSplitLine
    use Grid    , only : TllGridInfo, Init, Done, Interpol
    use tmm_info, only : TMeteoInfo, Init, AddHistory

    ! --- in/out --------------------------------

    type(TTmMeteo), intent(inout)     ::  tmm
    character(len=*), intent(in)      ::  archivekey
    type(TDate), intent(in)           ::  tday, t1, t2
    type(TllGridInfo), intent(in)     ::  lli
    real, intent(out)                 ::  ll(:,:)
    type(TMeteoInfo), intent(out)     ::  tmi
    integer, intent(out)              ::  status

    ! --- const ------------------------------------

    character(len=*), parameter ::  rname = mname//'/tmm_Read_SR_OLS'

    ! grid size
    integer, parameter  ::  nlon = 361, nlat = 181

    ! --- local -------------------------------

    character(len=10)         ::  sourcetype
    character(len=MAX_RCKEY_LEN)        ::  sourcename

    integer                   ::  month
    character(len=MAX_FILENAME_LEN)        ::  fname
    logical                   ::  exist, opened
    integer                   ::  fu
    type(TllGridInfo)         ::  lli_ols
    real                      ::  sr_ols(nlon,nlat)
    integer                   ::  i, j

    ! --- begin -------------------------------

    ! split source key in type and name:
    call goSplitLine( archivekey, sourcetype, ':', sourcename, status )
    IF_NOTOK_RETURN(status=1)

    ! input TMPP fields or raw prism fields ?
    select case ( sourcetype )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! standard
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'standard' )

        ! dummy values:
        ll = 0.0

        ! set history info
        call Init( tmi, 'srols', 'm', status )
        call AddHistory( tmi, 'standard', status )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! read directly from hdf file
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'tmpp', 'tm5-hdf', 'tm5-nc' )

        ! read from surf file:
        call ReadField( tmm, archivekey, 'srols', 'm', tday, t1, t2, &
                             lli, 'n', ll, tmi, status )
        IF_NOTOK_RETURN(status=1)

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! convert raw data:
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'ecmwf-tmpp', 'ecmwf-tm5', 'ncep-cdc', 'ncep-gfs' )

        !call wrtgol( 'tmm read  : srols ', t1, '  -  ', t2 ); call goPr

        ! info ...
        call Init( tmi, 'srols', 'm', status )

        ! month
        call Get( t1, month=month )

        ! write filename:
        write (fname,'(a,"/SR_OLSSON-SR_OLSSON_360_180_",i2.2,".d")') trim(tmm%input_dir), month

        ! info ...
        call AddHistory( tmi, 'file=='//trim(fname), status )

        ! exist ?
        inquire( file=fname, exist=exist )
        if ( .not. exist ) then
          write (gol,'("Olsson SR file not found:")'); call goErr
          write (gol,'("  ",a)') trim(fname); call goErr
          TRACEBACK; status=1; return
        end if

        ! select free file unit:
        fu = 1234
        do
          inquire( unit=fu, opened=opened )
          if ( .not. opened ) exit
          fu = fu + 1
        end do

        ! open data file:
        open( fu, file=trim(fname), form='formatted', status='old', iostat=status )
        if (status/=0) then
          write (gol,'("while opening olsson data file:")'); call goErr
          write (gol,'("  ",a)') trim(fname); call goErr
          TRACEBACK; status=1; return
        end if

        ! read field:
        read (fu,*,iostat=status) sr_ols
        if (status/=0) then
          write (gol,'("while reading from olsson data file:")'); call goErr
          write (gol,'("  ",a)') trim(fname); call goErr
          TRACEBACK; status=1; return
        end if

        ! close file:
        close( fu, iostat=status )
        if (status/=0) then
          write (gol,'("while closing olsson data file:")'); call goErr
          write (gol,'("  ",a)') trim(fname); call goErr
          TRACEBACK; status=1; return
        end if

        ! setup grid definition:
        !   lon  :  -180.0  -179.0  ..  180.0    ( 1 deg resolution; 360 points; date line twice)
        !   lat  :   -90.0   -89.0  ..   90.0    ( 1 deg resolution; 180 points; includes poles)
        call Init( lli_ols, -180.00, 360.0/(nlon-1), nlon, -90.0, 180.0/(nlat-1), nlat, status )
        IF_NOTOK_RETURN(status=1)

        ! interpol
        do j = 1, lli%nlat
          do i = 1, lli%nlon
            call Interpol( lli_ols, sr_ols, lli%lon_deg(i), lli%lat_deg(j), ll(i,j) )
          end do
        end do

        ! info ...
        call AddHistory( tmi, 'horizontal_interpolation', status )

        ! done
        call Done( lli_ols, status )
        IF_NOTOK_RETURN(status=1)

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! error ...
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case default

        write (gol,'("unsupported source type `",a,"`")') trim(sourcetype); call goErr
        TRACEBACK; status=1; return

    end select

    ! ok
    status = 0

  end subroutine tmm_Read_SR_OLS


!  ! ==================================================================
!  ! ===
!  ! ===  pv/theta -> eqv.lat.
!  ! ===
!  ! ==================================================================
!
!
!  subroutine tmm_ReadEqvLat( tmm, archivekey, &
!                                  tday, t1, t2, &
!                                  lli, levi, &
!                                  sp, pv, theta, eqvlatb, eqvinds, &
!                                  status )
!
!    use GO, only : TDate, operator(-), operator(+), operator(/), goSplitLine
!    use Grid, only : TllGridInfo, TLevelInfo
!
!    use tmm_mf  , only : ReadEqvLatStuff
!    use tmm_info, only : TMeteoInfo
!
!    ! --- in/out --------------------------------
!
!    type(TTmMeteo), intent(inout)           ::  tmm
!
!    character(len=*), intent(in)            ::  archivekey
!
!    type(TDate), intent(in)                 ::  tday, t1, t2
!
!    type(TllGridInfo), intent(in)           ::  lli
!    type(TLevelInfo), intent(in)            ::  levi
!
!    real, intent(out)                       ::  sp(:,:)       ! Pa
!    real, intent(out)                       ::  pv(:,:,:)
!    real, intent(out)                       ::  theta(:,:,:)
!    real, intent(out)                       ::  eqvlatb(:,:)
!    integer, intent(out)                    ::  eqvinds(:,:,:)
!
!    integer, intent(out)                    ::  status
!
!    ! --- const --------------------------------------
!
!    character(len=*), parameter  ::  rname = mname//'/tmm_ReadEqvLat'
!
!    ! --- local -------------------------------
!
!    character(len=10)        ::  sourcetype
!    character(len=MAX_RCKEY_LEN)       ::  sourcename
!
!    integer                  ::  imf
!
!    type(TMeteoInfo)         ::  tmi
!
!    ! --- begin -------------------------------
!
!    ! split source key in type and name:
!    call goSplitLine( archivekey, sourcetype, ':', sourcename )
!
!    ! input TMPP fields or raw prism fields ?
!    select case ( sourcetype )
!
!      case ( 'tmpp' )
!
!        ! pv valid for [t1,t2]
!        call ReadField( tmm, archivekey, 'PVo', tday, t1, t2, &
!                             lli, 'n', levi, 'n', sp, pv, tmi, status )
!        IF_NOTOK_RETURN(status=1)
!
!        ! theta valid for [t1,t2]
!        call ReadField( tmm, archivekey, 'theta', tday, t1, t2, &
!                             lli, 'n', levi, 'n', sp, theta, tmi, status )
!        IF_NOTOK_RETURN(status=1)
!
!        !
!        ! eqv lat bounds and indices
!        !
!        ! select (eventually retrieve first) the meteo file with this param:
!        call SelectMF( tmm, 'i', archivekey, 'eqvlatb', tday, t1, t2, imf, status )
!        IF_NOTOK_RETURN(status=1)
!        !
!        ! read from selected file
!        call ReadEqvLatStuff( tmm%mf(imf), t1, t2, eqvlatb, eqvinds, status )
!        IF_NOTOK_RETURN(status=1)
!
!      case default
!
!        write (gol,'("unsupported source type `",a,"`")') trim(sourcetype)
!        TRACEBACK; status=1; return
!
!    end select
!
!    ! ok
!    status = 0
!
!  end subroutine tmm_ReadEqvLat


  ! ##########################################################################################
  ! ###
  ! ### output
  ! ###
  ! ##########################################################################################


  !
  ! call WriteField( tmmd, 'od-fc-ml60-glb3x2', 'T', 'K', tday, t1, t2, &
  !                        lli, 'n', sp, status )
  !

  subroutine tmm_WriteField_2d( tmm, archivekey, &
                                     tmi, paramkey, unit, tday, t1, t2, &
                                     lli, nuv, ll, status )

    use GO      , only : goErr, gol, goPr
    use GO      , only : TDate, wrtgol
    use Grid    , only : TllGridInfo
    use tmm_mf  , only : WriteRecord
    use tmm_info, only : TMeteoInfo

    ! --- in/out --------------------------------

    type(TTmMeteo), intent(inout)           ::  tmm
    character(len=*), intent(in)            ::  archivekey
    type(TMeteoInfo), intent(in)            ::  tmi
    character(len=*), intent(in)            ::  paramkey
    character(len=*), intent(in)            ::  unit
    type(TDate), intent(in)                 ::  tday, t1, t2
    type(TllGridInfo), intent(in)           ::  lli
    character(len=1), intent(in)            ::  nuv
    real, intent(in)                        ::  ll(:,:)
    integer, intent(out)                    ::  status

    ! --- const ------------------------------------

    character(len=*), parameter ::  rname = mname//'/tmm_WriteField_2d'

    ! --- local ----------------------------------------

    integer          ::  imf

    ! --- begin ----------------------------------

    !write(gol, '("tmm write : ", a, 2x, a, " - ", a, 2x, i3, " x ", i3)') &
        !trim(paramkey), trim(Pretty(t1)), trim(Pretty(t2)), lli%nlat, lli%nlon ; call goPr
    call wrtgol( 'tmm write : '//trim(paramkey)//' ', t1, '  -  ', t2 ); call goPr

    ! check shape of grid:
    if ( ((nuv == 'n') .and. ((size(ll,1) /= lli%nlon  ) .or. (size(ll,2) /= lli%nlat  ))) .or. &
         ((nuv == 'u') .and. ((size(ll,1) /= lli%nlon+1) .or. (size(ll,2) /= lli%nlat  ))) .or. &
         ((nuv == 'v') .and. ((size(ll,1) /= lli%nlon  ) .or. (size(ll,2) /= lli%nlat+1))) ) then
      write (gol,'("2d array does not mach with grid definition:")'); call goErr
      write (gol,'("  param  : ",a          )') paramkey   ; call goErr
      write (gol,'("  lli    : ",i3," x ",i3)') lli%nlon, lli%nlat; call goErr
      write (gol,'("  nuv    : ",a          )') nuv; call goErr
      write (gol,'("  ll     : ",i3," x ",i3)') shape(ll); call goErr
      TRACEBACK; status=1; return
    end if

    ! select index of already open meteo file or setup access to new one;
    call SelectMF( tmm, 'o', archivekey, paramkey, tday, t1, t2, imf, status )
    IF_NOTOK_RETURN(status=1)

    ! write
    call WriteRecord( tmm%mf(imf), tmi, paramkey, unit, tday, t1, t2, &
                        lli, nuv, ll, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine tmm_WriteField_2d


  !
  ! call WriteField( tmmd, 'od-fc-ml60-glb3x2', 'T', 'K', tday, t1, t2, &
  !                        lli, 'n', levi, spm, temper, status )
  !

  subroutine tmm_WriteField_3d( tmm, archivekey, &
                                     tmi, spname, paramkey, unit, tday, t1, t2, &
                                     lli, nuv, levi, nw, sp, ll, status )!, &
                                     !nlev )

    use GO      , only : goErr, gol, goPr
    use GO      , only : gol, goPr
    use GO      , only : TDate, wrtgol
    use Grid    , only : TllGridInfo, TLevelInfo
    use tmm_mf  , only : WriteRecord
    use tmm_info, only : TMeteoInfo

    ! --- in/out --------------------------------

    type(TTmMeteo), intent(inout)           ::  tmm
    character(len=*), intent(in)            ::  archivekey
    type(TMeteoInfo), intent(in)            ::  tmi
    character(len=*), intent(in)            ::  spname
    character(len=*), intent(in)            ::  paramkey
    character(len=*), intent(in)            ::  unit
    type(TDate), intent(in)                 ::  tday, t1, t2
    type(TllGridInfo), intent(in)           ::  lli
    character(len=1), intent(in)            ::  nuv
    type(TLevelInfo), intent(in)            ::  levi
    character(len=1), intent(in)            ::  nw
    real, intent(in)                        ::  sp(:,:)        ! Pa
    real, intent(in)                        ::  ll(:,:,:)
    integer, intent(out)                    ::  status

    !integer, intent(in), optional           ::  nlev

    ! --- const ------------------------------------

    character(len=*), parameter ::  rname = mname//'/tmm_WriteField_3d'

    ! --- local ----------------------------------------

    integer          ::  imf

    ! --- begin ----------------------------------

    call wrtgol( 'tmm write : '//trim(paramkey)//' ', t1, '  -  ', t2 ); call goPr

    ! check shape of grid:
    if ( ((nuv == 'n') .and. ((size(ll,1) /= lli%nlon  ) .or. (size(ll,2) /= lli%nlat  ))) .or. &
         ((nuv == 'u') .and. ((size(ll,1) /= lli%nlon+1) .or. (size(ll,2) /= lli%nlat  ))) .or. &
         ((nuv == 'v') .and. ((size(ll,1) /= lli%nlon  ) .or. (size(ll,2) /= lli%nlat+1))) .or. &
         ((nuv == 'n') .and. ((size(sp,1) /= lli%nlon  ) .or. (size(sp,2) /= lli%nlat  ))) .or. &
         ((nuv == 'u') .and. ((size(sp,1) /= lli%nlon+1) .or. (size(sp,2) /= lli%nlat  ))) .or. &
         ((nuv == 'v') .and. ((size(sp,1) /= lli%nlon  ) .or. (size(sp,2) /= lli%nlat+1))) .or. &
         ((nw  == 'n') .and. (size(ll,3) > levi%nlev  )) .or. &
         ((nw  == 'w') .and. (size(ll,3) > levi%nlev+1)) ) then
      write (gol,'("3d arrays do not match with grid definition:")'); call goErr
      write (gol,'("  param  : ",a          )') paramkey; call goErr
      write (gol,'("  lli    : ",i3," x ",i3         )') lli%nlon, lli%nlat; call goErr
      write (gol,'("  nuv    : ",a                   )') nuv; call goErr
      write (gol,'("  levi   : ",i3                  )') levi%nlev; call goErr
      write (gol,'("  nw     : ",a                   )') nw; call goErr
      write (gol,'("  sp     : ",i3," x ",i3         )') shape(sp); call goErr
      write (gol,'("  ll     : ",i3," x ",i3," x ",i3)') shape(ll); call goErr
      TRACEBACK; status=1; return
    end if

    ! select index of already open meteo file or setup access to new one;
    call SelectMF( tmm, 'o', archivekey, paramkey, tday, t1, t2, imf, status )
    IF_NOTOK_RETURN(status=1)

    ! write
    call WriteRecord( tmm%mf(imf), tmi, spname, paramkey, unit, tday, t1, t2, &
                        lli, nuv, levi, nw, sp, ll, status )!, &
                        !nlev=nlev )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine tmm_WriteField_3d




end module TMM
