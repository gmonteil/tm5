!###############################################################################
!
! Input/output of meteofiles produced by TMPP.
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

module tmm_mf_tmpp

  use GO      , only : gol, goErr, goPr, goBug
  use GO      , only : TDate
  use file_hdf, only : THdfFile, TSds
  use os_specs, only : MAX_FILENAME_LEN, MAX_RCKEY_LEN

  implicit none

  ! --- in/out ----------------------------

  private

  public  ::  TMeteoFile_tmpp
  public  ::  Init, Done
  public  ::  Get
  public  ::  ReadRecord
  public  ::  WriteRecord

  ! --- const ------------------------------

  character(len=*), parameter  ::  mname = 'tmm_mf_tmpp'

  ! ~~~ output keys and defaults

  ! current format version
  character(len=*), parameter  ::  output_format = 'tmm-1.0'

  ! integer/real kinds
  integer, parameter   ::  iknd    = 4        ! integer
  integer, parameter   ::  rknd    = 8        ! real
  integer, parameter   ::  iknd_ds = 2        ! integer for big ds
  integer, parameter   ::  rknd_ds = 4        ! real    for big ds

  ! z compression
  character(len=*), parameter  ::  compression = 'deflate'
  integer, parameter           ::  deflate_level = 6


  !--- type ---------------------------------

  type TMeteoFile_tmpp
    ! input/output ?
    character(len=1)       ::  io
    !
    ! field collection
    !
    character(len=MAX_FILENAME_LEN) ::  fname
    type(THdfFile)                  ::  hdf
    character(len=MAX_RCKEY_LEN)    ::  paramkeys
    type(TDate)                     ::  trange(2)
    !
    ! file
    !
    type(TSds)             ::  sds
    logical                ::  selected
    character(len=16)      ::  paramkey
    !
    ! pw specials
    !
    logical               ::  mfw_redir
    !
    ! surface stress specials
    !
    logical                ::  sstr_to_ewss_nsss
    !
    ! surface pressure field
    !
    logical                ::  spm_load      ! load spm ?
    logical                ::  spm_selected  ! data set filled ?
    type(TSds)             ::  spm_sds       ! data set for sp
    logical                ::  spm_incl      ! record included in main file ?
    logical                ::  spm_extr      ! record in external spm file ?
    character(len=MAX_FILENAME_LEN) ::  spm_fname
    type(THdfFile)         ::  spm_hdf
    logical                ::  spm_n_to_uv
    !
    ! output
    !
    logical                ::  output_initialised
    integer                ::  output_nrec
    integer                ::  output_ntrec
    character(len=20)      ::  output_names(20)
    integer                ::  output_nname
    !
    ! adhoc ...
    integer                ::  fixyear
    logical                ::  qad
  end type TMeteoFile_tmpp


  ! --- interfaces -------------------------

  interface Init
    module procedure mf_Init
  end interface

  interface Done
    module procedure mf_Done
  end interface

  interface Get
    module procedure mf_Get
  end interface

  interface ReadRecord
    module procedure mf_ReadRecord
  end interface

!  interface ReadEqvLatStuff
!    module procedure mf_ReadEqvLatStuff
!  end interface

  interface WriteRecord
    module procedure mf_WriteRecord_2d
    module procedure mf_WriteRecord_3d
  end interface


contains


  ! ==============================================================


  subroutine mf_Init( mf, io, dir, archivekeys, paramkey, &
                               tref, t1, t2, status )

    use GO, only : TDate, Set, Get, NewDate, AnyDate, IsAnyDate
    use GO, only : rTotal, operator(-), operator(>=)
    use GO, only : goVarValue, goWriteKeyNum

    ! --- in/out ----------------------------

    type(TMeteoFile_tmpp), intent(out)  ::  mf
    character(len=1), intent(in)        ::  io
    character(len=*), intent(in)        ::  dir
    character(len=*), intent(in)        ::  archivekeys
    character(len=*), intent(in)        ::  paramkey
    type(TDate), intent(in)             ::  tref, t1, t2
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_Init'

    ! --- local --------------------------------

    character(len=64)    ::  mf_mdir, mf_tres
    character(len=64)    ::  mf_class, mf_type, mf_grid, mf_levs
    character(len=1)     ::  mf_psep, mf_nsep

    character(len=64)    ::  mf_filekey

    character(len=4)     ::  mf_fckey

    type(TDate)           ::  tfile
    integer               ::  ccyy, mm, dd, dh
    type(TDate)           ::  tc

    logical               ::  exist

    ! --- begin --------------------------------

    ! store i/o :
    mf%io = io

    ! default flags:
    mf%sstr_to_ewss_nsss = .false.
    mf%mfw_redir         = .false.

    !
    ! extract fields from archivekey :
    !   mdir=ec-fg_3h-ml60-glb3x2;tres=_21p06
    !
    mf%fixyear  = -1
    call goVarValue( archivekeys, ';', 'fixyear', '=', mf%fixyear, status )
    if (status>0) then; TRACEBACK; status=1; return; end if

    mf_tres  = 'no_tres'
    call goVarValue( archivekeys, ';', 'tres', '=', mf_tres, status )
    if (status>0) then; TRACEBACK; status=1; return; end if
    !
    mf_class = 'no_class'
    call goVarValue( archivekeys, ';', 'class', '=', mf_class, status )
    if (status>0) then; TRACEBACK; status=1; return; end if
    !
    mf_type  = 'no_type'
    call goVarValue( archivekeys, ';', 'type', '=', mf_type, status )
    if (status>0) then; TRACEBACK; status=1; return; end if
    !
    mf_grid  = 'no_grid'
    call goVarValue( archivekeys, ';', 'grid', '=', mf_grid, status )
    if (status>0) then; TRACEBACK; status=1; return; end if
    !
    mf_levs  = 'no_mlevs'
    call goVarValue( archivekeys, ';', 'levs', '=', mf_levs, status )
    if (status>0) then; TRACEBACK; status=1; return; end if
    !
    ! path seperation character:
    mf_psep  = '/'
    call goVarValue( archivekeys, ';', 'pathsep', '=', mf_psep, status )
    if (status>0) then; TRACEBACK; status=1; return; end if
    !
    ! name seperation character:
    mf_nsep  = '-'
    call goVarValue( archivekeys, ';', 'namesep', '=', mf_nsep, status )
    if (status>0) then; TRACEBACK; status=1; return; end if
    !
    ! quick and dirty time checks ?
    mf%qad = .false.
    !call goVarValue( archivekeys, ';', 'qad', '=', mf%qad, status )
    !if (status>0) then; TRACEBACK; status=1; return; end if
    !
    ! adhoc flag
    mf%sstr_to_ewss_nsss = .false.
    call goVarValue( archivekeys, ';', 'sstr', '=', mf%sstr_to_ewss_nsss, status )
    if (status>0) then; TRACEBACK; status=1; return; end if

    !
    ! main file
    !

    ! by default, no surface pressure stuff ...
    mf%spm_load    = .false.
    mf%spm_incl    = .false.   ! no spm fields included in hdf file
    mf%spm_extr    = .false.   ! no extra spm file
    mf%spm_n_to_uv = .false.   ! no interpolation from 'n' to 'u' or 'v'

    ! default time resolution for tmpp produced files:
    mf_tres = '_21p06'

    ! * set mf_filekey (uvsp,t,etc) and parmeters:
    select case ( paramkey )
      case ( 'sp', 'pu', 'pv', 'mfu', 'mfv' )
        mf_filekey = 'uvsp'
        mf%paramkeys = '-sp-pu-pv-mfu-mfv-'
        mf%spm_load  = .true.
      case ( 'pw', 'mfw' )
        mf_filekey = 'w'
        mf%paramkeys = '-pw-mfw-'
        mf%spm_load  = .true.
      case ( 'T' )
        mf_filekey   = 't'
        mf%paramkeys = '-T-'
        mf%spm_load  = .true.
      case ( 'Q' )
        mf_filekey   = 'q'
        mf%paramkeys = '-Q-'
        mf%spm_load  = .true.
      case ( 'CLWC', 'CIWC', 'CC', 'CCO', 'CCU', &
             'clwc', 'ciwc', 'cc', 'cco', 'ccu'  )
        mf_filekey   = 'cld'
        mf%paramkeys = '-CLWC-CIWC-CC-CCO-CCU-'
        mf%spm_load  = .true.
      case ( 'eu', 'ed', 'du', 'dd', 'cloud_base', 'cloud_top', 'cloud_lfs' )
        mf_filekey   = 'sub'
        mf%paramkeys = '-eu-ed-du-dd-cloud_base-cloud_top-cloud_lfs-'
        mf%spm_load  = .true.
      case ( 'oro', 'lsm', 'sr', 'srols', 'srmer', &
             'cvl', 'cvh', &
             'tv01', 'tv02', 'tv03', 'tv04', 'tv05', &
             'tv06', 'tv07', 'tv08', 'tv09', 'tv10', &
             'tv11', 'tv12', 'tv13', 'tv14', 'tv15', &
             'tv16', 'tv17', 'tv18', 'tv19', 'tv20', &
             'swvl1', &
             'albedo', 'lsrh', 'ci', 'g10m', 'u10m', 'v10m', 'sd', &
             'lsp', 'cp', 'sf', 'sshf', 'slhf', 'blh', &
             't2m', 'd2m', 'ssr', 'sstr', 'src', 'raero', 'ustar', &
             'sst', 'sps', &
             'ewss', 'nsss' )
        mf_levs = 'sfc'
        mf_grid = 'glb1x1'
        mf_filekey = 'surf'
        mf_tres = '_21p03'
        mf%paramkeys= '-'
        mf%paramkeys= trim(mf%paramkeys)//'oro-lsm-sr-srols-srmer-'
        mf%paramkeys= trim(mf%paramkeys)//'cvl-cvh-'
        mf%paramkeys= trim(mf%paramkeys)//'tv01-tv02-tv03-tv04-tv05-'
        mf%paramkeys= trim(mf%paramkeys)//'tv06-tv07-tv08-tv09-tv10-'
        mf%paramkeys= trim(mf%paramkeys)//'tv11-tv12-tv13-tv14-tv15-'
        mf%paramkeys= trim(mf%paramkeys)//'tv16-tv17-tv18-tv19-tv20-'
        mf%paramkeys= trim(mf%paramkeys)//'swvl1-'
        mf%paramkeys= trim(mf%paramkeys)//'albedo-lsrh-ci-10fg-u10m-v10m-sd-'
        mf%paramkeys= trim(mf%paramkeys)//'lsp-cp-sf-sshf-slhf-blh-'
        mf%paramkeys= trim(mf%paramkeys)//'t2m-d2m-ssr-sstr-src-raero-ustar-'
        mf%paramkeys= trim(mf%paramkeys)//'sst-sps-'
        mf%paramkeys= trim(mf%paramkeys)//'ewss-nsss-'
      case ( 'spm' )
        mf_levs = 'ml1'
        mf_filekey = 'spm'
        mf_tres = '_00p06'
        mf%paramkeys = '-spm-'
      case default
        write (gol,'("unsupported paramkey `",a,"`")') paramkey; call goErr
        TRACEBACK; status=1; return
    end select

    ! convert input times to file name times:
    call GetTime( mf_filekey, mf_tres, tref, t1, t2, status, &
                      tfile=tfile, trange=mf%trange )
    IF_NOTOK_RETURN(status=1)

    ! adhoc: fixed year ?
    if ( mf%fixyear > 0 ) call Set( tfile, year=mf%fixyear )

    ! extract time values:
    call Get( tfile, year=ccyy, month=mm, day=dd )

    ! main file:
    write (mf%fname,'(a,a,a,a,a,a,i4.4,a,i2.2,a,a,a,a,a,a,"_",i4.4,2i2.2,a,a)') &
              trim(dir), mf_psep, &
              trim(mf_class), mf_nsep, trim(mf_type), mf_nsep, &
              ccyy, mf_nsep, mm, mf_nsep, trim(mf_levs), mf_nsep, trim(mf_grid), mf_nsep, &
              trim(mf_filekey), ccyy, mm, dd, trim(mf_tres), '.hdf'

    ! pw specials
    mf%mfw_redir   = .true.

    ! load surface pressure with 3d field ?
    if ( mf%spm_load ) then
      ! external surface pressure file:
      mf%spm_extr = .true.
      write (mf%spm_fname,'(a,a,a,a,a,a,i4.4,a,i2.2,a,a,a,a,a,a,"_",i4.4,2i2.2,a,a)') &
                trim(dir), mf_psep, &
                trim(mf_class), mf_nsep, trim(mf_type), mf_nsep, &
                ccyy, mf_nsep, mm, mf_nsep, 'ml1', mf_nsep, trim(mf_grid), mf_nsep, &
                'spm', ccyy, mm, dd, '_00p06', '.hdf'
      ! only cell surface pressure in tmpp output ...
      mf%spm_n_to_uv = .true.
    end if

    ! ok
    status = 0

  end subroutine mf_Init


  ! ***


  subroutine mf_Done( mf, status )

    !use file_hdf, only : Done

    ! --- in/out ------------------------------------

    type(TMeteoFile_tmpp), intent(inout)::  mf
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_Done'

    ! --- begin -------------------------------------

    ! files have been closed in ReadRecord/WriteRecord

    ! ok
    status = 0

  end subroutine mf_Done



  ! ***


  subroutine mf_Get( mf, status, trange1, trange2, paramkeys )

    use GO, only : TDate

    ! --- in/out ----------------------------

    type(TMeteoFile_tmpp), intent(in)         ::  mf
    integer, intent(out)                      ::  status

    type(TDate), intent(out), optional        ::  trange1, trange2
    character(len=*), intent(out), optional   ::  paramkeys

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_Get'

    ! --- local --------------------------------

    ! --- begin --------------------------------

    ! time range:
    if ( present(trange1) ) trange1 = mf%trange(1)
    if ( present(trange2) ) trange2 = mf%trange(2)

    ! parameter names:
    if ( present(paramkeys) ) paramkeys = mf%paramkeys

    ! ok
    status = 0

  end subroutine mf_Get



  ! ******************************************************************
  ! ***
  ! *** time range, parameters, file names
  ! ***
  ! ******************************************************************


  !
  ! Return time parameters:
  !  o tfile   :  date in filename
  !  o trange  :  time interval covered by fields in file
  !

  subroutine GetTime( filekey, tres, tref, t1, t2, status, tfile, trange, nrec )

    use GO, only : TDate, NewDate, AnyDate, Get, Set, wrtgol, IncrDate, IsAnyDate
    use GO, only : operator(<), operator(+), operator(-), rTotal

    ! --- in/out --------------------------------

    character(len=*), intent(in)            ::  filekey
    character(len=*), intent(in)            ::  tres
    type(TDate), intent(in)                 ::  tref, t1, t2
    integer, intent(out)                    ::  status

    type(TDate), intent(out), optional      ::  tfile
    type(TDate), intent(out), optional      ::  trange(2)
    integer, intent(out), optional          ::  nrec

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/GetTime'

    ! --- local --------------------------------

    integer          ::  year, month
    integer          ::  hour1, time6(6)
    integer          ::  dd, hh, step
    logical          ::  interval
    real             ::  dhr

    ! --- begin --------------------------------

    ! set day shift, start hour, and step
    select case ( tres )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! tmpp [21,21]
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( '_21p06', '_21p03', '_av21' )

        ! routine is called with  tref,t1,t2:
        !    t1,t1,t2
        !    t1,any,any    (oro and other constant fields)
        ! thus use tref to construct the file times
        if ( IsAnyDate(t1) ) then
          call Get( tref, hour=hour1 )
          interval = .false.
        else
          call Get( t1, hour=hour1 )
          interval = t1 < t2
        end if

        ! file ccyymmdd contains fields for (21,21];
        ! only uvsp is valid for [21,21] since it contains surface pressure for 21:00
        if ( present(tfile) ) then
          tfile = tref
          call Set( tfile, hour=0, min=0, sec=0 )
          if ( (hour1 > 21) .or. ((interval .or. filekey=='uvsp') .and. hour1==21) ) then
            tfile = tfile + IncrDate(day=1)
          end if
       end if

        ! fields by default valid for (21,21];
        ! only uvsp is valid for [21,21] since it contains surface pressure for 21:00
        if ( present(trange) ) then
          trange(1) = tref
          call Set( trange(1), hour=0, min=0, sec=0 )  !  00:00 today
          if ( (hour1 > 21) .or. ((interval .or. filekey=='uvsp') .and. hour1==21) ) then
            trange(1) = trange(1) + IncrDate(day=1)
          end if
          trange(1) = trange(1) - IncrDate(hour=3)     !  previous 21:00
          trange(2) = trange(1) + IncrDate(day=1)      !  next     21:00
          ! boundary not included in most cases:
          if ( filekey /= 'uvsp' ) trange(1) = trange(1) + IncrDate(mili=1)
        end if

        ! number of records in file:
        if ( present(nrec) ) then
          select case ( tres )
            case ( '_21p06' ) ; nrec = 24/6
            case ( '_21p03' ) ; nrec = 24/3
            case ( '_av21'  ) ; nrec = 24/24
            case default
              write (gol,'("unsupported tres for setting nrec : ",a)') tres; call goErr
              TRACEBACK; status=1; return
          end select
        end if

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! tm5 constant
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'constant' )

        ! no date in filename ...
        if ( present(tfile) ) tfile = AnyDate()

        ! fields always valid ...
        if ( present(trange) ) then
          trange(1) = AnyDate()
          trange(2) = AnyDate()
        end if

        ! only one output record in constant file:
        if ( present(nrec) ) nrec = 1

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! tm5 monthly file
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'month' )

        ! file ccyymmdd contains fields for this month:
        if ( present(tfile) ) then
          call Get( t1, year=year, month=month )
          tfile = NewDate( year=year, month=month, day=1 )
        end if

        ! field valid from begin to end of month:
        if ( present(trange) ) then
          call Get( t1, year=year, month=month )
          trange(1) = NewDate( year=year, month=month, day=1, hour=00  )
          month = month + 1
          if ( month > 12 ) then
            year = year + 1
            month = 1
          end if
          trange(2) = NewDate( year=year, month=month, day=1, hour=00  )
        end if

        ! only one output record in month file:
        if ( present(nrec) ) nrec = 1

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! tm5 [00,24]
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( '_00p06', '_00p03', '_an0tr6', '_fg006up4tr3', '_fc012up2tr3', '00p01' )

        ! file ccyymmdd contains fields for [00,24) :
        if ( present(tfile) ) then
          tfile = t1
          call Set( tfile, hour=0, min=0, sec=0 )
        end if

        ! fields valid for [00,24) :
        if ( present(trange) ) then
          trange(1) = t1
          call Set( trange(1), hour=0, min=0, sec=0 )  !  00:00 today
          trange(2) = trange(1) + IncrDate(hour=24) - IncrDate(mili=1)
        end if

        ! number of records in file:
        if ( present(nrec) ) then
          select case ( tres )
            case ( '_00p06'  ) ; nrec = 24/6
            case ( '_an0tr6' ) ; nrec = 24/6
            case ( '_00p03', '_fg006up4tr3', '_fc012up2tr3' )
              ! by default: 3 hourly files
              ! for forecasts after 12+72, only 6 hourly available:
              !  f0  [ 00, 24)  : 00 03 06 09 12 15 18 21   :  nrec=8
              !  f1  [ 24, 48)  : 00 03 06 09 12 15 18 21   :  nrec=8
              !  f2  [ 48, 72)  : 00 03 06 09 12 15 18 21   :  nrec=8
              !  f3  [ 72, 96)  : 00 03 06 09 12    18      :  nrec=6
              !  f4  [ 96,120)  : 00    06    12    18      :  nrec=4
              !   :
              !  f9  [192,216)  : 00    06    12    18      :  nrec=4
              !  f10 [216,240)  : 00    06    12            :  nrec=3
              dhr = rTotal( t1 - tref, 'hour' )
              if ( dhr < 72.0 ) then
                nrec = 8
              else if ( dhr < 96 ) then
                nrec = 6
              else if ( dhr < 216 ) then
                nrec = 4
              else
                nrec = 3
              end if
            case ( '_00p01'       ) ; nrec = 24/24
            case default
              write (gol,'("unsupported tres for setting nrec : ",a)') tres; call goErr
              TRACEBACK; status=1; return
          end select
        end if

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! ???
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case default

        write (gol,'("unsupported time resolution key:")'); call goErr
        write (gol,'("  ",a)') trim(tres); call goErr
        TRACEBACK; status=1; return

    end select


    ! ok
    status = 0

  end subroutine GetTime


  ! ******************************************************************
  ! ***
  ! *** input
  ! ***
  ! ******************************************************************


  subroutine mf_SelectRecord( mf, paramkey, t1, t2, status )

    use GO, only : TDate, NewDate, IncrDate, wrtgol, Set, Get, rTotal
    use GO, only : operator(-), operator(+), operator(/), operator(==), operator(/=)
    use GO, only : operator(>=), operator(<=), operator(>)
    use file_hdf, only : GetInfo
    use file_hdf, only : TSds, Select, CheckInfo, CheckAttribute, ReadAttribute

    ! --- in/out --------------------------------

    type(TMeteoFile_tmpp), intent(inout)    ::  mf
    character(len=*), intent(in)            ::  paramkey
    type(TDate), intent(in)                 ::  t1, t2
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_SelectRecord'

    ! --- local -------------------------------

    integer            ::  nsds, isds
    integer            ::  time1(6), time2(6)
    type(TDate)        ::  tmid, treq, thelp
    type(TDate)        ::  t_spm
    integer            ::  hour, dhour

    integer            ::  status1, status2
    integer            ::  time1s(6), time2s(6)
    type(TDate)        ::  tim1s, tim2s

    ! --- begin -------------------------------

    mf%selected = .false.

    ! initial no record is found ...
    status = 1

    ! number of data sets:
    call GetInfo( mf%hdf, status, num_datasets=nsds )
    if (status/=0) then; TRACEBACK; return; end if

    ! loop over all data sets
    do isds = 1, nsds

      call Select( mf%sds, mf%hdf, isds-1, status )
      if (status/=0) then; TRACEBACK; return; end if

      ! correct param ?
      status=-1; call CheckInfo( mf%sds, status, name=paramkey )
      if (status>0) then; TRACEBACK; status=1; return; end if

      ! adhoc: some data sets have other name ...
      if ( status < 0 ) then
        select case ( paramkey )
          case ( 'sp', 'spm' )
            status=-1; call CheckInfo( mf%sds, status, name='ps' )
          case ( 'sps' )
            status=-1; call CheckInfo( mf%sds, status, name='sp' )
          case ( 'mfu' )
            status=-1; call CheckInfo( mf%sds, status, name='pu' )
          case ( 'mfv' )
            status=-1; call CheckInfo( mf%sds, status, name='pv' )
          case ( 'mfw' )
            status=-1; call CheckInfo( mf%sds, status, name='pw' )
          case ( 'PVo' )
            status=-1; call CheckInfo( mf%sds, status, name='PV' )
          case ( 'sr' )
            status=-1; call CheckInfo( mf%sds, status, name='sr_ecm' )
          case ( 'srols' )
            status=-1; call CheckInfo( mf%sds, status, name='sr_ols' )
          case ( 'srmer' )
            status=-1; call CheckInfo( mf%sds, status, name='sr_mer' )
          case ( 'albedo' )
            status=-1; call CheckInfo( mf%sds, status, name='al' )
          case ( 'ewss', 'nsss' )
            if ( mf%sstr_to_ewss_nsss ) then
              status=-1; call CheckInfo( mf%sds, status, name='sstr' )
            end if
        end select
        if (status>0) then; TRACEBACK; status=1; return; end if
        ! try other data set ?
        if ( status < 0 ) cycle
      end if

      ! correct time ?
      select case ( paramkey )
        case ( 'oro' )
          ! constant field, skip time check
        case ( 'srols' )
          ! monthly field:
          !   in TMPP files: daily samples, valid for [21,21]
          !   in TM5  files: monthly value
          ! read time2 :
          call ReadAttribute( mf%sds, 'time2', time2s, status )
          if (status/=0) then; TRACEBACK; return; end if
          tim2s = NewDate( time6=time2s )
          ! [t1,t2] covers month; time2 should be in this month:
          if ( (t1 >= tim2s) .and. (tim2s <= t2) ) then
            ! ok
            status = 0
          end if
        case default
          ! extract time arrays:
          call Get( t1, time6=time1 )
          call Get( t2, time6=time2 )
          !
          ! replace year ?
          if ( mf%fixyear > 0 ) then
            ! year is valid for time interval:
            !   ( fixyear-1/12/31 21:00 , fixyear/12/31 21:00 ]
            if ( t1 == t2 ) then
              ! instant time
              thelp = NewDate( year=time1(1), month=12, day=31, hour=21 )
              if ( t1 > thelp ) then
                time1(1) = mf%fixyear - 1
              else
                time1(1) = mf%fixyear
              end if
              time2 = time1
            else
              ! time interval
              thelp = NewDate( year=time1(1), month=12, day=31, hour=21 )
              if ( t1 >= thelp ) then
                time1(1) = mf%fixyear - 1
              else
                time1(1) = mf%fixyear
              end if
              thelp = NewDate( time6=time1 ) + (t2-t1)
              call Get( thelp, time6=time2 )
            endif
          end if
          !
          if ( mf%qad ) then
            ! test if times match exactly:
            status1=-1 ; call CheckAttribute( mf%sds, 'time1', time1, status1 )
            if (status1>0) then; TRACEBACK; status=1; return; end if
            status2=-1 ; call CheckAttribute( mf%sds, 'time2', time2, status2 )
            if (status2>0) then; TRACEBACK; status=1; return; end if
            status = status1 + status2
            ! times do not match ?
            if ( status < 0 ) then
              ! read time1 and time2 :
              call ReadAttribute( mf%sds, 'time1', time1s, status )
              if (status/=0) then; TRACEBACK; return; end if
              call ReadAttribute( mf%sds, 'time2', time2s, status )
              if (status/=0) then; TRACEBACK; return; end if
              ! hours should match ...
              if ( (time1s(4) /= time1(4)) .or. (time2s(4) /= time2(4)) ) then
                status = -1
                cycle
              end if
              ! warning ...
              write (gol,'("WARNING - weak time check passed:")'); call goPr
              write (gol,'("WARNING -   parameter : ",a)') trim(paramkey); call goPr
              call wrtgol( 'WARNING -   t1        : ', t1 ); call goPr
              call wrtgol( 'WARNING -   t2        : ', t2 ); call goPr
              write (gol,'("WARNING -   file      : ",a)') trim(mf%fname); call goPr
              ! ok
              status = 0
            end if
            ! weak time check passed
          else
            !! tmpp fields valid for [21,21] should be read as valid for [00,24] ...
            !if ( rTotal(t2-t1,'hour') == 24.0 ) then
            !  tim1s = t1 + IncrDate(hour=3)
            !  tim2s = t2 + IncrDate(hour=3)
            !  call Get( tim1s, time6=time1 )
            !  call Get( tim2s, time6=time2 )
            !end if
            ! first try wether times exactely match;
            status1=-1 ; call CheckAttribute( mf%sds, 'time1', time1, status1 )
            if (status1>0) then; TRACEBACK; status=1; return; end if
            status2=-1 ; call CheckAttribute( mf%sds, 'time2', time2, status2 )
            if (status2>0) then; TRACEBACK; status=1; return; end if
            status = status1 + status2
            ! try mid time ?
            if ( status < 0 ) then
              if ( t1 == t2 ) then
                ! instant time requested; check wether this is mid of times in file:
                call ReadAttribute( mf%sds, 'time1', time1s, status )
                IF_NOTOK_RETURN(status=1)
                tim1s = NewDate( time6=time1s )
                call ReadAttribute( mf%sds, 'time2', time2s, status )
                IF_NOTOK_RETURN(status=1)
                tim2s = NewDate( time6=time2s )
                ! mid of times in file:
                tmid = tim1s + (tim2s-tim1s)/2
                ! requested time:
                treq = t1
                if ( mf%fixyear > 0 ) call Set( treq, year=mf%fixyear )
                ! not ok ? then try again
                if ( treq /= tmid ) then
                  status = -1
                  cycle
                end if
              else
                ! interval requested; check wether time in file is mid time:
                tmid = t1 + (t2-t1)/2
                if ( mf%fixyear > 0 ) call Set( tmid, year=mf%fixyear )
                call Get( tmid, time6=time1 )
                status1=-1 ; call CheckAttribute( mf%sds, 'time1', time1, status1 )
                if (status1>0) then; TRACEBACK; status=1; return; end if
                status2=-1 ; call CheckAttribute( mf%sds, 'time2', time1, status2 )
                if (status2>0) then; TRACEBACK; status=1; return; end if
                status = status1 + status2
                if ( status < 0 ) cycle
              end if
            end if
            ! time check passed
          end if
      end select

      ! found!
      exit
    end do

    ! not found ?
    if ( status /= 0 ) then
      write (gol,'("Unable to locate field in hdf file:")'); call goErr
      write (gol,'("  parameter : ",a)') trim(paramkey); call goErr
      call wrtgol( '  t1        : ', t1 ); call goErr
      call wrtgol( '  t2        : ', t2 ); call goErr
      call wrtgol( '  mf%tr(1)  : ', mf%trange(1) ); call goErr
      call wrtgol( '  mf%tr(2)  : ', mf%trange(2) ); call goErr
      write (gol,'("  file      : ",a)') trim(mf%fname); call goErr
      TRACEBACK; status = 1; return
    end if

    ! store paramkey
    mf%paramkey = paramkey

    ! ok
    mf%selected = .true.


    ! ~~~ surface pressure ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ( mf%spm_load ) then

      mf%spm_selected = .false.

      if ( paramkey /= 'sp' ) then

        ! only hours 00, 06, 12, 18
        dhour = nint(rTotal(t2-t1,'hour'))
        select case ( dhour )
          case ( 0 )
            t_spm = t1
          case ( 3 )
            call Get( t1, hour=hour )
            if ( modulo(hour,6) == 0 ) then
              t_spm = t1
            else
              t_spm = t2
            end if
          case ( 6 )
            t_spm = t1 + IncrDate( hour=3 )
          case default
            write (gol,'("do not know how to form time for mid surface pressure:")'); call goErr
            call wrtgol( '  t1 : ', t1 ); call goErr
            call wrtgol( '  t2 : ', t2 ); call goErr
            TRACEBACK; status = 1; return
        end select

        ! replace year ?
        if ( mf%fixyear > 0 ) then
          call Set( t_spm, year=mf%fixyear )
        end if

        ! initial no record is found ...
        status = 1

        ! number of data sets:
        if ( mf%spm_extr ) then
          call GetInfo( mf%spm_hdf, status, num_datasets=nsds )
        else
          call GetInfo( mf%hdf, status, num_datasets=nsds )
        end if
        if (status/=0) then; TRACEBACK; return; end if

        ! loop over all data sets:
        do isds = 1, nsds

          ! select data set:
          if ( mf%spm_extr ) then
            call Select( mf%spm_sds, mf%spm_hdf, isds-1, status )
          else
            call Select( mf%spm_sds, mf%hdf, isds-1, status )
          end if
          if (status/=0) then; TRACEBACK; return; end if

          ! correct param ?
          status=-1; call CheckInfo( mf%spm_sds, status, name='ps' )
          if (status>0) then; TRACEBACK; status=1; return; end if

          ! not found ? also try other names:
          if ( status < 0) then
            select case ( paramkey )
              case ( 'mfu', 'pu' )
                status=-1; call CheckInfo( mf%spm_sds, status, name='spu' )
              case ( 'mfv', 'pv' )
                status=-1; call CheckInfo( mf%spm_sds, status, name='spv' )
              case default
                status=-1; call CheckInfo( mf%spm_sds, status, name='sp' )
            end select
            if (status>0) then; TRACEBACK; status=1; return; end if
            ! try next ?
            if ( status < 0 ) cycle
          end if

          ! correct time ?
          call Get( t1, time6=time1 )
          status=-1; call CheckAttribute( mf%spm_sds, 'time1', time1, status )
          if (status>0) then; TRACEBACK; status=1; return; end if
          if (status==0) then
            call Get( t2, time6=time1 )
            status=-1; call CheckAttribute( mf%spm_sds, 'time2', time2, status )
            if (status>0) then; TRACEBACK; status=1; return; end if
          else
            ! try special spm times
            call Get( t_spm, time6=time1 )
            status=-1; call CheckAttribute( mf%spm_sds, 'time1', time1, status )
            if (status<0) cycle
            if (status>0) then; TRACEBACK; status=1; return; end if
            ! try special spm times
            call Get( t_spm, time6=time2 )
            status=-1; call CheckAttribute( mf%spm_sds, 'time2', time2, status )
            if (status<0) cycle
            if (status>0) then; TRACEBACK; status=1; return; end if
          end if

          ! found!
          exit
        end do

        ! not found ?
        if ( status /= 0 ) then
          write (gol,'("Unable to locate surface pressure field in hdf file:")'); call goErr
          call wrtgol( '  t1 - t2 : ', t1, ' - ', t2 ); call goErr
          call wrtgol( '  t_spm   : ', t_spm ); call goErr
          write (gol,'("  file  : ",a)') trim(mf%spm_fname); call goErr
          TRACEBACK; status = 1; return
        end if

        ! ok
        mf%spm_selected = .true.

      end if

    end if

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! ok
    status = 0

  end subroutine mf_SelectRecord


  ! ***


  !
  ! initialiase grid info from sds
  !

  subroutine lli_Init_mf( lli, nuv, mf, status )

    use file_hdf, only : ReadAttribute, CheckAttribute
    use Grid, only : TllGridInfo, Init

    ! --- in/out ----------------------------------

    type(TllGridInfo), intent(out)       ::  lli
    character(len=1), intent(in)         ::  nuv
    type(TMeteoFile_tmpp), intent(in)    ::  mf
    integer, intent(out)                 ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/lli_Init_mf'

    ! --- local -----------------------------------

    real           ::  lon_deg, dlon_deg
    integer        ::  nlon
    real           ::  lat_deg, dlat_deg
    integer        ::  nlat
    integer        ::  stat

    ! --- begin ------------------------------------

    ! record not selected ? then return
    if ( .not. mf%selected ) then
      write (gol,'("no record selected ...")'); call goErr
      TRACEBACK; status = 1; return
    end if

    ! check ...
    call CheckAttribute( mf%sds, 'gridtype' , 'll', status )
    if ( status /= 0 ) then
      write (gol,'("sds does not seem to contain ll grid")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! extract grid position parameters from sds:
    call ReadAttribute( mf%sds, 'lon_first', lon_deg , status )
    if (status/=0) then; TRACEBACK; return; end if
    call ReadAttribute( mf%sds, 'lon_inc  ', dlon_deg, status )
    if (status/=0) then; TRACEBACK; return; end if
    call ReadAttribute( mf%sds, 'lon_n    ', nlon    , status )
    if (status/=0) then; TRACEBACK; return; end if
    call ReadAttribute( mf%sds, 'lat_first', lat_deg , status )
    if (status/=0) then; TRACEBACK; return; end if
    call ReadAttribute( mf%sds, 'lat_inc  ', dlat_deg, status )
    if (status/=0) then; TRACEBACK; return; end if
    call ReadAttribute( mf%sds, 'lat_n    ', nlat    , status )
    if (status/=0) then; TRACEBACK; return; end if

    ! the just read values might define points on the cell boundaries,
    ! while lli should define the cell centers;
    ! fill correct values using nuv :
    select case ( nuv )
      case ( 'n' )
        ! sds attributes defined for lat/lon of center
        call Init( lli, lon_deg, dlon_deg, nlon, &
                        lat_deg, dlat_deg, nlat, status )
        if (status/=0) then; TRACEBACK; return; end if
      case ( 'u' )
        ! sds attributes defined for lat/lon of east/west bound
        call Init( lli, lon_deg+0.5*dlon_deg, dlon_deg, nlon-1, &
                        lat_deg             , dlat_deg, nlat  , status )
        if (status/=0) then; TRACEBACK; return; end if
      case ( 'v' )
        ! sds attributes defined for lat/lon of south/north bound
        call Init( lli, lon_deg             , dlon_deg, nlon  , &
                        lat_deg+0.5*dlat_deg, dlat_deg, nlat-1, status )
        if (status/=0) then; TRACEBACK; return; end if
      case default
        write (gol,'("unsupported nuv `",a,"`")') nuv; call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine lli_Init_mf


  ! ***


  !
  ! initialiase level info from sds
  !

  subroutine levi_Init_mf( levi, mf, status )

    use file_hdf, only : GetInfo, ReadAttribute
    use Grid, only : TLevelInfo, Init

    ! --- in/out ----------------------------------

    type(TLevelInfo), intent(out)        ::  levi
    type(TMeteoFile_tmpp), intent(in)    ::  mf
    integer, intent(out)                 ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/levi_Init_mf'

    ! --- local -----------------------------------

    integer                ::  data_rank
    integer                ::  lm
    real, allocatable      ::  a(:), b(:)

    ! --- begin ------------------------------------

    ! record not selected ? then return
    if ( .not. mf%selected ) then
      write (gol,'("no record selected ...")'); call goErr
      TRACEBACK; status = 1; return
    end if

    ! extract size of data array:
    call GetInfo( mf%sds, status, data_rank=data_rank )
    if ( status/= 0 ) then; TRACEBACK; status=1; return; end if

    ! 2D or 3D
    select case ( data_rank )

      case ( 2 )

        ! set dummy values ...
        call Init( levi, 1, (/0.0,0.0/), (/0.0,0.0/), status )
        if ( status/= 0 ) then; TRACEBACK; status=1; return; end if

      case ( 3 )

        ! read number of levels and hybride parameters:
        call ReadAttribute( mf%sds, 'lm', lm, status )

        ! extract hybride coeff
        allocate( a(lm+1), b(lm+1) )
        call ReadAttribute( mf%sds, 'at', a, status )
        if (status/=0) then; TRACEBACK; return; end if
        call ReadAttribute( mf%sds, 'bt', b, status )
        if (status/=0) then; TRACEBACK; return; end if

        ! fill ...
        call Init( levi, lm, a, b, status )
        if ( status/= 0 ) then; TRACEBACK; status=1; return; end if

        ! done
        deallocate( a, b )

      case default

        write (gol,'("unsupported data rank : ",i6)') data_rank; call goErr
        TRACEBACK; status=1; return

    end select

    ! ok
    status = 0

  end subroutine levi_Init_mf


  ! ***


  subroutine mf_GetField( mf, status, gridtype, ll, spm )

    use PArray, only : pa_SetShape
    use file_hdf, only : GetInfo, ReadData, ReadAttribute

    ! --- in/out --------------------------------

    type(TMeteoFile_tmpp), intent(inout)      ::  mf
    integer, intent(out)                      ::  status

    character(len=*), intent(out), optional   ::  gridtype
    real, pointer, optional                   ::  ll(:,:,:)
    real, pointer, optional                   ::  spm(:,:)


    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_GetField'

    ! --- local -------------------------------

    integer            ::  data_rank
    integer            ::  data_dims(3)

    ! --- begin -------------------------------

    ! record not selected ? then return
    if ( .not. mf%selected ) then
      write (gol,'("no record selected ...")'); call goErr
      TRACEBACK; status = 1; return
    end if

    ! initial data is not extracted ...
    status = 1

    ! return grid type ?
    if ( present(gridtype) ) then
      call ReadAttribute( mf%sds, 'gridtype', gridtype, status )
      if (status/=0) then; TRACEBACK; return; end if
    end if

    ! return 3d data array ?
    if ( present(ll) ) then
      ! extract data rank and shape:
      call GetInfo( mf%sds, status, data_rank=data_rank, data_dims=data_dims )
      if (status/=0) then; TRACEBACK; return; end if
      ! extract data array:
      select case ( data_rank )
        case ( 2 )
          data_dims(3) = 1
          call pa_SetShape( ll, data_dims )
          call ReadData( mf%sds, ll(:,:,1), status )
          if (status/=0) then; TRACEBACK; return; end if
        case ( 3 )
          call pa_SetShape( ll, data_dims )
          call ReadData( mf%sds, ll, status )
          if (status/=0) then; TRACEBACK; return; end if
        case default
          write (gol,'("unsupported data rank:",i6)') data_rank; call goErr
          TRACEBACK; status=1; return
      end select
      ! ReadData breaks on error ...
      status = 0
    end if

    ! return surface pressure array ?
    if ( present(spm) ) then
      if ( .not. mf%spm_selected ) then
        write (gol,'("no spm record selected ...")'); call goErr
        TRACEBACK; status = 1; return
      end if
      ! extract data rank and shape:
      call GetInfo( mf%spm_sds, status, data_rank=data_rank, data_dims=data_dims )
      if (status/=0) then; TRACEBACK; return; end if
      ! extract data array:
      select case ( data_rank )
        case ( 2 )
          call pa_SetShape( spm, data_dims(1:2) )
          call ReadData( mf%spm_sds, spm, status )
          if (status/=0) then; TRACEBACK; return; end if
        case default
          write (gol,'("unsupported data rank for spm:",i6)') data_rank; call goErr
          TRACEBACK; status=1; return
      end select
      ! ReadData breaks on error ...
      status = 0
    end if

    ! no arguments processed ?
    ! probably subroutine is called incorrectly ...
    if ( status /= 0 ) then
      write (gol,'("no arguments processed; wrong call ?")'); call goErr
      TRACEBACK; status = 1; return
    end if

    ! ok
    status = 0

  end subroutine mf_GetField



  ! ***


  subroutine mf_ReadRecord( mf, paramkey, t1, t2, nuv, nw, &
                                gridtype, levi, &
                                lli, ll, sp_ll, &
                                status )

    use PArray, only : pa_Done
    use GO, only : TDate
    use Grid, only : TllGridInfo, TLevelInfo
    use file_hdf_base, only : Init, Done   !needed on aster

    ! --- in/out -------------------------------

    type(TMeteoFile_tmpp), intent(inout) ::  mf
    character(len=*), intent(in)         ::  paramkey
    type(TDate), intent(in)              ::  t1, t2
    character(len=1), intent(in)         ::  nuv
    character(len=2), intent(out)        ::  gridtype
    type(TLevelInfo), intent(out)        ::  levi
    character(len=1), intent(in)         ::  nw
    type(TllGridInfo), intent(inout)     ::  lli
    real, pointer                        ::  ll(:,:,:)
    real, pointer                        ::  sp_ll(:,:)
    integer, intent(out)                 ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_ReadRecord'

    ! --- local -------------------------------

    logical               ::  exist

    real, allocatable    ::  mfw_n(:,:,:)

    real, allocatable     ::  sp_ll_n(:,:)
    integer               ::  i, j

    ! --- begin ---------------------------------

    ! input ?
    if ( mf%io /= 'i' ) then
      write (gol,'("file should have been opened for input, but io=",a)') mf%io; call goErr
      TRACEBACK; status=1; return
    end if

    ! open for reading:
    ! if it is opened in Init, to many hdf files would be open:
    !
    inquire( file=trim(mf%fname), exist=exist )
    if ( .not. exist ) then
      write (gol,'("main file does not exist:")'); call goErr
      write (gol,'("  ",a)') trim(mf%fname); call goErr
      TRACEBACK; status=1; return
    end if
    !
    call Init( mf%hdf, trim(mf%fname), 'read', status )
    if (status/=0) then; TRACEBACK; return; end if
    !
    call Init( mf%sds, status )
    if (status/=0) then; TRACEBACK; return; end if
    !
    mf%selected = .false.

    ! open spm file if necessary;
    ! if it is opened in Init, to many hdf files would be open:
    if ( mf%spm_load ) then
      !
      if ( mf%spm_extr ) then
        !
        inquire( file=trim(mf%spm_fname), exist=exist )
        if ( .not. exist ) then
          write (gol,'("spm file does not exist:")'); call goErr
          write (gol,'("  ",a)') trim(mf%spm_fname); call goErr
          TRACEBACK; status=1; return
        end if
        !
        call Init( mf%spm_hdf, trim(mf%spm_fname), 'read', status )
        if (status/=0) then; TRACEBACK; return; end if
        !
      end if
      !
      call Init( mf%spm_sds, status )
      if (status/=0) then; TRACEBACK; return; end if
      !
      mf%spm_selected = .false.
    end if

    ! select records in hdf files:
    call mf_SelectRecord( mf, paramkey, t1, t2, status )
    if ( status/= 0 ) then; TRACEBACK; status=1; return; end if

    ! always regular lat/lon grid ..
    gridtype = 'll'

    ! setup grid definition:
    call lli_Init_mf( lli, nuv, mf, status )
    if ( status/= 0 ) then; TRACEBACK; status=1; return; end if

    ! setup level definition:
    call levi_Init_mf( levi, mf, status )
    if ( status/= 0 ) then; TRACEBACK; status=1; return; end if

    ! fill data array:
    call mf_GetField( mf, status, ll=ll )
    if ( status/= 0 ) then; TRACEBACK; status=1; return; end if

    ! ***

    ! special treatment of pw fields ?
    if ( (paramkey == 'mfw') .or. (paramkey == 'pw') ) then

      ! add extra top level ?
      if ( size(ll,3) == levi%nlev ) then
        ! copy current pw ('n' levels)
        allocate( mfw_n(lli%im,lli%jm,levi%nlev) )
        mfw_n = ll
        ! renew target array
        deallocate( ll )
        allocate( ll(lli%nlon,lli%nlat,levi%nlev+1) )
        ! store new pw:
        ll(:,:,1:levi%nlev) = mfw_n
        ll(:,:,levi%nlev+1) = 0.0   ! kg/s
        ! clear
        deallocate( mfw_n )
      end if

      ! change flux direction ?
      if ( mf%mfw_redir ) ll = - ll     ! upwards (increasing level)

    end if

    ! ***

    ! special treatment of surface stress fields ?

    if ( ((paramkey=='ewss') .or. (paramkey=='nsss')) .and. mf%sstr_to_ewss_nsss ) then

      !   sstr**2  =  (sstr**2)/2 + (sstr**2)/2  =  nsss**2 + ewss**2
      !  thus  nsss = sqrt( (sstr**2)/2.0 )

      ll = sqrt( (ll**2) / 2.0 )

    end if

    ! ***

    ! fill surface pressure field ?
    if ( mf%spm_load .and. mf%spm_selected ) then

      ! read field:
      call mf_GetField( mf, status, spm=sp_ll )
      if ( status/= 0 ) then; TRACEBACK; status=1; return; end if

      ! convert from 'n' to 'u'/'v' ?
      if ( mf%spm_n_to_uv ) then
        select case ( nuv )
          case ( 'n' )
            ! no interpolation
          case ( 'u' )
            ! copy current sp ('n')
            allocate( sp_ll_n(lli%im,lli%jm) )
            sp_ll_n = sp_ll
            ! renew size of output array:
            deallocate( sp_ll )
            allocate( sp_ll(0:lli%im,lli%jm) )
            ! interpol:
            sp_ll(0,:) = 1.5*sp_ll_n(1,:) - 0.5*sp_ll_n(2,:)
            do i = 1, lli%im-1
              sp_ll(i,:) = 1.5*sp_ll_n(i,:) + 0.5*sp_ll_n(i+1,:)
            end do
            sp_ll(lli%im,:) = -0.5*sp_ll_n(lli%im-1,:) + 1.5*sp_ll_n(lli%im,:)
            ! clear
            deallocate( sp_ll_n )
          case ( 'v' )
            ! copy current sp ('n')
            allocate( sp_ll_n(lli%im,lli%jm) )
            sp_ll_n = sp_ll
            ! renew size of output array:
            deallocate( sp_ll )
            allocate( sp_ll(lli%im,0:lli%jm) )
            ! interpol:
            sp_ll(:,0) = 1.5*sp_ll_n(:,1) - 0.5*sp_ll_n(:,2)
            do j = 1, lli%jm-1
              sp_ll(:,j) = 1.5*sp_ll_n(:,j) + 0.5*sp_ll_n(:,j+1)
            end do
            sp_ll(:,lli%jm) = -0.5*sp_ll_n(:,lli%jm-1) + 1.5*sp_ll_n(:,lli%jm)
            ! clear
            deallocate( sp_ll_n )
          case default
            write (gol,'("unsupported nuv `",a,"`")') nuv; call goErr
            TRACEBACK; status=1; return
        end select
      end if

    else
      call pa_Done( sp_ll )
    end if

    ! close surface pressure file ?
    if ( mf%spm_load ) then
      call Done( mf%spm_sds, status )
      if (status/=0) then; TRACEBACK; return; end if
      !
      if ( mf%spm_extr ) then
        call Done( mf%spm_hdf, status )
        if (status/=0) then; TRACEBACK; return; end if
      end if
    end if

    ! close field file:
    call Done( mf%sds, status )
    if (status/=0) then; TRACEBACK; return; end if
    call Done( mf%hdf, status )
    if (status/=0) then; TRACEBACK; return; end if

    ! ok
    status = 0

  end subroutine mf_ReadRecord


  ! ******************************************************************
  ! ***
  ! *** output
  ! ***
  ! ******************************************************************


  subroutine WriteHeader( mf, lli, status, levi )

    use Binas, only : grav, ae
    use Grid, only : TllGridInfo, TLevelInfo
    use file_hdf, only : WriteAttribute

    ! --- in/out -------------------------------

    type(TMeteoFile_tmpp), intent(inout)     ::  mf
    type(TllGridInfo), intent(in)            ::  lli
    integer, intent(out)                     ::  status

    type(TLevelInfo), intent(in), optional   ::  levi

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/WriteHeader'

    ! --- begin ---------------------------------

    ! write header:
    call WriteAttribute( mf%hdf, 'fname'   , trim(mf%fname), status )
    if (status/=0) then; TRACEBACK; return; end if
    call WriteAttribute( mf%hdf, 'format'  , output_format, status )
    if (status/=0) then; TRACEBACK; return; end if
    call WriteAttribute( mf%hdf, 'gridtype', 'll'         , status )
    if (status/=0) then; TRACEBACK; return; end if

    ! save first and last lon/lat (center) for use with HIPHOP
    call WriteAttribute( mf%hdf, 'lonmin', lli%lon_deg(1)       , status, knd=rknd )
    if (status/=0) then; TRACEBACK; return; end if
    call WriteAttribute( mf%hdf, 'lonmax', lli%lon_deg(lli%nlon), status, knd=rknd )
    if (status/=0) then; TRACEBACK; return; end if
    call WriteAttribute( mf%hdf, 'latmin', lli%lat_deg(1)       , status, knd=rknd )
    if (status/=0) then; TRACEBACK; return; end if
    call WriteAttribute( mf%hdf, 'latmax', lli%lat_deg(lli%nlat), status, knd=rknd )
    if (status/=0) then; TRACEBACK; return; end if

    ! other useful stuff ...
    call WriteAttribute( mf%hdf, 'grav'   , grav       , status, knd=rknd )
    if (status/=0) then; TRACEBACK; return; end if
    call WriteAttribute( mf%hdf, 'ae'     , ae         , status, knd=rknd )
    if (status/=0) then; TRACEBACK; return; end if
    call WriteAttribute( mf%hdf, 'area_m2', lli%area_m2, status, knd=rknd )
    if (status/=0) then; TRACEBACK; return; end if

    ! level stuff
    if ( present(levi) ) then
      call WriteAttribute( mf%hdf, 'lm', levi%nlev, status, knd=iknd )
      if (status/=0) then; TRACEBACK; return; end if
      call WriteAttribute( mf%hdf, 'at', levi%a   , status, knd=rknd )
      if (status/=0) then; TRACEBACK; return; end if
      call WriteAttribute( mf%hdf, 'bt', levi%b   , status, knd=rknd )
      if (status/=0) then; TRACEBACK; return; end if
    end if

    ! ok
    status = 0

  end subroutine WriteHeader


  ! ***


  subroutine WriteSdsHeader( sds, tmi, unit, tref, t1, t2, lli, nuv, status, &
                                  levi, nw, nlev )

    use Binas   , only : p_global
    use file_hdf, only : TSds
    use file_hdf, only : SetDim, WriteAttribute, Compress
    use GO      , only : TDate, Get, rTotal, operator(-), IsAnyDate
    use Grid    , only : TllGridInfo, TLevelInfo
    use tmm_info, only : TMeteoInfo

    ! --- in/out -----------------------------

    type(TSds), intent(inout)               ::  sds
    type(TMeteoInfo), intent(in)            ::  tmi
    character(len=*), intent(in)            ::  unit
    type(TDate), intent(in)                 ::  tref, t1, t2
    type(TllGridInfo), intent(in)           ::  lli
    character(len=1), intent(in)            ::  nuv
    integer, intent(out)                    ::  status
    type(TLevelInfo), intent(in), optional  ::  levi
    character(len=1), intent(in), optional  ::  nw
    integer, intent(in), optional           ::  nlev

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/WriteSdsHeader'

    ! --- local -----------------------------

    integer             ::  time6(6)
    integer             ::  dhour

    ! --- begin -----------------------------

    ! *** history

    ! write history of meteo field:
    if ( len_trim(tmi%history) < 1 ) then
      call WriteAttribute( sds, 'history', '-', status )
      IF_NOTOK_RETURN(status=1)
    else
      call WriteAttribute( sds, 'history', trim(tmi%history), status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! *** unit

    ! write unit attribute
    call WriteAttribute( sds, 'unit', unit, status )
    IF_NOTOK_RETURN(status=1)

    ! *** time

    ! time interval in hours
    if ( IsAnyDate(t1) .or. IsAnyDate(t2) ) then
      dhour = 0
    else
      dhour = nint(rTotal( t2 - t2, 'hour' ))
    end if

    ! write old time attributes
    call Get( t1, time6=time6 )
    call WriteAttribute( sds, 'idate', time6, status, knd=iknd )
    IF_NOTOK_RETURN(status=1)
    call WriteAttribute( sds, 'dthrs', dhour, status, knd=iknd )  ! hours
    IF_NOTOK_RETURN(status=1)

    ! write new time attributes
    call Get( tref, time6=time6 )
      call WriteAttribute( sds, 'tref', time6, status, knd=iknd )
      IF_NOTOK_RETURN(status=1)
    call Get( t1, time6=time6 )
      call WriteAttribute( sds, 'time1', time6, status, knd=iknd )
      IF_NOTOK_RETURN(status=1)
    call Get( t2, time6=time6 )
      call WriteAttribute( sds, 'time2', time6, status, knd=iknd )
      IF_NOTOK_RETURN(status=1)
    time6 = 0
    time6(4) = dhour
    call WriteAttribute( sds, 'dtime', time6, status, knd=iknd )
    IF_NOTOK_RETURN(status=1)

    ! *** grid

    select case ( nuv )
      case ( 'n' )
        ! sds dimensions:
        call SetDim( sds, 1-1, 'LON'   , 'deg', lli%lon_deg, status, knd=rknd )
        IF_NOTOK_RETURN(status=1)
        call SetDim( sds, 2-1, 'LAT'   , 'deg', lli%lat_deg, status, knd=rknd )
        IF_NOTOK_RETURN(status=1)
        ! sds attributes:
        call WriteAttribute( sds, 'gridtype' , 'll'          , status )
        IF_NOTOK_RETURN(status=1)
        call WriteAttribute( sds, 'lon_first', lli%lon_deg(1), status, knd=rknd )
        IF_NOTOK_RETURN(status=1)
        call WriteAttribute( sds, 'lon_inc  ', lli%dlon_deg  , status, knd=rknd )
        IF_NOTOK_RETURN(status=1)
        call WriteAttribute( sds, 'lon_n    ', lli%nlon      , status, knd=iknd )
        IF_NOTOK_RETURN(status=1)
        call WriteAttribute( sds, 'lat_first', lli%lat_deg(1), status, knd=rknd )
        IF_NOTOK_RETURN(status=1)
        call WriteAttribute( sds, 'lat_inc  ', lli%dlat_deg  , status, knd=rknd )
        IF_NOTOK_RETURN(status=1)
        call WriteAttribute( sds, 'lat_n    ', lli%nlat      , status, knd=iknd )
        IF_NOTOK_RETURN(status=1)
      case ( 'u' )
        ! sds dimensions:
        call SetDim( sds, 1-1, 'LONP1' , 'deg', lli%blon_deg, status, knd=rknd )
        IF_NOTOK_RETURN(status=1)
        call SetDim( sds, 2-1, 'LAT'   , 'deg', lli%lat_deg , status, knd=rknd )
        IF_NOTOK_RETURN(status=1)
        ! sds attributes:
        call WriteAttribute( sds, 'gridtype' , 'll', status )
        IF_NOTOK_RETURN(status=1)
        call WriteAttribute( sds, 'lon_first', lli%blon_deg(0), status, knd=rknd )
        IF_NOTOK_RETURN(status=1)
        call WriteAttribute( sds, 'lon_inc  ', lli%dlon_deg   , status, knd=rknd )
        IF_NOTOK_RETURN(status=1)
        call WriteAttribute( sds, 'lon_n    ', lli%im+1       , status, knd=iknd )
        IF_NOTOK_RETURN(status=1)
        call WriteAttribute( sds, 'lat_first', lli%lat_deg(1) , status, knd=rknd )
        IF_NOTOK_RETURN(status=1)
        call WriteAttribute( sds, 'lat_inc  ', lli%dlat_deg   , status, knd=rknd )
        IF_NOTOK_RETURN(status=1)
        call WriteAttribute( sds, 'lat_n    ', lli%nlat       , status, knd=iknd )
        IF_NOTOK_RETURN(status=1)
      case ( 'v' )
        ! sds dimensions:
        call SetDim( sds, 1-1, 'LON'   , 'deg', lli%lon_deg , status, knd=rknd )
        IF_NOTOK_RETURN(status=1)
        call SetDim( sds, 2-1, 'LATP1' , 'deg', lli%blat_deg, status, knd=rknd )
        IF_NOTOK_RETURN(status=1)
        ! sds attributes:
        call WriteAttribute( sds, 'gridtype' , 'll', status )
        call WriteAttribute( sds, 'lon_first', lli%lon_deg(1) , status, knd=rknd )
        IF_NOTOK_RETURN(status=1)
        call WriteAttribute( sds, 'lon_inc  ', lli%dlon_deg   , status, knd=rknd )
        IF_NOTOK_RETURN(status=1)
        call WriteAttribute( sds, 'lon_n    ', lli%nlon       , status, knd=iknd )
        IF_NOTOK_RETURN(status=1)
        call WriteAttribute( sds, 'lat_first', lli%blat_deg(0), status, knd=rknd )
        IF_NOTOK_RETURN(status=1)
        call WriteAttribute( sds, 'lat_inc  ', lli%dlat_deg   , status, knd=rknd )
        IF_NOTOK_RETURN(status=1)
        call WriteAttribute( sds, 'lat_n    ', lli%jm+1       , status, knd=iknd )
        IF_NOTOK_RETURN(status=1)
      case default
        write (gol,'("unsupported nuv `",a,"`")') nuv; call goErr
        TRACEBACK; status=1; return
    end select

    ! *** levels

    if ( present(levi) ) then
      if ( .not. present(nw) ) then
        write (gol,'("optional levi requires nw")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! sds dimensions:
      select case ( nw )
        case ( '*' )
          if ( present(nlev) ) then
            call SetDim( sds, 3-1, 'HYBRID_SELECTED', 'Pa' , levi%fp0(1:nlev), status, knd=rknd )
            IF_NOTOK_RETURN(status=1)
          end if
        case ( 'n' )
          call SetDim( sds, 3-1, 'HYBRID', 'Pa' , levi%fp0, status, knd=rknd )
          IF_NOTOK_RETURN(status=1)
        case ( 'w' )
          call SetDim( sds, 3-1, 'HYBRIDh', 'Pa', levi%p0, status, knd=rknd )
          IF_NOTOK_RETURN(status=1)
        case default
          write (gol,'("unsupported nw `",a,"` xxx")') nw; call goErr
          TRACEBACK; status=1; return
      end select
      ! sds attributes:
      call WriteAttribute( sds, 'lm', levi%nlev, status, knd=iknd )
      IF_NOTOK_RETURN(status=1)
      call WriteAttribute( sds, 'at', levi%a   , status, knd=rknd )
      IF_NOTOK_RETURN(status=1)
      call WriteAttribute( sds, 'bt', levi%b   , status, knd=rknd )
      IF_NOTOK_RETURN(status=1)
    end if

    ! *** data compression

    call Compress( sds, compression, status, deflate_level=deflate_level )
    IF_NOTOK_RETURN(status=1)

    ! *** end

    ! ok
    status = 0

  end subroutine WriteSdsHeader


  ! ***


  subroutine WriteStatus( mf, msg, status )

    ! --- in/out -------------------------------

    type(TMeteoFile_tmpp), intent(inout) ::  mf
    character(len=*), intent(in)         ::  msg
    integer, intent(out)                 ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/WriteStatus'

    ! --- local ------------------------------

    integer       ::  fu
    logical       ::  opened

    ! --- begin ---------------------------------

    ! select unused file unit:
    fu = 1234
    do
      inquire( unit=fu, opened=opened )
      if ( .not. opened ) exit
      fu = fu + 1
    end do

    ! open:
    open( fu, file=trim(mf%fname)//'.status', form='formatted', iostat=status )
    if (status/=0) then
      write (gol,'("opening status file:")'); call goErr
      write (gol,'("  file   : ",a)') trim(mf%fname); call goErr
      TRACEBACK; status=1; return
    end if

    ! write message:
    write (fu,'(a)',iostat=status) msg
    if (status/=0) then
      write (gol,'("writing status:")'); call goErr
      write (gol,'("  file   : ",a)') trim(mf%fname); call goErr
      write (gol,'("  msg    : ",a)') msg; call goErr
      TRACEBACK; status=1; return
    end if

    ! done:
    close( fu, iostat=status )
    if (status/=0) then
      write (gol,'("closing status file:")'); call goErr
      write (gol,'("  file   : ",a)') trim(mf%fname); call goErr
      TRACEBACK; status=1; return
    end if

    ! ok
    status = 0

  end subroutine WriteStatus


  ! ***


  subroutine AddOutputName( mf, name, status )

    ! --- in/out -------------------------------

    type(TMeteoFile_tmpp), intent(inout) ::  mf
    character(len=*), intent(in)         ::  name
    integer, intent(out)                 ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/AddOuputName'

    ! --- local ------------------------------

    integer       ::  iname

    ! --- begin ---------------------------------

    ! if not present yet, add name
    iname = 1
    do
      ! add ?
      if ( iname > mf%output_nname ) then
        ! place to store ?
        if ( iname > size(mf%output_names) ) then
          write (gol,'("length of mf%output_names array too small:")'); call goErr
          do iname = 1, size(mf%output_names)
            write (gol,'("  ",i3," ",a)') iname, trim(mf%output_names(iname)); call goErr
          end do
          TRACEBACK; status=1; return
        end if
        ! long enough for name ?
        if ( len(mf%output_names(iname)) < len(name) ) then
          write (gol,'("length of mf%output_names too small:")'); call goErr
          write (gol,'("  len(mf%output_names(i)) : ",i4)') len(mf%output_names(iname)); call goErr
          write (gol,'("  len(mname)              : ",i4)') len(name); call goErr
          TRACEBACK; status=1; return
        end if
        ! store:
        mf%output_names(iname) = name
        ! increase counter:
        mf%output_nname = iname
        ! leave:
        exit
      end if
      ! found ? then leave loop
      if ( mf%output_names(iname) == name ) exit
      ! next index
      iname = iname + 1
    end do

    ! ok
    status = 0

  end subroutine AddOutputName


  ! ***


  subroutine mf_WriteRecord_2d( mf, tmi, paramkey, unit, tref, t1, t2, &
                                lli, nuv, ll, status )

    use GO      , only : TDate
    use Grid    , only : TllGridInfo
    use file_hdf, only : TSds, WriteData
    use file_hdf_base, only :  Init, Done   ! needed on aster
    use tmm_info, only : TMeteoInfo

    ! --- in/out -------------------------------

    type(TMeteoFile_tmpp), intent(inout) ::  mf
    type(TMeteoInfo), intent(in)         ::  tmi
    character(len=*), intent(in)         ::  paramkey, unit
    type(TDate), intent(in)              ::  tref, t1, t2
    type(TllGridInfo), intent(in)        ::  lli
    character(len=1), intent(in)         ::  nuv
    real, intent(in)                     ::  ll(:,:)
    integer, intent(out)                 ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_WriteRecord_2d'

    ! --- local ------------------------------

    type(TSds)    ::  sds

    ! --- begin ---------------------------------

    ! output ?
    if ( mf%io /= 'o' ) then
      write (gol,'("file should have been opened for output, but io=",a)') mf%io; call goErr
      TRACEBACK; status=1; return
    end if

    ! new or existing ?
    if ( .not. mf%output_initialised ) then
      ! open new file, destroy old:
      call Init( mf%hdf, trim(mf%fname), 'create', status )
      IF_NOTOK_RETURN(status=1)
      ! write file header:
      call WriteHeader( mf, lli, status )
      IF_NOTOK_RETURN(status=1)
      ! status new
      call WriteStatus( mf, 'in-progress', status )
      IF_NOTOK_RETURN(status=1)
      ! no records written yet:
      mf%output_nrec = 0
      ! now the file is initialised
      mf%output_initialised = .true.
    else
      ! re-open file:
      call Init( mf%hdf, trim(mf%fname), 'write', status )
      IF_NOTOK_RETURN(status=1)
    endif


    ! *** data set

    ! add record name:
    call AddOutputName( mf, paramkey, status )
    IF_NOTOK_RETURN(status=1)

    ! init data set:
    call Init( sds, mf%hdf, paramkey, shape(ll), 'real', status, knd=rknd_ds )
    IF_NOTOK_RETURN(status=1)

    ! write unit, time, and grid info:
    call WriteSdsHeader( sds, tmi, unit, tref, t1, t2, lli, nuv, status )
    IF_NOTOK_RETURN(status=1)

    ! write grid
    call WriteData( sds, ll, status )
    IF_NOTOK_RETURN(status=1)

    ! done:
    call Done( sds, status )
    IF_NOTOK_RETURN(status=1)

    ! *** completed ?

    ! next record has been written:
    mf%output_nrec = mf%output_nrec + 1

    ! completed ? then re-write status file:
    if ( mf%output_nrec == mf%output_ntrec*mf%output_nname ) then
      call WriteStatus( mf, 'completed', status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! *** close

    call Done( mf%hdf, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine mf_WriteRecord_2d


  ! ***


  subroutine mf_WriteRecord_3d( mf, tmi, spname, paramkey, unit, tref, t1, t2, &
                                lli, nuv, levi, nw, ps, ll, status )

    use GO      , only : TDate
    use Grid    , only : TllGridInfo, TLevelInfo
    use file_hdf, only : TSds, WriteData
    use file_hdf_base, only : Init, Done  !needed on aster
    use tmm_info, only : TMeteoInfo

    ! --- in/out -------------------------------

    type(TMeteoFile_tmpp), intent(inout) ::  mf
    type(TMeteoInfo), intent(in)         ::  tmi
    character(len=*), intent(in)         ::  spname, paramkey, unit
    type(TDate), intent(in)              ::  tref, t1, t2
    type(TllGridInfo), intent(in)        ::  lli
    character(len=1), intent(in)         ::  nuv
    type(TLevelInfo), intent(in)         ::  levi
    character(len=1), intent(in)         ::  nw
    real, intent(in)                     ::  ps(:,:)
    real, intent(in)                     ::  ll(:,:,:)
    integer, intent(out)                 ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/mf_WriteRecord_3d'

    ! --- local ------------------------------

    type(TSds)    ::  sds
    integer       ::  iname

    ! --- begin ---------------------------------

    ! output ?
    if ( mf%io /= 'o' ) then
      write (gol,'("file should have been opened for output, but io=",a)') mf%io; call goErr
      TRACEBACK; status=1; return
    end if

    ! new or existing ?
    if ( .not. mf%output_initialised ) then
      ! open new file, destroy old:
      call Init( mf%hdf, trim(mf%fname), 'create', status )
      IF_NOTOK_RETURN(status=1)
      ! write file header:
      call WriteHeader( mf, lli, status, levi )
      IF_NOTOK_RETURN(status=1)
      ! status new
      call WriteStatus( mf, 'in-progress', status )
      IF_NOTOK_RETURN(status=1)
      ! no records written yet:
      mf%output_nrec = 0
      ! now the file is initialised
      mf%output_initialised = .true.
    else
      ! re-open file:
      call Init( mf%hdf, trim(mf%fname), 'write', status )
      IF_NOTOK_RETURN(status=1)
    endif


    ! *** data set

    ! add record name:
    call AddOutputName( mf, paramkey, status )
    IF_NOTOK_RETURN(status=1)

    ! init data set:
    call Init( sds, mf%hdf, paramkey, shape(ll), 'real', status, knd=rknd_ds )
    IF_NOTOK_RETURN(status=1)

    ! write unit, time, and grid info:
    call WriteSdsHeader( sds, tmi, unit, tref, t1, t2, lli, nuv, status, levi, nw )
    IF_NOTOK_RETURN(status=1)

    ! write grid
    call WriteData( sds, ll, status )
    IF_NOTOK_RETURN(status=1)

    ! done:
    call Done( sds, status )
    IF_NOTOK_RETURN(status=1)

    ! *** surface pressure

    ! init data set:
    call Init( sds, mf%hdf, spname, shape(ps), 'real', status, knd=rknd_ds )
    IF_NOTOK_RETURN(status=1)

    ! write unit, time, and grid info:
    call WriteSdsHeader( sds, tmi, 'Pa', tref, t1, t2, lli, nuv, status )
    IF_NOTOK_RETURN(status=1)

    ! write grid
    call WriteData( sds, ps, status )
    IF_NOTOK_RETURN(status=1)

    ! done:
    call Done( sds, status )
    IF_NOTOK_RETURN(status=1)

    ! *** completed ?

    ! next record has been written:
    mf%output_nrec = mf%output_nrec + 1

    ! completed ? then re-write status file:
    if ( mf%output_nrec == mf%output_ntrec*mf%output_nname ) then
      call WriteStatus( mf, 'completed', status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! *** close

    call Done( mf%hdf, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine mf_WriteRecord_3d


end module tmm_mf_tmpp
