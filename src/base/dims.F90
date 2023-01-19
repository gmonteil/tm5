!###############################################################################
!
! name of dimension file: dims_$proj_$nlevtm.f90
!   $proj    - the project name
!   $nlevtm  - the amount of vertical layers...
!
! version 1.0
!
!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module dims

  use binas, only : pi
  use os_specs, only : DUMMY_STR_LEN, MAX_FILENAME_LEN

  use dims_grid
  use dims_levels
  use global_types

  implicit none


  ! --- in/out -------------------------------------

  public

  ! these should remain private:
  private :: pi
  private :: mname


  ! --- const -------------------------------------

  character(len=*), parameter  ::  mname = 'module Dims'

  ! --- const --------------------------------------

  character(len=MAX_FILENAME_LEN) :: datadir  = ''

  real,parameter          :: gtor = pi/180.0
  real,parameter          :: one  = 1.0
  real,parameter          :: zero = 0.0
  integer,parameter       :: nlat180 = 180
  integer,parameter       :: nlon360 = 360
  real,dimension(nlat180) :: dxy11      ! surface area in 1x1

  real,parameter          :: dlat = dy*gtor
  real,parameter          :: dlon = dx*gtor

  ! extra grid:
  integer, parameter  ::  nregions_all = nregions_max + 1
  integer, parameter  ::  iglbsfc = nregions_max + 1

  ! define the parent child structure:

  integer,dimension(nregions,0:nregions) :: children

  ! children(r,0) = total number of children for region r (zero if childless)
  ! children(r,1) = number of first  child of region r,
  ! children(r,2) = number of second child of region r, etc.

  ! grid-coordinates of each region with respect its parent
  ! (will be calculated at start program)

  integer,dimension(1:nregions) :: ibeg,iend,jbeg,jend,lbeg,lend
  integer,dimension(nregions)   :: isr,ier,jsr,jer  ! scope of the region
  ! depends on NPtouch/SPtouch, Xcyc
  ! calculated in determine_children_etc
  ! _________________________________________________________________________
  !
  ! sequence of steps that is used to process the algorithm.
  ! nsplitsteps contans total steps
  ! n_operators the number of different operators
  ! The routine process_region calls the corresponding routines,
  ! depending on splitorder.
  ! x = x-advection
  ! y = y-advection
  ! z = z-advection
  ! v = vertical mixing SGS
  ! s = sources
  ! c = chemistry....

!  ! number of operators:
!  integer, parameter          ::  n_operators = 6
!  ! number of steps in operator splitting:
!  integer,parameter           ::  nsplitsteps = n_operators * 2
!  ! 1-character codes for processes in operator splitting:
!  character, parameter        ::  splitorder(nsplitsteps)  = &
!       (/'x','y','z','v','s','c','c','s','v','z','y','x'/)

  ! number of operators:
  integer, parameter          ::  n_operators = 7
  ! number of steps in operator splitting:
  integer,parameter           ::  nsplitsteps = n_operators * 2
  ! 1-character codes for processes in operator splitting:
  character, parameter        ::  splitorder(nsplitsteps) = &
       (/'x','y','z','v','d','s','c','c','s','d','v','z','y','x'/)

  ! splitorderzoom contains the fully expanded list of operations
  ! for the children
  character,dimension(nregions,maxref*nsplitsteps) :: splitorderzoom
  ! status keeps track of the operations per region
  integer,dimension(nregions)                      :: status

  integer,parameter                                :: zoom_mode = 1
  ! other modes not implemented...!
  logical,parameter                                :: zoom2D = .true.
  ! no zooming on vertical

  ! _________________________________________________________________________
  ! advection scheme:
#ifdef secmom
  character(5),parameter :: adv_scheme = '2nd_m'
#else
  character(5),parameter :: adv_scheme = 'slope'
#endif

  ! limits the slopes to physical values
  logical                :: limits       = .true.
  logical                :: limits_extra = .true.   ! used in secmom advection

  ! numbers of CFL violations and CFL numbers
  integer,dimension(nregions,3)    :: nxi
  integer,dimension(nregions,3)    :: nloop_max = 0
  real,dimension(nregions,3)       :: xi

  ! small number with respect to the altitude unit (used in the advectz part)
  real,parameter:: epsz=0.0001
  ! _________________________________________________________________________
  ! some timing variables to be used in chemistry applications
!  real                        :: sec_day,sec_month,sec_year
  real, parameter             :: sec_day = 86400.0, sec_normal_year = 31536000.0, sec_leap_year = 31622400.0
!  integer,dimension(12)       :: mlen       !length of the 12 month in days

  ! print debug output? default no:
  logical                     :: okdebug     = .false.
  logical                     :: okdebug_tmm = .false.

  logical                     :: splitsave=.true.     ! 1 file per region
  integer                     :: revert=1    ! if -1 reverses time and winds...
  integer                     :: istart
  integer                     :: ndiag,ninst,ncheck,ntrans
  integer                     :: itau,itaui,itaue,itaut,itau0,nwrite
  integer,dimension(nregions) :: itaur   ! itau count for the different regions

  ! read of 6-hourly fields is staggered with three hours.
  !integer                     :: nread = 6*3600
  !integer,parameter           :: staggered = 3*3600

  ! output frequency:
  integer                     :: nsec_read

  integer                     :: ndyn
  integer                     :: ndyn_max
  integer,dimension(6)        :: idate,idatei,idatee,idatet,idate0
  logical, target             :: newyr,newmonth,newday,cdebug,newsrun,newhour(nregions)
  integer                     :: julday0,iyear0
  integer                     :: icalendo
  integer                     :: ndiagp1,ndiagp2
  integer                     :: nstep,nstep0
  real                        :: cpu0,cpu1
  real,dimension(nregions)    :: areag

  !    main control variables, accessible through namelist 'inputz'
  !    all times (unless noted) are given in seconds.
  !    internally model time is kept in seconds (variables itau...) since
  !    1st-jan-iyear0, 00:00:00
  !   (iyear0 now defined as the actual year at start)
!---------------------------------------------------------------------------
  !    name    type      default   purpose
  !    ----    ----      -------   -------
  !
  !    ndyn    integer   1*3600    length of full advection step
  !    nconv   integer   1*3600    interval for convection calculation
  !    ndiff   integer   0         interval for horizontal diffusion calc
  !    nchem   integer   0         interval for chemistry calculations
  !    nsrce   integer   24*3600   interval for source calculation
  !    limits  logical   .true.    if set to .true. then
  !                                the slopes are limited
  !                                such that no negative tracer
  !                                masses should occur
  !    istart  integer   10        start/restart options:
  !                                1  coldstart with initial fields set to 0
  !                                2  coldstart with initial fields computed
  !                                   in sr trace1 in sources_sinks...
  !                                3  coldstart with initial
  !                                   fields read from model output (save file)
  !                                4  coldstart with initial
  !                                   fields read from model output stored
  !                                   in mixing ratio (no slopes).
  !    nread   integer   12*3600        interval for input of massfluxes and convection info
  !    nwrite  integer   0         interval for alternate output of restart
  !                                status on files save1.b and save2.b
  !    ninst   integer   0         interval for output of instantaneous
  !                                tracer mix ratio fields
  !    ncheck  integer   0         interval for output of tracer mix ratio
  !                                at checkpoints
  !    ndiag   integer   12*3600   interval for computing mean quantities
  !CMK ndiagp1 and ndiagp2 not implemented yet...
  !    ndiagp1 integer   -2        interval for output of
  !                                -1   daily
  !                                -2   monthly
  !                                -3   yearly
  !                                >=0  interval in seconds
  !    ndiagp2 integer   -2        interval for output of time averaged fields
  !                                -1   daily
  !                                -2   monthly
  !                                -3   yearly
  !                                >=0  interval in seconds
  !
  !    name      type     default  purpose
  !    ----      ----     -------  -------
  !
  !    icalendo  integer  2        calendar type
  !                                1 permanent 360 day year calendar
  !                                2 real calendar
  !                                3 permanent 365 day year calendar
  !                                4 permanent 366 day year calendar
  !    iyear0    integer  1980     base year for calendar calculations
  !                                (because of overflow problems this should
  !                                deviate on a 32 bit machine
  !                                not more than +-65 years from
  !                                any year actually used in the
  !                                model runs----> iyear0 now just the run year
  !
  !  date/times are expressed as yr,month,day,hour,min,sec
  !
  !    idatei(6) integer (1980 1 1 0 0 0) date/time for start of model run
  !    idatee(6) integer (1980 1 1 0 0 0) date/time for end of model run
  !    idatet(6) integer (1980 1 1 0 0 0) date/time after which instan-
  !                                taneous output is written (controlled
  !                                by 'ninst')
  !
  !    cdebug    logical  false    if true then output of debug info is
  !                                written on file debug.d
  !                                I old TM3 debugging mostly off!
  !
  !    okdebug   logical  true     TM5 debugging
  !    itau      integer           current model time
  !    idate(6)  integer           date corresponding to itau
  !    itaui     integer           start time (corresponds to idatei)
  !    itaue     integer           end time (corresponds to idatee)
  !    itaut     integer           time after which instantaneous output is
  !                                written (corresponds to idatet)
  !    itau0                       time/date when diagnostic arrays
  !    idate0(6) integer           were last reset
  !    julday0   integer           julian day of base time 1st-jan-iyear0, 0h
  !                                Needed only when icalendo == 2
  !    idacc(8)  integer           counters:
  !                                idacc(1)   no of times averaged tracer
  !                                mix ratio is calculated
  !                                others are not used at present
  !    newyr     logical .true.    if at beginning of a new year
  !    newmonth  logical .true.    if at beginning of a new month
  !    newday    logical .true.    if at beginning of a new day (i.e. at 00Z)
  !    newsrun   logical .true.    if at beginning of a new run or
  !                                   at beginning of a continuation run
  !    nstep     integer           advection step counter for current run
  !                                or continuation run
  !    nstep0    integer           not needed
  !    cpu0      real              process time at beginning of run (in sec)
  !    cpu1      real              process time at last reset time instant
  !    areag     real(nregions)    surface of globe and regions
  !    itaur     integer(nregions) time counter per region
  !
  !---------------------------------------------------------------------------

  character(len=DUMMY_STR_LEN) :: xlabel
  !
  !    variable      type      purpose
  !    --------      ----      -------
  !
  !    xlabel        char*160  run text label.
  !                            last 8 characters contain model version info
  !
  !----------------------------------------------------------------------------
  integer,dimension(nregions) :: unit_mix
  !
  integer,parameter :: kinput0=5
  !                        main control output
  integer,parameter :: kmain=6
  !                        secondary control input
  integer,parameter :: kdebug=9
  !                        temporary scratch files
  integer,parameter :: ktemp1=1

  !    czeta          real 1.        scaling factor for convection
  !    czetak         real 1.        scaling factor for vertical diffusion

  real              :: czeta,czetak


  ! levels not zoomed yet ...
  integer, parameter  ::  zbeg(nregions_max) = 0
  integer, parameter  ::  zend(nregions_max) = lm(1)



contains


  ! ==============================================


  subroutine CheckShape( shp1, shp2 )

    ! --- in/out -------------------------------------

    integer, intent(in)         ::  shp1(:)
    integer, intent(in)         ::  shp2(:)

    ! --- const -------------------------------------

    character(len=*), parameter  ::  rname = mname//', CheckShape'

    ! --- begin -------------------------------------

    if ( size(shp1) /= size(shp2) ) then
      write (*,'("ERROR - array shapes should have same length:")')
      write (*,'("ERROR -   shp1 : ",i4)') shp1
      write (*,'("ERROR -   shp2 : ",i4)') shp2
      write (*,'("ERROR in ",a)') rname; stop
    end if

    if ( any( shp1 /= shp2 ) ) then
      write (*,'("ERROR - array shapes are not equal:")')
      write (*,'("ERROR -   shp1 : ",i4)') shp1
      write (*,'("ERROR -   shp2 : ",i4)') shp2
      write (*,'("ERROR in ",a)') rname; stop
    end if

  end subroutine CheckShape


end module dims
