!#######################################################################
!
! GO_Fields - storage and operation on 3D fields
!
! USAGE
!
!   use GO_Field
!
!   ! dimensions:
!   integer, parameter     ::  nx = 4, ny = 3, nz = 2
!
!   ! data:
!   type(T_Instant_Field_Series)    ::  F
!   type(TDate)                     ::  t
!   integer                         ::  nreceive
!   integer                         ::  ireceive
!   real                            ::  receive_t
!   real                            ::  data(nx,ny,nz)
!   integer                         ::  status
!
!   ! Initialize time series of 3D field.
!   ! Define the content by a name and units.
!   ! Define the temporal interpolation with a step and the units,
!   ! at the moment only linear interpolation is supported.
!   call Instant_Field_Series_Init( F, 'temperature', 'K', &
!                   'interpolation=linear;step=3;units=hour',
!                   status )
!   IF_NOTOK_STOP
!
!   ! current time, F%data should became valid for this:
!   t = NewDate( 2011, 03, 31, 18, 58 )
!
!   ! setup for time t; return status -1 if new data should be put:
!   call Instant_Field_Series_Setup( F, t, status )
!   IF_ERROR_STOP
!   if ( status < 0 ) then
!     ! prepare to receive new data:
!     call Instant_Field_Series_Setup_Prepare( F, nreceive, status )
!     IF_NOTOK_STOP
!     ! loop over fields to be received:
!     do ireceive = 1, nreceive
!       ! get time value:
!       call Instant_Field_Series_Setup_InqTime( F, ireceive, receive_t, status )
!       IF_NOTOK_STOP
!       ! obtain data from somewhere:
!       call read_input_data( ..., receive_t, ..., data, status )
!       IF_NOTOK_STOP
!       ! info ...
!       call wrtgol( '    put field valid for ', receive_t ); call goPr
!       ! store; optionally define the lower and upper bounds.
!       call Instant_Field_Series_Setup_Put( F, ireceive, data, receive_t, status, &
!                                              lbo=(/1,1,1/), ubo=(/nx,ny,nz/) )
!       IF_NOTOK_STOP
!     end do
!   end if
!
!   ! clear:
!   call Instant_Field_Series_Done( F, status )
!   if (status/=0) stop
!
!
! INHERITENCE
!
!   T_Field
!     T_Instant_Field
!       T_Instant_Field_Series
!
!
!#######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (line",i5,")")') __FILE__, __LINE__; call goErr
!
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status >0) then; TRACEBACK; action; return; end if
!
!#######################################################################


module GO_Field

  use GO_Print, only : gol, goPr, goErr
  use GO_Date , only : TDate

  implicit none
  
  
  ! ---- in/out ----------------------------------
  
  private
  
  public  ::  T_Instant_Field_Series

  public  ::  Instant_Field_Series_Init, Instant_Field_Series_Done
  public  ::  Instant_Field_Series_Setup
  public  ::  Instant_Field_Series_Setup_Prepare
  public  ::  Instant_Field_Series_Setup_InqTime
  public  ::  Instant_Field_Series_Setup_Put

  
  ! ---- const -----------------------------------
  
  ! module name:
  character(len=*), parameter        ::  mname = 'GO_Fields'
  
  ! character lengths:
  integer, parameter      ::  LEN_NAME  = 512
  integer, parameter      ::  LEN_UNITS = 64
  
  ! key values:
  integer, parameter      ::  INSTF_BASE = 1
  integer, parameter      ::  INSTF_TEND = 2
  
  
  ! --- types ------------------------------------

  type T_Field
    ! description:
    character(len=LEN_NAME)   ::  name
    character(len=LEN_UNITS)  ::  units
    ! data shape and bounds:
    integer                   ::  shp(3)
    integer                   ::  lbo(3), ubo(3)
    ! data:
    real, pointer             ::  data(:,:,:)
    ! filled ?
    logical                   ::  with_data
  end type
  
  ! *
  
  type, extends(T_Field) :: T_Instant_Field
    ! time for which field is valid:
    type(TDate)               ::  t
  end type T_Instant_Field
  
  ! *
  
  type, extends(T_Instant_Field) :: T_Instant_Field_Series
    ! temporal interpolation:
    character(len=32)         ::  interpolation
    integer                   ::  interpolation_step
    character(len=32)         ::  interpolation_step_units
    ! tendency, interval for which it is valid:
    real, pointer             ::  ddata_dt(:,:,:)   ! units/s
    type(TDate)               ::  tt(2)
    logical                   ::  with_tend
    ! flags to know if and which data is to be receied:
    type(TDate)               ::  setup_t
    integer                   ::  nreceive
    type(TDate)               ::  receive_time(2)
    integer                   ::  receive_targ(2)
  end type T_Instant_Field_Series
  

contains


  ! ====================================================================
  ! ===
  ! === Field
  ! ===
  ! ====================================================================


  subroutine Field_Init( F, name, units, status )
    
    ! --- in/out ---------------------------------
    
    class(T_Field), intent(inout)                 ::  F
    character(len=*), intent(in)                  ::  name, units
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Field_Init'
      
    ! --- begin ----------------------------------
    
    ! store:
    F%name = trim(name)
    F%units = trim(units)
    
    ! no data yet:
    nullify( F%data )
    F%shp = -999
    F%lbo = -999
    F%ubo = -999

    ! no content yet:
    F%with_data = .false.

    ! ok
    status = 0
    
  end subroutine Field_Init    


  ! ***
  

  subroutine Field_Alloc( F, shp, status, lbo, ubo )
    
    ! --- in/out ---------------------------------
    
    class(T_Field), intent(inout)                 ::  F
    integer, intent(in)                           ::  shp(3)
    integer, intent(out)                          ::  status
    integer, intent(in), optional                 ::  lbo(3), ubo(3)
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Field_Alloc'
      
    ! --- begin ----------------------------------
    
    ! bounds provided?
    if ( present(lbo) .or. present(ubo) ) then
      ! check ...
      if ( (.not. present(lbo)) .or. (.not. present(ubo)) ) then
        write (gol,'("both lbo and ubo should be provided")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! store:
      F%lbo = lbo
      F%ubo = ubo
      ! set shape:
      F%shp = F%ubo - F%lbo + 1
      ! check ...
      if ( any(F%shp /= shp) ) then
        write (gol,'("shape does not match with bounds :")'); call goErr
        write (gol,'("  shp  : ",3i6)') shp; call goErr
        write (gol,'("  lbo  : ",3i6)') lbo; call goErr
        write (gol,'("  ubo  : ",3i6)') ubo; call goErr
        TRACEBACK; status=1; return
      end if
    else
      ! set default bounds:
      F%lbo = 1
      F%ubo = shp
      ! store shape:
      F%shp = shp
    end if

    ! storage:
    allocate( F%data(F%lbo(1):F%ubo(1),F%lbo(2):F%ubo(2),F%lbo(3):F%ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Field_Alloc
  

  ! ====================================================================
  ! ===
  ! === Instant_Field
  ! ===
  ! ====================================================================


  subroutine Instant_Field_Init( F, name, units, status )

    use GO_Date  , only : AnyDate
    
    ! --- in/out ---------------------------------
    
    class(T_Instant_Field), intent(out)         ::  F
    character(len=*), intent(in)                ::  name, units
    integer, intent(out)                        ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Instant_Field_Init'

    ! --- begin ----------------------------------
    
    ! default field initialization:
    call Field_Init( F, name, units, status )
    IF_NOTOK_RETURN(status=1)

    ! no time assigned yet:
    F%t = AnyDate()
    
    ! ok
    status = 0
    
  end subroutine Instant_Field_Init
  

  ! ====================================================================
  ! ===
  ! === Instant_Field_Series
  ! ===
  ! ====================================================================


  subroutine Instant_Field_Series_Init( F, name, units, interp, status )
  
    use GO_String, only : goVarValue
    use GO_Date  , only : AnyDate
    
    ! --- in/out ---------------------------------
    
    type(T_Instant_Field_Series), intent(out)   ::  F
    character(len=*), intent(in)                ::  name, units
    character(len=*), intent(in)                ::  interp
    integer, intent(out)                        ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Instant_Field_Series_Init'

    ! --- begin ----------------------------------
    
    ! initialize the instant field:
    call Field_Init( F, name, units, status )
    IF_NOTOK_RETURN(status=1)
    
    ! extract interpolation key:
    call goVarValue( interp, ';', 'interpolation', '=', F%interpolation, status )
    IF_NOTOK_RETURN(status=1)
    ! which type ?
    select case ( F%interpolation )
      !~ linear:
      case ( 'linear' )
        ! read time step:
        call goVarValue( interp, ';', 'step', '=', F%interpolation_step, status )
        IF_NOTOK_RETURN(status=1)
        ! read time step units:
        call goVarValue( interp, ';', 'units', '=', F%interpolation_step_units, status )
        IF_NOTOK_RETURN(status=1)
        ! no tendency data yet:    
        nullify( F%ddata_dt )
        F%tt(1) = AnyDate()
        F%tt(2) = AnyDate()
        F%with_tend = .false.
      !~ unknown
      case default
        write (gol,'("unsupported temporal interpolation `",a,"`")') trim(F%interpolation); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0
    
  end subroutine Instant_Field_Series_Init
  

  ! ***


  subroutine Instant_Field_Series_Done( F, status )
  
    use GO_Date, only : AnyDate
  
    ! --- in/out ---------------------------------
    
    type(T_Instant_Field_Series), intent(inout)    ::  F
    integer, intent(out)                    ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Instant_Field_Series_Done'
      
    ! --- begin ----------------------------------
    
    ! no data anymore:
    F%name = 'None'
    F%units = 'None'
    
    ! no current data:
    if ( F%with_data ) deallocate( F%data )
    F%t = AnyDate()
    
    ! no tendency data anymore:
    if ( F%with_tend ) deallocate( F%ddata_dt )
    F%tt(1) = AnyDate()
    F%tt(2) = AnyDate()
    
    ! ok
    status = 0
    
  end subroutine Instant_Field_Series_Done
  

  ! ***


  subroutine Instant_Field_Series_Setup( F, t, status )
  
    use GO_Date, only : operator(/=), operator(<=)
    use GO_Date, only : wrtgol
 
    ! --- in/out ---------------------------------
    
    type(T_Instant_Field_Series), intent(inout)    ::  F
    type(TDate), intent(in)                 ::  t
    integer, intent(out)                    ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Instant_Field_Series_Setup'

    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! store target time:
    F%setup_t = t
    
    ! no data present yet ?
    if ( .not. F%with_data ) then
      !! info ...
      !write (gol,'(a,": no data filled yet")') rname; call goPr
      ! new data required ...
      status = -1; return
    else
      ! data defined; not valid for this time ?
      if ( F%t /= t ) then
        !! info ...
        !write (gol,'(a,": data not valid for target time")') rname; call goPr
        ! tendency not defined yet ?
        if ( .not. F%with_tend ) then
          !! info ...
          !write (gol,'(a,": no tendency data yet; need to receive new data")') rname; call goPr
          ! new data required:
          status = -1; return
        else if ( (F%tt(1) <= t) .and. (t <= F%tt(2)) ) then
          !! info ...
          !call wrtgol( rname//': target time ', t, ' in interpolation interval ', F%tt(1), ' - ', F%tt(2) ); call goPr
          ! interpolate:
          call Instant_Field_Series_Interpol( F, t, status )
          IF_NOTOK_RETURN(status=1)
        else 
          !! info ...
          !write (gol,'(a,": target time outside interpolation interval; need to receive new data")') rname; call goPr
          ! new data required:
          status = -1; return
        end if
      end if
    end if
    
    ! ok
    status = 0
    
  end subroutine Instant_Field_Series_Setup
  

  ! ***


  ! Interpolate to t given available data and tendency.
  ! Tool for use within this module only.
  
  subroutine Instant_Field_Series_Interpol( F, t, status )
  
    use GO_Date, only : operator(/=), operator(-), rTotal
    use GO_Date, only : wrtgol
 
    ! --- in/out ---------------------------------
    
    type(T_Instant_Field_Series), intent(inout)    ::  F
    type(TDate), intent(in)                 ::  t
    integer, intent(out)                    ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Instant_Field_Series_Interpol'

    ! --- local ----------------------------------
    
    real      ::  dt_sec
    
    ! --- begin ----------------------------------
    
    !! info ...
    !call wrtgol( rname//': interpolate between [', F%tt(1), ',', F%tt(2), '] to ', t ); call goPr

    ! data defined by current base and tendency;
    ! interpolate:  data(t) = data(F%t) + ddata(F%t)/dt * (t-F%t)
    ! timestep in seconds:
    dt_sec = rTotal( t - F%t, 'sec' )
    ! update data:
    F%data = F%data + F%ddata_dt * dt_sec
    ! store new time:
    F%t = t
    
    ! ok
    status = 0
    
  end subroutine Instant_Field_Series_Interpol
  

  ! ***


  subroutine Instant_Field_Series_Setup_Prepare( F, nreceive, status )
  
    use GO_Date, only : TDate
    use GO_Date, only : operator(==), operator(/=), operator(<), operator(<=)
    use GO_Date, only : wrtgol
    use GO_Date, only : Get_Surrounding_Interval
 
    ! --- in/out ---------------------------------
    
    type(T_Instant_Field_Series), intent(inout)    ::  F
    integer, intent(out)                           ::  nreceive
    integer, intent(out)                           ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Instant_Field_Series_Setup_Prepare'
    
    ! --- local ----------------------------------
    
    type(TDate)   ::  tt(2)

    ! --- begin ----------------------------------
    
    !! info ...
    !call wrtgol( rname//': prepare setup for ', F%setup_t ); call goPr
    
    ! some fields need to be changed; 
    ! get interpolation points given temporal resolution:
    call Get_Surrounding_Interval( F%setup_t, &
                                   F%interpolation_step, F%interpolation_step_units, &
                                   tt, status )
    IF_NOTOK_RETURN(status=1)
    !! info ...
    !call wrtgol( rname//': surrounding interval ', tt(1), ' - ', tt(2) ); call goPr

    ! single interpolation point (at target time) or interval ?
    !~ single point
    if ( tt(1) == tt(2) ) then
      ! already data present ?
      if ( F%with_data ) then
        ! try if current data is already valid for target time;
        ! this routine should not have been called ...
        if ( F%t == F%setup_t ) then
          call wrtgol( 'no need to prepare setup, current data already valid for ', F%setup_t ); call goErr
          TRACEBACK; status=1; return
        end if
        ! idem for interpolation:
        if ( F%with_tend .and. (F%tt(1) <= F%setup_t) .and. (F%setup_t <= F%tt(2)) ) then
          call wrtgol( 'no need to prepare setup, interpolation interval includes ', F%setup_t ); call goErr
          TRACEBACK; status=1; return
        end if
      end if
      ! receive a single field valid for the target time, use it for the base:
      F%nreceive = 1
      F%receive_time(1) = tt(1)
      F%receive_targ(1) = INSTF_BASE
    !~ interval, interpolate to target time
    else
      ! already data present ?
      if ( F%with_data ) then
        ! check if is valid for one of the interval bounds:
        if ( F%t == tt(1) ) then
          ! receive the other field to compute the tendency:
          F%nreceive = 1
          F%receive_time(1) = tt(2)
          F%receive_targ(1) = INSTF_TEND
        else if ( F%t == tt(2) ) then
          ! receive the other field to compute the tendency:
          F%nreceive = 1
          F%receive_time(1) = tt(1)
          F%receive_targ(1) = INSTF_TEND
        else if ( F%with_tend .and. (F%tt(1) <= tt(1)) .and. (tt(1) <= F%tt(2)) ) then
          ! interpolate to this bound:
          call Instant_Field_Series_Interpol( F, tt(1), status )
          IF_NOTOK_RETURN(status=1)
          ! receive the other field to compute the tendency:
          F%nreceive = 1
          F%receive_time(1) = tt(2)
          F%receive_targ(1) = INSTF_TEND
        else if ( F%with_tend .and. (F%tt(1) <= tt(2)) .and. (tt(2) <= F%tt(2)) ) then
          ! interpolate to this bound:
          call Instant_Field_Series_Interpol( F, tt(2), status )
          IF_NOTOK_RETURN(status=1)
          ! receive the other field to compute the tendency:
          F%nreceive = 1
          F%receive_time(1) = tt(1)
          F%receive_targ(1) = INSTF_TEND
        else
          ! receive fields for both interval bounds:
          F%nreceive = 2
          F%receive_time(1) = tt(1)
          F%receive_targ(1) = INSTF_BASE
          F%receive_time(2) = tt(2)
          F%receive_targ(2) = INSTF_TEND
        end if
      else
        ! receive fields for both interval bounds:
        F%nreceive = 2
        F%receive_time(1) = tt(1)
        F%receive_targ(1) = INSTF_BASE
        F%receive_time(2) = tt(2)
        F%receive_targ(2) = INSTF_TEND
      end if  ! data present
    end if  ! single field or interval
    
    ! return value:
    nreceive = F%nreceive

    ! ok
    status = 0
    
  end subroutine Instant_Field_Series_Setup_Prepare
  

  ! ***


  subroutine Instant_Field_Series_Setup_InqTime( F, ip, t, status )
 
    ! --- in/out ---------------------------------
    
    type(T_Instant_Field_Series), intent(inout)    ::  F
    integer, intent(in)                     ::  ip
    type(TDate), intent(out)                ::  t
    integer, intent(out)                    ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Instant_Field_Series_Setup_InqTime'

    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (ip < 1) .or. (ip > F%nreceive) ) then
      write (gol,'("index ",i6," of interpolation point outside bounds 1 .. ",i6)') F%nreceive; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! extract:
    t = F%receive_time(ip)
    
    ! ok
    status = 0
    
  end subroutine Instant_Field_Series_Setup_InqTime
  

  ! ***


  subroutine Instant_Field_Series_Setup_Put( F, ip, data, t, status, lbo, ubo )
  
    use GO_Date, only : rTotal, operator(-), operator(==), operator(/=)
    use GO_Date, only : wrtgol
 
    ! --- in/out ---------------------------------
    
    type(T_Instant_Field_Series), intent(inout)     ::  F
    integer, intent(in)                             ::  ip
    real, intent(in)                                ::  data(:,:,:)
    type(TDate), intent(in)                         ::  t
    integer, intent(out)                            ::  status
    integer, intent(in), optional                   ::  lbo(3), ubo(3)

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Instant_Field_Series_Setup_Put'

    ! --- local ----------------------------------

    real        ::  dt_sec
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (ip < 1) .or. (ip > F%nreceive) ) then
      write (gol,'("index ",i6," of interpolation point outside bounds 1 .. ",i6)') ip, F%nreceive; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! need to allocate ?
    if ( .not. associated(F%data) ) then
      ! allocate data:
      call Field_Alloc( F, shape(data), status, lbo=lbo, ubo=ubo )
      IF_NOTOK_RETURN(status=1)
    else    
      ! check ...
      if ( any(shape(data) /= F%shp) ) then
        write (gol,'("shape of data does not match with definition:")'); call goErr
        write (gol,'("  data     : ",3i6)') shape(data); call goErr
        write (gol,'("  defined  : ",3i6)') F%shp; call goErr
        TRACEBACK; status=1; return
      end if
    end if

    ! check ...
    if ( t /= F%receive_time(ip) ) then
      write (gol,'("time of received field ",i6," does not match with setup:")') ip; call goErr
      call wrtgol( '  setup time    : ', F%receive_time(ip) ); call goErr
      call wrtgol( '  received time : ', t ); call goErr
      TRACEBACK; status=1; return
    end if

    ! what to do with this data ?
    select case ( F%receive_targ(ip) )
      !~~ store as base:
      case ( INSTF_BASE )
        ! store data:
        F%data = data
        ! store time:
        F%t = t
        ! set flag:
        F%with_data = .true.
        !! info ...
        !call wrtgol( rname//': stored base data for ', F%t ); call goPr
        !write (gol,'(a,": average value : ",f12.4)') rname, sum(F%data)/size(F%data); call goPr
      !~~ use for tendency then ?
      case ( INSTF_TEND )
        ! setup tendency if necessary:
        if ( .not. associated(F%ddata_dt) ) then
          ! storage:
          allocate( F%ddata_dt(F%lbo(1):F%ubo(1),F%lbo(2):F%ubo(2),F%lbo(3):F%ubo(3)), stat=status )
          IF_NOTOK_RETURN(status=1)
        end if
        ! check ...
        if ( t == F%t ) then
          call wrtgol( rname//': data already present for time ', t ); call goErr
          TRACEBACK; status=1; return
        end if
        ! timestep in seconds:
        dt_sec = rTotal( t - F%t, 'sec' )
        ! compute and store tendency:
        F%ddata_dt = ( data - F%data )/dt_sec
        ! store time range:
        F%tt(1) = F%t
        F%tt(2) = t
        ! set flag:
        F%with_tend = .true.
        !! info ...
        !call wrtgol( rname//': stored tendency data for ', F%tt(1), ' - ', F%tt(2) ); call goPr
        !write (gol,'(a,": average value end data : ",f12.4)') rname, sum(data)/size(data); call goPr
        !write (gol,'(a,": average value tendency : ",e12.4," over ",f12.6," sec")') rname, &
        !                     sum(F%ddata_dt)/size(F%ddata_dt), dt_sec; call goPr
        ! interpolate:
        call Instant_Field_Series_Interpol( F, F%setup_t, status )
        IF_NOTOK_RETURN(status=1)
      !~~ strange ...
      case default
        ! something went wrong; bug ?
        write (gol,'("unsupported destination code ",i6," for interpolation point ",i6)') F%receive_targ(ip), ip; call goErr
        TRACEBACK; status=1; return
    end select
    
    ! reset destination to dummy value:
    F%receive_targ(ip) = -1

    ! ok
    status = 0
    
  end subroutine Instant_Field_Series_Setup_Put
  

end module GO_Field



!!#######################################################################
!!###
!!### Test
!!###
!!#######################################################################
!!
!#define IF_NOTOK_STOP if (status/=0) then; TRACEBACK; stop; end if
!#define IF_ERROR_STOP if (status >0) then; TRACEBACK; stop; end if
!!
!! f90 -o test.x go_fu.F90 go_print.F90 go_string.F90 go_date.F90 go_field.F90
!
!program test
!
!  use GO_Print , only : gol, goPr, goErr
!  use GO_Date  , only : TDate, NewDate, wrtgol, IsAnyDate
!  use GO_Fields
!  
!  implicit none
!  
!  ! dimensions:
!  integer, parameter     ::  nx = 4, ny = 3, nz = 2
!
!  ! data:
!  type(T_Instant_Field_Series)       ::  F
!  integer                     ::  hour
!  type(TDate)                 ::  t
!  integer                     ::  nreceive
!  type(TDate)                 ::  receive_time
!  integer                     ::  ireceive
!  real                        ::  data_in(nx,ny,nz)
!  integer                     ::  status
!  
!  ! info ...
!  write (gol,'("test:")'); call goPr
!  write (gol,'("test: testing GO_Fields module")'); call goPr
!  write (gol,'("test:")'); call goPr
!
!  ! Initialize 3D field; provide name and units,
!  ! and description of temporal interpolation.
!  ! optionally provide shape or lower/upper bounds:
!  call Instant_Field_Series_Init( F, 'temperature', 'K', shape(data_in), &
!                 'interpolation=linear;step=3;units=hour', status )
!  IF_NOTOK_STOP
!  
!  ! loop over hours:
!  do hour = 0, 6
!  
!    write (gol,'("test:")'); call goPr
!    write (gol,'("test: >>> hour ",i2)') hour; call goPr
!    write (gol,'("test:")'); call goPr
!
!    ! current time:
!    t = NewDate( 2011, 03, 31, hour, 0 )
!
!    ! info ...
!    call wrtgol( 'test: target time : ', t ); call goPr
!
!    ! info ...
!    write (gol,'("test: setup ...")'); call goPr
!    ! setup for time t, or return status -1 if new data should be received:
!    call Instant_Field_Series_Setup( F, t, status )
!    IF_ERROR_STOP
!    if ( status < 0 ) then
!      ! info ...
!      write (gol,'("test: something to be done ...")'); call goPr
!      ! interpolate, or prepare to receive new data:
!      call Instant_Field_Series_Setup_Prepare( F, nreceive, status )
!      IF_NOTOK_STOP
!      ! info ...
!      write (gol,'("test: number of fields to receive : ",i6)') nreceive; call goPr
!      ! might need to receive new data ...
!      do ireceive = 1, nreceive
!        ! get time value:
!        call Instant_Field_Series_Setup_InqTime( F, ireceive, receive_time, status )
!        IF_NOTOK_STOP
!        ! obtain data from somewhere; here fill with hour plus decimal minutes ...
!        data_in = receive_time%hour + receive_time%min/60.0/10.0
!        ! info ...
!        call wrtgol( 'test: put field valid for ', receive_time ); call goPr
!        write (gol,'("test: average value : ",f12.4)') sum(data_in)/size(data_in); call goPr
!        ! store:
!        call Instant_Field_Series_Setup_Put( F, ireceive, data_in, receive_time, status )
!        IF_NOTOK_STOP
!      end do
!    end if
!    ! info ...
!    call wrtgol( 'test: data valid for : ', F%t ); call goPr
!    write (gol,'("test: average value : ",f12.4)') sum(F%data)/size(F%data); call goPr
!    
!  end do
!
!  ! clear:
!  call Instant_Field_Series_Done( F, status )
!  IF_NOTOK_STOP
!  
!  ! info ...
!  write (gol,'("test:")'); call goPr
!  write (gol,'("test: end.")'); call goPr
!  write (gol,'("test:")'); call goPr
!  
!end program test


