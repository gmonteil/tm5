!#################################################################
!
! Fortran module around UDUnits .
!
! USAGE
!
!   use UDUnits
!
!   integer                     ::  status
!   integer(UD_POINTER_KIND)    ::  unit, unit2
!   character(len=64)           ::  spec, spec2
!   real(8)                     ::  slope, offset
!
!   ! * module initialisation
!
!   call UDUnits_Init( status )
!   if (status/=UDUNITS_NOERR) then
!     print *, trim(UDUnits_StrError(status))
!     stop
!   end if
!
!   ! * high levell routines
!
!   spec = 'kg/s'
!   call UDUnits_Standard( spec, spec2, status )
!   write (*,'("standard name of `",a,"` is `",a,"`")') trim(spec), trim(spec2)
!
!   spec = 'gram/cm3' ; spec2 = 'kg/m3'
!   call UDUnits_ConversionFactor( spec, spec2, slope, status )
!   write (*,'("conversion factor from `",a,"` to `",a,"` is ",f12.4)') trim(spec), trim(spec2), slope
!
!   ! * low level routines
!
!   call UDUnits_Make( unit, status )
!   call UDUnits_Make( unit2, status )
!
!   call UDUnits_Decode( 'kg', unit, status )
!   call UDUnits_Encode( unit, spec, status )
!
!   call UDUnits_Convert( unit, unit2, slope, offset, status )
!
!   ! * done with module
!
!   call UDUnits_Done( status )
!
! HISTORY
!   2010 feb, Arjo Segers, JRC
!
!#################################################################


module UDUnits

  use UDUnits_Inc,  only : UD_POINTER_KIND
  use os_specs,     only : WRITE_STR_LEN, MAX_FILENAME_LEN

  implicit none


  ! --- in/out -----------------------------------

  private

  public    ::  UD_POINTER_KIND
  public    ::  UDUNITS_NOERR, UDUnits_StrError

  public    ::  UDUnits_Init, UDUnits_Done
  public    ::  UDUnits_Make
  public    ::  UDUnits_Decode, UDUnits_Encode
  public    ::  UDUnits_Convert

  public    ::  UDUnits_Standard
  public    ::  UDUnits_ConversionFactor


  ! --- const --------------------------------

  ! module name:
  character(len=*), parameter  ::  mname = 'UDUnits'

  ! name of environment variable with path to data file:
  character(len=*), parameter ::  env_var = 'UDUNITS_PATH'

!  ! unit should be of type :  integer(UD_POINTER_KIND)
!  integer, parameter  ::  UD_POINTER_KIND = 4

  ! no error:
  integer, parameter  ::  UDUNITS_NOERR = 0

!  ! error codes:
!  integer, parameter  ::  UT_EOF      =   1
!  integer, parameter  ::  UT_ENOFILE  =  -1
!  integer, parameter  ::  UT_ESYNTAX  =  -2
!  integer, parameter  ::  UT_EUNKNOWN =  -3
!  integer, parameter  ::  UT_EIO      =  -4
!  integer, parameter  ::  UT_EINVALID =  -5
!  integer, parameter  ::  UT_ENOINIT  =  -6
!  integer, parameter  ::  UT_ECONVERT =  -7
!  integer, parameter  ::  UT_EALLOC   =  -8
!  integer, parameter  ::  UT_ENOROOM  =  -9
!  integer, parameter  ::  UT_ENOTTIME = -10
!
  !integer, parameter  ::  UT_MAXNUM_BASE_QUANTITIES = 10

  ! storage for latest error:
  integer, parameter            ::  error_status  = -100
  character(len=WRITE_STR_LEN)  ::  error_message = ''

  ! maximum length of specifications:
  integer, parameter  ::  spec_len = 64



contains


  ! ====================================================================
  ! ===
  ! === module routines
  ! ===
  ! ====================================================================


  subroutine UDUnits_Init( status )

    !use UDUnits_Inc, only : udunits_inc_test

    ! --- in/out ---------------------------------

    integer, intent(out)      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/UDUnits_Init'

    ! --- external -------------------------------

    ! Initialize the units package:
    integer, external  :: UTOpen

    ! --- local ----------------------------------

    character(len=MAX_FILENAME_LEN) ::  UDUnits_path
    integer                         ::  length

    ! --- begin ----------------------------------

    !call udunits_inc_test()

    ! following the manuals, the path to the udunits data file is
    ! taken from the environment variable UDUnits_PATH if not specified;
    ! this does not seem to work properly however, and therefore
    ! the path is explicitly taken from the environment:
    call Get_Environment_Variable( env_var, UDUnits_path, length, status )
    if (status/=0) then
      write (error_message,'("could not get environment variable `",a,"`")') trim(env_var)
      status=error_status; return
    end if

    ! Initialize the units package:
    status = UTOpen( trim(UDUnits_path) )
    if (status/=0) write (error_message,'("could not initialize from data file `",a,"`")') trim(UDUnits_path)

  end subroutine UDUnits_Init


  ! ***


  subroutine UDUnits_Done( status )

    ! --- in/out ---------------------------------

    integer, intent(out)      ::  status

    ! --- external -------------------------------

    ! --- begin ----------------------------------

    ! function UTFree not available in Fortran interface ...

    ! ok
    status = 0

  end subroutine UDUnits_Done


  ! ====================================================================
  ! ===
  ! === error messages
  ! ===
  ! ====================================================================


  function UDUnits_StrError( status )

    use UDUnits_Inc, only : UT_EOF, UT_ENOFILE, UT_ESYNTAX, UT_EUNKNOWN, &
                            UT_EIO, UT_EINVALID, UT_ENOINIT, UT_ECONVERT, &
                            UT_EALLOC, UT_ENOROOM, UT_ENOTTIME

    ! --- in/out ---------------------------------

    character(len=WRITE_STR_LEN)    ::  UDUnits_StrError
    integer, intent(inout)          ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/UDUnits_StrError'

    ! --- begin ----------------------------------

    ! display message:
    select case ( status )
      ! supported:
      case ( UT_EOF        ) ; UDUnits_StrError = 'End of file'
      case ( UDUNITS_NOERR ) ; UDUnits_StrError = ''
      case ( UT_ENOFILE    ) ; UDUnits_StrError = 'Units file does not exist'
      case ( UT_ESYNTAX    ) ; UDUnits_StrError = 'Syntax error'
      case ( UT_EUNKNOWN   ) ; UDUnits_StrError = 'Unknown unit specification'
      case ( UT_EIO        ) ; UDUnits_StrError = 'I/O error while accessing the units file'
      case ( UT_EINVALID   ) ; UDUnits_StrError = 'Invalid value'
      case ( UT_ENOINIT    ) ; UDUnits_StrError = 'Package has not be initialized'
      case ( UT_ECONVERT   ) ; UDUnits_StrError = 'Conversion error'
      case ( UT_EALLOC     ) ; UDUnits_StrError = 'Memory allocation failure'
      case ( UT_ENOROOM    ) ; UDUnits_StrError = 'No room for result'
      case ( UT_ENOTTIME   ) ; UDUnits_StrError = 'No time value'
      ! other ...
      case ( error_status ) ; UDUnits_StrError = ''
      ! unknown:
      case default
        write (UDUnits_StrError,'("Unknown error status from UDUnits routine : ",i6)') status
    end select

    ! add error buffer:
    if ( status /= 0 ) then
      UDUnits_StrError = trim(UDUnits_StrError)//'; '//trim(error_message)
    end if

  end function UDUnits_StrError


  ! ====================================================================
  ! ===
  ! === low level routines
  ! ===
  ! ====================================================================


  subroutine UDUnits_Make( unit, status )

    ! --- in/out ---------------------------------

    integer(UD_POINTER_KIND), intent(out)   ::  unit
    integer, intent(out)            ::  status

    ! --- external -------------------------------

    ! set return status:
    integer(UD_POINTER_KIND), external  :: UTMake

    ! --- begin ----------------------------------

    ! Create a new unit:
    unit = UTMake()

    ! set return status:
    status = 0
    if ( unit < 0 ) status = int(unit)

  end subroutine UDUnits_Make


  ! ***


  subroutine UDUnits_Decode( spec, unit, status )

    ! --- in/out ---------------------------------

    character(len=*), intent(in)      ::  spec
    integer(UD_POINTER_KIND), intent(in)      ::  unit
    integer, intent(out)              ::  status

    ! --- external -------------------------------

    ! Decode a formatted specification into a unit:
    integer, external  :: UTDec

    ! --- begin ----------------------------------

    ! Decode a formatted specification into a unit:
    status = UTDec( spec, unit )
    if (status/=0) write (error_message,'("could not decode `",a,"`")') trim(spec)

  end subroutine UDUnits_Decode


  ! ***


  subroutine UDUnits_Encode( unit, spec, status )

    ! --- in/out ---------------------------------

    integer(UD_POINTER_KIND), intent(in)      ::  unit
    character(len=*), intent(out)     ::  spec
    integer, intent(out)              ::  status

    ! --- external -------------------------------

    ! Encode a unit into a formatted specification:
    integer, external  :: UTEnc

    ! --- begin ----------------------------------

    ! Encode a unit into a formatted specification:
    status = UTEnc( unit, spec )
    if (status/=0) write (error_message,'("could not encode from unit into formatted specification")')

  end subroutine UDUnits_Encode


  ! ***


  subroutine UDUnits_Convert( unit_from, unit_to, slope, intercept, status )

    ! --- in/out ---------------------------------

    integer(UD_POINTER_KIND), intent(in)      ::  unit_from
    integer(UD_POINTER_KIND), intent(in)      ::  unit_to
    real(8), intent(out)              ::  slope, intercept
    integer, intent(out)              ::  status

    ! --- external -------------------------------

    ! Convert from one unit to another:
    integer, external  :: UTCvt

    ! --- local ----------------------------------

    character(len=spec_len)   ::  spec_from, spec_to

    ! --- begin ----------------------------------

    ! Convert from one unit to another:
    status = UTCvt( unit_from, unit_to, slope, intercept )
    if (status/=0) then
      call UDUnits_Encode( unit_from, spec_from, status )
      if (status/=0) then
        write (error_message,'("could not convert units; failed to convert unit_from to specification")')
        status = error_status; return
      end if
      call UDUnits_Encode( unit_to, spec_to, status )
      if (status/=0) then
        write (error_message,'("could not convert from `",a,"`; failed to convert unit_to to specification")') trim(spec_from)
        status = error_status; return
      end if
      write (error_message,'("could not convert from `",a,"` to `",a,"`")') trim(spec_from), trim(spec_to)
      status = error_status; return
    end if

  end subroutine UDUnits_Convert


  ! ====================================================================
  ! ===
  ! === high level routines
  ! ===
  ! ====================================================================


  subroutine UDUnits_Standard( spec_from, spec_to, status )

    ! --- in/out ---------------------------------

    character(len=*), intent(in)      ::  spec_from
    character(len=*), intent(out)     ::  spec_to
    integer, intent(out)              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/UDUnits_Standard'

    ! --- external -------------------------------

    ! Convert from one unit to another:
    integer, external  :: UTCvt

    ! --- local ----------------------------------

    integer(UD_POINTER_KIND)      ::  unit_from

    ! --- begin ----------------------------------

    ! setup unit:
    call UDUnits_Make( unit_from, status )
    if (status/=0) return
    ! fill with secification:
    call UDUnits_Decode( spec_from, unit_from, status )
    if (status/=0) return
    ! extract standard name:
    call UDUnits_Encode( unit_from, spec_to, status )
    if (status/=0) return

    ! ok
    status = 0

  end subroutine UDUnits_Standard


  ! ***


  subroutine UDUnits_ConversionFactor( spec_from, spec_to, factor, status )

    ! --- in/out ---------------------------------

    character(len=*), intent(in)      ::  spec_from
    character(len=*), intent(in)      ::  spec_to
    real(8), intent(out)              ::  factor
    integer, intent(out)              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/UDUnits_ConversionFactor'

    ! --- external -------------------------------

    ! Convert from one unit to another:
    integer, external  :: UTCvt

    ! --- local ----------------------------------

    integer(UD_POINTER_KIND)      ::  unit_from
    integer(UD_POINTER_KIND)      ::  unit_to
    real(8)                    ::  offset

    ! --- begin ----------------------------------

    ! input unit:
    call UDUnits_Make( unit_from, status )
    if (status/=0) return
    call UDUnits_Decode( spec_from, unit_from, status )
    if (status/=0) return

    ! output unit:
    call UDUnits_Make( unit_to, status )
    if (status/=0) return
    call UDUnits_Decode( spec_to, unit_to, status )
    if (status/=0) return

    ! Convert from one unit to another:
    call UDUnits_Convert( unit_from, unit_to, factor, offset, status )
    if (status/=0) return

    ! check ...
    if ( offset /= 0.0d0 ) then
      write (error_message,*) 'found conversion offset unequal to zero : ', offset
      status=error_status; return
    end if

    ! ok
    status = 0

  end subroutine UDUnits_ConversionFactor


end module UDUnits


!! ######################################################################
!! ###
!! ### test
!! ###
!! ######################################################################
!
!! f90 -o test.x udunits_inc.F udunits.F90 -I${UDUNITS_HOME}/include -L${UDUNITS_HOME}/lib -ludunits  &&  ./test.x
!
!program test_udunits
!
!  use UDUnits
!
!  implicit none
!
!  integer                     ::  status
!  integer(UD_POINTER_KIND)    ::  unit, unit2
!  character(len=64)           ::  spec, spec2
!  real(8)                     ::  slope, offset
!
!  write (*,'("begin test_udunits")')
!
!  write (*,'("   UD_POINTER_KIND : ",i4)') UD_POINTER_KIND
!
!  ! * module initialisation
!
!  call UDUnits_Init( status )
!  if (status/=UDUNITS_NOERR) then
!    print *, trim(UDUnits_StrError(status))
!    stop
!  end if
!
!  ! * high levell routines
!
!  spec = 'kg/s'
!  call UDUnits_Standard( spec, spec2, status )
!  write (*,'("  standard name of `",a,"` is `",a,"`")') trim(spec), trim(spec2)
!
!  spec = 'gram/cm3' ; spec2 = 'kg/m3'
!  call UDUnits_ConversionFactor( spec, spec2, slope, status )
!  write (*,'("  conversion factor from `",a,"` to `",a,"` is ",f12.4)') trim(spec), trim(spec2), slope
!
!  ! * low level routines
!
!  call UDUnits_Make( unit, status )
!  call UDUnits_Make( unit2, status )
!
!  call UDUnits_Decode( 'kg', unit, status )
!  call UDUnits_Encode( unit, spec, status )
!
!  call UDUnits_Convert( unit, unit2, slope, offset, status )
!
!  ! * done with module
!
!  call UDUnits_Done( status )
!
!  ! *
!
!  write (*,'("end test_udunits")')
!
!end program test_udunits
