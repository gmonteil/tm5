!----------------------------------------------------------
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!----------------------------------------------------------

module em_file_info

  use GO, only : gol, goPr, goErr

  implicit none

  ! --- in/out -----------------------------

  private

  public  :: AnalyseLine, GetFileInfo

  ! --- const ------------------------------

  character(len=*), parameter  ::  mname = 'em_file_info'

  ! --- type -------------------------------

contains

  subroutine AnalyseLine( eml_line, species, source, category, year, &
       fscal, status )

    use GO, only : goReadFromLine
    use GO, only : goVarValue
    use os_specs, only : DUMMY_STR_LEN

    ! --- in/out -------------------------------------

    character(len=*), intent(in)             :: eml_line
    character(len=*), intent(out), optional  :: species
    character(len=*), intent(out), optional  :: source
    character(len=*), intent(out), optional  :: category
    integer, intent(out), optional           :: year
    real, intent(out), optional              :: fscal
    integer, intent(out)                     :: status

    ! --- const --------------------------------------

    character(len=*), parameter  :: rname = mname//'/AnalyseLine'

    ! --- local --------------------------------------

    character(len=DUMMY_STR_LEN) :: str, sp, ca, so
    integer                      :: ye
    real                         :: fs

    ! --- begin ---------------------------------

    ! Initial values
    str = eml_line
    ye = -1
    fs = 1.

    ! Split input line
    call goReadFromLine( str, sp, status, sep=' ' )
    IF_NOTOK_RETURN(status=1)
    str = adjustl( str )
    call goReadFromLine( str, ca, status, sep=' ' )
    IF_NOTOK_RETURN(status=1)
    str = adjustl( str )
    call goReadFromLine( str, ye, status, sep=' ' )
    IF_NOTOK_RETURN(status=1)
    str = adjustl( str )
    call goReadFromLine( str, so, status, sep=' ' )
    IF_NOTOK_RETURN(status=1)
    str = adjustl( str )

    ! extra fields:
    !    in old format the scale number :   1.23e4
    !    in new format                  :   scale=1.234e3   keywords=ant
    if ( trim(str) /= '' ) then
      ! test for one of the known keywords:
      if ( str(1:6) == 'scale=' ) then
        ! split line at white space, read value for scale:
        fs = 1.0
        call goVarValue( trim(str), ' ', 'scale', '=', fs, status )
        IF_ERROR_RETURN(status=1)
      else if ( str(1:9) == 'keywords=' ) then
        ! no need to process keywords here, only used for postprocessing
      else
        ! just a number, this is the scale factor:
        call goReadFromLine( str, fs, status, sep=' ' )
        IF_NOTOK_RETURN(status=1)
      end if
    end if

    if ( present(species) ) species = sp
    if ( present(source ) ) source = so
    if ( present(year)    ) year = ye
    if ( present(category)) category = ca
    if ( present(fscal)   ) fscal = fs

    status = 0

  end subroutine AnalyseLine


  ! ***


  !
  ! search in data base for corresponding file archive
  !

  subroutine GetFileInfo( req_source, req_species, req_category, req_year, status, &
                            filename, filetype, unit, nx, ny, nz, np, periods, keywords )

    ! --- modules ------------------------------------

    use GO, only : TTextFile, Init, Done, ReadLine
    use GO, only : goReadFromLine, goReplace
    use GO, only : Days_In_Year
    use os_specs, only : DUMMY_STR_LEN

    ! --- in/out -------------------------------------

    character(len=*), intent(in)            :: req_source
    character(len=*), intent(in)            :: req_species
    character(len=*), intent(in)            :: req_category
    integer, intent(in)                     :: req_year
    integer, intent(out)                    :: status
    character(len=*), intent(out), optional :: filename
    character(len=*), intent(out), optional :: filetype
    character(len=*), intent(out), optional :: unit
    integer, intent(out), optional          :: nx, ny, nz, np
    character(len=*), intent(out), optional :: periods
    character(len=*), intent(out), optional :: keywords

    ! --- const --------------------------------------

    character(len=*), parameter  :: rname = mname//'/GetFileInfo'
    character(len=*), parameter  :: info_fname = 'em_file_info.dat'

    ! --- local -------------------------------------

    character(len=DUMMY_STR_LEN) :: str
    integer           :: n
    character(len=DUMMY_STR_LEN):: fn_f
    character(len=32) :: ft_f
    character(len=32) :: so_f
    character(len=32) :: sp_f
    character(len=32) :: un_f
    character(len=64) :: ca_f
    character(len=32) :: year_descr
    integer           :: year_range(2)
    integer           :: nx_f, ny_f, nz_f, np_f
    character(len=32) :: periods_f
    type(TTextFile)   :: file

    ! --- begin ---------------------------------

    ! Open info file as commented text file
    call Init( file, info_fname, status, comment='!' )
    IF_NOTOK_RETURN(status=1)

    n = 0 ! counter for line number
    do
      n = n+1
      call ReadLine( file, str, status )
      if ( status < 0 ) then
         write (*,'("ERROR - no data file found with requested contents")')
         write (*,'("      - source:   ",a)') trim(req_source)
         write (*,'("      - species:  ",a)') trim(req_species)
         write (*,'("      - category: ",a)') trim(req_category)
         write (*,'("      - year:     ",i4.4)') req_year
         TRACEBACK; status=1; return
      end if
      if ( status > 0 ) then
         TRACEBACK; status=1; return
      end if

      ! space seperated lines with prescibed format:
      !read( str, '(a21,x,a1,x,a14,x,a5,x,a17,x,a14,5i5)', iostat=status ) &
      !     fn_f, ft_f, so_f, sp_f, un_f, ca_f, ye_f, nx_f, ny_f, nz_f, np_f
      !if ( status /= 0 ) then
      !   write (*,'("ERROR - problem reading `",a,"`, line",i7)') &
      !        trim(info_fname), n
      !   TRACEBACK; status=1; return
      !end if

      ! comma seperated fields:
      call goReadFromLine( str, fn_f, status )   ! file name
      IF_NOTOK_RETURN(status=1)
      call goReadFromLine( str, ft_f, status )   ! file type
      IF_NOTOK_RETURN(status=1)
      call goReadFromLine( str, so_f, status )   ! source
      IF_NOTOK_RETURN(status=1)
      call goReadFromLine( str, sp_f, status )   ! specie
      IF_NOTOK_RETURN(status=1)
      call goReadFromLine( str, un_f, status )   ! unit
      IF_NOTOK_RETURN(status=1)
      call goReadFromLine( str, ca_f, status )   ! category
      IF_NOTOK_RETURN(status=1)

      ! year description: '2000' or '2000-2005'
      call goReadFromLine( str, year_descr, status )
      IF_NOTOK_RETURN(status=1)
      ! read first value:
      call goReadFromLine( year_descr, year_range(1), status, sep='-' )
      IF_NOTOK_RETURN(status=1)
      ! value left ?
      if ( len_trim(year_descr) > 0 ) then
        ! read end of range:
        call goReadFromLine( year_descr, year_range(2), status )
        IF_NOTOK_RETURN(status=1)
      else
        ! single value, end of range is same as start:
        year_range(2) = year_range(1)
      end if

      ! grid:
      call goReadFromLine( str, nx_f, status )
      IF_NOTOK_RETURN(status=1)
      call goReadFromLine( str, ny_f, status )
      IF_NOTOK_RETURN(status=1)
      call goReadFromLine( str, nz_f, status )
      IF_NOTOK_RETURN(status=1)

      !call goReadFromLine( str, np_f, status )
      !IF_NOTOK_RETURN(status=1)
      ! flexible temporal resolution:
      call goReadFromLine( str, periods_f, status )
      IF_NOTOK_RETURN(status=1)
      select case ( trim(periods_f) )
        case ( 'hourly' )
          ! 24 hours per day ...
          np_f = Days_In_Year( req_year ) * 24   ! hours
        case ( 'daily' )
          ! one field per day:
          np_f = Days_In_Year( req_year ) ! days
        case default
          ! just a number:
          read (periods_f,*,iostat=status) np_f
          if ( status/=0 ) then
            write (gol,'("unsupported period string (or could not read number) : ",a)') trim(periods_f); call goErr
            TRACEBACK; status=1; return
          end if
      end select

      if ( req_source == so_f .and. req_species == sp_f .and. &
           req_category == ca_f .and. &
           (req_year >= year_range(1)) .and. (req_year <= year_range(2)) ) then
         if ( present(filename)  ) then
           ! initial value:
           filename   = fn_f
           ! replace some values if necessary:
           call goReplace( filename, '<yyyy>', '(i4.4)', req_year, status )
           IF_NOTOK_RETURN(status=1)
         end if
         if ( present(filetype)  ) filetype   = trim(ft_f)
         if ( present(unit)      ) unit       = trim(un_f)
         if ( present(nx)        ) nx         = nx_f
         if ( present(ny)        ) ny         = ny_f
         if ( present(nz)        ) nz         = nz_f
         if ( present(np)        ) np         = np_f
         if ( present(periods)   ) periods    = trim(periods_f)
         if ( present(keywords)  ) keywords   = trim(str)   ! remainder:  storage=sparse
         exit
      end if

    end do   ! lines in catalogue file

    call Done( file, status )
    IF_NOTOK_RETURN(status=1)

    status = 0

  end subroutine GetFileInfo

end module em_file_info
