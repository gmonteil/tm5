!
! Interface to NCEP Cray Binary file
!

module file_ncb

  use os_specs, only : MAX_FILENAME_LEN

  implicit none


  ! --- in/out ----------------------------------------------

  private

  public  ::  TNcepCrayBin
  public  ::  Init, Done
  public  ::  ReadRecord


  ! --- const -----------------------------------------------

  character(len=*), parameter ::  mname = 'file_ncb'

  ! data kinds
  integer, parameter      ::  ncb_ikind = 4
  integer, parameter      ::  ncb_rkind = 4

  ! header length (chars)
  integer, parameter      ::  ncb_lhead = 32

  ! length of rec2 in 4-bytes real:
  integer, parameter      ::  lrec2 = 250

  ! max levels
  integer, parameter      ::  maxlev = 100


  ! --- types --------------------------------------------

  type TNcepCrayBin
    ! file unit and name:
    integer                      ::  fu
    character(len=MAX_FILENAME_LEN) ::  fname
    ! header
    character(len=ncb_lhead)     ::  head
    ! time etc
    integer                      ::  fcst_hr
    integer                      ::  idate(4)
    ! levels
    integer                      ::  nlev
    real, pointer                ::  sigma_half(:)
    real, pointer                ::  sigma_full(:)
    ! spectral fields:
    integer                      ::  shT, shn
  end type TNcepCrayBin


  ! --- interfaces ----------------------------------------

  interface Init
    module procedure ncb_Init
  end interface

  interface Done
    module procedure ncb_Done
  end interface

  interface ReadRecord
    module procedure ncb_ReadRecord_2d
    module procedure ncb_ReadRecord_3d
  end interface


contains


  ! =======================================================================


  subroutine ncb_Init( ncb, fname, status )

    ! --- in/out ----------------------------------------

    type(TNcepCrayBin), intent(out)   ::  ncb
    character(len=*), intent(in)      ::  fname
    integer, intent(out)              ::  status

    ! --- const ------------------------------------------

    character(len=*), parameter ::  rname = mname//'/ncb_Init'

    ! --- local ------------------------------------------

    logical               ::  exist, opened
    integer               ::  i, sig0, ext0
    integer(ncb_ikind)    ::  idum
    real(ncb_rkind)       ::  rec2(lrec2)

    ! --- begin ------------------------------------------

    ! store file name:
    ncb%fname = fname

    ! file exist ?
    inquire( file=trim(ncb%fname), exist=exist )
    if ( .not. exist ) then
      write (*,'("ERROR - file not found :")')
      write (*,'("ERROR -   ",a)') trim(ncb%fname)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! select free file unit:
    ncb%fu = 1234
    do
      inquire( unit=ncb%fu, opened=opened )
      if ( .not. opened ) exit
      ncb%fu = ncb%fu + 1
    end do

    ! open file:
#ifdef __INTEL_COMPILER
    open( ncb%fu, file=trim(ncb%fname), status='old', action='read', &
                form='unformatted', convert='big_endian', iostat=status )
#else
    open( ncb%fu, file=trim(ncb%fname), status='old', action='read', &
                form='unformatted', iostat=status )
#endif
    if ( status /= 0 ) then
      write (*,'("ERROR - opening file :")')
      write (*,'("ERROR -   ",a)') trim(ncb%fname)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! read header:
    read (ncb%fu,iostat=status) ncb%head
    if ( status /= 0 ) then
      write (*,'("ERROR - reading header from file :")')
      write (*,'("ERROR -   ",a)') trim(ncb%fname)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! read record 2:
    read (ncb%fu,iostat=status) rec2
    if ( status /= 0 ) then
      write (*,'("ERROR - reading record 2 from file :")')
      write (*,'("ERROR -   ",a)') trim(ncb%fname)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! base indices of sigma info and extra info:
    sig0 = 5
    ext0 = sig0 + maxlev+1 + maxlev

    ! number of levels:
    ncb%nlev = int(rec2(ext0+2))

    ! extract times
    ncb%fcst_hr  = nint(rec2(1))
    ncb%idate(4) = transfer(rec2(2),idum)   ! hour
    ncb%idate(2) = transfer(rec2(3),idum)   ! month
    ncb%idate(3) = transfer(rec2(4),idum)   ! day
    ncb%idate(1) = transfer(rec2(5),idum)   ! year

    ! extract sigma half levels
    allocate( ncb%sigma_half(ncb%nlev+1) )
    ncb%sigma_half = rec2(sig0+1:sig0+ncb%nlev+1)

    ! extract sigma full levels:
    allocate( ncb%sigma_full(ncb%nlev) )
    ncb%sigma_full = rec2(sig0+(ncb%nlev+1)+1:sig0+(ncb%nlev+1)+ncb%nlev)

    ! extract spectral truncation; compute number of complex coeff:
    ncb%shT = int(rec2(ext0+1))
    ncb%shn = ( ncb%shT + 1 ) * ( ncb%shT + 2 ) / 2

    ! ok
    status = 0

  end subroutine ncb_Init


  ! ***


  subroutine ncb_Done( ncb, status )

    ! --- in/out ----------------------------------------

    type(TNcepCrayBin), intent(inout)   ::  ncb
    integer, intent(out)                ::  status

    ! --- const ------------------------------------------

    character(len=*), parameter ::  rname = mname//'/ncb_Done'

    ! --- begin ------------------------------------------

    ! clear
    deallocate( ncb%sigma_half )
    deallocate( ncb%sigma_full )

    ! close file:
    close( ncb%fu, iostat=status )
    if ( status /= 0 ) then
      write (*,'("ERROR - closing file :")')
      write (*,'("ERROR -   ",a)') trim(ncb%fname)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! ok
    status = 0

  end subroutine ncb_Done


  ! ***


  subroutine ncb_ReadRecord_2d( ncb, paramkey, sh, status )

    ! --- in/out ----------------------------------------

    type(TNcepCrayBin), intent(inout)   ::  ncb
    character(len=*), intent(in)        ::  paramkey
    complex, intent(out)                ::  sh(:)
    integer, intent(out)                ::  status

    ! --- const ------------------------------------------

    character(len=*), parameter ::  rname = mname//'/ncb_ReadRecord_2d'

    ! --- local ------------------------------------------

    integer                          ::  irec
    integer                          ::  i
    complex(ncb_rkind), allocatable  ::  rec_c(:)

    ! --- begin ------------------------------------------

    ! record 2 has been read

    ! setup record buffer:
    allocate( rec_c(ncb%shn  ) )

    ! determine record number:
    select case ( paramkey )
      case ( 'oro'  ) ; irec = 3
      case ( 'LNSP' ) ; irec = 4
      case default
        write (*,'("ERROR - do not know record number for param ",a)') paramkey
        write (*,'("ERROR in ",a)') rname; status=1; return
    end select

    ! loop over records, last read record is target:
    do i = 3, irec
      ! read coeff:
      read (ncb%fu,iostat=status) rec_c
      if ( status /= 0 ) then
        write (*,'("ERROR - reading record from file :")')
        write (*,'("ERROR -   file   : ",a )') trim(ncb%fname)
        write (*,'("ERROR -   record : ",i4)') i
        write (*,'("ERROR in ",a)') rname; status=1; return
      end if
    end do

    ! set output array (kind conversion ?)
    sh = rec_c

    ! clear
    deallocate( rec_c )

    ! ok
    status = 0

  end subroutine ncb_ReadRecord_2d



  ! ***


  subroutine ncb_ReadRecord_3d( ncb, paramkey, sh, status )

    ! --- in/out ----------------------------------------

    type(TNcepCrayBin), intent(inout)   ::  ncb
    character(len=*), intent(in)        ::  paramkey
    complex, intent(out)                ::  sh(:,:)
    integer, intent(out)                ::  status

    ! --- const ------------------------------------------

    character(len=*), parameter ::  rname = mname//'/ncb_ReadRecord_3d'

    ! --- local ------------------------------------------

    integer                          ::  irec
    integer                          ::  i
    complex(ncb_rkind), allocatable  ::  rec_c(:)

    ! --- begin ------------------------------------------

    ! record 2 has been read

    ! setup record buffer:
    allocate( rec_c(ncb%shn  ) )

    ! different level ordering:
    !   D, VO  :  D1, VO1, D2, VO2, ...
    !   other  :  T1, T2, ...
    !
    select case ( paramkey )

      case ( 'Tv'  )

        ! skip oro, lnsp
        irec = 2
        do i = 1, 2
          irec = irec + 1
          ! read coeff:
          read (ncb%fu,iostat=status) rec_c
          if ( status /= 0 ) then
            write (*,'("ERROR - reading record from file :")')
            write (*,'("ERROR -   file   : ",a )') trim(ncb%fname)
            write (*,'("ERROR -   record : ",i4)') irec
            write (*,'("ERROR in ",a)') rname; status=1; return
          end if
        end do

        ! read virtual temperature
        do i = 1, ncb%nlev
          ! read coeff:
          irec = irec + 1
          read (ncb%fu,iostat=status) rec_c
          if ( status /= 0 ) then
            write (*,'("ERROR - reading record from file :")')
            write (*,'("ERROR -   file   : ",a )') trim(ncb%fname)
            write (*,'("ERROR -   record : ",i4)') irec
            write (*,'("ERROR in ",a)') rname; status=1; return
          end if
          ! store:
          sh(:,i) = rec_c
        end do

      case ( 'D', 'VO'  )

        ! skip oro, lnsp, temperature
        irec = 2
        do i = 1, 2+ncb%nlev
          irec = irec + 1
          ! read coeff:
          read (ncb%fu,iostat=status) rec_c
          if ( status /= 0 ) then
            write (*,'("ERROR - reading record from file :")')
            write (*,'("ERROR -   file   : ",a )') trim(ncb%fname)
            write (*,'("ERROR -   record : ",i4)') irec
            write (*,'("ERROR in ",a)') rname; status=1; return
          end if
        end do

        ! read pairs D/VO
        do i = 1, ncb%nlev
          ! read coeff D:
          irec = irec + 1
          read (ncb%fu,iostat=status) rec_c
          if ( status /= 0 ) then
            write (*,'("ERROR - reading record from file :")')
            write (*,'("ERROR -   file   : ",a )') trim(ncb%fname)
            write (*,'("ERROR -   record : ",i4)') irec
            write (*,'("ERROR in ",a)') rname; status=1; return
          end if
          ! store ?
          if ( paramkey == 'D' ) sh(:,i) = rec_c
          ! read coeff VO:
          irec = irec + 1
          read (ncb%fu,iostat=status) rec_c
          if ( status /= 0 ) then
            write (*,'("ERROR - reading record from file :")')
            write (*,'("ERROR -   file   : ",a )') trim(ncb%fname)
            write (*,'("ERROR -   record : ",i4)') irec
            write (*,'("ERROR in ",a)') rname; status=1; return
          end if
          ! store ?
          if ( paramkey == 'VO' ) sh(:,i) = rec_c
        end do

      case ( 'Q'  )

        ! skip oro, lnsp, Tv, VO/D
        irec = 2
        do i = 1, 2 + 3*ncb%nlev
          irec = irec + 1
          ! read coeff:
          read (ncb%fu,iostat=status) rec_c
          if ( status /= 0 ) then
            write (*,'("ERROR - reading record from file :")')
            write (*,'("ERROR -   file   : ",a )') trim(ncb%fname)
            write (*,'("ERROR -   record : ",i4)') irec
            write (*,'("ERROR in ",a)') rname; status=1; return
          end if
        end do

        ! read humid:
        do i = 1, ncb%nlev
          ! read coeff:
          irec = irec + 1
          read (ncb%fu,iostat=status) rec_c
          if ( status /= 0 ) then
            write (*,'("ERROR - reading record from file :")')
            write (*,'("ERROR -   file   : ",a )') trim(ncb%fname)
            write (*,'("ERROR -   record : ",i4)') irec
            write (*,'("ERROR in ",a)') rname; status=1; return
          end if
          ! store:
          sh(:,i) = rec_c
        end do

      case default
        write (*,'("ERROR - unsupported param ",a)') paramkey
        write (*,'("ERROR in ",a)') rname; status=1; return
    end select

    ! clear
    deallocate( rec_c )

    ! ok
    status = 0

  end subroutine ncb_ReadRecord_3d



end module file_ncb


