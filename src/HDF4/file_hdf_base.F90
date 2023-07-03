module file_hdf_base

  use GO, only : gol, goPr, goErr
  use os_specs, only : MAX_FILENAME_LEN

  implicit none

  ! --- in/out --------------------------

  private

  public   ::  wpi

  public   ::  wp_int8, wp_int16, wp_int32, wp_int64
  public   ::  wp_float32, wp_float64

  public   ::  MAX_DATA_RANK
  public   ::  SD_UNLIMITED

  public   ::  THdfFile
  public   ::  TSds
  public   ::  TSdsDim

  public   ::  Init, Done

  public   ::  Defined
  public   ::  Select
  public   ::  GetInfo, CheckInfo
  public   ::  Compress

  public   ::  SetName

  public   ::  FindAttribute
  public   ::  GetAttributeInfo, CheckAttributeInfo

  public   ::  FindDataSet


  ! --- const ----------------------------

  character(len=*), parameter  ::  mname = 'file_hdf_base'

  ! ** hdf constants

  include "hdf.f90"

  ! ** working precision of hdf library

  integer, parameter    ::  wpi = 4

  ! ** working precisions of data

  integer, parameter    ::  wp_int8  = 1
  integer, parameter    ::  wp_int16 = 2
  integer, parameter    ::  wp_int32 = 4
  integer, parameter    ::  wp_int64 = 8

  integer, parameter    ::  wp_float32 = 4
  integer, parameter    ::  wp_float64 = 8

  ! ** maximum array ranks

  integer, parameter    ::  MAX_DATA_RANK = 32


  ! --- types ---------------------------

  ! ~~ scientific data set:

  type TSds
    ! internal id:
    integer(wpi)                  ::  id
    ! hdf file name:
    character(len=MAX_FILENAME_LEN) ::  hdfname
    ! name:
    character(len=64)             ::  name
    ! data specification:
    integer                       ::  dfnt
    character(len=3)              ::  typ
    integer                       ::  knd
    integer                       ::  rnk
    integer                       ::  shp(7)
  end type TSds

  ! ~~ dimension

  type TSdsDim
    ! internal id:
    integer(wpi)                  ::  id
  end type TSdsDim

  ! ~~ hdf file

  type THdfFile
    ! internal id:
    integer(wpi)                  ::  id
    ! file name
    character(len=MAX_FILENAME_LEN) ::  fname
  end type THdfFile


  ! --- interfaces ------------------------

  interface Init
    module procedure sds_Init
    module procedure sds_Init_select
    module procedure sds_Init_create
    module procedure dim_Init
    module procedure hdf_Init
  end interface

  interface Done
    module procedure sds_Done
    module procedure dim_Done
    module procedure hdf_Done
  end interface

  interface Defined
    module procedure sds_Defined
  end interface

  interface Select
    module procedure sds_Select_index
    module procedure dim_Select
  end interface

  interface GetInfo
    module procedure sds_GetInfo
    module procedure hdf_GetInfo
  end interface

  interface CheckInfo
    module procedure sds_CheckInfo
  end interface

  interface Compress
    module procedure sds_Compress
  end interface


  interface SetName
    module procedure dim_SetName
  end interface

  interface FindAttribute
    module procedure obj_FindAttribute
    module procedure sds_FindAttribute
    module procedure hdf_FindAttribute
  end interface

  interface GetAttributeInfo
    module procedure obj_GetAttributeInfo
    module procedure sds_GetAttributeInfo
    module procedure hdf_GetAttributeInfo
  end interface

  interface CheckAttributeInfo
    module procedure obj_CheckAttributeInfo
    module procedure sds_CheckAttributeInfo
    module procedure hdf_CheckAttributeInfo
  end interface

  interface FindDataSet
    module procedure hdf_FindDataSet
  end interface


contains


  ! ############################################################
  ! ###
  ! ### tools
  ! ###
  ! ############################################################


  !
  ! compare character strings case independent
  !

  logical function leq( s1, s2 )

    ! --- in/out ------------------------

    character(len=*), intent(in)  ::  s1, s2

    ! --- local -------------------------

    character(len=2)   ::  cc
    integer            ::  k

    ! --- begin -------------------------

    if ( len_trim(s1) /= len_trim(s2) ) then
      leq = .false.
      return
    end if

    do k = 1, len_trim(s1)

      select case ( s1(k:k) )
        case ( 'A', 'a' ); cc = 'Aa'
        case ( 'B', 'b' ); cc = 'Bb'
        case ( 'C', 'c' ); cc = 'Cc'
        case ( 'D', 'd' ); cc = 'Dd'
        case ( 'E', 'e' ); cc = 'Ee'
        case ( 'F', 'f' ); cc = 'Ff'
        case ( 'G', 'g' ); cc = 'Gg'
        case ( 'H', 'h' ); cc = 'Hh'
        case ( 'I', 'i' ); cc = 'Ii'
        case ( 'J', 'j' ); cc = 'Jj'
        case ( 'K', 'k' ); cc = 'Kk'
        case ( 'L', 'l' ); cc = 'Ll'
        case ( 'M', 'm' ); cc = 'Mm'
        case ( 'N', 'n' ); cc = 'Nn'
        case ( 'O', 'o' ); cc = 'Oo'
        case ( 'P', 'p' ); cc = 'Pp'
        case ( 'Q', 'q' ); cc = 'Qq'
        case ( 'R', 'r' ); cc = 'Rr'
        case ( 'S', 's' ); cc = 'Ss'
        case ( 'T', 't' ); cc = 'Tt'
        case ( 'U', 'u' ); cc = 'Uu'
        case ( 'V', 'v' ); cc = 'Vv'
        case ( 'W', 'w' ); cc = 'Ww'
        case ( 'X', 'x' ); cc = 'Xx'
        case ( 'Y', 'y' ); cc = 'Yy'
        case ( 'Z', 'z' ); cc = 'Zz'
        case default; cc = '**'
      end select

      if ( cc == '**' ) then
        if ( s2(k:k) /= s1(k:k) ) then
          leq = .false.
          return
        end if
      else
        if ( s2(k:k) /= cc(1:1) .and. s2(k:k) /= cc(2:2) ) then
          leq = .false.
          return
        end if
      end if

    end do

    leq = .true.

  end function leq


  ! ############################################################
  ! ###
  ! ### objects
  ! ###
  ! ############################################################


  subroutine obj_FindAttribute( obj_id, name, attr_index, status )

    ! --- in/out -------------------------

    integer(wpi), intent(in)                  ::  obj_id
    character(len=*), intent(in)              ::  name
    integer, intent(out)                      ::  attr_index
    integer, intent(inout)                    ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/obj_FindAttribute'

    ! --- local -------------------------------

    integer              ::  istat
    logical              ::  verbose

    ! --- external ----------------------------

    integer(wpi), external   ::  sfFAttr

    ! --- begin -------------------------------

    ! write error messages ?
    verbose = status == 0

    ! extract id of attribute:
    attr_index = sfFAttr( obj_id, name )
    if ( attr_index == FAIL ) then
      if ( verbose ) then
        write (gol,'("finding attribute `",a,"`")') trim(name); call goErr
        write (gol,'("in ",a)') rname; call goErr
      end if
      status=-1; return
    end if

    ! ok
    status = 0

  end subroutine obj_FindAttribute


  ! ***

  !
  ! argument attr_index   : 0,..,n-1
  !

  subroutine obj_GetAttributeInfo( obj_id, attr_index, status, &
                            name, data_type, data_type_descr, n_values )

    ! --- in/out -------------------------

    integer(wpi), intent(in)                  ::  obj_id
    integer, intent(in)                       ::  attr_index
    integer, intent(out)                      ::  status

    character(len=*), intent(out), optional   ::  name
    integer, intent(out), optional            ::  data_type
    character(len=1), intent(out), optional   ::  data_type_descr
    integer, intent(out), optional            ::  n_values

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/obj_GetAttributeInfo'

    ! --- local -------------------------------

    character(len=64)    ::  attr_name
    integer              ::  attr_data_type
    integer              ::  attr_n_values

    ! --- external ----------------------------

    integer(wpi), external   ::  sfGAInfo

    ! --- begin -------------------------------

    ! extract info:
    status = sfGAInfo( obj_id, attr_index, attr_name, attr_data_type, attr_n_values )
    if ( status /= SUCCEED ) then
      write (gol,'("getting attribute info")') trim(name); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    ! return values:
    if ( present(name)      ) name      = attr_name

    if ( present(data_type) ) data_type = attr_data_type
    if ( present(data_type_descr) ) then
      select case ( attr_data_type )
        case ( DFNT_INT8, DFNT_INT16, DFNT_INT32,  DFNT_INT64 )
          data_type_descr = 'i'
        case ( DFNT_FLOAT32,  DFNT_FLOAT64 )
          data_type_descr = 'r'
        case ( DFNT_CHAR )
          data_type_descr = 's'
        case default
          write (gol,'("do not know the data type description")'); call goErr
          write (gol,'("  attribute name      : ",a)') trim(attr_name); call goErr
          write (gol,'("  attribute data type : ",i6)') attr_data_type; call goErr
          write (gol,'("in ",a)') rname; call goErr; status=1; return
      end select
    end if

    if ( present(n_values)  ) n_values  = attr_n_values

    ! ok
    status = 0

  end subroutine obj_GetAttributeInfo


  ! ***


  !
  ! argument attr_index   : 0,..,n-1
  !

  subroutine obj_CheckAttributeInfo( obj_id, attr_index, status, &
                                       name, data_type, n_values )

    ! --- in/out -------------------------

    integer(wpi), intent(in)                  ::  obj_id
    integer, intent(in)                       ::  attr_index
    integer, intent(inout)                    ::  status

    character(len=*), intent(in), optional    ::  name
    integer, intent(in), optional             ::  data_type
    integer, intent(in), optional             ::  n_values

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/obj_CheckAttributeInfo'

    ! --- local -------------------------------

    logical              ::  verbose
    character(len=64)    ::  attr_name
    integer              ::  attr_data_type
    integer              ::  attr_n_values

    ! --- begin -------------------------------

    ! write error messages ?
    verbose = status == 0

    ! check name
    if ( present(name) ) then
      call GetAttributeInfo( obj_id, attr_index, status, name=attr_name )
      if (status/=0) then; write (gol,'("in ",a)') rname; call goErr; status=1; return; end if
      if ( .not. leq(attr_name,name) ) then
        if ( verbose ) then
          write (gol,'("found different attribute name :")'); call goErr
          write (gol,'("  requested : ",a)') trim(name); call goErr
          write (gol,'("  found     : ",a)') trim(attr_name); call goErr
        end if
        status=-1; return
      end if
    end if

    ! check data type
    if ( present(data_type) ) then
      call GetAttributeInfo( obj_id, attr_index, status, &
                               data_type=attr_data_type, name=attr_name )
      if (status/=0) then; write (gol,'("in ",a)') rname; call goErr; status=1; return; end if
      if ( attr_data_type /= data_type ) then
        if ( verbose ) then
          write (gol,'("found different data type :")'); call goErr
          write (gol,'("  requested : ",i6)') data_type; call goErr
          write (gol,'("  found     : ",i6)') attr_data_type; call goErr
          write (gol,'("  attribute :")') trim(attr_name); call goErr
        end if
        status=-1; return
      end if
    end if

    ! check number of values:
    if ( present(n_values) ) then
      call GetAttributeInfo( obj_id, attr_index, status, &
                               n_values=attr_n_values, name=attr_name )
      if (status/=0) then; write (gol,'("in ",a)') rname; call goErr; status=1; return; end if
      if ( attr_n_values /= n_values ) then
        if ( verbose ) then
          write (gol,'("found different number of values :")'); call goErr
          write (gol,'("  requested : ")') n_values; call goErr
          write (gol,'("  found     : ")') attr_n_values; call goErr
          write (gol,'("  attribute : ")') trim(attr_name); call goErr
        end if
        status=-1; return
      end if
    end if

    ! ok
    status = 0

  end subroutine obj_CheckAttributeInfo



  ! ############################################################
  ! ###
  ! ### scientific data sets
  ! ###
  ! ############################################################


  ! ================================================================
  ! init, done
  ! ================================================================


  subroutine sds_Init( sds, status )

    ! --- in/out -----------------------------

    type(Tsds), intent(out)           ::  sds
    integer, intent(out)              ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/sds_Init'

    ! --- begin ------------------------------

    ! dummy ...
    sds%hdfname = 'unknown-hdf-file'
    sds%typ = 'xxx'

    ! no id yet
    sds%id = -1

    ! ok
    status = 0

  end subroutine sds_Init


  ! ***


  subroutine sds_Done( sds, status )

    ! --- in/out -----------------------------

    type(Tsds), intent(inout)         ::  sds
    integer, intent(out)              ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/sds_Done'

    ! --- external ----------------------------

    integer(wpi), external   ::  sfEndAcc

    ! --- begin ------------------------------

    if ( sds%id /= -1 ) then
      status = sfEndAcc( sds%id )
      if ( status == FAIL ) then
        write (gol,'("ending scientific data set ",i6)') sds%id; call goErr
        write (gol,'("  hdf file name : ",a)') trim(sds%hdfname); call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
      end if
    end if

    ! ok
    status = 0

  end subroutine sds_Done


  ! ***


  logical function sds_Defined( sds )

    ! --- in/out ------------------------------

    type(TSds), intent(in)    ::  sds

    ! --- begin ------------------------------

    sds_Defined = sds%id /= -1

  end function sds_Defined


  ! ================================================================
  ! === select sds
  ! ================================================================


  subroutine sds_Select_index( sds, hdf, ind, status )

    ! --- in/out -------------------------

    type(TSds), intent(out)            ::  sds
    type(THdfFile), intent(in)         ::  hdf
    integer, intent(in)                ::  ind
    integer, intent(out)               ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/sds_Select_index'

    ! --- external ------------------------

    integer(wpi), external  ::  sfSelect

    ! --- begin ---------------------------

    sds%id = sfSelect( hdf%id, ind )   ! <-- 0,..,n-1
    if ( sds%id == FAIL ) then
      write (gol,'("unable to locate data set with index ",i6)') ind; call goErr
      write (gol,'("  hdf file name : ",a)') trim(hdf%fname); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    ! ok
    status = 0

  end subroutine sds_Select_index


  ! ***


  subroutine sds_Init_select( sds, hdf, name, status )

    ! --- in/out -------------------------

    type(Tsds), intent(out)            ::  sds
    type(THdfFile), intent(inout)      ::  hdf
    character(len=*), intent(in)       ::  name
    integer, intent(inout)             ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/sds_Init_select'

    ! --- local -------------------------------

    integer              ::  sds_index

    ! --- external ------------------------

    integer(wpi), external  ::  sfN2Index

    ! --- begin -------------------------------

    ! default init
    call Init( sds, status )
    if (status/=0) then; write (gol,'("in ",a)') rname; call goErr; status=1; return; end if

    ! fill hdf and sds names:
    sds%hdfname = hdf%fname
    sds%name = name

    ! search for the record
    sds_index = sfN2Index( hdf%id, name )
    if ( sds_index == FAIL ) then
      write (gol,'("converting sds name to index :")'); call goErr
      write (gol,'("  sds name : ",a)') trim(sds%name); call goErr
      write (gol,'("  hdf name : ",a)') trim(sds%hdfname); call goErr
      write (gol,'("in ",a)') rname; call goErr
      status=1; return
    end if

    ! select sds id:
    call Select( sds, hdf, sds_index, status )
    if (status/=0) then; write (gol,'("in ",a)') rname; call goErr; status=1; return; end if

    ! ok
    status = 0

  end subroutine sds_Init_select


  ! =============================================================
  ! === sds info
  ! =============================================================


  subroutine sds_GetInfo( sds, status, &
                          name, data_rank, data_dims, data_type, num_attrs )

    ! --- in/out -------------------------

    type(TSds), intent(in)                    ::  sds
    integer, intent(out)                      ::  status

    character(len=*), intent(out), optional   ::  name
    integer, intent(out), optional            ::  data_rank
    integer, intent(out), optional            ::  data_type
    integer, intent(out), optional            ::  data_dims(:)
    integer, intent(out), optional            ::  num_attrs

    ! --- local -------------------------------

    integer                ::  istat
    character(len=64)      ::  sds_name
    integer                ::  sds_data_rank, sds_data_type
    integer                ::  sds_data_dims(MAX_DATA_RANK)
    integer                ::  sds_num_attrs

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/sds_GetInfo'

    ! --- external ----------------------------

    integer(wpi), external   ::  sfGInfo

    ! --- begin -------------------------------

    ! extract info about record:
    istat = sfGInfo( sds%id, sds_name, sds_data_rank, sds_data_dims, sds_data_type, sds_num_attrs )
    if ( istat /= SUCCEED ) then
      write (gol,'("error getting info")'); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    ! return values:
    if ( present(name)      ) name      = sds_name
    if ( present(data_rank) ) data_rank = sds_data_rank
    if ( present(data_type) ) data_type = sds_data_type
    if ( present(data_dims) ) data_dims = sds_data_dims(1:size(data_dims))
    if ( present(num_attrs) ) num_attrs = sds_num_attrs

    ! ok
    status = 0

  end subroutine sds_GetInfo


  ! ***


  subroutine sds_CheckInfo( sds, status, &
                            name, data_rank, data_dims, data_type, num_attrs )

    ! --- in/out -------------------------

    type(TSds), intent(in)                     ::  sds
    integer, intent(inout)                     ::  status

    character(len=*), intent(in), optional     ::  name
    integer, intent(in), optional              ::  data_rank
    integer, intent(in), optional              ::  data_type
    integer, intent(in), optional              ::  data_dims(:)
    integer, intent(in), optional              ::  num_attrs

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/sds_CheckInfo'

    ! --- local -------------------------------

    logical                ::  verbose
    character(len=64)      ::  sds_name
    integer                ::  sds_data_rank, sds_data_type
    integer, allocatable   ::  sds_data_dims(:)
    integer                ::  sds_num_attrs

    ! --- begin -------------------------------

    ! write error messages ?
    verbose = status == 0

    ! check name
    if ( present(name) ) then
      call GetInfo( sds, status, name=sds_name )
      if (status/=0) then; write (gol,'("in ",a)') rname; call goErr; status=1; return; end if
      if ( .not. leq(sds_name,name) ) then
        if ( verbose ) then
          write (gol,'("found different name :")'); call goErr
          write (gol,'("  requested : ",a)') trim(name); call goErr
          write (gol,'("  found     : ",a)') trim(sds_name); call goErr
          write (gol,'("in ",a)') rname; call goErr
        end if
        status=-1; return
      end if
    end if

    ! check data rank
    if ( present(data_rank) ) then
      call GetInfo( sds, status, data_rank=sds_data_rank )
      if (status/=0) then; write (gol,'("in ",a)') rname; call goErr; status=1; return; end if
      if ( sds_data_rank /= data_rank ) then
        if ( verbose ) then
          write (gol,'("found different data rank :")'); call goErr
          write (gol,'("  requested : ",i6)') data_rank; call goErr
          write (gol,'("  found     : ",i6)') sds_data_rank; call goErr
          write (gol,'("in ",a)') rname; call goErr
        end if
        status=-1; return
      end if
    end if

    ! check data type
    if ( present(data_type) ) then
      call GetInfo( sds, status, data_type=sds_data_type )
      if (status/=0) then; write (gol,'("in ",a)') rname; call goErr; status=1; return; end if
      if ( sds_data_type /= data_type ) then
        if ( verbose ) then
          write (gol,'("found different data type :")'); call goErr
          write (gol,'("  requested : ",i6)') data_type; call goErr
          write (gol,'("  found     : ",i6)') sds_data_type; call goErr
          write (gol,'("in ",a)') rname; call goErr
        end if
        status=-1; return
      end if
    end if

    ! check data dimensions
    if ( present(data_dims) ) then
      allocate( sds_data_dims(size(data_dims)) )
      call GetInfo( sds, status, data_dims=sds_data_dims )
      if (status/=0) then; write (gol,'("in ",a)') rname; call goErr; status=1; return; end if
      if ( any( sds_data_dims /= data_dims ) ) then
        if ( verbose ) then
          write (gol,'("different data dims :")'); call goErr
          write (gol,'("  requested : ",7i4)') data_dims; call goErr
          write (gol,'("  found     : ",7i4)') sds_data_dims; call goErr
          write (gol,'("in ",a)') rname; call goErr
        end if
        status=-1; return
      end if
      deallocate( sds_data_dims )
    end if

    ! check number of attributes:
    if ( present(num_attrs) ) then
      call GetInfo( sds, status, num_attrs=sds_num_attrs )
      if (status/=0) then; write (gol,'("in ",a)') rname; call goErr; status=1; return; end if
      if ( sds_num_attrs /= num_attrs ) then
        if ( verbose ) then
          write (gol,'("different data num_attrs :")'); call goErr
          write (gol,'("  requested : ")') num_attrs; call goErr
          write (gol,'("  found     : ")') sds_num_attrs; call goErr
          write (gol,'("in ",a)') rname; call goErr
        end if
        status=-1; return
      end if
    end if

    ! ok
    status = 0

  end subroutine sds_CheckInfo


  ! =============================================================
  ! === create sds data
  ! =============================================================

  !
  !  'int8'       'integer(1)'    'int'|'integer', bits=8 |knd=1
  !  'int16'      'integer(2)'    'int'|'integer', bits=16|knd=2
  !  'int32'      'integer(4)'    'int'|'integer', bits=32|knd=4
  !  'int64'      'integer(8)'    'int'|'integer', bits=64|knd=8
  !
  !  'float32'    'real(4)'       'float'|'real', bits=32|knd=4
  !  'float64'    'real(8)'       'float'|'real', bits=64|knd=8
  !
  !  'char'
  !

  subroutine sds_Init_create( sds, hdf, name, shp, typekey, status, &
                                   knd, bits )

    ! --- in/out -------------------------

    type(TSds), intent(out)                   ::  sds
    type(THdfFile), intent(inout)             ::  hdf
    character(len=*), intent(in)              ::  name
    integer, intent(in)                       ::  shp(:)
    character(len=*), intent(in)              ::  typekey
    integer, intent(out)                      ::  status

    integer, intent(in), optional             ::  knd
    integer, intent(in), optional             ::  bits

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/sds_Init_create'

    ! --- local -------------------------------

    integer               ::  dfnt
    character(len=3)      ::  dtyp
    integer               ::  dbits, dknd

    ! --- external ----------------------------

    integer(wpi), external  ::  sfCreate

    ! --- begin -------------------------------

    ! default initialisation:
    call Init( sds, status )
    if (status/=0) then; write (gol,'("in ",a)') rname; call goErr; status=1; return; end if

    ! fill hdf file name
    sds%hdfname = hdf%fname

    ! determine data type:
    select case ( typekey )
      case ( 'int8', 'integer(1)' )
        dfnt = DFNT_INT8
        dtyp = 'int'
        dknd = 1
      case ( 'int16', 'integer(2)' )
        dfnt = DFNT_INT16
        dtyp = 'int'
        dknd = 2
      case ( 'int32', 'integer(4)' )
        dfnt = DFNT_INT32
        dtyp = 'int'
        dknd = 4
      case ( 'int64', 'integer(8)' )
        dfnt = DFNT_INT64
        dtyp = 'int'
        dknd = 8
      case ( 'int', 'integer' )
        if ( present(bits) ) then
          dbits = bits
        else if ( present(knd) ) then
          dbits = knd * 8
        else
          dbits = kind(1) * 8
        end if
        select case ( dbits )
          case ( 8 )
            dfnt = DFNT_INT8
            dknd = 1
          case ( 16 )
            dfnt = DFNT_INT16
            dknd = 2
          case ( 32 )
            dfnt = DFNT_INT32
            dknd = 4
          case ( 64 )
            dfnt = DFNT_INT64
            dknd = 8
          case default
            write (gol,'("integer data not implemented for dbits=",i6)') dbits; call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
        end select
        dtyp = 'int'
      case ( 'float32', 'real(4)' )
        dfnt = DFNT_FLOAT32
        dtyp = 'flt'
        dknd = 4
      case ( 'float64', 'real(8)' )
        dfnt = DFNT_FLOAT64
        dtyp = 'flt'
        dknd = 8
      case ( 'float', 'real' )
        if ( present(bits) ) then
          dbits = bits
        else if ( present(knd) ) then
          dbits = knd * 8
        else
          dbits = kind(1) * 8
        end if
        select case ( dbits )
          case ( 32 )
            dfnt = DFNT_FLOAT32
            dknd = 4
          case ( 64 )
            dfnt = DFNT_FLOAT64
            dknd = 8
          case default
            write (gol,'("real data not implemented for dbits=",i6)') dbits; call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
        end select
        dtyp = 'flt'
      case ( 'char' )
        dfnt = DFNT_CHAR
        dtyp = 'chr'
        dknd = 1
      case default
        write (gol,'("typekey not implemented: ",a)') trim(typekey); call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
    end select

    ! store type and kind
    sds%dfnt = dfnt
    sds%typ  = dtyp
    sds%knd  = dknd

    ! store rank
    sds%rnk  = size(shp)
    if ( sds%rnk < 1 .or. sds%rnk > 7 ) then
      write (gol,'("invalid rank : ",i4)') sds%rnk; call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    ! store shape
    sds%shp(1:sds%rnk) = shp

    ! start new record:
    sds%id = sfCreate( hdf%id, name, sds%dfnt, sds%rnk, sds%shp(1:sds%rnk) )
    if ( sds%id == FAIL ) then
      write (gol,'("from sfCreate :")'); call goErr
      write (gol,'("  name     : ",a)') trim(name); call goErr
      write (gol,'("  hdf file : ",a)') trim(sds%hdfname); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    ! ok
    status = 0

  end subroutine sds_Init_create


  ! ================================================================
  ! compression
  ! ================================================================


  subroutine sds_Compress( sds, compression, status, &
                                skip_size, deflate_level )

    ! --- in/out -------------------------

    type(Tsds), intent(inout)           ::  sds
    character(len=*), intent(in)        ::  compression
    integer, intent(out)                ::  status

    integer, intent(in), optional       ::  skip_size
    integer, intent(in), optional       ::  deflate_level

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/sds_Compress'

    ! --- local -------------------------------

    integer                ::  comp_type
    integer                ::  comp_prm(1)

    ! --- external ---------------------------

    integer(wpi), external   :: sfsCompress

    ! --- begin -------------------------------

    ! default compression parameters:
    comp_prm = (/ 0 /)

    ! set compression type and parameters given key:
    select case ( compression )
      case ( 'none'    )
        comp_type = COMP_CODE_NONE
      case ( 'rle'     )   ! run-length encoding
        comp_type = COMP_CODE_RLE
      case ( 'skphuff' )   ! skipping Huffman
        comp_type = COMP_CODE_SKPHUFF
        comp_prm = (/ 1 /)
        if ( present(skip_size) ) comp_prm(1) = skip_size
      case ( 'deflate' )   ! gzip
        comp_type = COMP_CODE_DEFLATE
        comp_prm = (/ 6 /)
        if ( present(deflate_level) ) comp_prm(1) = deflate_level
      case default
        write (gol,'("unknown compression type : ",a)') trim(compression); call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
    end select

    ! call HDF routine:
    status = sfsCompress( sds%id, comp_type, comp_prm )
    if ( status == FAIL ) then
      write (gol,'("from sfsCompress : ")'); call goErr
      write (gol,'("  compression    : ",a )') trim(compression); call goErr
      write (gol,'("  compress type  : ",i6)') comp_type; call goErr
      write (gol,'("  compress param : ",i6)') comp_prm; call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    ! ok
    status = 0

  end subroutine sds_Compress


  ! =============================================================
  ! === sds attributes
  ! =============================================================


  subroutine sds_FindAttribute( sds, name, attr_index, status )

    ! --- in/out -------------------------

    type(TSds), intent(in)                    ::  sds
    character(len=*), intent(in)              ::  name
    integer, intent(out)                      ::  attr_index
    integer, intent(inout)                    ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/sds_FindAttribute'

    ! --- local -------------------------------

    logical              ::  verbose

    ! --- begin -------------------------------

    ! write error messages ?
    verbose = status == 0

    call FindAttribute( sds%id, name, attr_index, status )
    if (status<0) then
      ! not found ..
      if (verbose) then; write (gol,'("in ",a)') rname; call goErr; end if
      status=-1; return
    else if ( status == 0 ) then
      ! ok
      status=0; return
    else
      ! error
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if
  end subroutine sds_FindAttribute


  ! ***


  subroutine sds_GetAttributeInfo( sds, attr_index, status, &
                                     name, data_type, data_type_descr, &
                                     n_values )

    ! --- in/out -------------------------

    type(TSds), intent(in)                    ::  sds
    integer, intent(in)                       ::  attr_index
    integer, intent(out)                      ::  status

    character(len=*), intent(out), optional   ::  name
    integer, intent(out), optional            ::  data_type
    character(len=1), intent(out), optional   ::  data_type_descr
    integer, intent(out), optional            ::  n_values

     ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/sds_GetAttributeInfo'

   ! --- begin -------------------------------

    call GetAttributeInfo( sds%id, attr_index, status, &
                            name=name, &
                            data_type=data_type, data_type_descr=data_type_descr, &
                            n_values=n_values )
    if (status/=0) then; write (gol,'("in ",a)') rname; call goErr; status=1; return; end if

    ! ok
    status = 0

  end subroutine sds_GetAttributeInfo


  ! ***


  subroutine sds_CheckAttributeInfo( sds, attr_index, status, &
                                       name, data_type, n_values )

    ! --- in/out -------------------------

    type(TSds), intent(in)                    ::  sds
    integer, intent(in)                       ::  attr_index
    integer, intent(inout)                    ::  status

    character(len=*), intent(in), optional    ::  name
    integer, intent(in), optional             ::  data_type
    integer, intent(in), optional             ::  n_values

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/sds_CheckAttributeInfo'

    ! --- local -------------------------------

    logical              ::  verbose

    ! --- begin -------------------------------

    ! write error messages ?
    verbose = status == 0

    call CheckAttributeInfo( sds%id, attr_index, &
                              name=name, data_type=data_type, n_values=n_values, &
                              status=status )
    if ( status < 0 ) then
      ! error
      if (verbose) then; write (gol,'("in ",a)') rname; call goErr; end if
      status=-1; return
    else if ( status == 0 ) then
      ! ok
      status=0; return
    else
      ! error
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

  end subroutine sds_CheckAttributeInfo



  ! ############################################################
  ! ###
  ! ### dimensions
  ! ###
  ! ############################################################


  ! ================================================================
  ! init, done
  ! ================================================================


  subroutine dim_Init( sdim, status )

    ! --- in/out ------------------------------

    type(TSdsDim), intent(out)         ::  sdim
    integer, intent(out)               ::  status

    ! --- begin -------------------------------

    sdim%id = FAIL

    ! ok
    status = 0

  end subroutine dim_Init


  ! *


  subroutine dim_Done( sdim, status )

    ! --- in/out ------------------------------

    type(TSdsDim), intent(inout)       ::  sdim
    integer, intent(out)               ::  status

    ! --- begin -------------------------------

    ! nothing to be done

    ! ok
    status = 0

  end subroutine dim_Done


  ! ================================================================
  ! select
  ! ================================================================

  !
  ! argument ind         : 0,..,n-1
  !

  subroutine dim_Select( sdim, sds, ind, status )

    ! --- in/out -------------------------

    type(TSdsDim), intent(out)         ::  sdim
    type(TSds), intent(in)             ::  sds
    integer, intent(in)                ::  ind
    integer, intent(out)               ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/dim_Select'

    ! --- external ------------------------

    integer(wpi), external  ::  sfDimID

    ! --- begin ---------------------------

    sdim%id = sfDimID( sds%id, ind )   ! <-- 0,..,n-1
    if ( sdim%id == FAIL ) then
      write (gol,'("error selecting dimension :")'); call goErr
      write (gol,'("  index    : ",i6)') ind; call goErr
      write (gol,'("  sds name : ",a)') trim(sds%name); call goErr
      write (gol,'("  hdf name : ",a)') trim(sds%hdfname); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    ! ok
    status = 0

  end subroutine dim_Select


  ! ================================================================
  ! set dimension name
  ! ================================================================


  subroutine dim_SetName( sdim, name, status )

    ! --- in/out -------------------------

    type(TSdsDim), intent(inout)              ::  sdim
    character(len=*), intent(in)              ::  name
    integer, intent(out)                      ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/dim_SetName'

    ! --- external ---------------------------

    integer(wpi), external   :: sfSDmName

    ! --- begin ---------------------------

    ! set dimension name
    status = sfSDmName( sdim%id, name )
    if ( status == FAIL ) then
      write (gol,'("setting dimension name :")'); call goErr
      write (gol,'("  dim name : ",a)') name; call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    ! ok
    status = 0

  end subroutine dim_SetName




  ! ############################################################
  ! ###
  ! ### hdf files
  ! ###
  ! ############################################################


  ! ================================================================
  ! init, done
  ! ================================================================


  subroutine hdf_Init( hdf, fname, key, status )

    ! --- in/out ------------------------------

    ! !ARGUMENTS:
    type(THdfFile), intent(out)        ::  hdf
    character(len=*), intent(in)       ::  fname
    character(len=*), intent(in)       ::  key
    integer, intent(out)               ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/hdf_Init'

    ! --- local -------------------------------

    integer          ::  dfacc

    ! --- external ----------------------------

    integer(wpi), external  ::  sfStart

    ! --- begin -------------------------------

    ! code to open file:
    select case ( key )
      case ( 'read' )
        dfacc = DFACC_READ
      case ( 'write' )
        dfacc = DFACC_WRITE
       case ( 'create' )
        dfacc = DFACC_CREATE
     case default
        write (gol,'("do not know what how to access hdf for `",a,"`:")') key; call goErr
        write (gol,'("  file name : ",a)') trim(fname); call goErr
        write (gol,'("in ",a)') rname; call goErr; status=1; return
    end select

    ! open file:
    hdf%id = sfStart( fname, dfacc )
    if ( hdf%id == FAIL ) then
      write (gol,'("from starting access to hdf file:")'); call goErr
      write (gol,'("  file name  : ",a)') trim(fname); call goErr
      write (gol,'("  access key : ",a)') trim(key); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    ! save file name:
    hdf%fname = fname

    ! ok
    status = 0

  end subroutine hdf_Init


  ! ***


  subroutine hdf_Done( hdf, status )

    ! --- in/out ------------------------------

    ! !ARGUMENTS:
    type(THdfFile), intent(out)        ::  hdf
    integer, intent(out)               ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/hdf_Done'

    ! --- external ----------------------------

    integer(wpi), external  ::  sfEnd

    ! --- begin -------------------------------

    ! close file:
    status = sfEnd( hdf%id )
    if ( status == FAIL ) then
      write (gol,'("while closing HDF file:")'); call goErr
      write (gol,'("  file name : ",a)') trim(hdf%fname); call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    ! ok
    status = 0

  end subroutine hdf_Done


  ! ================================================================
  !  info
  ! ================================================================


  subroutine hdf_GetInfo( hdf, status, num_datasets, num_global_attrs )

    ! --- in/out -------------------------

    type(THdfFile), intent(inout)      ::  hdf
    integer, intent(out)               ::  status

    integer, intent(out), optional     ::  num_datasets
    integer, intent(out), optional     ::  num_global_attrs

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/hdf_GetInfo'

    ! --- local -------------------------------

    integer              ::  istat
    integer              ::  f_num_datasets, f_num_global_attrs

    ! --- external ----------------------------

    integer(wpi), external  ::  sfFInfo

    ! --- begin -------------------------------

    ! extract info
    istat = sfFInfo( hdf%id, f_num_datasets, f_num_global_attrs )
    if ( istat == FAIL ) then
      write (gol,'("from sfFInfo :")'); call goErr
      write (gol,'("  hdf file : ",a)') hdf%fname; call goErr
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    ! return result
    if ( present(num_datasets)     ) num_datasets     = f_num_datasets
    if ( present(num_global_attrs) ) num_global_attrs = f_num_global_attrs

    ! ok
    status = 0

  end subroutine hdf_GetInfo


  ! =============================================================
  ! === hdf data sets
  ! =============================================================


  subroutine hdf_FindDataSet( hdf, name, sds_index, status )

    ! --- in/out -------------------------

    type(THdfFile), intent(in)                ::  hdf
    character(len=*), intent(in)              ::  name
    integer, intent(out)                      ::  sds_index
    integer, intent(inout)                    ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/hdf_FindDataSet'

    ! --- external -------------------------------

    integer(wpi), external  ::  sfN2Index

    ! --- local -------------------------------

    logical              ::  verbose

    ! --- begin -------------------------------

    ! write error messages ?
    verbose = status == 0

    ! find index from name:
    sds_index = sfN2Index( hdf%id, name )
    if ( status < 0 ) then
      ! not found ...
      if (verbose) then
        write (gol,'("data set not found ")'); call goErr
        write (gol,'("  name       : ",a)') trim(name)
        write (gol,'("  file name  : ",a)') trim(hdf%fname); call goErr
        write (gol,'("in ",a)') rname; call goErr; status=-1; return
      end if
    else if ( status == 0 ) then
      ! ok
      status=0; return
    else
      ! error ...
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

  end subroutine hdf_FindDataSet


  ! =============================================================
  ! === hdf global attributes
  ! =============================================================


  subroutine hdf_FindAttribute( hdf, name, attr_index, status )

    ! --- in/out -------------------------

    type(THdfFile), intent(in)                ::  hdf
    character(len=*), intent(in)              ::  name
    integer, intent(out)                      ::  attr_index
    integer, intent(inout)                    ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/hdf_FindAttribute'

    ! --- local -------------------------------

    logical              ::  verbose

    ! --- begin -------------------------------

    ! write error messages ?
    verbose = status == 0

    ! find attribute index from name:
    call FindAttribute( hdf%id, name, attr_index, status )
    if ( status < 0 ) then
      ! not found ...
      if (verbose) then; write (gol,'("in ",a)') rname; call goErr; end if
      status=-1; return
    else if ( status == 0 ) then
      ! ok
      status=0; return
    else
      ! error ...
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

  end subroutine hdf_FindAttribute


  ! ***


  subroutine hdf_GetAttributeInfo( hdf, attr_index, status, &
                                     name, &
                                     data_type, data_type_descr, &
                                     n_values )

    ! --- in/out ----------------------------

    type(THdfFile), intent(in)                ::  hdf
    integer, intent(in)                       ::  attr_index
    integer, intent(inout)                    ::  status

    character(len=*), intent(out), optional   ::  name
    integer, intent(out), optional            ::  data_type
    character(len=1), intent(out), optional   ::  data_type_descr
    integer, intent(out), optional            ::  n_values

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/hdf_GetAttributeInfo'

    ! --- begin -------------------------------

    call GetAttributeInfo( hdf%id, attr_index, status, &
                            name=name, &
                            data_type=data_type, data_type_descr=data_type_descr, &
                            n_values=n_values )
    if (status/=0) then; write (gol,'("in ",a)') rname; call goErr; status=1; return; end if

    ! ok
    status = 0

  end subroutine hdf_GetAttributeInfo


  ! ***


  subroutine hdf_CheckAttributeInfo( hdf, attr_index, status, &
                                          name, data_type, n_values )

    ! --- in/out -------------------------

    type(THdfFile), intent(in)                ::  hdf
    integer, intent(in)                       ::  attr_index
    integer, intent(inout)                    ::  status

    character(len=*), intent(in), optional    ::  name
    integer, intent(in), optional             ::  data_type
    integer, intent(in), optional             ::  n_values

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/hdf_CheckAttributeInfo'

    ! --- local -------------------------------

    logical              ::  verbose

    ! --- begin -------------------------------

    ! write error messages ?
    verbose = status == 0

    call CheckAttributeInfo( hdf%id, attr_index, status, &
                              name=name, data_type=data_type, n_values=n_values )
    if ( status < 0 ) then
      ! check failed ...
      if (verbose) then; write (gol,'("in ",a)') rname; call goErr; end if
      status=-1; return
    else if ( status == 0 ) then
      ! ok
      status = 0; return
    else
      ! error ...
      write (gol,'("in ",a)') rname; call goErr; status=1; return
    end if

    ! ok
    status = 0

  end subroutine hdf_CheckAttributeInfo


end module file_hdf_base
