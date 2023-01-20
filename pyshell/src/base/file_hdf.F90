!
!ProTeX: 1.14-AJS
!
!BOI
!
! !TITLE:        File_HDF - Interface to HDF4 files
! !AUTHORS:      Arjo Segers, Roeland van Oss, Gijs van Soest, Jos van Geffen
! !AFFILIATION:  KNMI
! !DATE:         \today
!
! !INTRODUCTION: Introduction
!
!   The routines facilitate reading/writing of integer or real arrays.
!   The following type/kinds are possible:
!
!   \begin{tabular}{cc}
!     integer(1)  &           \\
!     integer(2)  &           \\
!     integer(4)  &  real(4)  \\
!     integer(8)  &  real(8)  \\
!   \end{tabular}
!
!   Subroutines convert where possible between types and/or kinds:
!   \begin{itemize}
!     \item
!       integer/real data from a hdf file is converted to the kind of
!       the output array provided to reading routines;
!     \item
!       integer data might be read into a real array;
!     \item
!       integer/real arrays might be written in an other kind than that
!       of the original array.
!   \end{itemize}
!
!
! !INTRODUCTION: Compile and link
!
!   The module includes the file "hdf.f90" from the HDF4 distribution.
!   The include directory with this file should be specified in an
!   argument of the 'configure' script:
!   \bv
!     ./configure --includes=/usr/local/include
!   \ev
!
!   The library should be linked together with the default hdf libraries.
!
!
! !INTRODUCTION: Usage
!
! \bv
!
!   use file_hdf
!
!   HDF files
!   ---------
!
!     type(THdfFile)    ::  hdf
!
!     call Init( hdf, '/data/test.hdf', 'read'|'create'|'write', status )
!
!     call GetInfo( hdf, status [,num_datasets=] [,num_global_attrs=] )
!
!     call Done( hdf, status )
!
!
!   Scientific Data Sets
!   --------------------
!
!     type(TSds)        ::  sds
!
!
!     !
!     ! *** basic access
!     !
!
!     call Init( sds, status )
!
!     !
!     ! select sds in hdf file;
!     !   index  : sets numbered 0,..,num_datasets-1
!     !
!     call Select( sds, hdf, index, status )
!
!     !
!     ! extract properties
!     !
!     call   GetInfo( sds, status [,name=] [,data_rank=] [,data_dims=] [,data_type=] [,num_attrs=]           )
!
!     !
!     ! check sds params:
!     !   name  :  case independent
!     !
!     !   status : on input  : if non-zero, do not print error messages
!     !            on output : <0 check failed, =0 check success, >0 error
!     !
!     call CheckInfo( sds, status [,name=] [,data_rank=] [,data_dims=] [,data_type=] [,num_attrs=] )
!
!     call Done( sds, status )
!
!
!     !
!     ! *** init for reading
!     !
!
!     call Init( sds, hdf, 'name', status )
!
!     call ReadData( sds, x(:,:,..), status [,start=(/0,0,../)] [,stride=(/1,1,../)] )
!
!     call Done( sds, status )
!
!
!     !
!     ! *** init for creation
!     !
!
!     call Init( sds, hdf, 'name', shape(x), <typekey>, status [,knd=] [,bits=] )
!
!     ! The last value of the shape array might be the
!     ! hdf constant 'SD_UNLIMITED' to specify an unlimited dimension.
!     !
!     ! The <typekey> with optional knd or bits specify
!     ! how the data is storred:
!     !
!     !   'int8'     |  'integer(1)'  |  ('int'  |'integer'), (bits=8 |knd=1)
!     !   'int16'    |  'integer(2)'  |  ('int'  |'integer'), (bits=16|knd=2)
!     !   'int32'    |  'integer(4)'  |  ('int'  |'integer'), (bits=32|knd=4)
!     !   'int64'    |  'integer(8)'  |  ('int'  |'integer'), (bits=64|knd=8)
!     !   'float32'  |  'real(4)'     |  ('float'|'real'   ), (bits=32|knd=4)
!     !   'float64'  |  'real(8)'     |  ('float'|'real'   ), (bits=64|knd=8)
!     !
!     !   'char'
!     !
!
!     call Compress( sds, 'none'   , status )
!     call Compress( sds, 'rle'    , status )
!     call Compress( sds, 'skphuff', status [,skip_size=1..] )
!     call Compress( sds, 'deflate', status [,deflate_level=0..9] )
!
!     call WriteData( sds, x, status [,start=(/0,0,../)] [,stride=(/1,1,../)] )
!
!     call Done( sds, status )
!
!
!     !
!     ! *** init for creation with unlimited dimension
!     !
!
!     call Init( sds, hdf, 'name', (/2,3,..,SD_UNLIMITED/), &
!                                          <typekey>, status [,knd=] [,bits=] )
!
!     call WriteData( sds, x(:,:,..,1), status, start=(/0,0,..,0/) [,stride=(/1,1,../)] )
!     call WriteData( sds, x(:,:,..,1), status, start=(/0,0,..,1/) [,stride=(/1,1,../)] )
!                                                      :
!
!     call Done( sds, status )
!
!
!   SDS dimensions
!   --------------
!
!     Simple usage:
!
!       call SetDim( sds, dim_index, 'name', 'unit', (/../), status [,knd=] )
!
!     or using the expert routines:
!
!       type(TSdsDim)  ::  dim
!       call Init( dim, status )
!       call Select( dim, sds, dim_index, status )   ! 0,..,n-1
!       call SetName( dim, 'name'. status )
!       call SetScale( dim, (/../), status [,knd=] )
!       call Done( dim, status )
!
!
!   Attributes
!   ----------
!
!     !
!     ! select attribute index of certain attribute;
!     ! status /= 0 if not found
!     !
!     call FindAttribute( hdf|sds|dim, 'name', attr_index, status )
!
!     !
!     ! get attribute params:
!     !   attr_index       :  sets numbered 0,..,num_attrs-1
!     !   data_type        :  internal hdf number, see hdf manual
!     !   data_type_descr  :  'i', 'r', 's'
!     !                       (no 'l', since logicals are storred as integers)
!     !
!     call   GetAttributeInfo( hdf|sds|dim, attr_index, status &
!                                [,name=] &
!                                [,data_type=] [,data_type_descr=c] &
!                                [,n_values=] )
!
!     !
!     ! check attribute params:
!     !   name        :  case independent
!     !   attr_index  :  sets numbered 0,..,num_attrs-1
!     !
!     !   status : on input  : if zero, print error messages;
!     !            on output : <0 check failed, =0 check success, >0 error
!     !
!     call CheckAttributeInfo( hdf|sds|dim, attr_index, status [,name=] [,data_type=] [,n_values=] )
!
!     call  ReadAttribute( hdf|sds|dim, 'name', i, status )
!     call  ReadAttribute( hdf|sds|dim, 'name', r, status )
!     call  ReadAttribute( hdf|sds|dim, 'name', l, status )
!     call  ReadAttribute( hdf|sds|dim, 'name', s, status )
!
!     !
!     !   status : on input  : if non-zero, do not print error messages
!     !            on output : <0 check failed, =0 check success, >0 error
!     !
!     call CheckAttribute( hdf|sds|dim, 'name', 1     ,status )
!     call CheckAttribute( hdf|sds|dim, 'name', 1.0   ,status )
!     call CheckAttribute( hdf|sds|dim, 'name', .true.,status )
!     call CheckAttribute( hdf|sds|dim, 'name', 's'   ,status )
!
!     call WriteAttribute( hdf|sds|dim, 'name', 1     ,status [,knd=] )
!     call WriteAttribute( hdf|sds|dim, 'name', 1.0   ,status [,knd=] )
!     call WriteAttribute( hdf|sds|dim, 'name', .true.,status         )
!     call WriteAttribute( hdf|sds|dim, 'name', 's'   ,status         )
!
!
!   Time series
!   -----------
!
!     Write time series of n-D arrays to a hdf file.
!
!     type(TTimeSeriesHDF)   :: F
!
!     ! Open file, specify maximum number of data sets:
!     call Init( F, 'test.hdf', 40, status )
!
!     ! Add new record for field.
!     ! If field does not exist yet, it is created with the specified:
!     !  o unit as attribute
!     !  o type key : 'real(4)', 'integer(2)', etc
!     !  o maximum shape
!     ! For arrays, the maximum size should be provided.
!
!     ! first calls initialise the data sets and fill first temporal record:
!     call AddRecord( F, 'psurf' , 'comment', 'hPa', 'real(4)',            1000.0, status )
!     call AddRecord( F, 'press' , 'comment', 'hPa', 'real(4)', (/60/)   , pp(:) , status )
!     call AddRecord( F, 'kernel', 'comment', 'c/c', 'real(4)', (/20,20/), A(:,:), status )
!     call AddRecord( F, 'key', 5, 'ab123', status )
!
!     ! extra calls add records:
!     call AddRecord( F, 'press' , 'comment', 'hPa', 'real(4)', (/60/)   , pp(:) , status )
!     call AddRecord( F, 'press' , 'comment', 'hPa', 'real(4)', (/60/)   , pp(:) , status )
!
!     ! Close file:
!     call Done( F, status )
!
! \ev
!
!
! !INTRODUCTION: Examples
!
! \bv
!
!   ! Extract array from existing hdf file.
!   ! Sds is identified by name 'aaa' .
!   ! Attribute 'version' should be 'v1.23' .
!   ! Error messages are issued in case of any error.
!   type(THdfFile)      ::  hdf
!   type(TSds)          ::  sds
!   real                ::  x(10,20)
!   call Init( hdf, 'test.hdf', 'read', status )
!   call Init( sds, hdf, 'aaa', status )
!   call CheckAttribute( sds, 'version', 'v1.23', status )
!   call ReadData( sds, x, status )
!   call Done( sds, status )
!   call Done( hdf, status )
!
!   ! Extract array from existing hdf file.
!   ! Sds is identified by name 'aaa' and a time attribure.
!   ! Error messages are issued in case of any error.
!   type(THdfFile)      ::  hdf
!   type(TSds)          ::  sds
!   real                ::  x(10,20)
!   integer             ::  n, k
!   integer             ::  status
!   call Init( hdf, 'test.hdf', 'read', status )
!   call Init( sds, status )
!   call GetInfo( hdf, status, num_datasets=n )
!   k = 0
!   do
!     k = k + 1
!     if ( k > n ) stop 'sds not found'
!     call Select( sds, hdf, k-1, status )     ! <-- sets numbered 0..n-1 !
!     call CheckInfo( sds, name='aaa', status )
!     if ( status /= 0 ) cycle
!     call CheckAttribute( sds, 'time', (/2000,1,1,0,0,0/), status )
!     if ( status /= 0 ) cycle
!     call ReadData( sds, x, status )
!   end do
!   call Done( sds, status )
!   call Done( hdf, status )
!
!   ! Write array to a hdf file.
!   type(THdfFile)      ::  hdf
!   type(TSds)          ::  sds
!   real                ::  x(10,20)
!   x = 0.0
!   call Init( hdf, 'test.hdf', 'create', status )
!   call Init( sds, hdf, 'aaa', shape(x), 'real(4)', status )
!   call SetDim( sds, 0, 'lon', 'deg', lons, status )
!   call SetDim( sds, 1, 'lat', 'deg', lats, status )
!   call Compress( sds, 'deflate', status, deflate_level=6 )
!   call WriteData( sds, x, status )
!   call Done( sds, status )
!   call Done( hdf, status )
!
!
! \ev
!
! !INTRODUCTION: Source files
!
! \bv
!
!   configure            : script to generate source files and Makefile
!
!   file_hdf_base.f90    : definition of types, basic routines
!
!   file_hdf_iwp.F90.in  : template for routines dealing with integer
!                          variables of different kinds in files named:
!                            file_hdf_i4.f90
!                            file_hdf_i8.f90
!                          etc; the template contains keys '<wp>' that are
!                          replaced by '4', '8' etc via a sed command in
!                          the makefile
!
!   file_hdf_rwp.F90.in  : template for routines dealing with real
!                          variables of different kinds
!
!   file_hdf_l.f90       : routines dealing with logical variables
!
!   file_hdf_s.f90       : routines dealing with charcter string variables
!
!   file_hdf.f90.in      : Template for the main module that collects all other;
!                          comments '!i<wp>!' etc are removed if the
!                          corresponding kind specific modules are generated
!                          Also contains types and routines for time series.
!
!
! \ev
!
! !INTRODUCTION: See also
!
! \begin{itemize}
!   \item \htmladdnormallink{HDF Documentation}{http://hdf.ncsa.uiuc.edu/hdf4.html}
!  \end{itemize}
!
!
! \ev
!
!EOI

module file_hdf

  use os_specs, only : MAX_FILENAME_LEN

  use file_hdf_base, only : THdfFile, TSds, TSdsDim
  use file_hdf_base, only : SD_UNLIMITED
  use file_hdf_base, only : Init, Done
  use file_hdf_base, only : Select
  use file_hdf_base, only : GetInfo, CheckInfo
  use file_hdf_base, only : Compress
  use file_hdf_base, only : SetName
  use file_hdf_base, only : FindDataSet
  use file_hdf_base, only : FindAttribute, GetAttributeInfo, CheckAttributeInfo
  use file_hdf_i1, only : ReadAttribute, CheckAttribute, WriteAttribute, ReadData, WriteData, SetScale, SetDim
  use file_hdf_i2, only : ReadAttribute, CheckAttribute, WriteAttribute, ReadData, WriteData, SetScale, SetDim
  use file_hdf_i4, only : ReadAttribute, CheckAttribute, WriteAttribute, ReadData, WriteData, SetScale, SetDim
  use file_hdf_i8, only : ReadAttribute, CheckAttribute, WriteAttribute, ReadData, WriteData, SetScale, SetDim
  use file_hdf_r4, only : ReadAttribute, CheckAttribute, WriteAttribute, ReadData, WriteData, SetScale, SetDim
  use file_hdf_r8, only : ReadAttribute, CheckAttribute, WriteAttribute, ReadData, WriteData, SetScale, SetDim
  use file_hdf_l, only : ReadAttribute, CheckAttribute, WriteAttribute
  use file_hdf_s, only : ReadAttribute, CheckAttribute, WriteAttribute, ReadData, WriteData

  implicit none

  ! --- in/out --------------------------

  private

  public   ::  SD_UNLIMITED

  public   ::  THdfFile
  public   ::  TSds
  public   ::  TSdsDim

  public   ::  Init, Done
  public   ::  GetInfo, CheckInfo

  public   ::  Select
  public   ::  FindDataSet
  public   ::  ReadData
  public   ::  Compress, WriteData

  public   ::  SetName, SetScale, SetDim

  public   ::  FindAttribute
  public   ::  GetAttributeInfo, CheckAttributeInfo
  public   ::  ReadAttribute, CheckAttribute
  public   ::  WriteAttribute

  public   ::  TTimeSeriesHDF, TTimeSeriesSDS
  public   ::  AddRecord

  ! --- const ------------------------------

  ! module name
  character(len=*), parameter  ::  mname = 'file_hdf'

  ! default compression of hdf data sets
  character(len=*), parameter  ::  compression = 'none'
  !character(len=*), parameter  ::  compression = 'deflate'   ! <--- errors ...

  ! --- types -------------------------------

  type TTimeSeriesSDS
    type(TSds)             ::  sds
    character(len=30)      ::  name
    integer, pointer       ::  shp(:)
    integer                ::  istart
  end type TTimeSeriesSDS

  type TTimeSeriesHDF
    character(len=MAX_FILENAME_LEN) ::  fname
    type(THdfFile)                  ::  hdf
    integer                         ::  nfield
    integer                         ::  mfield
    type(TTimeSeriesSDS), pointer   ::  tss(:)
  end type TTimeSeriesHDF


  ! --- interfaces -------------------------

  interface Init
    module procedure tsh_Init
  end interface

  interface Done
    module procedure tsh_Done
  end interface

  interface AddRecord
    module procedure tsh_AddRecord_i
    module procedure tsh_AddRecord_i1
    module procedure tsh_AddRecord_i2
    module procedure tsh_AddRecord_i3
    module procedure tsh_AddRecord_r
    module procedure tsh_AddRecord_r1
    module procedure tsh_AddRecord_r2
    module procedure tsh_AddRecord_r3
    module procedure tsh_AddRecord_r4
    module procedure tsh_AddRecord_s
  end interface


contains


  ! ===========================================================


  subroutine tsh_Init( F, fname, mfield, status )

    ! --- in/out ------------------------------------

    type(TTimeSeriesHDF), intent(out)        ::  F
    character(len=*), intent(in)             ::  fname
    integer, intent(in)                      ::  mfield
    integer, intent(out)                     ::  status

    ! --- const ------------------------------------

    character(len=*), parameter  ::  rname = mname//'/tsh_Init'

    ! --- local ------------------------------------

    integer          ::  i

    ! --- begin -------------------------------------

    ! store file name
    F%fname = fname

    ! open new hdf file
    call Init( F%hdf, trim(F%fname), 'create', status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! set maximum nuber of fields
    if ( mfield <= 0 ) then
      write (*,'("ERROR - strange argument for maximum number of fields:")')
      write (*,'("ERROR -   mfield  : ",i6)') mfield
      write (*,'("ERROR -   file    : ",a)') trim(F%fname)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if
    F%mfield = mfield

    ! init fields
    allocate( F%tss(F%mfield) )
    do i = 1, F%mfield
      nullify( F%tss(i)%shp )
    end do

    ! no fields defined yet
    F%nfield = 0

    ! ok
    status = 0

  end subroutine tsh_Init


  ! ***


  subroutine tsh_Done( F, status )

    ! --- in/out ------------------------------------

    type(TTimeSeriesHDF), intent(inout)  ::  F
    integer, intent(out)                 ::  status

    ! --- const ------------------------------------

    character(len=*), parameter  ::  rname = mname//'/tsh_Done'

    ! --- local ------------------------------------

    integer              ::  i

    ! --- begin -------------------------------------

    ! close fields
    do i = 1, F%nfield
      ! close set
      call Done( F%tss(i)%sds, status )
      if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
      ! destroy arrays
      deallocate( F%tss(i)%shp )
    end do

    ! deallocate
    deallocate( F%tss )

    ! close hdf file
    call Done( F%hdf, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! ok
    status = 0

  end subroutine tsh_Done


  ! ***


  subroutine SearchField( F, name, comment, unit, shp, shp_x, typekey, k, status )

    ! --- in/out -----------------------

    type(TTimeSeriesHDF), intent(inout)   ::  F
    character(len=*), intent(in)          ::  name
    character(len=*), intent(in)          ::  comment
    character(len=*), intent(in)          ::  unit
    integer, intent(in)                   ::  shp(:)
    integer, intent(in)                   ::  shp_x(:)
    character(len=*), intent(in)          ::  typekey
    integer, intent(out)                  ::  k
    integer, intent(out)                  ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/SearchField'

    ! --- local --------------------------------------

    integer         ::  i

    ! --- local ---------------------------------------

    ! not found by default ...
    k = -1

    ! loop over current fields
    do i = 1, F%nfield
      ! this field ? then leave
      if ( F%tss(i)%name == name ) then
        k = i
        exit
      end if
    end do

    ! name ok ? then check shape
    if ( k > 0 ) then
      if ( size(shp_x) /= size(F%tss(k)%shp) ) then
        write (*,'("ERROR - rank of record is not the same as rank of open data set:")')
        write (*,'("ERROR -   data set name   : ",a)') F%tss(k)%name
        write (*,'("ERROR -   data set shape  : ",7i4)') F%tss(k)%shp
        write (*,'("ERROR -   requested shape : ",7i4)') shp
        write (*,'("ERROR in ",a)') rname; status=1; return
      end if
      if ( any( shp_x > F%tss(k)%shp ) ) then
        write (*,'("ERROR - max shape of record exceeds shape of input array:")')
        write (*,'("ERROR -   data set name   : ",a)') F%tss(k)%name
        write (*,'("ERROR -   data set shape  : ",7i4)') F%tss(k)%shp
        write (*,'("ERROR -   array shape     : ",7i4)') shp_x
        write (*,'("ERROR in ",a)') rname; status=1; return
      end if
      ! found and ok
      status=0; return
    end if

    ! not found, thus new number
    if ( F%nfield == F%mfield ) then
      write (*,'("ERROR - unable to add new field after ",i4)') F%mfield
      write (*,'("ERROR -   please increase argument `mfield` to Init procedure.")')
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if
    F%nfield = F%nfield + 1
    k = F%nfield

    ! init new field
    ! * store name
    F%tss(k)%name = name
    ! * init data set
    call Init( F%tss(k)%sds, F%hdf, name, shp, typekey, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    !call Compress( F%tss(k)%sds, compression=compression )
    ! * comment
    if ( len_trim(comment) > 0 ) then
      call WriteAttribute( F%tss(k)%sds, 'comment', comment, status )
      if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    end if
    ! * store unit ?
    if ( len_trim(unit) > 0 ) then
      call WriteAttribute( F%tss(k)%sds, 'unit', unit, status )
      if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    end if
    ! * store shape
    allocate( F%tss(k)%shp(size(shp)) )
    F%tss(k)%shp = shp
    ! * init record counter
    F%tss(k)%istart = 0

    ! ok
    status = 0

  end subroutine SearchField


  ! ****************************************************************
  ! *** integer
  ! ****************************************************************


  subroutine tsh_AddRecord_i( F, name, comment, unit, typekey, x, status )

    ! --- in/out -----------------------

    type(TTimeSeriesHDF), intent(inout)   ::  F
    character(len=*), intent(in)          ::  name
    character(len=*), intent(in)          ::  comment
    character(len=*), intent(in)          ::  unit
    character(len=*), intent(in)          ::  typekey
    integer, intent(in)                   ::  x
    integer, intent(out)                  ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/tsh_AddRecord_i'

    ! --- local --------------------------------------

    integer         ::  k
    integer         ::  start(1)

    ! --- local ---------------------------------------

    ! search index
    call SearchField( F, name, comment, unit, (/SD_UNLIMITED/), (/SD_UNLIMITED/), typekey, k, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! write record
    start = (/F%tss(k)%istart/)
    call WriteData( F%tss(k)%sds, (/x/), status, start=start )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! increase record base:
    F%tss(k)%istart = F%tss(k)%istart + 1

    ! ok
    status = 0

  end subroutine tsh_AddRecord_i


  ! ***


  subroutine tsh_AddRecord_i1( F, name, comment, unit, typekey, shp, x, status )

    use parray, only : pa_Init, pa_SetShape, pa_Done

    ! --- in/out -----------------------

    type(TTimeSeriesHDF), intent(inout)   ::  F
    character(len=*), intent(in)          ::  name
    character(len=*), intent(in)          ::  comment
    character(len=*), intent(in)          ::  unit
    character(len=*), intent(in)          ::  typekey
    integer, intent(in)                   ::  shp(:)
    integer, intent(in)                   ::  x(:)
    integer, intent(out)                  ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/tsh_AddRecord_i1'

    ! --- local --------------------------------------

    integer, pointer  ::  xa(:,:)
    integer           ::  k
    integer           ::  start(2)

    ! --- local ---------------------------------------

    ! search index
    call SearchField( F, name, comment, unit, &
                      (/shp,SD_UNLIMITED/), (/shape(x),SD_UNLIMITED/), &
                      typekey, k, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! write record
    call pa_Init( xa )
    call pa_SetShape( xa, (/shp,1/) )
    xa = 0
    xa(1:size(x),1) = x
    start = (/0,F%tss(k)%istart/)
    call WriteData( F%tss(k)%sds, xa, status, start=start )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    call pa_Done( xa )

    ! increase record base:
    F%tss(k)%istart = F%tss(k)%istart + 1

    ! ok
    status = 0

  end subroutine tsh_AddRecord_i1


  ! ***


  subroutine tsh_AddRecord_i2( F, name, comment, unit, typekey, shp, x, status )

    use parray, only : pa_Init, pa_SetShape, pa_Done

    ! --- in/out -----------------------

    type(TTimeSeriesHDF), intent(inout)   ::  F
    character(len=*), intent(in)          ::  name
    character(len=*), intent(in)          ::  comment
    character(len=*), intent(in)          ::  unit
    character(len=*), intent(in)          ::  typekey
    integer, intent(in)                   ::  shp(:)
    integer, intent(in)                   ::  x(:,:)
    integer, intent(out)                  ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/tsh_AddRecord_i2'

    ! --- local --------------------------------------

    integer, pointer  ::  xa(:,:,:)
    integer           ::  k
    integer           ::  start(3)

    ! --- local ---------------------------------------

    ! search index
    call SearchField( F, name, comment, unit, &
                      (/shp,SD_UNLIMITED/), (/shape(x),SD_UNLIMITED/), &
                      typekey, k, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! write record
    call pa_Init( xa )
    call pa_SetShape( xa, (/shp,1/) )
    xa = 0
    xa(1:size(x,1),1:size(x,2),1) = x
    start = (/0,0,F%tss(k)%istart/)
    call WriteData( F%tss(k)%sds, xa, status, start=start )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    call pa_Done( xa )

    ! increase record base:
    F%tss(k)%istart = F%tss(k)%istart + 1

    ! ok
    status = 0

  end subroutine tsh_AddRecord_i2


  ! ***


  subroutine tsh_AddRecord_i3( F, name, comment, unit, typekey, shp, x, status )

    use parray, only : pa_Init, pa_SetShape, pa_Done

    ! --- in/out -----------------------

    type(TTimeSeriesHDF), intent(inout)   ::  F
    character(len=*), intent(in)          ::  name
    character(len=*), intent(in)          ::  comment
    character(len=*), intent(in)          ::  unit
    character(len=*), intent(in)          ::  typekey
    integer, intent(in)                   ::  shp(:)
    integer, intent(in)                   ::  x(:,:,:)
    integer, intent(out)                  ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/tsh_AddRecord_i3'

    ! --- local --------------------------------------

    integer, pointer  ::  xa(:,:,:,:)
    integer           ::  k
    integer           ::  start(4)

    ! --- local ---------------------------------------

    ! search index
    call SearchField( F, name, comment, unit, &
                      (/shp,SD_UNLIMITED/), (/shape(x),SD_UNLIMITED/), &
                      typekey, k, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! write record
    call pa_Init( xa )
    call pa_SetShape( xa, (/shp,1/) )
    xa = 0
    xa(1:size(x,1),1:size(x,2),1:size(x,3),1) = x
    start = (/0,0,0,F%tss(k)%istart/)
    call WriteData( F%tss(k)%sds, xa, status, start=start )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    call pa_Done( xa )

    ! increase record base:
    F%tss(k)%istart = F%tss(k)%istart + 1

    ! ok
    status = 0

  end subroutine tsh_AddRecord_i3


  ! ****************************************************************
  ! *** real
  ! ****************************************************************


  subroutine tsh_AddRecord_r( F, name, comment, unit, typekey, x, status )

    ! --- in/out -----------------------

    type(TTimeSeriesHDF), intent(inout)   ::  F
    character(len=*), intent(in)          ::  name
    character(len=*), intent(in)          ::  comment
    character(len=*), intent(in)          ::  unit
    character(len=*), intent(in)          ::  typekey
    real, intent(in)                      ::  x
    integer, intent(out)                  ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/tsh_AddRecord_r'

    ! --- local --------------------------------------

    integer         ::  k
    integer         ::  start(1)

    ! --- local ---------------------------------------

    ! search index
    call SearchField( F, name, comment, unit, &
                         (/SD_UNLIMITED/), (/SD_UNLIMITED/), &
                         typekey, k, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! write record
    start = (/F%tss(k)%istart/)
    call WriteData( F%tss(k)%sds, (/x/), status, start=start )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! increase record base:
    F%tss(k)%istart = F%tss(k)%istart + 1

    ! ok
    status = 0

  end subroutine tsh_AddRecord_r


  ! ***


  subroutine tsh_AddRecord_r1( F, name, comment, unit, typekey, shp, x, status )

    use parray, only : pa_Init, pa_SetShape, pa_Done

    ! --- in/out -----------------------

    type(TTimeSeriesHDF), intent(inout)   ::  F
    character(len=*), intent(in)          ::  name
    character(len=*), intent(in)          ::  comment
    character(len=*), intent(in)          ::  unit
    character(len=*), intent(in)          ::  typekey
    integer, intent(in)                   ::  shp(:)
    real, intent(in)                      ::  x(:)
    integer, intent(out)                  ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/tsh_AddRecord_r1'

    ! --- local --------------------------------------

    real, pointer   ::  xa(:,:)
    integer         ::  k
    integer         ::  start(2)

    ! --- local ---------------------------------------

    ! search index
    call SearchField( F, name, comment, unit, &
                         (/shp,SD_UNLIMITED/), (/shape(x),SD_UNLIMITED/), &
                         typekey, k, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! write record
    call pa_Init( xa )
    call pa_SetShape( xa, (/shp,1/) )
    xa = 0.0
    xa(1:size(x),1) = x
    start = (/0,F%tss(k)%istart/)
    call WriteData( F%tss(k)%sds, xa, status, start=start )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    call pa_Done( xa )

    ! increase record base:
    F%tss(k)%istart = F%tss(k)%istart + 1

    ! ok
    status = 0

  end subroutine tsh_AddRecord_r1


  ! ***


  subroutine tsh_AddRecord_r2( F, name, comment, unit, typekey, shp, x, status )

    use parray, only : pa_Init, pa_SetShape, pa_Done

    ! --- in/out -----------------------

    type(TTimeSeriesHDF), intent(inout)   ::  F
    character(len=*), intent(in)          ::  name
    character(len=*), intent(in)          ::  comment
    character(len=*), intent(in)          ::  unit
    character(len=*), intent(in)          ::  typekey
    integer, intent(in)                   ::  shp(:)
    real, intent(in)                      ::  x(:,:)
    integer, intent(out)                  ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/tsh_AddRecord_r2'

    ! --- local --------------------------------------

    real, pointer   ::  xa(:,:,:)
    integer         ::  k
    integer         ::  start(3)

    ! --- local ---------------------------------------

    ! search index
    call SearchField( F, name, comment, unit, &
                         (/shp,SD_UNLIMITED/), (/shape(x),SD_UNLIMITED/), &
                         typekey, k, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! write record
    call pa_Init( xa )
    call pa_SetShape( xa, (/shp,1/) )
    xa = 0.0
    xa(1:size(x,1),1:size(x,2),1) = x
    start = (/0,0,F%tss(k)%istart/)
    call WriteData( F%tss(k)%sds, xa, status, start=start )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    call pa_Done( xa )

    ! increase record base:
    F%tss(k)%istart = F%tss(k)%istart + 1

    ! ok
    status = 0

  end subroutine tsh_AddRecord_r2


  ! ***


  subroutine tsh_AddRecord_r3( F, name, comment, unit, typekey, shp, x, status )

    use parray, only : pa_Init, pa_SetShape, pa_Done

    ! --- in/out -----------------------

    type(TTimeSeriesHDF), intent(inout)   ::  F
    character(len=*), intent(in)          ::  name
    character(len=*), intent(in)          ::  comment
    character(len=*), intent(in)          ::  unit
    character(len=*), intent(in)          ::  typekey
    integer, intent(in)                   ::  shp(:)
    real, intent(in)                      ::  x(:,:,:)
    integer, intent(out)                  ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/tsh_AddRecord_r3'

    ! --- local --------------------------------------

    real, pointer   ::  xa(:,:,:,:)
    integer         ::  k
    integer         ::  start(4)

    ! --- local ---------------------------------------

    ! search index
    call SearchField( F, name, comment, unit, &
                         (/shp,SD_UNLIMITED/), (/shape(x),SD_UNLIMITED/), &
                         typekey, k, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! write record
    call pa_Init( xa )
    call pa_SetShape( xa, (/shp,1/) )
    xa = 0.0
    xa(1:size(x,1),1:size(x,2),1:size(x,3),1) = x
    start = (/0,0,0,F%tss(k)%istart/)
    call WriteData( F%tss(k)%sds, xa, status, start=start )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    call pa_Done( xa )

    ! increase record base:
    F%tss(k)%istart = F%tss(k)%istart + 1

    ! ok
    status = 0

  end subroutine tsh_AddRecord_r3


  ! ***


  subroutine tsh_AddRecord_r4( F, name, comment, unit, typekey, shp, x, status )

    use parray, only : pa_Init, pa_SetShape, pa_Done

    ! --- in/out -----------------------

    type(TTimeSeriesHDF), intent(inout)   ::  F
    character(len=*), intent(in)          ::  name
    character(len=*), intent(in)          ::  comment
    character(len=*), intent(in)          ::  unit
    character(len=*), intent(in)          ::  typekey
    integer, intent(in)                   ::  shp(:)
    real, intent(in)                      ::  x(:,:,:,:)
    integer, intent(out)                  ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/tsh_AddRecord_r4'

    ! --- local --------------------------------------

    real, pointer   ::  xa(:,:,:,:,:)
    integer         ::  k
    integer         ::  start(5)

    ! --- local ---------------------------------------

    ! search index
    call SearchField( F, name, comment, unit, &
                         (/shp,SD_UNLIMITED/), (/shape(x),SD_UNLIMITED/), &
                         typekey, k, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! write record
    call pa_Init( xa )
    call pa_SetShape( xa, (/shp,1/) )
    xa = 0.0
    xa(1:size(x,1),1:size(x,2),1:size(x,3),1:size(x,4),1) = x
    start = (/0,0,0,0,F%tss(k)%istart/)
    call WriteData( F%tss(k)%sds, xa, status, start=start )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    call pa_Done( xa )

    ! increase record base:
    F%tss(k)%istart = F%tss(k)%istart + 1

    ! ok
    status = 0

  end subroutine tsh_AddRecord_r4


  ! ****************************************************************
  ! *** char
  ! ****************************************************************


  subroutine tsh_AddRecord_s( F, name, comment, slen, s, status )

    ! --- in/out -----------------------

    type(TTimeSeriesHDF), intent(inout)   ::  F
    character(len=*), intent(in)          ::  name
    character(len=*), intent(in)          ::  comment
    integer, intent(in)                   ::  slen
    character(len=*), intent(in)          ::  s
    integer, intent(out)                  ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/tsh_AddRecord_s'

    ! --- local --------------------------------------

    integer         ::  k
    integer         ::  start(2)

    ! --- local ---------------------------------------

    ! search index
    call SearchField( F, name, comment, 'none', &
                         (/slen,SD_UNLIMITED/), (/len(s),SD_UNLIMITED/), &
                         'char', k, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! write record
    start = (/0,F%tss(k)%istart/)
    call WriteData( F%tss(k)%sds, (/s/), status, start=start )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! increase record base:
    F%tss(k)%istart = F%tss(k)%istart + 1

    ! ok
    status = 0

  end subroutine tsh_AddRecord_s


end module file_hdf



! ########################################################################
! ###
! ### test program
! ###
! ########################################################################
!
!program test
!
!  use file_TimeSeriesHDF
!
!  type(TTimeSeriesHDF)   ::  F
!  integer                ::  status
!  integer                ::  i
!  real                   ::  a(3,4)
!  real                   ::  p(4)
!
!  print *, 'test: begin'
!
!  print *, 'test: open file ...'
!  call Init( F, 'test.hdf', 30, status )
!  if (status/=0) stop 'ERROR'
!
!  do i = 1, 10
!    call AddRecord( F, 'lon', 'longitude', 'deg', 'real(4)', 5.0*i, status )
!    if (status/=0) stop 'ERROR'
!    call AddRecord( F, 'lat', 'latitude', 'deg', 'real(4)', (/2/), (/2.0*i,3.0*i/), status )
!    if (status/=0) stop 'ERROR'
!  end do
!
!  !call AddRecord( F, 'lat', 'deg', 'real(4)', (/2.0*i,3.0*i,4.0/), status )
!  !if (status/=0) stop 'ERROR'
!
!  do i = 1, 3
!    call AddRecord( F, 'kernel', 'averaging kernel', 'c/c', 'real(4)', (/3,4/), a*0.0+i, status )
!    if (status/=0) stop 'ERROR'
!  end do
!
!  do i = 1, 3
!    call AddRecord( F, 'p', 'pressure', 'hPa', 'real(4)', (/6/), p*0.0+i, status )
!    if (status/=0) stop 'ERROR'
!  end do
!
!  do i = 1, 4
!    call AddRecord( F, 'key', 'record key', 5, 'abc', status )
!    if (status/=0) stop 'ERROR'
!  end do
!
!  print *, 'test: close file ...'
!  call Done( F, status )
!  if (status/=0) stop 'ERROR'
!
!  print *, 'test: end'
!
!end program test

