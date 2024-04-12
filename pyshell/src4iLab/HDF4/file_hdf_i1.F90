!
! Template for file_hdf routines with integer arguments.
!
! To generate kind specific versions, use:
!
!    /bin/sed -e 's/1/1/g' file_hdf_iwp.F90.in > file_hdf_i1.F90
!    /bin/sed -e 's/1/2/g' file_hdf_iwp.F90.in > file_hdf_i2.F90
!    /bin/sed -e 's/1/4/g' file_hdf_iwp.F90.in > file_hdf_i4.F90
!    /bin/sed -e 's/1/8/g' file_hdf_iwp.F90.in > file_hdf_i8.F90
!

module file_hdf_i1

  implicit none
  
  ! --- in/out --------------------------
  
  private
  
  public   ::  ReadData
  public   ::  WriteData
  
  public   ::  SetScale, SetDim

  public   ::  ReadAttribute, CheckAttribute
  public   ::  WriteAttribute
  
  ! --- const ----------------------------
  
  include "hdf.f90"
  
  character(len=*), parameter ::  mname = 'file_hdf_i1'
    

  ! --- interfaces ------------------------
  
  interface ReadData
    module procedure sds_ReadData_i1_1d
    module procedure sds_ReadData_i1_2d
    module procedure sds_ReadData_i1_3d
    module procedure sds_ReadData_i1_4d
    module procedure sds_ReadData_i1_5d
    module procedure sds_ReadData_i1_6d
    module procedure sds_ReadData_i1_7d
  end interface

  interface WriteData
    module procedure sds_WriteData_i1_1d
    module procedure sds_WriteData_i1_2d
    module procedure sds_WriteData_i1_3d
    module procedure sds_WriteData_i1_4d
    module procedure sds_WriteData_i1_5d
    module procedure sds_WriteData_i1_6d
    module procedure sds_WriteData_i1_7d
  end interface
  
  interface SetScale
    module procedure dim_SetScale_i1
  end interface
  
  interface SetDim
    module procedure sds_SetDim_i1
  end interface
  
  interface ReadAttribute
    module procedure obj_ReadAttribute_i1_0d
    module procedure obj_ReadAttribute_i1_1d
    !
    module procedure sds_ReadAttribute_i1_0d
    module procedure sds_ReadAttribute_i1_1d
    !
    module procedure dim_ReadAttribute_i1_0d
    module procedure dim_ReadAttribute_i1_1d
    !
    module procedure hdf_ReadAttribute_i1_0d
    module procedure hdf_ReadAttribute_i1_1d
  end interface
  
  interface CheckAttribute
    module procedure obj_CheckAttribute_i1_0d
    module procedure obj_CheckAttribute_i1_1d
    !
    module procedure sds_CheckAttribute_i1_0d
    module procedure sds_CheckAttribute_i1_1d
    !
    module procedure dim_CheckAttribute_i1_0d
    module procedure dim_CheckAttribute_i1_1d
    !
    module procedure hdf_CheckAttribute_i1_0d
    module procedure hdf_CheckAttribute_i1_1d
  end interface
  
  interface WriteAttribute
    module procedure obj_WriteAttribute_i1_0d
    module procedure obj_WriteAttribute_i1_1d
    !
    module procedure sds_WriteAttribute_i1_0d
    module procedure sds_WriteAttribute_i1_1d
    !
    module procedure dim_WriteAttribute_i1_0d
    module procedure dim_WriteAttribute_i1_1d
    !
    module procedure hdf_WriteAttribute_i1_0d
    module procedure hdf_WriteAttribute_i1_1d
  end interface

  
contains


  ! ############################################################
  ! ###
  ! ### objects
  ! ###
  ! ############################################################
  

  
  ! ================================================================
  ! ===
  ! === read attributes
  ! ===
  ! ================================================================
  
  
  subroutine obj_ReadAttribute_i1_0d( obj_id, name, i, status  )
  
    use file_hdf_base, only : wpi
    use file_hdf_base, only : wp_int8, wp_int16, wp_int32, wp_int64
    use file_hdf_base, only : FindAttribute, CheckAttributeInfo, GetAttributeInfo
    
    ! --- in/out -------------------------
    
    integer(wpi), intent(in)            ::  obj_id
    character(len=*), intent(in)        ::  name
    integer(1), intent(out)          ::  i
    integer, intent(out)                ::  status
    
    ! --- const --------------------------
    
    character(len=*), parameter ::  rname = mname//'/obj_ReadAttribute_i1_0d'

    ! --- local -------------------------------
    
    integer              ::  attr_index, data_type
    integer(wp_int8 )    ::  int8
    integer(wp_int16)    ::  int16
    integer(wp_int32)    ::  int32
    integer(wp_int64)    ::  int64
    
    ! --- external ----------------------------

    integer(wpi), external   ::  sfRNAtt

    ! --- begin -------------------------------
    
    call FindAttribute( obj_id, name, attr_index, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    call CheckAttributeInfo( obj_id, attr_index, status, n_values=1 )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    
    ! extract value:
    call GetAttributeInfo( obj_id, attr_index, status, data_type=data_type )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    select case ( data_type )
      case ( DFNT_INT8 )
        status = sfRNAtt( obj_id, attr_index, int8  ); i = int(int8 ,kind=1)
      case ( DFNT_INT16 )
        status = sfRNAtt( obj_id, attr_index, int16 ); i = int(int16,kind=1)
      case ( DFNT_INT32 )
        status = sfRNAtt( obj_id, attr_index, int32 ); i = int(int32,kind=1)
      case ( DFNT_INT64 )
        status = sfRNAtt( obj_id, attr_index, int64 ); i = int(int64,kind=1)
      case default
        write (*,'("ERROR - not implemented for data type ",i6)') data_type
        write (*,'("ERROR in ",a)') rname; status=1; return
    end select
    if ( status /= SUCCEED ) then
      write (*,'("ERROR - reading attribute : ",a)') trim(name)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine obj_ReadAttribute_i1_0d
  

  ! ***
      
  
  subroutine obj_ReadAttribute_i1_1d( obj_id, name, i, status )
  
    use file_hdf_base, only : wpi
    use file_hdf_base, only : wp_int8, wp_int16, wp_int32, wp_int64
    use file_hdf_base, only : FindAttribute, CheckAttributeInfo, GetAttributeInfo
    
    ! --- in/out -------------------------
    
    integer(wpi), intent(in)            ::  obj_id
    character(len=*), intent(in)        ::  name
    integer(1), intent(out)          ::  i(:)
    integer, intent(out)                ::  status
    
    ! --- const --------------------------
    
    character(len=*), parameter ::  rname = mname//'/obj_ReadAttribute_i1_1d'

    ! --- local -------------------------------
    
    integer                           ::  attr_index, data_type
    integer                           ::  n
    integer(wp_int8 ), allocatable    ::  int8 (:)
    integer(wp_int16), allocatable    ::  int16(:)
    integer(wp_int32), allocatable    ::  int32(:)
    integer(wp_int64), allocatable    ::  int64(:)
    
    ! --- external ----------------------------

    integer(wpi), external   ::  sfRNAtt

    ! --- begin -------------------------------
    
    n = size(i)
    
    call FindAttribute( obj_id, name, attr_index, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    call CheckAttributeInfo( obj_id, attr_index, status, n_values=n )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    
    ! extract value:
    call GetAttributeInfo( obj_id, attr_index, status, data_type=data_type )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    select case ( data_type )
      case ( DFNT_INT8 )
        allocate( int8 (n) )
        status = sfRNAtt( obj_id, attr_index, int8  )
        i = int(int8 ,kind=1)
        deallocate( int8  )
      case ( DFNT_INT16 )
        allocate( int16(n) )
        status = sfRNAtt( obj_id, attr_index, int16 )
        i = int(int16,kind=1)
        deallocate( int16 )
      case ( DFNT_INT32 )
        allocate( int32(n) )
        status = sfRNAtt( obj_id, attr_index, int32 )
        i = int(int32,kind=1)
        deallocate( int32 )
      case ( DFNT_INT64 )
        allocate( int64(n) )
        status = sfRNAtt( obj_id, attr_index, int64 )
        i = int(int64,kind=1)
        deallocate( int64 )
      case default
        write (*,'("ERROR - not implemented for data type ",i6)') data_type
        write (*,'("ERROR in ",a)') rname; status=1; return
    end select
    if ( status /= SUCCEED ) then
      write (*,'("ERROR - reading attribute : ",a)') trim(name)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine obj_ReadAttribute_i1_1d
  

  

  ! ================================================================
  ! ===
  ! === check attributes
  ! ===
  ! ================================================================
  
  
  subroutine obj_CheckAttribute_i1_0d( obj_id, name, i, status )
  
    use file_hdf_base, only : wpi

    ! --- in/out -------------------------
    
    integer(wpi), intent(in)            ::  obj_id
    character(len=*), intent(in)        ::  name
    integer(1), intent(in)           ::  i
    integer, intent(inout)              ::  status
    
    ! --- const --------------------------
    
    character(len=*), parameter ::  rname = mname//'/obj_CheckAttribute_i1_0d'

    ! --- local -------------------------------
    
    logical                    ::  verbose
    integer(1)              ::  attr_i
    
    ! --- begin -------------------------------
    
    ! write error messages ?
    verbose = status == 0

    call ReadAttribute( obj_id, name, attr_i, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    
    ! check:
    if ( attr_i /= i ) then
      if (verbose) then
        write (*,'("ERROR - foud different attribute values:")')
        write (*,'("ERROR -   attr name : ",a)') trim(name)
        write (*,'("ERROR -   requested : ",i6)') i
        write (*,'("ERROR -   found     : ",i6)') attr_i
        write (*,'("ERROR in ",a)') rname
      end if
      status=-1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine obj_CheckAttribute_i1_0d
  

  ! ***
      
  
  subroutine obj_CheckAttribute_i1_1d( obj_id, name, i, status )
  
    use file_hdf_base, only : wpi

    ! --- in/out -------------------------
    
    integer(wpi), intent(in)            ::  obj_id
    character(len=*), intent(in)        ::  name
    integer(1), intent(in)           ::  i(:)
    integer, intent(inout)              ::  status
   
    ! --- const --------------------------
    
    character(len=*), parameter ::  rname = mname//'/obj_CheckAttribute_i1_1d'

    ! --- local -------------------------------
    
    logical                      ::  verbose
    integer                      ::  n
    integer(1), allocatable   ::  attr_i(:)
    
    ! --- begin -------------------------------
    
    ! write error messages ?
    verbose = status == 0

    n = size(i)
    allocate( attr_i(n) )
    
    call ReadAttribute( obj_id, name, attr_i, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    
    ! not ok ?
    if ( any( attr_i /= i ) ) then
      if (verbose) then
        write (*,'("ERROR - foud different attribute values:")')
        write (*,'("ERROR -   attr name : ",a)') trim(name)
        write (*,'("ERROR -   requested : ",i6)') i
        write (*,'("ERROR -   found     : ",i6)') attr_i
        write (*,'("ERROR in ",a)') rname
      end if
      deallocate( attr_i )
      status=-1; return
    end if
    
    ! clear
    deallocate( attr_i )
    
    ! ok
    status = 0
    
  end subroutine obj_CheckAttribute_i1_1d

  

  ! ================================================================
  ! ===
  ! === write attributes
  ! ===
  ! ================================================================
  
  
  subroutine obj_WriteAttribute_i1_0d( obj_id, name, i, status, knd )
  
    use file_hdf_base, only : wpi
    use file_hdf_base, only : wp_int8, wp_int16, wp_int32, wp_int64

    ! --- in/out -------------------------
    
    integer(wpi), intent(in)            ::  obj_id
    character(len=*), intent(in)        ::  name
    integer(1), intent(in)           ::  i
    integer, intent(inout)              ::  status

    integer, intent(in), optional       ::  knd
    
    ! --- const --------------------------
    
    character(len=*), parameter ::  rname = mname//'/obj_WriteAttribute_i1_0d'

    ! --- local -------------------------------
    
    integer              ::  ikind

    ! --- external ----------------------------

    integer(wpi), external   ::  sfSNAtt

    ! --- begin -------------------------------

    ikind = kind(i)
    if ( present(knd) ) ikind = knd

    select case ( ikind )
      case ( 1 )
        status = sfSNAtt( obj_id, name, DFNT_INT8 , 1, int(i,kind=wp_int8 ) )
      case ( 2 )
        status = sfSNAtt( obj_id, name, DFNT_INT16, 1, int(i,kind=wp_int16) )
      case ( 4 )
        status = sfSNAtt( obj_id, name, DFNT_INT32, 1, int(i,kind=wp_int32) )
      case ( 8 )
        status = sfSNAtt( obj_id, name, DFNT_INT64, 1, int(i,kind=wp_int64) )
      case default
        write (*,'("ERROR - no implementation for writing with kind ",i2)') ikind
        write (*,'("ERROR in ",a)') rname; status=-1; return
    end select
    if ( status /= SUCCEED ) then
      write (*,'("ERROR - while writing attribute:")')
      write (*,'("ERROR -   attr name   : ",a)') name
      write (*,'("ERROR -   input kind  : ",i2)') kind(i)
      write (*,'("ERROR -   output kind : ",i2)') ikind
      write (*,'("ERROR in ",a)') rname; status=-1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine obj_WriteAttribute_i1_0d
    
  
  ! ***
  
  
  subroutine obj_WriteAttribute_i1_1d( obj_id, name, i, status, knd )
  
    use file_hdf_base, only : wpi
    use file_hdf_base, only : wp_int8, wp_int16, wp_int32, wp_int64

    ! --- in/out -------------------------
    
    integer(wpi), intent(in)            ::  obj_id
    character(len=*), intent(in)        ::  name
    integer(1), intent(in)           ::  i(:)
    integer, intent(inout)              ::  status

    integer, intent(in), optional       ::  knd
    
    ! --- const --------------------------
    
    character(len=*), parameter ::  rname = mname//'/obj_WriteAttribute_i1_1d'

    ! --- local -------------------------------
    
    integer              ::  ikind

    ! --- external ----------------------------

    integer(wpi), external   ::  sfSNAtt

    ! --- begin -------------------------------

    ikind = kind(i)
    if ( present(knd) ) ikind = knd

    select case ( ikind )
      case ( 1 )
        status = sfSNAtt( obj_id, name, DFNT_INT8 , size(i), int(i,kind=wp_int8 ) )
      case ( 2 )
        status = sfSNAtt( obj_id, name, DFNT_INT16, size(i), int(i,kind=wp_int16) )
      case ( 4 )
        status = sfSNAtt( obj_id, name, DFNT_INT32, size(i), int(i,kind=wp_int32) )
      case ( 8 )
        status = sfSNAtt( obj_id, name, DFNT_INT64, size(i), int(i,kind=wp_int64) )
      case default
        write (*,'("ERROR - no implementation for writing with kind ",i2)') ikind
        write (*,'("ERROR in ",a)') rname; status=-1; return
    end select
    if ( status /= SUCCEED ) then
      write (*,'("ERROR - while writing attribute:")')
      write (*,'("ERROR -   attr name   : ",a)') name
      write (*,'("ERROR -   input kind  : ",i2)') kind(i)
      write (*,'("ERROR -   output kind : ",i2)') ikind
      write (*,'("ERROR in ",a)') rname; status=-1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine obj_WriteAttribute_i1_1d
  
  
  
  ! ############################################################
  ! ###
  ! ### scientific data sets
  ! ###
  ! ############################################################
  

  ! ================================================================
  ! read attributes
  ! ================================================================
  
  
  subroutine sds_ReadAttribute_i1_0d( sds, name, i, status )
  
    use file_hdf_base, only : TSds

    ! --- in/out -------------------------
    
    type(Tsds), intent(in)             ::  sds
    character(len=*), intent(in)       ::  name
    integer(1), intent(out)         ::  i
    integer, intent(out)               ::  status
    
    ! --- const --------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sds_ReadAttribute_i1_0d'

    ! --- begin -------------------------------
    
    call ReadAttribute( sds%id, name, i, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
        
    ! ok
    status = 0
    
  end subroutine sds_ReadAttribute_i1_0d
  
  
  ! ***

  
  subroutine sds_ReadAttribute_i1_1d( sds, name, i, status )
  
    use file_hdf_base, only : TSds

    ! --- in/out -------------------------
    
    type(Tsds), intent(in)             ::  sds
    character(len=*), intent(in)       ::  name
    integer(1), intent(out)         ::  i(:)
    integer, intent(out)               ::  status
    
    ! --- const --------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sds_ReadAttribute_i1_1d'

    ! --- begin -------------------------------
    
    call ReadAttribute( sds%id, name, i, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    
    ! ok
    status = 0
    
  end subroutine sds_ReadAttribute_i1_1d
  
  
  ! =============================================================
  ! === check attributes
  ! =============================================================

  
  subroutine sds_CheckAttribute_i1_0d( sds, name, i, status )
  
    use file_hdf_base, only : TSds

    ! --- in/out -------------------------
    
    type(TSds), intent(in)             ::  sds
    character(len=*), intent(in)       ::  name
    integer(1), intent(in)          ::  i
    integer, intent(inout)             ::  status
    
    ! --- const --------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sds_CheckAttribute_i1_0d'

    ! --- local ------------------------------
    
    logical          ::  verbose
    
    ! --- begin ---------------------------
    
    ! write error messages ?
    verbose = status == 0
    
    call CheckAttribute( sds%id, name, i, status )
    if (status>0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    if (status<0) then
      if (verbose) write (*,'("ERROR in ",a)') rname
      status=-1; return
    end if

    ! ok
    status = 0
    
  end subroutine sds_CheckAttribute_i1_0d
  
  
  ! ***
  
  
  subroutine sds_CheckAttribute_i1_1d( sds, name, i, status )
  
    use file_hdf_base, only : TSds

    ! --- in/out -------------------------
    
    type(TSds), intent(in)             ::  sds
    character(len=*), intent(in)       ::  name
    integer(1), intent(in)          ::  i(:)
    integer, intent(inout)             ::  status
    
    ! --- const --------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sds_CheckAttribute_i1_1d'

    ! --- local ------------------------------
    
    logical          ::  verbose
    
    ! --- begin ---------------------------
    
    ! write error messages ?
    verbose = status == 0
    
    call CheckAttribute( sds%id, name, i, status )
    if (status>0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    if (status<0) then
      if (verbose) write (*,'("ERROR in ",a)') rname
      status=-1; return
    end if

    ! ok
    status = 0
    
  end subroutine sds_CheckAttribute_i1_1d
  
  

  ! ================================================================
  ! write attributes
  ! ================================================================
  
  
  subroutine sds_WriteAttribute_i1_0d( sds, name, i, status, knd )
  
    use file_hdf_base, only : TSds

    ! --- in/out -------------------------
    
    type(Tsds), intent(in)             ::  sds
    character(len=*), intent(in)       ::  name
    integer(1), intent(in)          ::  i
    integer, intent(out)               ::  status
    
    integer, intent(in), optional      ::  knd
    
    ! --- const --------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sds_WriteAttribute_i1_0d'

    ! --- begin -------------------------------
    
    call WriteAttribute( sds%id, name, i, status, knd )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
        
    ! ok
    status = 0
    
  end subroutine sds_WriteAttribute_i1_0d
  
  
  ! ***

  
  subroutine sds_WriteAttribute_i1_1d( sds, name, i, status, knd )
  
    use file_hdf_base, only : TSds

    ! --- in/out -------------------------
    
    type(Tsds), intent(in)             ::  sds
    character(len=*), intent(in)       ::  name
    integer(1), intent(in)          ::  i(:)
    integer, intent(out)               ::  status

    integer, intent(in), optional      ::  knd
    
    ! --- const --------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sds_WriteAttribute_i1_1d'

    ! --- begin -------------------------------
    
    call WriteAttribute( sds%id, name, i, status, knd )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    
    ! ok
    status = 0
    
  end subroutine sds_WriteAttribute_i1_1d
  
  
  ! =============================================================
  ! === read data
  ! =============================================================

  
  subroutine sds_ReadData_i1_1d( sds, data, status, start, stride )
  
    use parray, only : pa_Init, pa_SetShape, pa_Done
    use file_hdf_base, only : wpi
    use file_hdf_base, only : TSds
    use file_hdf_base, only : wp_int8, wp_int16, wp_int32, wp_int64
    use file_hdf_base, only : CheckInfo, GetInfo

    ! --- const ------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sds_ReadData_i1_1d'

    integer, parameter  ::  rank = 1
    
    ! --- in/out ----------------------------
    
    type(TSds), intent(in)              ::  sds
    integer(1), intent(out)          ::  data(:)
    integer, intent(out)                ::  status

    integer, intent(in), optional       ::  start(rank)
    integer, intent(in), optional       ::  stride(rank)
    
    ! --- local -------------------------------
    
    integer                             ::  data_type
    integer                             ::  the_start(rank)
    integer                             ::  the_stride(rank)
    integer(wp_int8 ), pointer          ::  int8 (:)
    integer(wp_int16), pointer          ::  int16(:)
    integer(wp_int32), pointer          ::  int32(:)
    integer(wp_int64), pointer          ::  int64(:)
    
    ! --- external ----------------------------
    
    integer(wpi), external  ::  sfRData

    ! --- begin -------------------------------
    
    ! check data rank and shape:
    !call CheckInfo( sds, data_rank=rank, data_dims=shape(data) )

    ! read data of specified kind:
    the_start  = 0; if ( present(start ) ) the_start  = start
    the_stride = 1; if ( present(stride) ) the_stride = stride
    call GetInfo( sds, status, data_type=data_type )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    select case ( data_type )
      case ( DFNT_INT8 )
        call pa_Init( int8  )
        call pa_SetShape( int8 , shape(data) )
        status = sfRData( sds%id, the_start, the_stride, shape(int8 ), int8  )
        data = int(int8,kind=1)
        call pa_Done( int8  )
      case ( DFNT_INT16 )
        call pa_Init( int16 )
        call pa_SetShape( int16, shape(data) )
        status = sfRData( sds%id, the_start, the_stride, shape(int16), int16 )
        data = int(int16,kind=1)
        call pa_Done( int16 )
      case ( DFNT_INT32 )
        call pa_Init( int32 )
        call pa_SetShape( int32, shape(data) )
        status = sfRData( sds%id, the_start, the_stride, shape(int32), int32 )
        data = int(int32,kind=1)
        call pa_Done( int32 )
      case ( DFNT_INT64 )
        call pa_Init( int64 )
        call pa_SetShape( int64, shape(data) )
        status = sfRData( sds%id, the_start, the_stride, shape(int64), int64 )
        data = int(int64,kind=1)
        call pa_Done( int64 )
      case default
        write (*,'("ERROR - not implemented for data type ",i6)') data_type
        write (*,'("ERROR in ",a)') rname; status=1; return
    end select
    if ( status == FAIL ) then
      write (*,'("ERROR - reading data `",a,"`")') trim(sds%name)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine sds_ReadData_i1_1d
  

  ! ***  


  subroutine sds_ReadData_i1_2d( sds, data, status, start, stride )
  
    use parray, only : pa_Init, pa_SetShape, pa_Done
    use file_hdf_base, only : wpi
    use file_hdf_base, only : TSds
    use file_hdf_base, only : wp_int8, wp_int16, wp_int32, wp_int64
    use file_hdf_base, only : CheckInfo, GetInfo

    ! --- const ------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sds_ReadData_i1_2d'

    integer, parameter  ::  rank = 2
    
    ! --- in/out ----------------------------
    
    type(TSds), intent(in)              ::  sds
    integer(1), intent(out)          ::  data(:,:)
    integer, intent(out)                ::  status

    integer, intent(in), optional       ::  start(rank)
    integer, intent(in), optional       ::  stride(rank)
    
    ! --- local -------------------------------
    
    integer                             ::  data_type
    integer                             ::  the_start(rank)
    integer                             ::  the_stride(rank)
    integer(wp_int8 ), pointer          ::  int8 (:,:)
    integer(wp_int16), pointer          ::  int16(:,:)
    integer(wp_int32), pointer          ::  int32(:,:)
    integer(wp_int64), pointer          ::  int64(:,:)
    
    ! --- external ----------------------------
    
    integer(wpi), external  ::  sfRData

    ! --- begin -------------------------------
    
    ! check data rank and shape:
    !call CheckInfo( sds, data_rank=rank, data_dims=shape(data) )

    ! read data of specified kind:
    the_start  = 0; if ( present(start ) ) the_start  = start
    the_stride = 1; if ( present(stride) ) the_stride = stride
    call GetInfo( sds, status, data_type=data_type )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    select case ( data_type )
      case ( DFNT_INT8 )
        call pa_Init( int8  )
        call pa_SetShape( int8 , shape(data) )
        status = sfRData( sds%id, the_start, the_stride, shape(int8 ), int8  )
        data = int(int8,kind=1)
        call pa_Done( int8  )
      case ( DFNT_INT16 )
        call pa_Init( int16 )
        call pa_SetShape( int16, shape(data) )
        status = sfRData( sds%id, the_start, the_stride, shape(int16), int16 )
        data = int(int16,kind=1)
        call pa_Done( int16 )
      case ( DFNT_INT32 )
        call pa_Init( int32 )
        call pa_SetShape( int32, shape(data) )
        status = sfRData( sds%id, the_start, the_stride, shape(int32), int32 )
        data = int(int32,kind=1)
        call pa_Done( int32 )
      case ( DFNT_INT64 )
        call pa_Init( int64 )
        call pa_SetShape( int64, shape(data) )
        status = sfRData( sds%id, the_start, the_stride, shape(int64), int64 )
        data = int(int64,kind=1)
        call pa_Done( int64 )
      case default
        write (*,'("ERROR - not implemented for data type ",i6)') data_type
        write (*,'("ERROR in ",a)') rname; status=1; return
    end select
    if ( status == FAIL ) then
      write (*,'("ERROR - reading data `",a,"`")') trim(sds%name)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine sds_ReadData_i1_2d
  

  ! ***  


  subroutine sds_ReadData_i1_3d( sds, data, status, start, stride )
  
    use parray, only : pa_Init, pa_SetShape, pa_Done
    use file_hdf_base, only : wpi
    use file_hdf_base, only : TSds
    use file_hdf_base, only : wp_int8, wp_int16, wp_int32, wp_int64
    use file_hdf_base, only : CheckInfo, GetInfo

    ! --- const ------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sds_ReadData_i1_3d'

    integer, parameter  ::  rank = 3
    
    ! --- in/out ----------------------------
    
    type(TSds), intent(in)              ::  sds
    integer(1), intent(out)          ::  data(:,:,:)
    integer, intent(out)                ::  status

    integer, intent(in), optional       ::  start(rank)
    integer, intent(in), optional       ::  stride(rank)
    
    ! --- local -------------------------------
    
    integer                             ::  data_type
    integer                             ::  the_start(rank)
    integer                             ::  the_stride(rank)
    integer(wp_int8 ), pointer          ::  int8 (:,:,:)
    integer(wp_int16), pointer          ::  int16(:,:,:)
    integer(wp_int32), pointer          ::  int32(:,:,:)
    integer(wp_int64), pointer          ::  int64(:,:,:)
    
    ! --- external ----------------------------
    
    integer(wpi), external  ::  sfRData

    ! --- begin -------------------------------
    
    ! check data rank and shape:
    !call CheckInfo( sds, data_rank=rank, data_dims=shape(data) )

    ! read data of specified kind:
    the_start  = 0; if ( present(start ) ) the_start  = start
    the_stride = 1; if ( present(stride) ) the_stride = stride
    call GetInfo( sds, status, data_type=data_type )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    select case ( data_type )
      case ( DFNT_INT8 )
        call pa_Init( int8  )
        call pa_SetShape( int8 , shape(data) )
        status = sfRData( sds%id, the_start, the_stride, shape(int8 ), int8  )
        data = int(int8,kind=1)
        call pa_Done( int8  )
      case ( DFNT_INT16 )
        call pa_Init( int16 )
        call pa_SetShape( int16, shape(data) )
        status = sfRData( sds%id, the_start, the_stride, shape(int16), int16 )
        data = int(int16,kind=1)
        call pa_Done( int16 )
      case ( DFNT_INT32 )
        call pa_Init( int32 )
        call pa_SetShape( int32, shape(data) )
        status = sfRData( sds%id, the_start, the_stride, shape(int32), int32 )
        data = int(int32,kind=1)
        call pa_Done( int32 )
      case ( DFNT_INT64 )
        call pa_Init( int64 )
        call pa_SetShape( int64, shape(data) )
        status = sfRData( sds%id, the_start, the_stride, shape(int64), int64 )
        data = int(int64,kind=1)
        call pa_Done( int64 )
      case default
        write (*,'("ERROR - not implemented for data type ",i6)') data_type
        write (*,'("ERROR in ",a)') rname; status=1; return
    end select
    if ( status == FAIL ) then
      write (*,'("ERROR - reading data `",a,"`")') trim(sds%name)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine sds_ReadData_i1_3d
  

  ! ***  


  subroutine sds_ReadData_i1_4d( sds, data, status, start, stride )
  
    use parray, only : pa_Init, pa_SetShape, pa_Done
    use file_hdf_base, only : wpi
    use file_hdf_base, only : TSds
    use file_hdf_base, only : wp_int8, wp_int16, wp_int32, wp_int64
    use file_hdf_base, only : CheckInfo, GetInfo

    ! --- const ------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sds_ReadData_i1_4d'

    integer, parameter  ::  rank = 4
    
    ! --- in/out ----------------------------
    
    type(TSds), intent(in)              ::  sds
    integer(1), intent(out)          ::  data(:,:,:,:)
    integer, intent(out)                ::  status

    integer, intent(in), optional       ::  start(rank)
    integer, intent(in), optional       ::  stride(rank)
    
    ! --- local -------------------------------
    
    integer                             ::  data_type
    integer                             ::  the_start(rank)
    integer                             ::  the_stride(rank)
    integer(wp_int8 ), pointer          ::  int8 (:,:,:,:)
    integer(wp_int16), pointer          ::  int16(:,:,:,:)
    integer(wp_int32), pointer          ::  int32(:,:,:,:)
    integer(wp_int64), pointer          ::  int64(:,:,:,:)
    
    ! --- external ----------------------------
    
    integer(wpi), external  ::  sfRData

    ! --- begin -------------------------------
    
    ! check data rank and shape:
    !call CheckInfo( sds, data_rank=rank, data_dims=shape(data) )

    ! read data of specified kind:
    the_start  = 0; if ( present(start ) ) the_start  = start
    the_stride = 1; if ( present(stride) ) the_stride = stride
    call GetInfo( sds, status, data_type=data_type )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    select case ( data_type )
      case ( DFNT_INT8 )
        call pa_Init( int8  )
        call pa_SetShape( int8 , shape(data) )
        status = sfRData( sds%id, the_start, the_stride, shape(int8 ), int8  )
        data = int(int8,kind=1)
        call pa_Done( int8  )
      case ( DFNT_INT16 )
        call pa_Init( int16 )
        call pa_SetShape( int16, shape(data) )
        status = sfRData( sds%id, the_start, the_stride, shape(int16), int16 )
        data = int(int16,kind=1)
        call pa_Done( int16 )
      case ( DFNT_INT32 )
        call pa_Init( int32 )
        call pa_SetShape( int32, shape(data) )
        status = sfRData( sds%id, the_start, the_stride, shape(int32), int32 )
        data = int(int32,kind=1)
        call pa_Done( int32 )
      case ( DFNT_INT64 )
        call pa_Init( int64 )
        call pa_SetShape( int64, shape(data) )
        status = sfRData( sds%id, the_start, the_stride, shape(int64), int64 )
        data = int(int64,kind=1)
        call pa_Done( int64 )
      case default
        write (*,'("ERROR - not implemented for data type ",i6)') data_type
        write (*,'("ERROR in ",a)') rname; status=1; return
    end select
    if ( status == FAIL ) then
      write (*,'("ERROR - reading data `",a,"`")') trim(sds%name)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine sds_ReadData_i1_4d
  

  ! ***  


  subroutine sds_ReadData_i1_5d( sds, data, status, start, stride )
  
    use parray, only : pa_Init, pa_SetShape, pa_Done
    use file_hdf_base, only : wpi
    use file_hdf_base, only : TSds
    use file_hdf_base, only : wp_int8, wp_int16, wp_int32, wp_int64
    use file_hdf_base, only : CheckInfo, GetInfo

    ! --- const ------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sds_ReadData_i1_5d'

    integer, parameter  ::  rank = 5
    
    ! --- in/out ----------------------------
    
    type(TSds), intent(in)              ::  sds
    integer(1), intent(out)          ::  data(:,:,:,:,:)
    integer, intent(out)                ::  status

    integer, intent(in), optional       ::  start(rank)
    integer, intent(in), optional       ::  stride(rank)
    
    ! --- local -------------------------------
    
    integer                             ::  data_type
    integer                             ::  the_start(rank)
    integer                             ::  the_stride(rank)
    integer(wp_int8 ), pointer          ::  int8 (:,:,:,:,:)
    integer(wp_int16), pointer          ::  int16(:,:,:,:,:)
    integer(wp_int32), pointer          ::  int32(:,:,:,:,:)
    integer(wp_int64), pointer          ::  int64(:,:,:,:,:)
    
    ! --- external ----------------------------
    
    integer(wpi), external  ::  sfRData

    ! --- begin -------------------------------
    
    ! check data rank and shape:
    !call CheckInfo( sds, data_rank=rank, data_dims=shape(data) )

    ! read data of specified kind:
    the_start  = 0; if ( present(start ) ) the_start  = start
    the_stride = 1; if ( present(stride) ) the_stride = stride
    call GetInfo( sds, status, data_type=data_type )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    select case ( data_type )
      case ( DFNT_INT8 )
        call pa_Init( int8  )
        call pa_SetShape( int8 , shape(data) )
        status = sfRData( sds%id, the_start, the_stride, shape(int8 ), int8  )
        data = int(int8,kind=1)
        call pa_Done( int8  )
      case ( DFNT_INT16 )
        call pa_Init( int16 )
        call pa_SetShape( int16, shape(data) )
        status = sfRData( sds%id, the_start, the_stride, shape(int16), int16 )
        data = int(int16,kind=1)
        call pa_Done( int16 )
      case ( DFNT_INT32 )
        call pa_Init( int32 )
        call pa_SetShape( int32, shape(data) )
        status = sfRData( sds%id, the_start, the_stride, shape(int32), int32 )
        data = int(int32,kind=1)
        call pa_Done( int32 )
      case ( DFNT_INT64 )
        call pa_Init( int64 )
        call pa_SetShape( int64, shape(data) )
        status = sfRData( sds%id, the_start, the_stride, shape(int64), int64 )
        data = int(int64,kind=1)
        call pa_Done( int64 )
      case default
        write (*,'("ERROR - not implemented for data type ",i6)') data_type
        write (*,'("ERROR in ",a)') rname; status=1; return
    end select
    if ( status == FAIL ) then
      write (*,'("ERROR - reading data `",a,"`")') trim(sds%name)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine sds_ReadData_i1_5d
  

  ! ***  


  subroutine sds_ReadData_i1_6d( sds, data, status, start, stride )
  
    use parray, only : pa_Init, pa_SetShape, pa_Done
    use file_hdf_base, only : wpi
    use file_hdf_base, only : TSds
    use file_hdf_base, only : wp_int8, wp_int16, wp_int32, wp_int64
    use file_hdf_base, only : CheckInfo, GetInfo

    ! --- const ------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sds_ReadData_i1_6d'

    integer, parameter  ::  rank = 6
    
    ! --- in/out ----------------------------
    
    type(TSds), intent(in)              ::  sds
    integer(1), intent(out)          ::  data(:,:,:,:,:,:)
    integer, intent(out)                ::  status

    integer, intent(in), optional       ::  start(rank)
    integer, intent(in), optional       ::  stride(rank)
    
    ! --- local -------------------------------
    
    integer                             ::  data_type
    integer                             ::  the_start(rank)
    integer                             ::  the_stride(rank)
    integer(wp_int8 ), pointer          ::  int8 (:,:,:,:,:,:)
    integer(wp_int16), pointer          ::  int16(:,:,:,:,:,:)
    integer(wp_int32), pointer          ::  int32(:,:,:,:,:,:)
    integer(wp_int64), pointer          ::  int64(:,:,:,:,:,:)
    
    ! --- external ----------------------------
    
    integer(wpi), external  ::  sfRData

    ! --- begin -------------------------------
    
    ! check data rank and shape:
    !call CheckInfo( sds, data_rank=rank, data_dims=shape(data) )

    ! read data of specified kind:
    the_start  = 0; if ( present(start ) ) the_start  = start
    the_stride = 1; if ( present(stride) ) the_stride = stride
    call GetInfo( sds, status, data_type=data_type )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    select case ( data_type )
      case ( DFNT_INT8 )
        call pa_Init( int8  )
        call pa_SetShape( int8 , shape(data) )
        status = sfRData( sds%id, the_start, the_stride, shape(int8 ), int8  )
        data = int(int8,kind=1)
        call pa_Done( int8  )
      case ( DFNT_INT16 )
        call pa_Init( int16 )
        call pa_SetShape( int16, shape(data) )
        status = sfRData( sds%id, the_start, the_stride, shape(int16), int16 )
        data = int(int16,kind=1)
        call pa_Done( int16 )
      case ( DFNT_INT32 )
        call pa_Init( int32 )
        call pa_SetShape( int32, shape(data) )
        status = sfRData( sds%id, the_start, the_stride, shape(int32), int32 )
        data = int(int32,kind=1)
        call pa_Done( int32 )
      case ( DFNT_INT64 )
        call pa_Init( int64 )
        call pa_SetShape( int64, shape(data) )
        status = sfRData( sds%id, the_start, the_stride, shape(int64), int64 )
        data = int(int64,kind=1)
        call pa_Done( int64 )
      case default
        write (*,'("ERROR - not implemented for data type ",i6)') data_type
        write (*,'("ERROR in ",a)') rname; status=1; return
    end select
    if ( status == FAIL ) then
      write (*,'("ERROR - reading data `",a,"`")') trim(sds%name)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine sds_ReadData_i1_6d
  

  ! ***  


  subroutine sds_ReadData_i1_7d( sds, data, status, start, stride )
  
    use parray, only : pa_Init, pa_SetShape, pa_Done
    use file_hdf_base, only : wpi
    use file_hdf_base, only : TSds
    use file_hdf_base, only : wp_int8, wp_int16, wp_int32, wp_int64
    use file_hdf_base, only : CheckInfo, GetInfo

    ! --- const ------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sds_ReadData_i1_7d'

    integer, parameter  ::  rank = 7
    
    ! --- in/out ----------------------------
    
    type(TSds), intent(in)              ::  sds
    integer(1), intent(out)          ::  data(:,:,:,:,:,:,:)
    integer, intent(out)                ::  status

    integer, intent(in), optional       ::  start(rank)
    integer, intent(in), optional       ::  stride(rank)
    
    ! --- local -------------------------------
    
    integer                             ::  data_type
    integer                             ::  the_start(rank)
    integer                             ::  the_stride(rank)
    integer(wp_int8 ), pointer          ::  int8 (:,:,:,:,:,:,:)
    integer(wp_int16), pointer          ::  int16(:,:,:,:,:,:,:)
    integer(wp_int32), pointer          ::  int32(:,:,:,:,:,:,:)
    integer(wp_int64), pointer          ::  int64(:,:,:,:,:,:,:)
    
    ! --- external ----------------------------
    
    integer(wpi), external  ::  sfRData

    ! --- begin -------------------------------
    
    ! check data rank and shape:
    !call CheckInfo( sds, data_rank=rank, data_dims=shape(data) )

    ! read data of specified kind:
    the_start  = 0; if ( present(start ) ) the_start  = start
    the_stride = 1; if ( present(stride) ) the_stride = stride
    call GetInfo( sds, status, data_type=data_type )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    select case ( data_type )
      case ( DFNT_INT8 )
        call pa_Init( int8  )
        call pa_SetShape( int8 , shape(data) )
        status = sfRData( sds%id, the_start, the_stride, shape(int8 ), int8  )
        data = int(int8,kind=1)
        call pa_Done( int8  )
      case ( DFNT_INT16 )
        call pa_Init( int16 )
        call pa_SetShape( int16, shape(data) )
        status = sfRData( sds%id, the_start, the_stride, shape(int16), int16 )
        data = int(int16,kind=1)
        call pa_Done( int16 )
      case ( DFNT_INT32 )
        call pa_Init( int32 )
        call pa_SetShape( int32, shape(data) )
        status = sfRData( sds%id, the_start, the_stride, shape(int32), int32 )
        data = int(int32,kind=1)
        call pa_Done( int32 )
      case ( DFNT_INT64 )
        call pa_Init( int64 )
        call pa_SetShape( int64, shape(data) )
        status = sfRData( sds%id, the_start, the_stride, shape(int64), int64 )
        data = int(int64,kind=1)
        call pa_Done( int64 )
      case default
        write (*,'("ERROR - not implemented for data type ",i6)') data_type
        write (*,'("ERROR in ",a)') rname; status=1; return
    end select
    if ( status == FAIL ) then
      write (*,'("ERROR - reading data `",a,"`")') trim(sds%name)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine sds_ReadData_i1_7d
  

  ! =============================================================
  ! === Write data
  ! =============================================================


  subroutine sds_WriteData_i1_1d( sds, data, status, start, stride )
  
    use parray, only : pa_Init, pa_SetShape, pa_Done
    use file_hdf_base, only : wpi
    use file_hdf_base, only : TSds
    use file_hdf_base, only : wp_int8, wp_int16, wp_int32, wp_int64
    use file_hdf_base, only : wp_float32, wp_float64

    ! --- const --------------------------
    
    character(len=*), parameter ::  rname = mname//'/sds_WriteData_i1_1d'

    integer, parameter   ::  rank = 1
  
    ! --- in/out -------------------------
    
    type(TSds), intent(in)                    ::  sds
    integer(1), intent(in)                 ::  data(:)
    integer, intent(out)                      ::  status
    
    integer, intent(in), optional             ::  start(rank)
    integer, intent(in), optional             ::  stride(rank)
    
    ! --- local -------------------------------
    
    integer               ::  strt(7), strd(7)
    integer               ::  shap(7), ns
    
    integer(wp_int8 ), pointer          ::  int8 (:)
    integer(wp_int16), pointer          ::  int16(:)
    integer(wp_int32), pointer          ::  int32(:)
    integer(wp_int64), pointer          ::  int64(:)
    real(wp_float32), pointer           ::  float32(:)
    real(wp_float64), pointer           ::  float64(:)
    
    ! --- external ----------------------------
    
    integer(wpi), external  ::  sfWData
    
    ! --- begin -------------------------------
    
    !! check shape
    !if ( any( shape(data) /= sds%shp(1:sds%rnk) ) ) then
    !  print *, 'Shape of data does not match shape specified during creation:'
    !  print *, '  data        : ', shape(data)
    !  print *, '  created for : ', sds%shp(1:sds%rnk)
    !  stop 'FATAL ERROR IN sds_WriteData_r1_1d'
    !end if

    ! set shape of data, extend with dimensions 1
    shap = 1
    shap(1:rank) = shape(data)
    ns = rank
    if ( present(start ) ) ns = size(start)

    ! write data:
    strt = 0; if ( present(start ) ) strt(1:ns) = start
    strd = 1; if ( present(stride) ) strd(1:ns) = stride
    select case ( sds%typ )
      case ( 'int' )
        select case ( sds%knd )
          case ( 1 )
            call pa_Init( int8 )
            call pa_SetShape( int8 , shape(data) )
            int8 = int(data,kind=1)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), int8 )
            call pa_Done( int8  )
          case ( 2 )
            call pa_Init( int16 )
            call pa_SetShape( int16 , shape(data) )
            int16 = int(data,kind=2)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), int16 )
            call pa_Done( int16  )
          case ( 4 )
            call pa_Init( int32 )
            call pa_SetShape( int32 , shape(data) )
            int32 = int(data,kind=4)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), int32 )
            call pa_Done( int32  )
          case ( 8 )
            call pa_Init( int64 )
            call pa_SetShape( int64 , shape(data) )
            int64 = int(data,kind=8)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), int64 )
            call pa_Done( int64  )
          case default
            write (*,'("ERROR - unsupported integer kind : ",i4)') sds%knd
            write (*,'("ERROR in ",a)') rname; status=1; return
        end select
      case ( 'flt' )
        select case ( sds%knd )
          case ( 4 )
            call pa_Init( float32 )
            call pa_SetShape( float32 , shape(data) )
            float32 = real(data,kind=4)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), float32 )
            call pa_Done( float32  )
          case ( 8 )
            call pa_Init( float64 )
            call pa_SetShape( float64 , shape(data) )
            float64 = real(data,kind=8)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), float64 )
            call pa_Done( float64  )
          case default
            write (*,'("ERROR - unsupported real kind : ",i4)') sds%knd
            write (*,'("ERROR in ",a)') rname; status=1; return
        end select
      case default
        write (*,'("ERROR - unknown sds%typ : ",a)') trim(sds%typ)
        write (*,'("ERROR in ",a)') rname; status=1; return
    end select
    if ( status == FAIL ) then
      write (*,'("ERROR - writing data set:")')
      write (*,'("ERROR -   data set : ",a)') trim(sds%name)
      write (*,'("ERROR -   hdf file : ",a)') trim(sds%hdfname)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine sds_WriteData_i1_1d
  
  
  ! ***

  
  subroutine sds_WriteData_i1_2d( sds, data, status, start, stride )
  
    use parray, only : pa_Init, pa_SetShape, pa_Done
    use file_hdf_base, only : wpi
    use file_hdf_base, only : TSds
    use file_hdf_base, only : wp_int8, wp_int16, wp_int32, wp_int64
    use file_hdf_base, only : wp_float32, wp_float64

    ! --- const --------------------------
    
    character(len=*), parameter ::  rname = mname//'/sds_WriteData_i1_2d'

    integer, parameter   ::  rank = 2
  
    ! --- in/out -------------------------
    
    type(TSds), intent(in)                    ::  sds
    integer(1), intent(in)                 ::  data(:,:)
    integer, intent(out)                      ::  status
    
    integer, intent(in), optional             ::  start(rank)
    integer, intent(in), optional             ::  stride(rank)
    
    ! --- local -------------------------------
    
    integer               ::  strt(7), strd(7)
    integer               ::  shap(7), ns
    
    integer(wp_int8 ), pointer          ::  int8 (:,:)
    integer(wp_int16), pointer          ::  int16(:,:)
    integer(wp_int32), pointer          ::  int32(:,:)
    integer(wp_int64), pointer          ::  int64(:,:)
    real(wp_float32), pointer           ::  float32(:,:)
    real(wp_float64), pointer           ::  float64(:,:)
    
    ! --- external ----------------------------
    
    integer(wpi), external  ::  sfWData
    
    ! --- begin -------------------------------
    
    !! check shape
    !if ( any( shape(data) /= sds%shp(1:sds%rnk) ) ) then
    !  print *, 'Shape of data does not match shape specified during creation:'
    !  print *, '  data        : ', shape(data)
    !  print *, '  created for : ', sds%shp(1:sds%rnk)
    !  stop 'FATAL ERROR IN sds_WriteData_r1_1d'
    !end if

    ! set shape of data, extend with dimensions 1
    shap = 1
    shap(1:rank) = shape(data)
    ns = rank
    if ( present(start ) ) ns = size(start)

    ! write data:
    strt = 0; if ( present(start ) ) strt(1:ns) = start
    strd = 1; if ( present(stride) ) strd(1:ns) = stride
    select case ( sds%typ )
      case ( 'int' )
        select case ( sds%knd )
          case ( 1 )
            call pa_Init( int8 )
            call pa_SetShape( int8 , shape(data) )
            int8 = int(data,kind=1)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), int8 )
            call pa_Done( int8  )
          case ( 2 )
            call pa_Init( int16 )
            call pa_SetShape( int16 , shape(data) )
            int16 = int(data,kind=2)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), int16 )
            call pa_Done( int16  )
          case ( 4 )
            call pa_Init( int32 )
            call pa_SetShape( int32 , shape(data) )
            int32 = int(data,kind=4)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), int32 )
            call pa_Done( int32  )
          case ( 8 )
            call pa_Init( int64 )
            call pa_SetShape( int64 , shape(data) )
            int64 = int(data,kind=8)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), int64 )
            call pa_Done( int64  )
          case default
            write (*,'("ERROR - unsupported integer kind : ",i4)') sds%knd
            write (*,'("ERROR in ",a)') rname; status=1; return
        end select
      case ( 'flt' )
        select case ( sds%knd )
          case ( 4 )
            call pa_Init( float32 )
            call pa_SetShape( float32 , shape(data) )
            float32 = real(data,kind=4)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), float32 )
            call pa_Done( float32  )
          case ( 8 )
            call pa_Init( float64 )
            call pa_SetShape( float64 , shape(data) )
            float64 = real(data,kind=8)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), float64 )
            call pa_Done( float64  )
          case default
            write (*,'("ERROR - unsupported real kind : ",i4)') sds%knd
            write (*,'("ERROR in ",a)') rname; status=1; return
        end select
      case default
        write (*,'("ERROR - unknown sds%typ : ",a)') trim(sds%typ)
        write (*,'("ERROR in ",a)') rname; status=1; return
    end select
    if ( status == FAIL ) then
      write (*,'("ERROR - writing data set:")')
      write (*,'("ERROR -   data set : ",a)') trim(sds%name)
      write (*,'("ERROR -   hdf file : ",a)') trim(sds%hdfname)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine sds_WriteData_i1_2d
  
  
  ! ***

  
  subroutine sds_WriteData_i1_3d( sds, data, status, start, stride )
  
    use parray, only : pa_Init, pa_SetShape, pa_Done
    use file_hdf_base, only : wpi
    use file_hdf_base, only : TSds
    use file_hdf_base, only : wp_int8, wp_int16, wp_int32, wp_int64
    use file_hdf_base, only : wp_float32, wp_float64

    ! --- const --------------------------
    
    character(len=*), parameter ::  rname = mname//'/sds_WriteData_i1_3d'

    integer, parameter   ::  rank = 3
  
    ! --- in/out -------------------------
    
    type(TSds), intent(in)                    ::  sds
    integer(1), intent(in)                 ::  data(:,:,:)
    integer, intent(out)                      ::  status
    
    integer, intent(in), optional             ::  start(rank)
    integer, intent(in), optional             ::  stride(rank)
    
    ! --- local -------------------------------
    
    integer               ::  strt(7), strd(7)
    integer               ::  shap(7), ns
    
    integer(wp_int8 ), pointer          ::  int8 (:,:,:)
    integer(wp_int16), pointer          ::  int16(:,:,:)
    integer(wp_int32), pointer          ::  int32(:,:,:)
    integer(wp_int64), pointer          ::  int64(:,:,:)
    real(wp_float32), pointer           ::  float32(:,:,:)
    real(wp_float64), pointer           ::  float64(:,:,:)
    
    ! --- external ----------------------------
    
    integer(wpi), external  ::  sfWData
    
    ! --- begin -------------------------------
    
    !! check shape
    !if ( any( shape(data) /= sds%shp(1:sds%rnk) ) ) then
    !  print *, 'Shape of data does not match shape specified during creation:'
    !  print *, '  data        : ', shape(data)
    !  print *, '  created for : ', sds%shp(1:sds%rnk)
    !  stop 'FATAL ERROR IN sds_WriteData_r1_1d'
    !end if

    ! set shape of data, extend with dimensions 1
    shap = 1
    shap(1:rank) = shape(data)
    ns = rank
    if ( present(start ) ) ns = size(start)

    ! write data:
    strt = 0; if ( present(start ) ) strt(1:ns) = start
    strd = 1; if ( present(stride) ) strd(1:ns) = stride
    select case ( sds%typ )
      case ( 'int' )
        select case ( sds%knd )
          case ( 1 )
            call pa_Init( int8 )
            call pa_SetShape( int8 , shape(data) )
            int8 = int(data,kind=1)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), int8 )
            call pa_Done( int8  )
          case ( 2 )
            call pa_Init( int16 )
            call pa_SetShape( int16 , shape(data) )
            int16 = int(data,kind=2)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), int16 )
            call pa_Done( int16  )
          case ( 4 )
            call pa_Init( int32 )
            call pa_SetShape( int32 , shape(data) )
            int32 = int(data,kind=4)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), int32 )
            call pa_Done( int32  )
          case ( 8 )
            call pa_Init( int64 )
            call pa_SetShape( int64 , shape(data) )
            int64 = int(data,kind=8)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), int64 )
            call pa_Done( int64  )
          case default
            write (*,'("ERROR - unsupported integer kind : ",i4)') sds%knd
            write (*,'("ERROR in ",a)') rname; status=1; return
        end select
      case ( 'flt' )
        select case ( sds%knd )
          case ( 4 )
            call pa_Init( float32 )
            call pa_SetShape( float32 , shape(data) )
            float32 = real(data,kind=4)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), float32 )
            call pa_Done( float32  )
          case ( 8 )
            call pa_Init( float64 )
            call pa_SetShape( float64 , shape(data) )
            float64 = real(data,kind=8)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), float64 )
            call pa_Done( float64  )
          case default
            write (*,'("ERROR - unsupported real kind : ",i4)') sds%knd
            write (*,'("ERROR in ",a)') rname; status=1; return
        end select
      case default
        write (*,'("ERROR - unknown sds%typ : ",a)') trim(sds%typ)
        write (*,'("ERROR in ",a)') rname; status=1; return
    end select
    if ( status == FAIL ) then
      write (*,'("ERROR - writing data set:")')
      write (*,'("ERROR -   data set : ",a)') trim(sds%name)
      write (*,'("ERROR -   hdf file : ",a)') trim(sds%hdfname)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine sds_WriteData_i1_3d
  
  
  ! ***

  
  subroutine sds_WriteData_i1_4d( sds, data, status, start, stride )
  
    use parray, only : pa_Init, pa_SetShape, pa_Done
    use file_hdf_base, only : wpi
    use file_hdf_base, only : TSds
    use file_hdf_base, only : wp_int8, wp_int16, wp_int32, wp_int64
    use file_hdf_base, only : wp_float32, wp_float64

    ! --- const --------------------------
    
    character(len=*), parameter ::  rname = mname//'/sds_WriteData_i1_4d'

    integer, parameter   ::  rank = 4
  
    ! --- in/out -------------------------
    
    type(TSds), intent(in)                    ::  sds
    integer(1), intent(in)                 ::  data(:,:,:,:)
    integer, intent(out)                      ::  status
    
    integer, intent(in), optional             ::  start(rank)
    integer, intent(in), optional             ::  stride(rank)
    
    ! --- local -------------------------------
    
    integer               ::  strt(7), strd(7)
    integer               ::  shap(7), ns
    
    integer(wp_int8 ), pointer          ::  int8 (:,:,:,:)
    integer(wp_int16), pointer          ::  int16(:,:,:,:)
    integer(wp_int32), pointer          ::  int32(:,:,:,:)
    integer(wp_int64), pointer          ::  int64(:,:,:,:)
    real(wp_float32), pointer           ::  float32(:,:,:,:)
    real(wp_float64), pointer           ::  float64(:,:,:,:)
    
    ! --- external ----------------------------
    
    integer(wpi), external  ::  sfWData
    
    ! --- begin -------------------------------
    
    !! check shape
    !if ( any( shape(data) /= sds%shp(1:sds%rnk) ) ) then
    !  print *, 'Shape of data does not match shape specified during creation:'
    !  print *, '  data        : ', shape(data)
    !  print *, '  created for : ', sds%shp(1:sds%rnk)
    !  stop 'FATAL ERROR IN sds_WriteData_r1_1d'
    !end if

    ! set shape of data, extend with dimensions 1
    shap = 1
    shap(1:rank) = shape(data)
    ns = rank
    if ( present(start ) ) ns = size(start)

    ! write data:
    strt = 0; if ( present(start ) ) strt(1:ns) = start
    strd = 1; if ( present(stride) ) strd(1:ns) = stride
    select case ( sds%typ )
      case ( 'int' )
        select case ( sds%knd )
          case ( 1 )
            call pa_Init( int8 )
            call pa_SetShape( int8 , shape(data) )
            int8 = int(data,kind=1)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), int8 )
            call pa_Done( int8  )
          case ( 2 )
            call pa_Init( int16 )
            call pa_SetShape( int16 , shape(data) )
            int16 = int(data,kind=2)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), int16 )
            call pa_Done( int16  )
          case ( 4 )
            call pa_Init( int32 )
            call pa_SetShape( int32 , shape(data) )
            int32 = int(data,kind=4)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), int32 )
            call pa_Done( int32  )
          case ( 8 )
            call pa_Init( int64 )
            call pa_SetShape( int64 , shape(data) )
            int64 = int(data,kind=8)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), int64 )
            call pa_Done( int64  )
          case default
            write (*,'("ERROR - unsupported integer kind : ",i4)') sds%knd
            write (*,'("ERROR in ",a)') rname; status=1; return
        end select
      case ( 'flt' )
        select case ( sds%knd )
          case ( 4 )
            call pa_Init( float32 )
            call pa_SetShape( float32 , shape(data) )
            float32 = real(data,kind=4)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), float32 )
            call pa_Done( float32  )
          case ( 8 )
            call pa_Init( float64 )
            call pa_SetShape( float64 , shape(data) )
            float64 = real(data,kind=8)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), float64 )
            call pa_Done( float64  )
          case default
            write (*,'("ERROR - unsupported real kind : ",i4)') sds%knd
            write (*,'("ERROR in ",a)') rname; status=1; return
        end select
      case default
        write (*,'("ERROR - unknown sds%typ : ",a)') trim(sds%typ)
        write (*,'("ERROR in ",a)') rname; status=1; return
    end select
    if ( status == FAIL ) then
      write (*,'("ERROR - writing data set:")')
      write (*,'("ERROR -   data set : ",a)') trim(sds%name)
      write (*,'("ERROR -   hdf file : ",a)') trim(sds%hdfname)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine sds_WriteData_i1_4d
  
  
  ! ***

  
  subroutine sds_WriteData_i1_5d( sds, data, status, start, stride )
  
    use parray, only : pa_Init, pa_SetShape, pa_Done
    use file_hdf_base, only : wpi
    use file_hdf_base, only : TSds
    use file_hdf_base, only : wp_int8, wp_int16, wp_int32, wp_int64
    use file_hdf_base, only : wp_float32, wp_float64

    ! --- const --------------------------
    
    character(len=*), parameter ::  rname = mname//'/sds_WriteData_i1_5d'

    integer, parameter   ::  rank = 5
  
    ! --- in/out -------------------------
    
    type(TSds), intent(in)                    ::  sds
    integer(1), intent(in)                 ::  data(:,:,:,:,:)
    integer, intent(out)                      ::  status
    
    integer, intent(in), optional             ::  start(rank)
    integer, intent(in), optional             ::  stride(rank)
    
    ! --- local -------------------------------
    
    integer               ::  strt(7), strd(7)
    integer               ::  shap(7), ns
    
    integer(wp_int8 ), pointer          ::  int8 (:,:,:,:,:)
    integer(wp_int16), pointer          ::  int16(:,:,:,:,:)
    integer(wp_int32), pointer          ::  int32(:,:,:,:,:)
    integer(wp_int64), pointer          ::  int64(:,:,:,:,:)
    real(wp_float32), pointer           ::  float32(:,:,:,:,:)
    real(wp_float64), pointer           ::  float64(:,:,:,:,:)
    
    ! --- external ----------------------------
    
    integer(wpi), external  ::  sfWData
    
    ! --- begin -------------------------------
    
    !! check shape
    !if ( any( shape(data) /= sds%shp(1:sds%rnk) ) ) then
    !  print *, 'Shape of data does not match shape specified during creation:'
    !  print *, '  data        : ', shape(data)
    !  print *, '  created for : ', sds%shp(1:sds%rnk)
    !  stop 'FATAL ERROR IN sds_WriteData_r1_1d'
    !end if

    ! set shape of data, extend with dimensions 1
    shap = 1
    shap(1:rank) = shape(data)
    ns = rank
    if ( present(start ) ) ns = size(start)

    ! write data:
    strt = 0; if ( present(start ) ) strt(1:ns) = start
    strd = 1; if ( present(stride) ) strd(1:ns) = stride
    select case ( sds%typ )
      case ( 'int' )
        select case ( sds%knd )
          case ( 1 )
            call pa_Init( int8 )
            call pa_SetShape( int8 , shape(data) )
            int8 = int(data,kind=1)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), int8 )
            call pa_Done( int8  )
          case ( 2 )
            call pa_Init( int16 )
            call pa_SetShape( int16 , shape(data) )
            int16 = int(data,kind=2)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), int16 )
            call pa_Done( int16  )
          case ( 4 )
            call pa_Init( int32 )
            call pa_SetShape( int32 , shape(data) )
            int32 = int(data,kind=4)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), int32 )
            call pa_Done( int32  )
          case ( 8 )
            call pa_Init( int64 )
            call pa_SetShape( int64 , shape(data) )
            int64 = int(data,kind=8)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), int64 )
            call pa_Done( int64  )
          case default
            write (*,'("ERROR - unsupported integer kind : ",i4)') sds%knd
            write (*,'("ERROR in ",a)') rname; status=1; return
        end select
      case ( 'flt' )
        select case ( sds%knd )
          case ( 4 )
            call pa_Init( float32 )
            call pa_SetShape( float32 , shape(data) )
            float32 = real(data,kind=4)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), float32 )
            call pa_Done( float32  )
          case ( 8 )
            call pa_Init( float64 )
            call pa_SetShape( float64 , shape(data) )
            float64 = real(data,kind=8)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), float64 )
            call pa_Done( float64  )
          case default
            write (*,'("ERROR - unsupported real kind : ",i4)') sds%knd
            write (*,'("ERROR in ",a)') rname; status=1; return
        end select
      case default
        write (*,'("ERROR - unknown sds%typ : ",a)') trim(sds%typ)
        write (*,'("ERROR in ",a)') rname; status=1; return
    end select
    if ( status == FAIL ) then
      write (*,'("ERROR - writing data set:")')
      write (*,'("ERROR -   data set : ",a)') trim(sds%name)
      write (*,'("ERROR -   hdf file : ",a)') trim(sds%hdfname)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine sds_WriteData_i1_5d
  
  
  ! ***

  
  subroutine sds_WriteData_i1_6d( sds, data, status, start, stride )
  
    use parray, only : pa_Init, pa_SetShape, pa_Done
    use file_hdf_base, only : wpi
    use file_hdf_base, only : TSds
    use file_hdf_base, only : wp_int8, wp_int16, wp_int32, wp_int64
    use file_hdf_base, only : wp_float32, wp_float64

    ! --- const --------------------------
    
    character(len=*), parameter ::  rname = mname//'/sds_WriteData_i1_6d'

    integer, parameter   ::  rank = 6
  
    ! --- in/out -------------------------
    
    type(TSds), intent(in)                    ::  sds
    integer(1), intent(in)                 ::  data(:,:,:,:,:,:)
    integer, intent(out)                      ::  status
    
    integer, intent(in), optional             ::  start(rank)
    integer, intent(in), optional             ::  stride(rank)
    
    ! --- local -------------------------------
    
    integer               ::  strt(7), strd(7)
    integer               ::  shap(7), ns
    
    integer(wp_int8 ), pointer          ::  int8 (:,:,:,:,:,:)
    integer(wp_int16), pointer          ::  int16(:,:,:,:,:,:)
    integer(wp_int32), pointer          ::  int32(:,:,:,:,:,:)
    integer(wp_int64), pointer          ::  int64(:,:,:,:,:,:)
    real(wp_float32), pointer           ::  float32(:,:,:,:,:,:)
    real(wp_float64), pointer           ::  float64(:,:,:,:,:,:)
    
    ! --- external ----------------------------
    
    integer(wpi), external  ::  sfWData
    
    ! --- begin -------------------------------
    
    !! check shape
    !if ( any( shape(data) /= sds%shp(1:sds%rnk) ) ) then
    !  print *, 'Shape of data does not match shape specified during creation:'
    !  print *, '  data        : ', shape(data)
    !  print *, '  created for : ', sds%shp(1:sds%rnk)
    !  stop 'FATAL ERROR IN sds_WriteData_r1_1d'
    !end if

    ! set shape of data, extend with dimensions 1
    shap = 1
    shap(1:rank) = shape(data)
    ns = rank
    if ( present(start ) ) ns = size(start)

    ! write data:
    strt = 0; if ( present(start ) ) strt(1:ns) = start
    strd = 1; if ( present(stride) ) strd(1:ns) = stride
    select case ( sds%typ )
      case ( 'int' )
        select case ( sds%knd )
          case ( 1 )
            call pa_Init( int8 )
            call pa_SetShape( int8 , shape(data) )
            int8 = int(data,kind=1)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), int8 )
            call pa_Done( int8  )
          case ( 2 )
            call pa_Init( int16 )
            call pa_SetShape( int16 , shape(data) )
            int16 = int(data,kind=2)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), int16 )
            call pa_Done( int16  )
          case ( 4 )
            call pa_Init( int32 )
            call pa_SetShape( int32 , shape(data) )
            int32 = int(data,kind=4)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), int32 )
            call pa_Done( int32  )
          case ( 8 )
            call pa_Init( int64 )
            call pa_SetShape( int64 , shape(data) )
            int64 = int(data,kind=8)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), int64 )
            call pa_Done( int64  )
          case default
            write (*,'("ERROR - unsupported integer kind : ",i4)') sds%knd
            write (*,'("ERROR in ",a)') rname; status=1; return
        end select
      case ( 'flt' )
        select case ( sds%knd )
          case ( 4 )
            call pa_Init( float32 )
            call pa_SetShape( float32 , shape(data) )
            float32 = real(data,kind=4)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), float32 )
            call pa_Done( float32  )
          case ( 8 )
            call pa_Init( float64 )
            call pa_SetShape( float64 , shape(data) )
            float64 = real(data,kind=8)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), float64 )
            call pa_Done( float64  )
          case default
            write (*,'("ERROR - unsupported real kind : ",i4)') sds%knd
            write (*,'("ERROR in ",a)') rname; status=1; return
        end select
      case default
        write (*,'("ERROR - unknown sds%typ : ",a)') trim(sds%typ)
        write (*,'("ERROR in ",a)') rname; status=1; return
    end select
    if ( status == FAIL ) then
      write (*,'("ERROR - writing data set:")')
      write (*,'("ERROR -   data set : ",a)') trim(sds%name)
      write (*,'("ERROR -   hdf file : ",a)') trim(sds%hdfname)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine sds_WriteData_i1_6d
  
  
  ! ***

  
  subroutine sds_WriteData_i1_7d( sds, data, status, start, stride )
  
    use parray, only : pa_Init, pa_SetShape, pa_Done
    use file_hdf_base, only : wpi
    use file_hdf_base, only : TSds
    use file_hdf_base, only : wp_int8, wp_int16, wp_int32, wp_int64
    use file_hdf_base, only : wp_float32, wp_float64

    ! --- const --------------------------
    
    character(len=*), parameter ::  rname = mname//'/sds_WriteData_i1_7d'

    integer, parameter   ::  rank = 7
  
    ! --- in/out -------------------------
    
    type(TSds), intent(in)                    ::  sds
    integer(1), intent(in)                 ::  data(:,:,:,:,:,:,:)
    integer, intent(out)                      ::  status
    
    integer, intent(in), optional             ::  start(rank)
    integer, intent(in), optional             ::  stride(rank)
    
    ! --- local -------------------------------
    
    integer               ::  strt(7), strd(7)
    integer               ::  shap(7), ns
    
    integer(wp_int8 ), pointer          ::  int8 (:,:,:,:,:,:,:)
    integer(wp_int16), pointer          ::  int16(:,:,:,:,:,:,:)
    integer(wp_int32), pointer          ::  int32(:,:,:,:,:,:,:)
    integer(wp_int64), pointer          ::  int64(:,:,:,:,:,:,:)
    real(wp_float32), pointer           ::  float32(:,:,:,:,:,:,:)
    real(wp_float64), pointer           ::  float64(:,:,:,:,:,:,:)
    
    ! --- external ----------------------------
    
    integer(wpi), external  ::  sfWData
    
    ! --- begin -------------------------------
    
    !! check shape
    !if ( any( shape(data) /= sds%shp(1:sds%rnk) ) ) then
    !  print *, 'Shape of data does not match shape specified during creation:'
    !  print *, '  data        : ', shape(data)
    !  print *, '  created for : ', sds%shp(1:sds%rnk)
    !  stop 'FATAL ERROR IN sds_WriteData_r1_1d'
    !end if

    ! set shape of data, extend with dimensions 1
    shap = 1
    shap(1:rank) = shape(data)
    ns = rank
    if ( present(start ) ) ns = size(start)

    ! write data:
    strt = 0; if ( present(start ) ) strt(1:ns) = start
    strd = 1; if ( present(stride) ) strd(1:ns) = stride
    select case ( sds%typ )
      case ( 'int' )
        select case ( sds%knd )
          case ( 1 )
            call pa_Init( int8 )
            call pa_SetShape( int8 , shape(data) )
            int8 = int(data,kind=1)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), int8 )
            call pa_Done( int8  )
          case ( 2 )
            call pa_Init( int16 )
            call pa_SetShape( int16 , shape(data) )
            int16 = int(data,kind=2)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), int16 )
            call pa_Done( int16  )
          case ( 4 )
            call pa_Init( int32 )
            call pa_SetShape( int32 , shape(data) )
            int32 = int(data,kind=4)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), int32 )
            call pa_Done( int32  )
          case ( 8 )
            call pa_Init( int64 )
            call pa_SetShape( int64 , shape(data) )
            int64 = int(data,kind=8)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), int64 )
            call pa_Done( int64  )
          case default
            write (*,'("ERROR - unsupported integer kind : ",i4)') sds%knd
            write (*,'("ERROR in ",a)') rname; status=1; return
        end select
      case ( 'flt' )
        select case ( sds%knd )
          case ( 4 )
            call pa_Init( float32 )
            call pa_SetShape( float32 , shape(data) )
            float32 = real(data,kind=4)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), float32 )
            call pa_Done( float32  )
          case ( 8 )
            call pa_Init( float64 )
            call pa_SetShape( float64 , shape(data) )
            float64 = real(data,kind=8)
            status = sfWData( sds%id, strt(1:ns), strd(1:ns), shap(1:ns), float64 )
            call pa_Done( float64  )
          case default
            write (*,'("ERROR - unsupported real kind : ",i4)') sds%knd
            write (*,'("ERROR in ",a)') rname; status=1; return
        end select
      case default
        write (*,'("ERROR - unknown sds%typ : ",a)') trim(sds%typ)
        write (*,'("ERROR in ",a)') rname; status=1; return
    end select
    if ( status == FAIL ) then
      write (*,'("ERROR - writing data set:")')
      write (*,'("ERROR -   data set : ",a)') trim(sds%name)
      write (*,'("ERROR -   hdf file : ",a)') trim(sds%hdfname)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine sds_WriteData_i1_7d
  
  
  
  ! ############################################################
  ! ###
  ! ### dimensions
  ! ###
  ! ############################################################
  

  ! ================================================================
  ! set dimension scale
  ! ================================================================
  
  
  subroutine dim_SetScale_i1( sdim, scale, status, knd )
  
    use file_hdf_base, only : wpi
    use file_hdf_base, only : TSdsDim
    use file_hdf_base, only : wp_int8, wp_int16, wp_int32, wp_int64

    ! --- in/out -------------------------
    
    type(TSdsDim), intent(in)           ::  sdim
    integer(1), intent(in)           ::  scale(:)
    integer, intent(out)                ::  status

    integer, intent(in), optional       ::  knd
    
    ! --- const --------------------------
    
    character(len=*), parameter ::  rname = mname//'/dim_SetScale_i1'

    ! --- local ---------------------------
    
    integer              ::  ikind

    ! --- external ---------------------------
    
    integer(wpi), external   :: sfSDScale
  
    ! --- begin ---------------------------
      
    ikind = kind(scale)
    if ( present(knd) ) ikind = knd

    select case ( ikind )
      case ( 1 )
        status = sfSDScale( sdim%id, size(scale), DFNT_INT8 , int(scale,kind=wp_int8 ) )
      case ( 2 )
        status = sfSDScale( sdim%id, size(scale), DFNT_INT16, int(scale,kind=wp_int16) )
      case ( 4 )
        status = sfSDScale( sdim%id, size(scale), DFNT_INT32, int(scale,kind=wp_int32) )
      case ( 8 )
        status = sfSDScale( sdim%id, size(scale), DFNT_INT64, int(scale,kind=wp_int64) )
      case default
        write (*,'("ERROR - unsupported integer kind : ",i4)') ikind
        write (*,'("ERROR in ",a)') rname; status=1; return
    end select
    if ( status /= SUCCEED ) then
      write (*,'("ERROR - writing scale")')
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! ok
    status = 0
    
  end subroutine dim_SetScale_i1
  


  ! ================================================================
  ! set dimension stuff
  ! ================================================================
  
  
  subroutine sds_SetDim_i1( sds, dim_index, name, unit, scale, status, knd )
  
    use file_hdf_base, only : TSds, TSdsDim, Init, Done
    use file_hdf_base, only : MAX_DATA_RANK
    use file_hdf_base, only : GetInfo
    use file_hdf_base, only : Select, SetName
    use file_hdf_s, only : WriteAttribute
    
    ! --- in/out -------------------------
    
    type(TSds), intent(in)              ::  sds
    integer, intent(in)                 ::  dim_index
    character(len=*), intent(in)        ::  name
    character(len=*), intent(in)        ::  unit
    integer(1), intent(in)           ::  scale(:)
    integer, intent(out)                ::  status

    integer, intent(in), optional       ::  knd
    
    ! --- const --------------------------
    
    character(len=*), parameter ::  rname = mname//'/sds_SetDim_i1'

    ! --- local ---------------------------
    
    integer              ::  data_rank, data_dims(0:MAX_DATA_RANK-1)
    type(TSdsDim)        ::  sdim

    ! --- begin ---------------------------
    
    call GetInfo( sds, status, data_rank=data_rank, data_dims=data_dims )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    if ( dim_index < 0 .or. dim_index > data_rank-1 ) then
      write (*,'("ERROR - wrong dimension index : ",i4)') dim_index
      write (*,'("ERROR - expecting range 0, .., ",i4)') data_rank-1
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if
    
    call Init( sdim, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    call Select( sdim, sds, dim_index, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    call SetName( sdim, name, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    call WriteAttribute( sdim, 'unit', unit, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    if ( size(scale) /= data_dims(dim_index) ) then
      write (*,'("ERROR - wrong scale length : ",i4)') size(scale)
      write (*,'("ERROR - expecting length ",i4," for dim index ",i4)') data_dims(dim_index), dim_index
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    call SetScale( sdim, scale, status, knd )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    call Done( sdim, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
      
    ! ok
    status = 0
    
  end subroutine sds_SetDim_i1
  


  ! ================================================================
  ! read attributes
  ! ================================================================
  
  
  subroutine dim_ReadAttribute_i1_0d( sdim, name, i, status )
  
    use file_hdf_base, only : TSdsDim

    ! --- in/out -------------------------
    
    type(TSdsDim), intent(in)           ::  sdim
    character(len=*), intent(in)        ::  name
    integer(1), intent(out)          ::  i
    integer, intent(out)                ::  status

    ! --- const --------------------------
    
    character(len=*), parameter ::  rname = mname//'/dim_ReadAttribute_i1_0d'

    ! --- begin -------------------------------
    
    call ReadAttribute( sdim%id, name, i, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
        
    ! ok
    status = 0
    
  end subroutine dim_ReadAttribute_i1_0d
  
  
  ! ***

  
  subroutine dim_ReadAttribute_i1_1d( sdim, name, i, status )
  
    use file_hdf_base, only : TSdsDim

    ! --- in/out -------------------------
    
    type(TSdsDim), intent(in)           ::  sdim
    character(len=*), intent(in)        ::  name
    integer(1), intent(out)          ::  i(:)
    integer, intent(out)                ::  status
    
    ! --- const --------------------------
    
    character(len=*), parameter ::  rname = mname//'/dim_ReadAttribute_i1_1d'

    ! --- begin -------------------------------
    
    call ReadAttribute( sdim%id, name, i, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    
    ! ok
    status = 0
    
  end subroutine dim_ReadAttribute_i1_1d
  
  
  ! =============================================================
  ! === check attributes
  ! =============================================================

  
  subroutine dim_CheckAttribute_i1_0d( sdim, name, i, status )
  
    use file_hdf_base, only : TSdsDim
  
    ! --- in/out -------------------------
    
    type(TSdsDim), intent(in)           ::  sdim
    character(len=*), intent(in)        ::  name
    integer(1), intent(in)           ::  i
    integer, intent(inout)              ::  status
    
    ! --- const --------------------------
    
    character(len=*), parameter ::  rname = mname//'/dim_CheckAttribute_i1_0d'

    ! --- local ------------------------------
    
    logical          ::  verbose
    
    ! --- begin ---------------------------
    
    ! write error messages ?
    verbose = status == 0
    
    call CheckAttribute( sdim%id, name, i, status )
    if (status>0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    if (status<0) then
      if (verbose) write (*,'("ERROR in ",a)') rname
      status=-1; return
    end if

    ! ok
    status = 0
    
  end subroutine dim_CheckAttribute_i1_0d
  
  
  ! ***
  
  
  subroutine dim_CheckAttribute_i1_1d( sdim, name, i, status )
  
    use file_hdf_base, only : TSdsDim
  
    ! --- in/out -------------------------
    
    type(TSdsDim), intent(in)          ::  sdim
    character(len=*), intent(in)       ::  name
    integer(1), intent(in)          ::  i(:)
    integer, intent(inout)             ::  status
    
    ! --- const --------------------------
    
    character(len=*), parameter ::  rname = mname//'/dim_CheckAttribute_i1_1d'

    ! --- local ------------------------------
    
    logical          ::  verbose
    
    ! --- begin ---------------------------
    
    ! write error messages ?
    verbose = status == 0
    
    call CheckAttribute( sdim%id, name, i, status )
    if (status>0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    if (status<0) then
      if (verbose) write (*,'("ERROR in ",a)') rname
      status=-1; return
    end if

    ! ok
    status = 0
    
  end subroutine dim_CheckAttribute_i1_1d
  
  

  ! ================================================================
  ! write attributes
  ! ================================================================
  
  
  subroutine dim_WriteAttribute_i1_0d( sdim, name, i, status, knd )
  
    use file_hdf_base, only : TSdsDim
  
    ! --- in/out -------------------------
    
    type(TSdsDim), intent(in)          ::  sdim
    character(len=*), intent(in)       ::  name
    integer(1), intent(in)          ::  i
    integer, intent(out)               ::  status

    integer, intent(in), optional      ::  knd
    
    ! --- const --------------------------
    
    character(len=*), parameter ::  rname = mname//'/dim_WriteAttribute_i1_0d'

    ! --- begin -------------------------------
    
    call WriteAttribute( sdim%id, name, i, status, knd )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
        
    ! ok
    status = 0
    
  end subroutine dim_WriteAttribute_i1_0d
  
  
  ! ***

  
  subroutine dim_WriteAttribute_i1_1d( sdim, name, i, status, knd )
  
    use file_hdf_base, only : TSdsDim
  
    ! --- in/out -------------------------
    
    type(TSdsDim), intent(in)          ::  sdim
    character(len=*), intent(in)       ::  name
    integer(1), intent(in)          ::  i(:)
    integer, intent(out)               ::  status

    integer, intent(in), optional      ::  knd
    
    ! --- const --------------------------
    
    character(len=*), parameter ::  rname = mname//'/dim_WriteAttribute_i1_1d'

    ! --- begin -------------------------------
    
    call WriteAttribute( sdim%id, name, i, status, knd )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    
    ! ok
    status = 0
    
  end subroutine dim_WriteAttribute_i1_1d
  
  
  ! ############################################################
  ! ###
  ! ### hdf files
  ! ###
  ! ############################################################
  

  ! ================================================================
  ! read attributes
  ! ================================================================
  
  
  subroutine hdf_ReadAttribute_i1_0d( hdf, name, i, status )
  
    use file_hdf_base, only : THdfFile
  
    ! --- in/out -------------------------
    
    type(THdfFile), intent(in)         ::  hdf
    character(len=*), intent(in)       ::  name
    integer(1), intent(out)         ::  i
    integer, intent(out)               ::  status
    
    ! --- const --------------------------
    
    character(len=*), parameter ::  rname = mname//'/hdf_ReadAttribute_i1_0d'

    ! --- begin -------------------------------
    
    call ReadAttribute( hdf%id, name, i, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
        
    ! ok
    status = 0
    
  end subroutine hdf_ReadAttribute_i1_0d
  
  
  ! ***

  
  subroutine hdf_ReadAttribute_i1_1d( hdf, name, i, status )
  
    use file_hdf_base, only : THdfFile
  
    ! --- in/out -------------------------
    
    type(THdfFile), intent(in)         ::  hdf
    character(len=*), intent(in)       ::  name
    integer(1), intent(out)         ::  i(:)
    integer, intent(out)               ::  status
    
    ! --- const --------------------------
    
    character(len=*), parameter ::  rname = mname//'/hdf_ReadAttribute_i1_1d'

    ! --- begin -------------------------------
    
    call ReadAttribute( hdf%id, name, i, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    
    ! ok
    status = 0
    
  end subroutine hdf_ReadAttribute_i1_1d
  
  
  ! =============================================================
  ! === check attributes
  ! =============================================================

  
  subroutine hdf_CheckAttribute_i1_0d( hdf, name, i, status )
  
    use file_hdf_base, only : THdfFile
  
    ! --- in/out -------------------------
    
    type(THdfFile), intent(in)         ::  hdf
    character(len=*), intent(in)       ::  name
    integer(1), intent(in)          ::  i
    integer, intent(inout)             ::  status
    
    ! --- const --------------------------
    
    character(len=*), parameter ::  rname = mname//'/hdf_CheckAttribute_i1_0d'

    ! --- local ------------------------------
    
    logical          ::  verbose
    
    ! --- begin ---------------------------
    
    ! write error messages ?
    verbose = status == 0
    
    call CheckAttribute( hdf%id, name, i, status )
    if (status>0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    if (status<0) then
      if (verbose) write (*,'("ERROR in ",a)') rname
      status=-1; return
    end if

    ! ok
    status = 0
    
  end subroutine hdf_CheckAttribute_i1_0d
  
  
  ! ***
  
  
  subroutine hdf_CheckAttribute_i1_1d( hdf, name, i, status )
  
    use file_hdf_base, only : THdfFile
  
    ! --- in/out -------------------------
    
    type(THdfFile), intent(in)         ::  hdf
    character(len=*), intent(in)       ::  name
    integer(1), intent(in)          ::  i(:)
    integer, intent(inout)             ::  status
    
    ! --- const --------------------------
    
    character(len=*), parameter ::  rname = mname//'/hdf_CheckAttribute_i1_1d'

    ! --- local ------------------------------
    
    logical          ::  verbose
    
    ! --- begin ---------------------------
    
    ! write error messages ?
    verbose = status == 0
    
    call CheckAttribute( hdf%id, name, i, status )
    if (status>0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    if (status<0) then
      if (verbose) write (*,'("ERROR in ",a)') rname
      status=-1; return
    end if

    ! ok
    status = 0
    
  end subroutine hdf_CheckAttribute_i1_1d
  
  

  ! ================================================================
  ! write attributes
  ! ================================================================
  
  
  subroutine hdf_WriteAttribute_i1_0d( hdf, name, i, status, knd )
  
    use file_hdf_base, only : THdfFile
  
    ! --- in/out -------------------------
    
    type(THdfFile), intent(in)         ::  hdf
    character(len=*), intent(in)       ::  name
    integer(1), intent(in)          ::  i
    integer, intent(out)               ::  status

    integer, intent(in), optional      ::  knd
    
    ! --- const --------------------------
    
    character(len=*), parameter ::  rname = mname//'/hdf_WriteAttribute_i1_0d'

    ! --- begin -------------------------------
    
    call WriteAttribute( hdf%id, name, i, status, knd )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
        
    ! ok
    status = 0
    
  end subroutine hdf_WriteAttribute_i1_0d
  
  
  ! ***

  
  subroutine hdf_WriteAttribute_i1_1d( hdf, name, i, status, knd )
  
    use file_hdf_base, only : THdfFile
  
    ! --- in/out -------------------------
    
    type(THdfFile), intent(in)         ::  hdf
    character(len=*), intent(in)       ::  name
    integer(1), intent(in)          ::  i(:)
    integer, intent(out)               ::  status

    integer, intent(in), optional      ::  knd
    
    ! --- const --------------------------
    
    character(len=*), parameter ::  rname = mname//'/hdf_WriteAttribute_i1_1d'

    ! --- begin -------------------------------
    
    call WriteAttribute( hdf%id, name, i, status, knd )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    
    ! ok
    status = 0
    
  end subroutine hdf_WriteAttribute_i1_1d
  
  
end module file_hdf_i1

