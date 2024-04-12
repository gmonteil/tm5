module file_hdf_s

  implicit none
  
  ! --- in/out --------------------------
  
  private
  
  public   ::  ReadData, WriteData
  public   ::  ReadAttribute, CheckAttribute
  public   ::  WriteAttribute

  ! --- const ----------------------------
  
  include "hdf.f90"
  
  character(len=*), parameter  ::  mname = 'file_hdf_s'
    

  ! --- interfaces ------------------------
  
  interface ReadData
    module procedure sds_ReadData_s_1d
    module procedure sds_ReadData_s_2d
  end interface

  interface WriteData
    module procedure sds_WriteData_s_1d
    module procedure sds_WriteData_s_2d
    module procedure sds_WriteData_s_3d
    module procedure sds_WriteData_s_4d
  end interface
  
  interface ReadAttribute
    module procedure obj_ReadAttribute_s_0d
    module procedure obj_ReadAttribute_s_1d
    !
    module procedure sds_ReadAttribute_s_0d
    module procedure sds_ReadAttribute_s_1d
    !
    module procedure dim_ReadAttribute_s_0d
    module procedure dim_ReadAttribute_s_1d
    !
    module procedure hdf_ReadAttribute_s_0d
    module procedure hdf_ReadAttribute_s_1d
  end interface
  
  interface CheckAttribute
    module procedure obj_CheckAttribute_s_0d
    module procedure obj_CheckAttribute_s_1d
    !
    module procedure sds_CheckAttribute_s_0d
    module procedure sds_CheckAttribute_s_1d
    !
    module procedure dim_CheckAttribute_s_0d
    module procedure dim_CheckAttribute_s_1d
    !
    module procedure hdf_CheckAttribute_s_0d
    module procedure hdf_CheckAttribute_s_1d
  end interface
  
  interface WriteAttribute
    module procedure obj_WriteAttribute_s_0d
    module procedure obj_WriteAttribute_s_1d
    !
    module procedure sds_WriteAttribute_s_0d
    module procedure sds_WriteAttribute_s_1d
    !
    module procedure dim_WriteAttribute_s_0d
    module procedure dim_WriteAttribute_s_1d
    !
    module procedure hdf_WriteAttribute_s_0d
    module procedure hdf_WriteAttribute_s_1d
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
  
  
  subroutine obj_ReadAttribute_s_0d( obj_id, name, s, status )
  
    use file_hdf_base, only : wpi
    use file_hdf_base, only : FindAttribute, CheckAttributeInfo, GetAttributeInfo

    ! --- in/out -------------------------
    
    integer(wpi), intent(in)            ::  obj_id
    character(len=*), intent(in)        ::  name
    character(len=*), intent(inout)     ::  s
    integer, intent(out)                ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/obj_ReadAttribute_s_0d'
    
    ! --- local -------------------------------
    
    integer              ::  attr_index, data_type
    integer              ::  n_values, n
    
    character(len=len(s)) ::  sdum
    
    ! --- external ----------------------------

    integer(wpi), external   ::  sfRCAtt

    ! --- begin -------------------------------
    
    call FindAttribute( obj_id, name, attr_index, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    call GetAttributeInfo( obj_id, attr_index, status, n_values=n_values )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    n = min( n_values, len(s) )
    
    ! extract value:
    call GetAttributeInfo( obj_id, attr_index, status, data_type=data_type )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    
    ! read:
    select case ( data_type )
      case ( DFNT_CHAR )
        status = sfRCAtt( obj_id, attr_index, sdum )
        s = trim(sdum(1:n))
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
    
  end subroutine obj_ReadAttribute_s_0d


  ! ***
      
  
  subroutine obj_ReadAttribute_s_1d( obj_id, name, s, status )
  
    use file_hdf_base, only : wpi
    use file_hdf_base, only : FindAttribute, CheckAttributeInfo, GetAttributeInfo

    ! --- in/out -------------------------
    
    integer(wpi), intent(in)            ::  obj_id
    character(len=*), intent(in)        ::  name
    character(len=*), intent(out)       ::  s(:)
    integer, intent(out)                ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/obj_ReadAttribute_s_1d'
    
    ! --- local -------------------------------
    
    integer              ::  attr_index, data_type
    
    ! --- external ----------------------------

    integer(wpi), external   ::  sfRCAtt

    ! --- begin -------------------------------
    
    call FindAttribute( obj_id, name, attr_index, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    !call CheckAttributeInfo( obj_id, attr_index, n_values=len(s)*size(s) )
    !if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! extract value:
    call GetAttributeInfo( obj_id, attr_index, status, data_type=data_type )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! read:
    select case ( data_type )
      case ( DFNT_CHAR )
        status = sfRCAtt( obj_id, attr_index, s )
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
    
  end subroutine obj_ReadAttribute_s_1d




  ! ================================================================
  ! ===
  ! === check attributes
  ! ===
  ! ================================================================
  
  
  subroutine obj_CheckAttribute_s_0d( obj_id, name, s, status )
  
    use file_hdf_base, only : wpi
  
    ! --- in/out -------------------------
    
    integer(wpi), intent(in)           ::  obj_id
    character(len=*), intent(in)       ::  name
    character(len=*), intent(in)       ::  s
    integer, intent(inout)             ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/obj_CheckAttribute_s_0d'
    
    ! --- local -------------------------------
    
    logical                      ::  verbose
    character(len=len(s))        ::  attr_s
    
    ! --- begin -------------------------------
    
    ! write error messages ?
    verbose = status == 0

    call ReadAttribute( obj_id, name, attr_s, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    
    ! check:
    if ( attr_s /= s ) then
      if (verbose) then
        write (*,'("ERROR - foud different attribute values:")')
        write (*,'("ERROR -   attr name : ",a)') trim(name)
        write (*,'("ERROR -   requested : ",a)') trim(s)
        write (*,'("ERROR -   found     : ",a)') trim(attr_s)
        write (*,'("ERROR in ",a)') rname
      end if
      status=-1; return
    end if

    ! ok
    status = 0
    
  end subroutine obj_CheckAttribute_s_0d


  ! ***
      
  
  subroutine obj_CheckAttribute_s_1d( obj_id, name, s, status )
  
    use file_hdf_base, only : wpi
  
    ! --- in/out -------------------------
    
    integer(wpi), intent(in)            ::  obj_id
    character(len=*), intent(in)        ::  name
    character(len=*), intent(in)        ::  s(:)
    integer, intent(inout)              ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/obj_CheckAttribute_s_1d'
    
    ! --- local -------------------------------
    
    logical                      ::  verbose
    character(len=len(s))        ::  attr_s(size(s))
    
    ! --- begin -------------------------------
    
    ! write error messages ?
    verbose = status == 0

    call ReadAttribute( obj_id, name, attr_s, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    
    ! check
    if ( any( attr_s /= s ) ) then
      if (verbose) then
        write (*,'("ERROR - foud different attribute values:")')
        write (*,'("ERROR -   attr name : ",a)') trim(name)
        write (*,'("ERROR -   requested : ",a," ...")') trim(s(1))
        write (*,'("ERROR -   found     : ",a," ...")') trim(attr_s(1))
        write (*,'("ERROR in ",a)') rname
      end if
      status=-1; return
    end if

    ! ok
    status = 0
    
  end subroutine obj_CheckAttribute_s_1d



  ! ================================================================
  ! ===
  ! === write attributes
  ! ===
  ! ================================================================
  
  
  subroutine obj_WriteAttribute_s_0d( obj_id, name, s, status )
  
    use file_hdf_base, only : wpi

    ! --- in/out -------------------------
    
    integer(wpi), intent(in)            ::  obj_id
    character(len=*), intent(in)        ::  name
    character(len=*), intent(in)        ::  s
    integer, intent(out)                ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/obj_WriteAttribute_s_0d'
    
    ! --- external ----------------------------

    integer(wpi), external   ::  sfSCAtt

    ! --- begin -------------------------------
    
    status = sfSCAtt( obj_id, name, DFNT_CHAR, len(s), s )
    if ( status /= SUCCEED ) then
      write (*,'("ERROR - error writing attribute ",a)') trim(name)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! ok
    status = 0
    
  end subroutine obj_WriteAttribute_s_0d
    
  
  ! ***
  
  
  subroutine obj_WriteAttribute_s_1d( obj_id, name, s, status )
  
    use file_hdf_base, only : wpi

    ! --- in/out -------------------------
    
    integer(wpi), intent(in)            ::  obj_id
    character(len=*), intent(in)        ::  name
    character(len=*), intent(in)        ::  s(:)
    integer, intent(out)                ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/obj_WriteAttribute_s_1d'
    
    ! --- external ----------------------------

    integer(wpi), external   ::  sfSCAtt

    ! --- begin -------------------------------
    
    status = sfSCAtt( obj_id, name, DFNT_CHAR, len(s)*size(s), s )
    if ( status /= SUCCEED ) then
      write (*,'("ERROR - error writing attribute ",a)') trim(name)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! ok
    status = 0
    
  end subroutine obj_WriteAttribute_s_1d
    
  
  ! ############################################################
  ! ###
  ! ### scientific data sets
  ! ###
  ! ############################################################
  

  ! ================================================================
  ! get attributes
  ! ================================================================
   
  
  subroutine sds_ReadAttribute_s_0d( sds, name, s, status )
  
    use file_hdf_base, only : TSds

    ! --- in/out -------------------------
    
    type(Tsds), intent(in)              ::  sds
    character(len=*), intent(in)        ::  name
    character(len=*), intent(inout)     ::  s
    integer, intent(out)                ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sds_ReadAttribute_s_0d'
    
    ! --- begin -------------------------------
    
    call ReadAttribute( sds%id, name, s, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! ok
    status = 0
    
  end subroutine sds_ReadAttribute_s_0d
    
  
  ! ***
  
  
  subroutine sds_ReadAttribute_s_1d( sds, name, s, status )
  
    use file_hdf_base, only : TSds

    ! --- in/out -------------------------
    
    type(Tsds), intent(in)              ::  sds
    character(len=*), intent(in)        ::  name
    character(len=*), intent(out)       ::  s(:)
    integer, intent(inout)              ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sds_ReadAttribute_s_1d'
    
    ! --- begin -------------------------------
    
    call ReadAttribute( sds%id, name, s, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! ok
    status = 0
    
  end subroutine sds_ReadAttribute_s_1d

    
  
  ! =============================================================
  ! === check attributes
  ! =============================================================

  
  subroutine sds_CheckAttribute_s_0d( sds, name, s, status )
  
    use file_hdf_base, only : TSds

    ! --- in/out -------------------------
    
    type(TSds), intent(in)             ::  sds
    character(len=*), intent(in)       ::  name
    character(len=*), intent(in)       ::  s
    integer, intent(inout)             ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sds_CheckAttribute_s_0d'
    
    ! --- local ------------------------------
    
    logical          ::  verbose
    
    ! --- begin ---------------------------
    
    ! write error messages ?
    verbose = status == 0
    
    call CheckAttribute( sds%id, name, s, status )
    if (status>0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    if (status<0) then
      if (verbose) write (*,'("ERROR in ",a)') rname
      status=-1; return
    end if

    ! ok
    status = 0
    
  end subroutine sds_CheckAttribute_s_0d
  
  
  ! ***
  
  
  subroutine sds_CheckAttribute_s_1d( sds, name, s, status )
  
    use file_hdf_base, only : TSds

    ! --- in/out -------------------------
    
    type(TSds), intent(in)             ::  sds
    character(len=*), intent(in)       ::  name
    character(len=*), intent(in)       ::  s(:)
    integer, intent(inout)             ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sds_CheckAttribute_s_1d'
    
    ! --- local ------------------------------
    
    logical          ::  verbose
    
    ! --- begin ---------------------------
    
    ! write error messages ?
    verbose = status == 0
    
    call CheckAttribute( sds%id, name, s, status )
    if (status>0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    if (status<0) then
      if (verbose) write (*,'("ERROR in ",a)') rname
      status=-1; return
    end if

    ! ok
    status = 0
    
  end subroutine sds_CheckAttribute_s_1d
  
  
  ! ================================================================
  ! write attributes
  ! ================================================================
   
  
  subroutine sds_WriteAttribute_s_0d( sds, name, s, status )
  
    use file_hdf_base, only : TSds

    ! --- in/out -------------------------
    
    type(Tsds), intent(in)             ::  sds
    character(len=*), intent(in)       ::  name
    character(len=*), intent(in)       ::  s
    integer, intent(out)               ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sds_WriteAttribute_s_0d'
    
    ! --- begin -------------------------------
    
    call WriteAttribute( sds%id, name, s, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! ok
    status = 0
    
  end subroutine sds_WriteAttribute_s_0d
    
  
  ! ***
  
  
  subroutine sds_WriteAttribute_s_1d( sds, name, s, status )
  
    use file_hdf_base, only : TSds

    ! --- in/out -------------------------
    
    type(Tsds), intent(in)              ::  sds
    character(len=*), intent(in)        ::  name
    character(len=*), intent(in)        ::  s(:)
    integer, intent(out)                ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sds_WriteAttribute_s_1d'
    
    ! --- begin -------------------------------
    
    call WriteAttribute( sds%id, name, s, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! ok
    status = 0
    
  end subroutine sds_WriteAttribute_s_1d
    
  

  ! =============================================================
  ! === read data
  ! =============================================================

  
  subroutine sds_ReadData_s_1d( sds, data, status, start, stride )
  
    use file_hdf_base, only : wpi
    use file_hdf_base, only : TSds
    use file_hdf_base, only : CheckInfo, GetInfo

    ! --- const ------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sds_ReadData_s_1d'
    
    integer, parameter  ::  rank = 1
    
    ! --- in/out ----------------------------
    
    type(TSds), intent(in)              ::  sds
    character(len=*), intent(out)       ::  data
    integer, intent(out)                ::  status

    integer, intent(in), optional       ::  start(rank)
    integer, intent(in), optional       ::  stride(rank)
    
    ! --- local -------------------------------
    
    integer                             ::  data_type
    integer                             ::  the_start(rank)
    integer                             ::  the_stride(rank)
    integer                             ::  data_dims(rank)
    
    ! --- external ----------------------------
    
    integer(wpi), external  ::  sfRData

    ! --- begin -------------------------------
    
    ! check data rank and shape:
    !call CheckInfo( sds, data_rank=rank, data_dims=shape(data) )

    ! set dims etc
    the_start  = 0; if ( present(start ) ) the_start  = start
    the_stride = 1; if ( present(stride) ) the_stride = stride
    data_dims = (/len(data)/)

    ! read data of specified kind:
    call GetInfo( sds, status, data_type=data_type )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! read
    select case ( data_type )
      case ( DFNT_CHAR )
        status = sfRData( sds%id, the_start, the_stride, data_dims, data )
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
    
  end subroutine sds_ReadData_s_1d
  

  ! ***  


  subroutine sds_ReadData_s_2d( sds, data, status, start, stride )
  
    use file_hdf_base, only : wpi
    use file_hdf_base, only : TSds
    use file_hdf_base, only : wp_int8, wp_int16, wp_int32, wp_int64
    use file_hdf_base, only : CheckInfo, GetInfo

    ! --- const ------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sds_ReadData_s_2d'
    
    integer, parameter  ::  rank = 2
    
    ! --- in/out ----------------------------
    
    type(TSds), intent(in)              ::  sds
    character(len=*), intent(out)       ::  data(:)
    integer, intent(out)                ::  status

    integer, intent(in), optional       ::  start(rank)
    integer, intent(in), optional       ::  stride(rank)
    
    ! --- local -------------------------------
    
    integer                             ::  data_type
    integer                             ::  the_start(rank)
    integer                             ::  the_stride(rank)
    integer                             ::  data_dims(rank)
    
    ! --- external ----------------------------
    
    integer(wpi), external  ::  sfRData

    ! --- begin -------------------------------
    
    ! check data rank and shape:
    !call CheckInfo( sds, data_rank=rank, data_dims=shape(data) )

    ! set dims etc
    the_start  = 0; if ( present(start ) ) the_start  = start
    the_stride = 1; if ( present(stride) ) the_stride = stride
    data_dims = (/len(data),shape(data)/)

    ! read data of specified kind:
    call GetInfo( sds, status, data_type=data_type )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! read
    select case ( data_type )
      case ( DFNT_CHAR )
        status = sfRData( sds%id, the_start, the_stride, data_dims, data )
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
    
  end subroutine sds_ReadData_s_2d
  


  ! =============================================================
  ! === Write data
  ! =============================================================


  subroutine sds_WriteData_s_1d( sds, data, status, start, stride )
  
    use file_hdf_base, only : wpi
    use file_hdf_base, only : TSds

    ! --- const --------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sds_WriteData_s_1d'
    
    integer, parameter   ::  rank = 1
  
    ! --- in/out -------------------------
    
    type(TSds), intent(in)                    ::  sds
    character(len=*), intent(in)              ::  data
    integer, intent(out)                      ::  status

    integer, intent(in), optional             ::  start(rank)
    integer, intent(in), optional             ::  stride(rank)
    
    ! --- local -------------------------------
    
    integer               ::  the_start(rank), the_stride(rank)
    integer               ::  data_dims(rank)
    
    ! --- external ----------------------------
    
    integer(wpi), external  ::  sfWData
    
    ! --- begin -------------------------------
    
    !! check shape
    !if ( any( shape(data) /= sds%shp(1:sds%rnk) ) ) then
    !  print *, 'Shape of data does not match shape specified during creation:'
    !  print *, '  data        : ', shape(data)
    !  print *, '  created for : ', sds%shp(1:sds%rnk)
    !  stop 'FATAL ERROR IN sds_WriteData_i<wp>_1d'
    !end if

    ! write data:
    the_start  = 0; if ( present(start ) ) the_start  = start
    the_stride = 1; if ( present(stride) ) the_stride = stride
    data_dims = (/len(data)/)
    select case ( sds%typ )
      case ( 'chr' )
        status = sfWData( sds%id, the_start, the_stride, data_dims, data )
      case default
        write (*,'("ERROR - unknown sds%typ `",a,"`")') sds%typ
        write (*,'("ERROR in ",a)') rname; status=1; return
    end select
    if ( status == FAIL ) then
      write (*,'("ERROR - error writing data `",a,"`")') trim(sds%name)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! ok
    status = 0
    
  end subroutine sds_WriteData_s_1d
  
  
  ! ***

  
  subroutine sds_WriteData_s_2d( sds, data, status, start, stride )
  
    use file_hdf_base, only : wpi
    use file_hdf_base, only : TSds

    ! --- const --------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sds_WriteData_s_2d'

    integer, parameter   ::  rank = 2
  
    ! --- in/out -------------------------
    
    type(TSds), intent(in)                    ::  sds
    character(len=*), intent(in)              ::  data(:)
    integer, intent(out)                      ::  status

    integer, intent(in), optional             ::  start(rank)
    integer, intent(in), optional             ::  stride(rank)
    
    ! --- local -------------------------------
    
    integer               ::  the_start(rank), the_stride(rank)
    integer               ::  data_dims(rank)
    
    ! --- external ----------------------------
    
    integer(wpi), external  ::  sfWData
    
    ! --- begin -------------------------------
    
    !! check shape
    !if ( any( shape(data) /= sds%shp(1:sds%rnk) ) ) then
    !  print *, 'Shape of data does not match shape specified during creation:'
    !  print *, '  data        : ', shape(data)
    !  print *, '  created for : ', sds%shp(1:sds%rnk)
    !  stop 'FATAL ERROR IN sds_WriteData_i<wp>_2d'
    !end if

    ! write data:
    the_start  = 0; if ( present(start ) ) the_start  = start
    the_stride = 1; if ( present(stride) ) the_stride = stride
    data_dims = (/len(data),shape(data)/)
    select case ( sds%typ )
      case ( 'chr' )
        status = sfWData( sds%id, the_start, the_stride, data_dims, data )
      case default
        write (*,'("ERROR - unknown sds%typ `",a,"`")') sds%typ
        write (*,'("ERROR in ",a)') rname; status=1; return
    end select
    if ( status == FAIL ) then
      write (*,'("ERROR - error writing data `",a,"`")') trim(sds%name)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! ok
    status = 0
    
  end subroutine sds_WriteData_s_2d
  
  
  ! ***

  
  subroutine sds_WriteData_s_3d( sds, data, status, start, stride )
  
    use file_hdf_base, only : wpi
    use file_hdf_base, only : TSds

    ! --- const --------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sds_WriteData_s_3d'

    integer, parameter   ::  rank = 3
  
    ! --- in/out -------------------------
    
    type(TSds), intent(in)                    ::  sds
    character(len=*), intent(in)              ::  data(:,:)
    integer, intent(out)                      ::  status

    integer, intent(in), optional             ::  start(rank)
    integer, intent(in), optional             ::  stride(rank)
    
    ! --- local -------------------------------
    
    integer               ::  the_start(rank), the_stride(rank)
    integer               ::  data_dims(rank)
    
    ! --- external ----------------------------
    
    integer(wpi), external  ::  sfWData
    
    ! --- begin -------------------------------
    
    !! check shape
    !if ( any( shape(data) /= sds%shp(1:sds%rnk) ) ) then
    !  print *, 'Shape of data does not match shape specified during creation:'
    !  print *, '  data        : ', shape(data)
    !  print *, '  created for : ', sds%shp(1:sds%rnk)
    !  stop 'FATAL ERROR IN sds_WriteData_i<wp>_2d'
    !end if

    ! write data:
    the_start  = 0; if ( present(start ) ) the_start  = start
    the_stride = 1; if ( present(stride) ) the_stride = stride
    data_dims = (/len(data),shape(data)/)
    select case ( sds%typ )
      case ( 'chr' )
        status = sfWData( sds%id, the_start, the_stride, data_dims, data )
      case default
        write (*,'("ERROR - unknown sds%typ `",a,"`")') sds%typ
        write (*,'("ERROR in ",a)') rname; status=1; return
    end select
    if ( status == FAIL ) then
      write (*,'("ERROR - error writing data `",a,"`")') trim(sds%name)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! ok
    status = 0
    
  end subroutine sds_WriteData_s_3d
  
  
  ! ***

  
  subroutine sds_WriteData_s_4d( sds, data, status, start, stride )
  
    use file_hdf_base, only : wpi
    use file_hdf_base, only : TSds

    ! --- const --------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sds_WriteData_s_4d'

    integer, parameter   ::  rank = 4
  
    ! --- in/out -------------------------
    
    type(TSds), intent(in)                    ::  sds
    character(len=*), intent(in)              ::  data(:,:,:)
    integer, intent(out)                      ::  status

    integer, intent(in), optional             ::  start(rank)
    integer, intent(in), optional             ::  stride(rank)
    
    ! --- local -------------------------------
    
    integer               ::  the_start(rank), the_stride(rank)
    integer               ::  data_dims(rank)
    
    ! --- external ----------------------------
    
    integer(wpi), external  ::  sfWData
    
    ! --- begin -------------------------------
    
    !! check shape
    !if ( any( shape(data) /= sds%shp(1:sds%rnk) ) ) then
    !  print *, 'Shape of data does not match shape specified during creation:'
    !  print *, '  data        : ', shape(data)
    !  print *, '  created for : ', sds%shp(1:sds%rnk)
    !  stop 'FATAL ERROR IN sds_WriteData_i<wp>_2d'
    !end if

    ! write data:
    the_start  = 0; if ( present(start ) ) the_start  = start
    the_stride = 1; if ( present(stride) ) the_stride = stride
    data_dims = (/len(data),shape(data)/)
    select case ( sds%typ )
      case ( 'chr' )
        status = sfWData( sds%id, the_start, the_stride, data_dims, data )
      case default
        write (*,'("ERROR - unknown sds%typ `",a,"`")') sds%typ
        write (*,'("ERROR in ",a)') rname; status=1; return
    end select
    if ( status == FAIL ) then
      write (*,'("ERROR - error writing data `",a,"`")') trim(sds%name)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! ok
    status = 0
    
  end subroutine sds_WriteData_s_4d
  
  

  ! ############################################################
  ! ###
  ! ### dimensions
  ! ###
  ! ############################################################
  

  ! ================================================================
  ! get attributes
  ! ================================================================
   
  
  subroutine dim_ReadAttribute_s_0d( sdim, name, s, status )
  
    use file_hdf_base, only : TSdsDim

    ! --- in/out -------------------------
    
    type(TSdsDim), intent(in)           ::  sdim
    character(len=*), intent(in)        ::  name
    character(len=*), intent(inout)     ::  s
    integer, intent(out)                ::  status

    ! --- const --------------------------
    
    character(len=*), parameter  ::  rname = mname//'/dim_ReadAttribute_s_0d'

    ! --- begin -------------------------------
    
    call ReadAttribute( sdim%id, name, s, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! ok
    status = 0
    
  end subroutine dim_ReadAttribute_s_0d
    
  
  ! ***
  
  
  subroutine dim_ReadAttribute_s_1d( sdim, name, s, status )
  
    use file_hdf_base, only : TSdsDim

    ! --- in/out -------------------------
    
    type(TSdsDim), intent(in)           ::  sdim
    character(len=*), intent(in)        ::  name
    character(len=*), intent(out)       ::  s(:)
    integer, intent(out)                ::  status
    
    ! --- const --------------------------
    
    character(len=*), parameter  ::  rname = mname//'/dim_ReadAttribute_s_1d'

    ! --- begin -------------------------------
    
    call ReadAttribute( sdim%id, name, s, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! ok
    status = 0
    
  end subroutine dim_ReadAttribute_s_1d

    
  
  ! =============================================================
  ! === check attributes
  ! =============================================================

  
  subroutine dim_CheckAttribute_s_0d( sdim, name, s, status )
  
    use file_hdf_base, only : TSdsDim

    ! --- in/out -------------------------
    
    type(TSdsDim), intent(in)           ::  sdim
    character(len=*), intent(in)        ::  name
    character(len=*), intent(in)        ::  s
    integer, intent(inout)              ::  status
    
    ! --- const --------------------------
    
    character(len=*), parameter  ::  rname = mname//'/dim_CheckAttribute_s_0d'

    ! --- local ------------------------------
    
    logical          ::  verbose
    
    ! --- begin ---------------------------
    
    ! write error messages ?
    verbose = status == 0
    
    call CheckAttribute( sdim%id, name, s, status )
    if (status>0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    if (status<0) then
      if (verbose) write (*,'("ERROR in ",a)') rname
      status=-1; return
    end if

    ! ok
    status = 0
    
  end subroutine dim_CheckAttribute_s_0d
  
  
  ! ***
  
  
  subroutine dim_CheckAttribute_s_1d( sdim, name, s, status )
  
    use file_hdf_base, only : TSdsDim

    ! --- in/out -------------------------
    
    type(TSdsDim), intent(in)           ::  sdim
    character(len=*), intent(in)        ::  name
    character(len=*), intent(in)        ::  s(:)
    integer, intent(inout)              ::  status
    
    ! --- const --------------------------
    
    character(len=*), parameter  ::  rname = mname//'/dim_CheckAttribute_s_1d'

    ! --- local ------------------------------
    
    logical          ::  verbose
    
    ! --- begin ---------------------------
    
    ! write error messages ?
    verbose = status == 0
    
    call CheckAttribute( sdim%id, name, s, status )
    if (status>0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    if (status<0) then
      if (verbose) write (*,'("ERROR in ",a)') rname
      status=-1; return
    end if

    ! ok
    status = 0
    
  end subroutine dim_CheckAttribute_s_1d
  
  
  ! ================================================================
  ! write attributes
  ! ================================================================
   
  
  subroutine dim_WriteAttribute_s_0d( sdim, name, s, status )
  
    use file_hdf_base, only : TSdsDim

    ! --- in/out -------------------------
    
    type(TSdsDim), intent(in)           ::  sdim
    character(len=*), intent(in)        ::  name
    character(len=*), intent(in)        ::  s
    integer, intent(out)                ::  status
    
    ! --- const --------------------------
    
    character(len=*), parameter  ::  rname = mname//'/dim_WriteAttribute_s_0d'

    ! --- begin -------------------------------
    
    call WriteAttribute( sdim%id, name, s, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! ok
    status = 0
    
  end subroutine dim_WriteAttribute_s_0d
    
  
  ! ***
  
  
  subroutine dim_WriteAttribute_s_1d( sdim, name, s, status )
  
    use file_hdf_base, only : TSdsDim

    ! --- in/out -------------------------
    
    type(TSdsDim), intent(in)           ::  sdim
    character(len=*), intent(in)        ::  name
    character(len=*), intent(in)        ::  s(:)
    integer, intent(out)                ::  status
    
    ! --- const --------------------------
    
    character(len=*), parameter  ::  rname = mname//'/dim_WriteAttribute_s_1d'

    ! --- begin -------------------------------
    
    call WriteAttribute( sdim%id, name, s, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! ok
    status = 0
    
  end subroutine dim_WriteAttribute_s_1d
    

  
  ! ############################################################
  ! ###
  ! ### hdf files
  ! ###
  ! ############################################################
  

  ! ================================================================
  ! get attributes
  ! ================================================================
   
  
  subroutine hdf_ReadAttribute_s_0d( hdf, name, s, status )
  
    use file_hdf_base, only : THdfFile
  
    ! --- in/out -------------------------
    
    type(THdfFile), intent(in)          ::  hdf
    character(len=*), intent(in)        ::  name
    character(len=*), intent(inout)     ::  s
    integer, intent(out)                ::  status
    
    ! --- const --------------------------
    
    character(len=*), parameter  ::  rname = mname//'/hdf_ReadAttribute_s_0d'

    ! --- begin -------------------------------
    
    call ReadAttribute( hdf%id, name, s, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! ok
    status = 0
    
  end subroutine hdf_ReadAttribute_s_0d
    
  
  ! ***
  
  
  subroutine hdf_ReadAttribute_s_1d( hdf, name, s, status )
  
    use file_hdf_base, only : THdfFile
  
    ! --- in/out -------------------------
    
    type(THdfFile), intent(in)          ::  hdf
    character(len=*), intent(in)        ::  name
    character(len=*), intent(out)       ::  s(:)
    integer, intent(out)                ::  status
    
    ! --- const --------------------------
    
    character(len=*), parameter  ::  rname = mname//'/hdf_ReadAttribute_s_1d'

    ! --- begin -------------------------------
    
    call ReadAttribute( hdf%id, name, s, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! ok
    status = 0
    
  end subroutine hdf_ReadAttribute_s_1d

    
  
  ! =============================================================
  ! === check attributes
  ! =============================================================

  
  subroutine hdf_CheckAttribute_s_0d( hdf, name, s, status )
  
    use file_hdf_base, only : THdfFile
  
    ! --- in/out -------------------------
    
    type(THdfFile), intent(in)          ::  hdf
    character(len=*), intent(in)        ::  name
    character(len=*), intent(in)        ::  s
    integer, intent(inout)              ::  status
    
    ! --- const --------------------------
    
    character(len=*), parameter  ::  rname = mname//'/hdf_CheckAttribute_s_0d'

    ! --- local ------------------------------
    
    logical          ::  verbose
    
    ! --- begin ---------------------------
    
    ! write error messages ?
    verbose = status == 0
    
    call CheckAttribute( hdf%id, name, s, status )
    if (status>0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    if (status<0) then
      if (verbose) write (*,'("ERROR in ",a)') rname
      status=-1; return
    end if

    ! ok
    status = 0
    
  end subroutine hdf_CheckAttribute_s_0d
  
  
  ! ***
  
  
  subroutine hdf_CheckAttribute_s_1d( hdf, name, s, status )
  
    use file_hdf_base, only : THdfFile
  
    ! --- in/out -------------------------
    
    type(THdfFile), intent(in)          ::  hdf
    character(len=*), intent(in)        ::  name
    character(len=*), intent(in)        ::  s(:)
    integer, intent(inout)              ::  status
    
    ! --- const --------------------------
    
    character(len=*), parameter  ::  rname = mname//'/hdf_CheckAttribute_s_1d'

    ! --- local ------------------------------
    
    logical          ::  verbose
    
    ! --- begin ---------------------------
    
    ! write error messages ?
    verbose = status == 0
    
    call CheckAttribute( hdf%id, name, s, status )
    if (status>0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    if (status<0) then
      if (verbose) write (*,'("ERROR in ",a)') rname
      status=-1; return
    end if

    ! ok
    status = 0
    
  end subroutine hdf_CheckAttribute_s_1d
  
  
  ! ================================================================
  ! write attributes
  ! ================================================================
   
  
  subroutine hdf_WriteAttribute_s_0d( hdf, name, s, status )
  
    use file_hdf_base, only : THdfFile
  
    ! --- in/out -------------------------
    
    type(THdfFile), intent(in)          ::  hdf
    character(len=*), intent(in)        ::  name
    character(len=*), intent(in)        ::  s
    integer, intent(out)                ::  status
    
    ! --- const --------------------------
    
    character(len=*), parameter  ::  rname = mname//'/hdf_WriteAttribute_s_0d'

    ! --- begin -------------------------------
    
    call WriteAttribute( hdf%id, name, s, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! ok
    status = 0
    
  end subroutine hdf_WriteAttribute_s_0d
    
  
  ! ***
  
  
  subroutine hdf_WriteAttribute_s_1d( hdf, name, s, status )
  
    use file_hdf_base, only : THdfFile
  
    ! --- in/out -------------------------
    
    type(THdfFile), intent(in)          ::  hdf
    character(len=*), intent(in)        ::  name
    character(len=*), intent(in)        ::  s(:)
    integer, intent(out)                ::  status
    
    ! --- const --------------------------
    
    character(len=*), parameter  ::  rname = mname//'/hdf_WriteAttribute_s_1d'

    ! --- begin -------------------------------
    
    call WriteAttribute( hdf%id, name, s, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! ok
    status = 0
    
  end subroutine hdf_WriteAttribute_s_1d
    
  

end module file_hdf_s

 
