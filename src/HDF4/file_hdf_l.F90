module file_hdf_l

  implicit none
  
  ! --- in/out --------------------------
  
  private
  
  public   ::  ReadAttribute, CheckAttribute
  public   ::  WriteAttribute
  

  ! --- const ----------------------------
  
  include "hdf.f90"
  
  character(len=*), parameter  ::  mname = 'file_hdf_l'
    

  ! --- interfaces ------------------------
  
  interface ReadAttribute
    module procedure obj_ReadAttribute_l_0d
    module procedure obj_ReadAttribute_l_1d
    !
    module procedure sds_ReadAttribute_l_0d
    module procedure sds_ReadAttribute_l_1d
    !
    module procedure dim_ReadAttribute_l_0d
    module procedure dim_ReadAttribute_l_1d
    !
    module procedure hdf_ReadAttribute_l_0d
    module procedure hdf_ReadAttribute_l_1d
  end interface
  
  interface CheckAttribute
    module procedure obj_CheckAttribute_l_0d
    module procedure obj_CheckAttribute_l_1d
    !
    module procedure sds_CheckAttribute_l_0d
    module procedure sds_CheckAttribute_l_1d
    !
    module procedure dim_CheckAttribute_l_0d
    module procedure dim_CheckAttribute_l_1d
    !
    module procedure hdf_CheckAttribute_l_0d
    module procedure hdf_CheckAttribute_l_1d
  end interface
  
  interface WriteAttribute
    module procedure obj_WriteAttribute_l_0d
    module procedure obj_WriteAttribute_l_1d
    !
    module procedure sds_WriteAttribute_l_0d
    module procedure sds_WriteAttribute_l_1d
    !
    module procedure dim_WriteAttribute_l_0d
    module procedure dim_WriteAttribute_l_1d
    !
    module procedure hdf_WriteAttribute_l_0d
    module procedure hdf_WriteAttribute_l_1d
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
  
  
  subroutine obj_ReadAttribute_l_0d( obj_id, name, l, status )
  
    use file_hdf_base, only : wpi
    use file_hdf_base, only : FindAttribute, CheckAttributeInfo, GetAttributeInfo

    ! --- in/out -------------------------
    
    integer(wpi), intent(in)            ::  obj_id
    character(len=*), intent(in)        ::  name
    logical, intent(out)                ::  l
    integer, intent(out)                ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/obj_ReadAttribute_l_0d'
    
    ! --- local -------------------------------
    
    integer              ::  attr_index, data_type
    
    ! --- external ----------------------------

    integer(wpi), external   ::  sfRNAtt

    ! --- begin -------------------------------
    
    ! get index:
    call FindAttribute( obj_id, name, attr_index, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! check number of values:
    call CheckAttributeInfo( obj_id, attr_index, status, n_values=1 )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    
    ! extract value:
    call GetAttributeInfo( obj_id, attr_index, status, data_type=data_type )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! read:
    select case ( data_type )
      case ( DFNT_INT32 )
        status = sfRNAtt( obj_id, attr_index, 1, l  )
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
    
  end subroutine obj_ReadAttribute_l_0d


  ! ***
      
  
  subroutine obj_ReadAttribute_l_1d( obj_id, name, l, status )
  
    use file_hdf_base, only : wpi
    use file_hdf_base, only : FindAttribute, CheckAttributeInfo, GetAttributeInfo

    ! --- in/out -------------------------
    
    integer(wpi), intent(in)            ::  obj_id
    character(len=*), intent(in)        ::  name
    logical, intent(out)                ::  l(:)
    integer, intent(out)                ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/obj_ReadAttribute_l_1d'
    
    ! --- local -------------------------------
    
    integer              ::  attr_index, data_type
    
    ! --- external ----------------------------

    integer(wpi), external   ::  sfRNAtt

    ! --- begin -------------------------------
    
    ! get index:
    call FindAttribute( obj_id, name, attr_index, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! check number of values:
    call CheckAttributeInfo( obj_id, attr_index, status, n_values=size(l) )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    
    ! extract value:
    call GetAttributeInfo( obj_id, attr_index, status, data_type=data_type )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! read:
    select case ( data_type )
      case ( DFNT_INT32 )
        status = sfRNAtt( obj_id, attr_index, size(l), l )
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
    
  end subroutine obj_ReadAttribute_l_1d




  ! ================================================================
  ! ===
  ! === check attributes
  ! ===
  ! ================================================================
  
  
  subroutine obj_CheckAttribute_l_0d( obj_id, name, l, status )
  
    use file_hdf_base, only : wpi
  
    ! --- in/out -------------------------
    
    integer(wpi), intent(in)           ::  obj_id
    character(len=*), intent(in)       ::  name
    logical, intent(in)                ::  l
    integer, intent(inout)             ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/obj_CheckAttribute_l_0d'
    
    ! --- local -------------------------------
    
    logical        ::  verbose
    logical        ::  attr_l
    
    ! --- begin -------------------------------
    
    ! write error messages ?
    verbose = status == 0

    ! read data:
    call ReadAttribute( obj_id, name, attr_l, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    
    ! check
    if ( attr_l .neqv. l ) then
      if (verbose) then
        write (*,'("ERROR - foud different attribute values:")')
        write (*,'("ERROR -   attr name : ",a)') trim(name)
        write (*,'("ERROR -   requested : ",l2)') l
        write (*,'("ERROR -   found     : ",l2)') attr_l
        write (*,'("ERROR in ",a)') rname
      end if
      status=-1; return
    end if

    ! ok
    status = 0
    
  end subroutine obj_CheckAttribute_l_0d


  ! ***
      
  
  subroutine obj_CheckAttribute_l_1d( obj_id, name, l, status )
  
    use file_hdf_base, only : wpi
  
    ! --- in/out -------------------------
    
    integer(wpi), intent(in)            ::  obj_id
    character(len=*), intent(in)        ::  name
    logical, intent(in)                 ::  l(:)
    integer, intent(inout)              ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/obj_CheckAttribute_l_1d'
    
    ! --- local -------------------------------
    
    logical        ::  verbose
    logical        ::  attr_l(size(l))
    
    ! --- begin -------------------------------
    
    ! write error messages ?
    verbose = status == 0

    ! read data:
    call ReadAttribute( obj_id, name, attr_l, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    
    ! check
    if ( any( attr_l .neqv. l ) ) then
      if (verbose) then
        write (*,'("ERROR - foud different attribute values:")')
        write (*,'("ERROR -   attr name : ",a)') trim(name)
        write (*,'("ERROR -   requested : ",10l2)') l
        write (*,'("ERROR -   found     : ",10l2)') attr_l
        write (*,'("ERROR in ",a)') rname
      end if
      status=-1; return
    end if

    ! ok
    status = 0
    
  end subroutine obj_CheckAttribute_l_1d



  ! ================================================================
  ! ===
  ! === write attributes
  ! ===
  ! ================================================================
  
  
  subroutine obj_WriteAttribute_l_0d( obj_id, name, l, status )
  
    use file_hdf_base, only : wpi

    ! --- in/out -------------------------
    
    integer(wpi), intent(in)            ::  obj_id
    character(len=*), intent(in)        ::  name
    logical, intent(in)                 ::  l
    integer, intent(out)                ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/obj_WriteAttribute_l_0d'
    
    ! --- external ----------------------------

    integer(wpi), external   ::  sfSNAtt

    ! --- begin -------------------------------
    
    ! write attribute:
    status = sfSNAtt( obj_id, name, DFNT_INT32, 1, l )
    if ( status /= SUCCEED ) then
      write (*,'("ERROR - error writing attribute ",a)') trim(name)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine obj_WriteAttribute_l_0d
  
  
  ! ***
  
  
  subroutine obj_WriteAttribute_l_1d( obj_id, name, l, status )
  
    use file_hdf_base, only : wpi

    ! --- in/out -------------------------
    
    integer(wpi), intent(in)            ::  obj_id
    character(len=*), intent(in)        ::  name
    logical, dimension(:), intent(in)   ::  l
    integer, intent(out)                ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/obj_WriteAttribute_l_1d'
    
    ! --- external ----------------------------

    integer(wpi), external   ::  sfSNAtt

    ! --- begin -------------------------------
    
    ! write attribute:
    status = sfSNAtt( obj_id, name, DFNT_INT32, size(l), l )
    if ( status /= SUCCEED ) then
      write (*,'("ERROR - error writing attribute ",a)') trim(name)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine obj_WriteAttribute_l_1d
  
  
    
  ! ############################################################
  ! ###
  ! ### scientific data sets
  ! ###
  ! ############################################################
  

  ! ================================================================
  ! get attributes
  ! ================================================================
  
    
  subroutine sds_ReadAttribute_l_0d( sds, name, l, status )
  
    use file_hdf_base, only : TSds

    ! --- in/out -------------------------
    
    type(Tsds), intent(in)              ::  sds
    character(len=*), intent(in)        ::  name
    logical, intent(out)                ::  l
    integer, intent(out)                ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sds_ReadAttribute_l_0d'
    
    ! --- begin -------------------------------
    
    call ReadAttribute( sds%id, name, l, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    
    ! ok
    status = 0
    
  end subroutine sds_ReadAttribute_l_0d
  
  
  ! ***
  
  
  subroutine sds_ReadAttribute_l_1d( sds, name, l, status )
  
    use file_hdf_base, only : TSds

    ! --- in/out -------------------------
    
    type(Tsds), intent(in)              ::  sds
    character(len=*), intent(in)        ::  name
    logical, intent(out)                ::  l(:)
    integer, intent(out)                ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sds_ReadAttribute_l_1d'
    
    ! --- begin -------------------------------
    
    call ReadAttribute( sds%id, name, l, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    
    ! ok
    status = 0
    
  end subroutine sds_ReadAttribute_l_1d
  
  
  ! =============================================================
  ! === check attributes
  ! =============================================================

  
  subroutine sds_CheckAttribute_l_0d( sds, name, l, status )
  
    use file_hdf_base, only : TSds

    ! --- in/out -------------------------
    
    type(TSds), intent(in)             ::  sds
    character(len=*), intent(in)       ::  name
    logical, intent(in)                ::  l
    integer, intent(inout)             ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sds_CheckAttribute_l_0d'
    
    ! --- local ------------------------------
    
    logical          ::  verbose
    
    ! --- begin ---------------------------
    
    ! write error messages ?
    verbose = status == 0
    
    call CheckAttribute( sds%id, name, l, status )
    if (status>0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    if (status<0) then
      if (verbose) write (*,'("ERROR in ",a)') rname
      status=-1; return
    end if

    ! ok
    status = 0
    
  end subroutine sds_CheckAttribute_l_0d
  
  
  ! ***
  
  
  subroutine sds_CheckAttribute_l_1d( sds, name, l, status )
  
    use file_hdf_base, only : TSds

    ! --- in/out -------------------------
    
    type(TSds), intent(in)             ::  sds
    character(len=*), intent(in)       ::  name
    logical, intent(in)                ::  l(:)
    integer, intent(inout)             ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sds_CheckAttribute_l_1d'
    
    ! --- local ------------------------------
    
    logical          ::  verbose
    
    ! --- begin ---------------------------
    
    ! write error messages ?
    verbose = status == 0
    
    call CheckAttribute( sds%id, name, l, status )
    if (status>0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    if (status<0) then
      if (verbose) write (*,'("ERROR in ",a)') rname
      status=-1; return
    end if

    ! ok
    status = 0
    
  end subroutine sds_CheckAttribute_l_1d
  
  
  
  ! ================================================================
  ! write attributes
  ! ================================================================
  
    
  subroutine sds_WriteAttribute_l_0d( sds, name, l, status )
  
    use file_hdf_base, only : TSds

    ! --- in/out -------------------------
    
    type(Tsds), intent(in)              ::  sds
    character(len=*), intent(in)        ::  name
    logical, intent(in)                 ::  l
    integer, intent(out)                ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sds_WriteAttribute_l_0d'
    
    ! --- begin -------------------------------
    
    call WriteAttribute( sds%id, name, l, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    
    ! ok
    status = 0
    
  end subroutine sds_WriteAttribute_l_0d
  
  
  ! ***
  
  
  subroutine sds_WriteAttribute_l_1d( sds, name, l, status )
  
    use file_hdf_base, only : TSds

    ! --- in/out -------------------------
    
    type(Tsds), intent(in)              ::  sds
    character(len=*), intent(in)        ::  name
    logical, intent(in)                 ::  l(:)
    integer, intent(out)                ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sds_WriteAttribute_l_1d'
    
    ! --- begin -------------------------------
    
    call WriteAttribute( sds%id, name, l, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    
    ! ok
    status = 0
    
  end subroutine sds_WriteAttribute_l_1d
  
  
  ! ############################################################
  ! ###
  ! ### dimensions
  ! ###
  ! ############################################################
  

  ! ================================================================
  ! get attributes
  ! ================================================================
  
    
  subroutine dim_ReadAttribute_l_0d( sdim, name, l, status )
  
    use file_hdf_base, only : TSdsDim
  
    ! --- in/out -------------------------
    
    type(TSdsDim), intent(in)           ::  sdim
    character(len=*), intent(in)        ::  name
    logical, intent(out)                ::  l
    integer, intent(out)                ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/dim_ReadAttribute_l_0d'
    
    ! --- begin -------------------------------
    
    call ReadAttribute( sdim%id, name, l, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    
    ! ok
    status = 0
    
  end subroutine dim_ReadAttribute_l_0d
  
  
  ! ***
  
  
  subroutine dim_ReadAttribute_l_1d( sdim, name, l, status )
  
    use file_hdf_base, only : TSdsDim
  
    ! --- in/out -------------------------
    
    type(TSdsDim), intent(in)           ::  sdim
    character(len=*), intent(in)        ::  name
    logical, intent(out)                ::  l(:)
    integer, intent(out)                ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/dim_ReadAttribute_l_1d'
    
    ! --- begin -------------------------------
    
    call ReadAttribute( sdim%id, name, l, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    
    ! ok
    status = 0
    
  end subroutine dim_ReadAttribute_l_1d
  
  
  ! =============================================================
  ! === check attributes
  ! =============================================================

  
  subroutine dim_CheckAttribute_l_0d( sdim, name, l, status )
  
    use file_hdf_base, only : TSdsDim
  
    ! --- in/out -------------------------
    
    type(TSdsDim), intent(in)          ::  sdim
    character(len=*), intent(in)       ::  name
    logical, intent(in)                ::  l
    integer, intent(inout)             ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/dim_CheckAttribute_l_0d'
    
    ! --- local ------------------------------
    
    logical          ::  verbose
    
    ! --- begin ---------------------------
    
    ! write error messages ?
    verbose = status == 0
    
    call CheckAttribute( sdim%id, name, l, status )
    if (status>0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    if (status<0) then
      if (verbose) write (*,'("ERROR in ",a)') rname
      status=-1; return
    end if

    ! ok
    status = 0
    
  end subroutine dim_CheckAttribute_l_0d
  
  
  ! ***
  
  
  subroutine dim_CheckAttribute_l_1d( sdim, name, l, status )
  
    use file_hdf_base, only : TSdsDim
  
    ! --- in/out -------------------------
    
    type(TSdsDim), intent(in)          ::  sdim
    character(len=*), intent(in)       ::  name
    logical, intent(in)                ::  l(:)
    integer, intent(inout)             ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/dim_CheckAttribute_l_1d'
    
    ! --- local ------------------------------
    
    logical          ::  verbose
    
    ! --- begin ---------------------------
    
    ! write error messages ?
    verbose = status == 0
    
    call CheckAttribute( sdim%id, name, l, status )
    if (status>0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    if (status<0) then
      if (verbose) write (*,'("ERROR in ",a)') rname
      status=-1; return
    end if

    ! ok
    status = 0
    
  end subroutine dim_CheckAttribute_l_1d
  
  
  
  ! ================================================================
  ! write attributes
  ! ================================================================
  
    
  subroutine dim_WriteAttribute_l_0d( sdim, name, l, status )
  
    use file_hdf_base, only : TSdsDim
  
    ! --- in/out -------------------------
    
    type(TSdsDim), intent(in)          ::  sdim
    character(len=*), intent(in)       ::  name
    logical, intent(in)                ::  l
    integer, intent(out)                ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/dim_WriteAttribute_l_0d'
    
    ! --- begin -------------------------------
    
    call WriteAttribute( sdim%id, name, l, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    
    ! ok
    status = 0
    
  end subroutine dim_WriteAttribute_l_0d
  
  
  ! ***
  
  
  subroutine dim_WriteAttribute_l_1d( sdim, name, l, status )
  
    use file_hdf_base, only : TSdsDim
  
    ! --- in/out -------------------------
    
    type(TSdsDim), intent(in)           ::  sdim
    character(len=*), intent(in)        ::  name
    logical, intent(in)                 ::  l(:)
    integer, intent(out)                ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/dim_WriteAttribute_l_1d'
    
    ! --- begin -------------------------------
    
    call WriteAttribute( sdim%id, name, l, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    
    ! ok
    status = 0
    
  end subroutine dim_WriteAttribute_l_1d
  

  ! ############################################################
  ! ###
  ! ### hdf files
  ! ###
  ! ############################################################
  

  ! ================================================================
  ! get attributes
  ! ================================================================
  
    
  subroutine hdf_ReadAttribute_l_0d( hdf, name, l, status )
  
    use file_hdf_base, only : THdfFile
  
    ! --- in/out -------------------------
    
    type(THdfFile), intent(in)          ::  hdf
    character(len=*), intent(in)        ::  name
    logical, intent(out)                ::  l
    integer, intent(out)                ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/hdf_ReadAttribute_l_0d'
    
    ! --- begin -------------------------------
    
    call ReadAttribute( hdf%id, name, l, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    
  end subroutine hdf_ReadAttribute_l_0d
  
  
  ! ***
  
  
  subroutine hdf_ReadAttribute_l_1d( hdf, name, l, status )
  
    use file_hdf_base, only : THdfFile
  
    ! --- in/out -------------------------
    
    type(THdfFile), intent(in)          ::  hdf
    character(len=*), intent(in)        ::  name
    logical, intent(out)                ::  l(:)
    integer, intent(out)                ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/hdf_ReadAttribute_l_1d'
    
    ! --- begin -------------------------------
    
    call ReadAttribute( hdf%id, name, l, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    
    ! ok
    status = 0
    
  end subroutine hdf_ReadAttribute_l_1d
  
  
  ! =============================================================
  ! === check attributes
  ! =============================================================

  
  subroutine hdf_CheckAttribute_l_0d( hdf, name, l, status )
  
    use file_hdf_base, only : THdfFile
  
    ! --- in/out -------------------------
    
    type(THdfFile), intent(in)         ::  hdf
    character(len=*), intent(in)       ::  name
    logical, intent(in)                ::  l
    integer, intent(inout)             ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/hdf_CheckAttribute_l_0d'
    
    ! --- local ------------------------------
    
    logical          ::  verbose
    
    ! --- begin ---------------------------
    
    ! write error messages ?
    verbose = status == 0
    
    call CheckAttribute( hdf%id, name, l, status )
    if (status>0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    if (status<0) then
      if (verbose) write (*,'("ERROR in ",a)') rname
      status=-1; return
    end if

    ! ok
    status = 0
    
  end subroutine hdf_CheckAttribute_l_0d
  
  
  ! ***
  
  
  subroutine hdf_CheckAttribute_l_1d( hdf, name, l, status )
  
    use file_hdf_base, only : THdfFile
  
    ! --- in/out -------------------------
    
    type(THdfFile), intent(in)         ::  hdf
    character(len=*), intent(in)       ::  name
    logical, intent(in)                ::  l(:)
    integer, intent(inout)             ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/hdf_CheckAttribute_l_1d'
    
    ! --- local ------------------------------
    
    logical          ::  verbose
    
    ! --- begin ---------------------------
    
    ! write error messages ?
    verbose = status == 0
    
    call CheckAttribute( hdf%id, name, l, status )
    if (status>0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    if (status<0) then
      if (verbose) write (*,'("ERROR in ",a)') rname
      status=-1; return
    end if

    ! ok
    status = 0
    
  end subroutine hdf_CheckAttribute_l_1d
  
  
  
  ! ================================================================
  ! write attributes
  ! ================================================================
  
    
  subroutine hdf_WriteAttribute_l_0d( hdf, name, l, status )
  
    use file_hdf_base, only : THdfFile
  
    ! --- in/out -------------------------
    
    type(THdfFile), intent(in)         ::  hdf
    character(len=*), intent(in)       ::  name
    logical, intent(in)                ::  l
    integer, intent(inout)             ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/hdf_WriteAttribute_l_0d'
    
    ! --- begin -------------------------------
    
    call WriteAttribute( hdf%id, name, l, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    
    ! ok
    status = 0
    
  end subroutine hdf_WriteAttribute_l_0d
  
  
  ! ***
  
  
  subroutine hdf_WriteAttribute_l_1d( hdf, name, l, status )
  
    use file_hdf_base, only : THdfFile
  
    ! --- in/out -------------------------
    
    type(THdfFile), intent(in)          ::  hdf
    character(len=*), intent(in)        ::  name
    logical, intent(in)                 ::  l(:)
    integer, intent(out)                ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/hdf_WriteAttribute_l_1d'
    
    ! --- begin -------------------------------
    
    call WriteAttribute( hdf%id, name, l, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    
    ! ok
    status = 0
    
  end subroutine hdf_WriteAttribute_l_1d
  
  

end module file_hdf_l

 
