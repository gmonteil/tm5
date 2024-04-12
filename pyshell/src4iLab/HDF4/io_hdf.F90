!###############################################################################
!
! NAME
!   io_hdf
!
! Subroutines to read arrays from / write arrays to HDF 4 files
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

module io_hdf

  implicit none

  private

  public :: io_write
  public :: io_write3d_32d,  io_write2d_32d, io_write3d_32dr, io_write2d_32dr, io_write3d_32, io_write4d_32

  public :: io_read3d_32,    io_read4d_32
  public :: io_read2d_i16g,  io_read2d_i16d
  public :: io_read3d_32d,   io_read2d_32g,  io_read3d_32g
  public :: io_read4d_32g,   io_read4d_32d,  io_read3d_32dr
  public :: io_read2d_32d,   io_read2d_32dr

  public :: DFNT_INT8, DFNT_INT16, DFNT_INT32, DFNT_INT64
  public :: DFNT_FLOAT32, DFNT_FLOAT64
  public :: DFACC_READ, DFACC_CREATE, DFACC_WRITE, DFNT_CHAR
  public :: SUCCEED, FAIL

  !

  interface operator (+)
     module procedure upper_case
  end interface

  interface io_write
     module procedure io_write4d_32
     module procedure io_write3d_32
     module procedure io_write2d_32
     module procedure io_write1d_32
     module procedure io_write2d_i32
     module procedure io_write2d_i16
     module procedure io_write1d_i16
     module procedure io_write4d_i32
  end interface

  include 'hdf.f90'

  logical,parameter               :: okdebug = .false.
  integer,parameter               :: comp_type = 1 
  integer,parameter,dimension(1)  :: comp_prm = (/ 1 /)

contains



  subroutine io_write4d_32( sd_id, im, labim, jm, labjm, &
       lm, lablm, nt, labnt, &
       data, name, idate)
    implicit none

    ! in/out:
    integer,intent(in) :: im
    integer,intent(in) :: jm
    integer,intent(in) :: lm
    integer,intent(in) :: nt
    integer,intent(in) :: sd_id
    character(len=*),intent(in) :: name
    character(len=*),intent(in) :: labim
    character(len=*),intent(in) :: labjm
    character(len=*),intent(in) :: lablm
    character(len=*),intent(in) :: labnt
    integer,dimension(6),optional   :: idate
    real,dimension(im,jm,lm,nt),intent(in) :: data

    ! local
    integer            :: sfcreate, sfscompress, sfwdata, sfendacc
    integer            :: sfdimid, sfsdmname, sfsnatt
    integer            :: rank = 4    ,istat
    integer,dimension(4) :: start = (/ 0,0,0,0 /)
    integer,dimension(4) :: stride= (/1,1,1,1/)
    integer            :: sds_id, dimid0, dimid1, dimid2, dimid3

    sds_id = sfcreate(sd_id,name, DFNT_FLOAT32, rank, (/im,jm,lm,nt/) )
    if(present(idate)) istat =  sfsnatt(sds_id,'idate',  DFNT_INT32, 6, idate)
    dimid2 = sfdimid(sds_id, 0)
    istat = sfsdmname(dimid2,labim)
    dimid1 = sfdimid(sds_id, 1)
    istat = sfsdmname(dimid1,labjm)
    dimid0 = sfdimid(sds_id, 2)
    istat = sfsdmname(dimid0, lablm)
    dimid3 = sfdimid(sds_id, 3)
    istat = sfsdmname(dimid3, labnt)
    !istat = sfscompress(sds_id,comp_type,comp_prm)
    istat  = sfwdata( sds_id, (/0,0,0,0/), (/1,1,1,1/), (/im,jm,lm,nt/), &
         real(data,kind=4) )   
    istat = sfendacc(sds_id)

  end subroutine io_write4d_32



  subroutine io_write4d_i32( sd_id, im, labim, jm, labjm, &
       lm, lablm, nt, labnt, &
       data, name, idate )
    implicit none

    ! in/out:
    integer,intent(in)            :: im,jm,lm,nt
    integer,intent(in)            :: sd_id
    character(len=*),intent(in)   :: name
    character(len=*),intent(in)   :: labim,labjm,lablm,labnt
    integer,dimension(6),optional :: idate
    integer,dimension(im,jm,lm,nt),intent(in) :: data

    ! local
    integer              :: sfcreate, sfscompress, sfwdata, sfendacc
    integer              :: sfdimid, sfsdmname, sfsnatt
    integer              :: rank = 4    ,istat
    integer,dimension(4) :: start = (/ 0,0,0,0 /)
    integer,dimension(4) :: stride= (/1,1,1,1/)
    integer              :: sds_id, dimid0, dimid1, dimid2, dimid3

    sds_id = sfcreate(sd_id,name, DFNT_INT32, rank, (/im,jm,lm,nt/) )
    if(present(idate)) istat =  sfsnatt(sds_id,'idate',  DFNT_INT32, 6, idate)
    dimid2 = sfdimid(sds_id, 0)
    istat = sfsdmname(dimid2,labim)
    dimid1 = sfdimid(sds_id, 1)
    istat = sfsdmname(dimid1,labjm)
    dimid0 = sfdimid(sds_id, 2)
    istat = sfsdmname(dimid0, lablm)
    dimid3 = sfdimid(sds_id, 3)
    istat = sfsdmname(dimid3, labnt)
    !istat = sfscompress(sds_id,comp_type,comp_prm)
    istat  = sfwdata(sds_id, (/0,0,0,0/),(/1,1,1,1/) ,(/im,jm,lm,nt/) ,data)   
    istat = sfendacc(sds_id)

  end subroutine io_write4d_I32



  subroutine io_write3d_32(sd_id,im,labim,jm,labjm,lm,lablm,data,name,idate)
    implicit none

    ! in/out:
    integer,intent(in)                  :: im,jm,lm
    integer,intent(in)                  :: sd_id
    character(len=*),intent(in)         :: name
    character(len=*),intent(in)         :: labim,labjm,lablm
    real,dimension(im,jm,lm),intent(in) :: data
    integer,dimension(6),optional       :: idate

    ! local
    integer            :: sfcreate, sfscompress, sfwdata, sfendacc
    integer            :: sfdimid, sfsdmname, sfsnatt
    integer            :: rank = 3    ,istat
    integer,dimension(3) :: start = (/ 0,0,0 /), stride= (/1,1,1/)
    integer            :: sds_id, dimid0, dimid1, dimid2

    sds_id = sfcreate(sd_id,name, DFNT_FLOAT32, rank, (/im,jm,lm/) )
    if(present(idate)) istat =  sfsnatt(sds_id,'idate',  DFNT_INT32, 6, idate)
    dimid2 = sfdimid(sds_id, 0)
    istat = sfsdmname(dimid2,labim)
    dimid1 = sfdimid(sds_id, 1)
    istat = sfsdmname(dimid1,labjm)
    dimid0 = sfdimid(sds_id, 2)
    istat = sfsdmname(dimid0, lablm)
    !istat = sfscompress(sds_id,comp_type,comp_prm)
    istat  = sfwdata( sds_id, (/0,0,0/) ,(/1,1,1/) ,(/im,jm,lm/), &
         real(data,kind=4) )   
    istat = sfendacc(sds_id)

  end subroutine io_write3d_32



  subroutine io_read3d_32(sd_id,im,jm,lm,data,name,ifail)
    implicit none

    ! in/out:
    integer, intent(in)                  :: im,jm,lm
    integer, intent(in)                  :: sd_id
    character(len=*),intent(in)          :: name
    integer,intent(out)                  :: ifail
    real,dimension(im,jm,lm),intent(out) :: data

    ! local
    integer, parameter  :: MAX_VAR_DIMS = 32
    character(len=64)   :: xname
    integer             :: index, rank, sds_id
    integer             :: istat, attributes, num_type
    integer             :: sffinfo, sfselect, sfginfo
    integer             :: sfendacc, sfend, sfrnatt, sfrcatt
    integer             :: sffattr, sfrdata, sfn2index

    integer,dimension(MAX_VAR_DIMS) :: dim_sizes

    real(kind=4),dimension(:,:,:),allocatable :: hdfr

! introduce workaround to handle FLOAT64
! (was required because some mmix files were written in FLOAT64)
!PBi
    real(kind=8),dimension(:,:,:),allocatable :: hdfr_8
!PBf

    ifail = 1

    index = sfn2index(sd_id,name)
    if ( index == -1 ) return

    sds_id = sfselect(sd_id, index)
    istat = sfginfo(sds_id, xname, rank, dim_sizes, num_type, attributes)

!PBi
    if (num_type == DFNT_FLOAT32) then
      print *, 'FLOAT32'
    else if (num_type == DFNT_FLOAT64) then
      print *, 'FLOAT64'
    else
      print *, 'io_read3d_32: invalid data type', num_type
      return
    endif

!PBf


    write(*,*) 'io_read3d_32: name = ', name
    write(*,*) 'io_read3d_32: rank = ', rank
    write(*,*) 'io_read3d_32: dims = ', dim_sizes(1:rank)

    if(  rank == 3                .and. &
         dim_sizes(1) == im       .and. &
         dim_sizes(2) == jm       .and. &
         dim_sizes(3) == lm              ) then
!PBi
    if (num_type == DFNT_FLOAT32) then
!PBf
       allocate(hdfr(im,jm,lm))
       istat =    sfrdata(sds_id, (/0,0,0/),(/1,1,1/),(/im,jm,lm/),hdfr)
!PBi
    else if (num_type == DFNT_FLOAT64) then
       allocate(hdfr_8(im,jm,lm))
       istat =    sfrdata(sds_id, (/0,0,0/),(/1,1,1/),(/im,jm,lm/),hdfr_8)
    endif
!PBf
       if ( istat == SUCCEED ) then
          if(okdebug) print*, 'io_read3d_32: Successfully retrieved ' &
               //name//' from file'
       else
          print*, 'io_read3d_32: Failed to read '//name//' from file'
          return
       end if

!PBi
    if (num_type == DFNT_FLOAT32) then
!PBf
       data = hdfr
       istat = sfendacc(sds_id)
       deallocate(hdfr)
       ifail = 0
!PBi
    else if (num_type == DFNT_FLOAT64) then
       data = hdfr_8
       istat = sfendacc(sds_id)
       deallocate(hdfr_8)
       ifail = 0
    endif
!PBf

    end if

  end subroutine io_read3d_32



  subroutine io_read4d_32(sd_id,im,jm,lm,nt,data,name,ifail)
    implicit none

    ! in/out:
    integer, intent(in)            :: im,jm,lm,nt
    integer, intent(in)            :: sd_id
    character(len=*),intent(in)    :: name
    real,dimension(im,jm,lm,nt),intent(out)  :: data
    integer,intent(out)            :: ifail

    ! local
    integer, parameter   :: MAX_VAR_DIMS = 32
    character(len=64)    :: xname
    integer              :: index, rank, sds_id
    integer              :: istat, attributes, num_type
    integer              :: sffinfo, sfselect, sfginfo
    integer              :: sfendacc, sfend, sfrnatt, sfrcatt
    integer              :: sffattr, sfrdata, sfn2index
    integer,dimension(MAX_VAR_DIMS) :: dim_sizes

    real(kind=4),dimension(:,:,:,:),allocatable :: hdfr

    ifail = 1

    index = sfn2index(sd_id,name)
    if (index == -1) return

    sds_id = sfselect(sd_id, index)
    istat = sfginfo(sds_id, xname, rank, dim_sizes, num_type, attributes)

    write(*,*) 'io_read4d_32: name = ', name
    write(*,*) 'io_read4d_32: rank = ', rank
    write(*,*) 'io_read4d_32: dims = ', dim_sizes(1:rank)

    if(          rank == 4        .and. &
         dim_sizes(1) == im       .and. &
         dim_sizes(2) == jm       .and. &
         dim_sizes(3) == lm       .and. &
         dim_sizes(4) == nt      ) then
       allocate(hdfr(im,jm,lm,nt))
       istat =    sfrdata(sds_id, (/0,0,0,0/),(/1,1,1,1/),(/im,jm,lm,nt/),hdfr)
       if ( istat == SUCCEED ) then
          if ( okdebug ) print*, 'io_read4d_32: Successfully retrieved '// &
               name//' from file'
       else
          print*, 'io_read4d_32: Failed to read '//name//' from file'
          return
       end if

       data = hdfr
       istat = sfendacc(sds_id)
       deallocate(hdfr)
       ifail = 0
    end if

  end subroutine io_read4d_32


  subroutine io_write1d_i16(sd_id,im,labim,data,name,idate)
    implicit none

    ! in/out:
    integer,intent(in)            :: im
    integer,intent(in)            :: sd_id
    character(len=*),intent(in)   :: name, labim
    integer,dimension(6),optional :: idate
    integer,dimension(im),intent(in) :: data

    ! local
    integer            :: sfcreate, sfscompress, sfwdata, sfendacc
    integer            :: sfdimid, sfsdmname, sfsnatt
    integer            :: rank = 1    ,istat
    integer,dimension(1) :: start = (/ 0 /), stride= (/1/)
    integer            :: sds_id, dimid1, dimid2

    sds_id = sfcreate(sd_id,name, DFNT_INT16, rank, (/im/) )
    if(present(idate)) istat =  sfsnatt(sds_id,'idate',  DFNT_INT32, 6, idate)
    dimid2 = sfdimid(sds_id, 0)
    istat = sfsdmname(dimid2,labim)
    !istat = sfscompress(sds_id,comp_type,comp_prm)
    istat  = sfwdata(sds_id, start, stride ,(/im/) , int(data,kind=2))   
    istat = sfendacc(sds_id)

  end subroutine io_write1d_i16



  subroutine io_write1d_32(sd_id,im,labim,data,name,idate)
    implicit none

    ! in/out:
    integer,intent(in)            :: im
    integer,intent(in)            :: sd_id
    character(len=*),intent(in)   :: name, labim
    integer,dimension(6),optional :: idate
    real,dimension(im),intent(in) :: data

    ! local
    integer            :: sfcreate, sfscompress, sfwdata, sfendacc
    integer            :: sfdimid, sfsdmname, sfsnatt
    integer            :: rank = 1    ,istat
    integer,dimension(1) :: start = (/ 0 /), stride= (/1/)
    integer            :: sds_id, dimid1, dimid2

    sds_id = sfcreate(sd_id,name, DFNT_FLOAT32, rank, (/im/) )
    if(present(idate)) istat =  sfsnatt(sds_id,'idate',  DFNT_INT32, 6, idate)
    dimid2 = sfdimid(sds_id, 0)
    istat = sfsdmname(dimid2,labim)
    !istat = sfscompress(sds_id,comp_type,comp_prm)
    istat  = sfwdata(sds_id, start, stride ,(/im/) , real(data,kind=4))   
    istat = sfendacc(sds_id)

  end subroutine io_write1d_32



  subroutine io_write2d_32(sd_id,im,labim,jm,labjm,data,name,idate)
    implicit none

    ! in/out:
    integer,intent(in)               :: im,jm
    integer,intent(in)               :: sd_id
    character(len=*),intent(in)      :: name,labim,labjm
    integer,dimension(6),optional    :: idate
    real,dimension(im,jm),intent(in) :: data

    ! local
    integer              :: sfcreate, sfscompress, sfwdata, sfendacc
    integer              :: sfdimid, sfsdmname, sfsnatt
    integer              :: rank = 2    ,istat
    integer,dimension(2) :: start = (/ 0,0 /), stride= (/1,1/)
    integer              :: sds_id, dimid1, dimid2

    sds_id = sfcreate(sd_id,name, DFNT_FLOAT32, rank, (/im,jm/) )
    if(present(idate)) istat =  sfsnatt(sds_id,'idate',  DFNT_INT32, 6, idate)
    dimid2 = sfdimid(sds_id, 0)
    istat = sfsdmname(dimid2,labim)
    dimid1 = sfdimid(sds_id, 1)
    istat = sfsdmname(dimid1,labjm)
    !istat = sfscompress(sds_id,comp_type,comp_prm)
    istat  = sfwdata(sds_id, start, stride ,(/im,jm/) , real(data,kind=4))   
    istat = sfendacc(sds_id)

  end subroutine io_write2d_32



  subroutine io_read2d_I16g(sd_id,im,jm,data,name,ifail,index,idate)
    implicit none

    ! in/out:
    integer, intent(in)                  :: im,jm
    integer, intent(in)                  :: sd_id
    character(len=*),intent(in)          :: name
    integer,dimension(im,jm),intent(out) :: data
    integer,intent(out)                  :: ifail
    integer,dimension(6), optional       :: idate

    ! local
    integer,dimension(6) :: idate_file

    integer, parameter   :: MAX_VAR_DIMS = 32
    character(len=64)    :: xname,attr_name
    integer,optional     :: index
    integer(kind=4)      :: rank
    integer              :: sds_id
    integer              :: istat, attributes, num_type, n_values, data_type
    integer              :: sffinfo, sfselect, sfginfo, sfrattr
    integer              :: sfendacc, sfend, sfrnatt, sfrcatt
    integer              :: sffattr, sfrdata, sfn2index
    integer              :: sfgainfo,attr_index,lname
    integer              :: num_ds,num_at,ids
    integer              :: idebug = 0
    integer(kind=4),dimension(MAX_VAR_DIMS):: dim_sizes

    ifail = 1

    lname = len_trim(name)
    istat = sffinfo(sd_id,num_ds,num_at)

    if( present(index) ) then
       sds_id = sfselect(sd_id, index)   
       istat = sfginfo(sds_id, xname, rank, dim_sizes, num_type, attributes)
       if ( idebug >= 100 ) then
          write(*,*) 'io_read2d_I16g: name = ', name,xname(1:lname)
          write(*,*) 'io_read2d_I16g: rank = ', rank
          write(*,*) 'io_read2d_I16g: dims = ', dim_sizes(1:rank)
          write(*,*) 'io_read2d_I16g: attr = ', attributes
       end if
       ! check rank and name....
       names: if(+xname(1:lname) == +name(1:lname) .and.   &
            rank == 2        .and. &
            dim_sizes(1) == im       .and. &
            dim_sizes(2) == jm      ) then

          if( present(idate) ) then
             attr_index = sffattr(sds_id, 'idate')
             istat = sfgainfo(sds_id, attr_index, attr_name, &
                  data_type,  n_values)
             istat = sfrattr(sds_id, attr_index, idate_file)
             if ( idebug >= 100 ) then
                write(*,*) 'io_read2d_I16g: idate = ', idate,idate_file
             end if
             if( sum(abs(idate-idate_file)) == 0 ) then
                call read_it   !everything OK....proceed
                !cmk index=index+1   !set index to next position 
                return         !and return
             end if
          else
             call read_it   !everything OK....proceed
             !cmk index=index+1   !set index to next position
             return         !and return
          end if
       end if names
       istat = sfendacc(sds_id)   !close 'wrong' ds
    end if

    dsloop: do ids=0,num_ds-1
       sds_id = sfselect(sd_id, ids)   
       istat = sfginfo(sds_id, xname, rank, dim_sizes, num_type, attributes)
       if ( idebug >= 100 ) then
          write(*,*) 'io_read2d_I16g: name = ', name,xname(1:lname)
          write(*,*) 'io_read2d_I16g: rank = ', rank
          write(*,*) 'io_read2d_I16g: dims = ', dim_sizes(1:rank)
          write(*,*) 'io_read2d_I16g: attr = ', attributes
       end if
       if(+xname(1:lname) /= +name(1:lname) .or.   &   !check rank and name....
            rank /= 2        .or. &
            dim_sizes(1) /= im       .or. &
            dim_sizes(2) /= jm      ) then

          istat = sfendacc(sds_id)   !close 'wrong' ds
          cycle dsloop
       end if

       if(present (idate)) then
          attr_index = sffattr(sds_id, 'idate')
          istat = sfgainfo(sds_id, attr_index, attr_name, data_type,  n_values)
          istat = sfrattr(sds_id, attr_index, idate_file)
          if ( idebug >= 100 ) then
             write(*,*) 'io_read2d_I16g: idate = ', idate,idate_file
          end if
          if( sum(abs(idate-idate_file)) /= 0 ) then 
             istat = sfendacc(sds_id)   !close 'wrong' ds
             cycle dsloop
          end if
       end if
       call read_it     !everything OK
       ! set index to next position if present...
       if ( present(index) ) index=ids+1
       return
    end do dsloop

    print*, 'io_read2d_I16g: Could not find '//name//' from file'
    if ( present(idate) ) print *,'io_read2d_I16g: With date..',idate

  contains 

    subroutine read_it
      integer :: sfrdata
      integer(kind=2),dimension(:,:),allocatable :: hdfr

      allocate(hdfr(im,jm))
      istat =    sfrdata(sds_id, (/0,0/),(/1,1/),(/im,jm/),hdfr)
      if ( istat == SUCCEED ) then
         if ( idebug >= 100 ) &
              print*, 'io_read2d_I16g: Successfully retrieved '// &
              name//' from file'
      else
         print*, 'io_read2d_I16g: Failed to read '//name//' from file'
         return
      end if
      data = hdfr
      istat = sfendacc(sds_id)
      deallocate(hdfr)
      ifail = 0
    end subroutine read_it

  end subroutine io_read2d_I16g



  subroutine io_write2d_I16(sd_id,im,jm,data,name,idate)

    implicit none

    ! in/out:
    integer,intent(in)                          :: im,jm
    integer,intent(in)                          :: sd_id
    character(len=*),intent(in)                 :: name
    integer(kind=2),dimension(im,jm),intent(in) :: data
    integer,dimension(6),optional               :: idate

    ! local
    integer              :: sfcreate, sfscompress, sfwdata, sfendacc, sfsnatt
    integer              :: rank = 2    ,istat
    integer,dimension(2) :: start = (/ 0,0 /), stride= (/1,1/)
    integer              :: sds_id

    sds_id = sfcreate(sd_id, name, DFNT_INT16, rank, (/im,jm/) )
    if( present(idate) ) istat =  sfsnatt(sds_id,'idate',  DFNT_INT32, 6, idate)
    !print*,'sfcreate returns dataset id:',sds_id
    !istat = sfscompress(sds_id,comp_type,comp_prm)
    !print*,'sfscompress returns:',istat
    istat  = sfwdata(sds_id, start, stride, (/im,jm/) , data)
    !print*,'sfwdata returns:',istat
    istat = sfendacc(sds_id)
    !print*, 'sfendacc returns: ', istat

  end subroutine io_write2d_I16



  subroutine io_write2d_I32(sd_id,im,labim,jm,labjm,data,name,idate)
    implicit none

    ! in/out:
    integer,intent(in)                  :: im,jm
    integer,intent(in)                  :: sd_id
    character(len=*),intent(in)         :: name,labim,labjm
    integer,dimension(im,jm),intent(in) :: data
    integer,dimension(6),optional       :: idate

    ! local
    integer            :: sfcreate, sfscompress, sfwdata, sfendacc
    integer            :: sfdimid, sfsdmname, sfsnatt
    integer            :: rank = 2    ,istat
    integer,dimension(2) :: start = (/ 0,0 /), stride= (/1,1/)
    integer            :: sds_id, dimid1, dimid2

    sds_id = sfcreate(sd_id,name, DFNT_INT32, rank, (/im,jm/) )
    if(present(idate)) istat =  sfsnatt(sds_id,'idate',  DFNT_INT32, 6, idate)
    dimid2 = sfdimid(sds_id, 0)
    istat = sfsdmname(dimid2,labim)
    dimid1 = sfdimid(sds_id, 1)
    istat = sfsdmname(dimid1,labjm)
    !istat = sfscompress(sds_id,comp_type,comp_prm)
    istat  = sfwdata(sds_id, start, stride ,(/im,jm/) , data)   
    istat = sfendacc(sds_id)

  end subroutine io_write2d_I32



  subroutine io_read2d_I16D(sd_id,im,jm,data,name,index,ifail,idate)
    implicit none

    ! in/out:
    integer,intent(in)              :: im,jm
    integer,intent(in)              :: sd_id
    character(len=*),intent(in)     :: name
    integer,dimension(im,jm),intent(out) :: data
    integer,intent(inout)           :: index
    integer,intent(out)             :: ifail
    integer,dimension(6),intent(in) :: idate

    ! local
    integer, parameter              :: MAX_VAR_DIMS = 32
    integer,dimension(6)            :: idate_file
    integer                         :: sffinfo, sfselect, sfginfo, sfrdata
    integer                         :: sfendacc, sffattr, sfgainfo, sfrattr
    integer                         :: rank = 2    ,istat
    integer,dimension(2)            :: start = (/ 0,0 /), stride= (/1,1/)
    integer(kind=2),dimension(:,:),allocatable :: hdfr
    integer                         :: sds_id, ndatasets, nglobat, i
    integer                         :: xrank, xtype, natt, j
    integer,dimension(MAX_VAR_DIMS) :: dim_sizes
    character(len=64)     :: xname,attr_name
    integer                         :: num_type, n_values, data_type
    integer                         :: attributes, attr_index

    ifail = 1

    sds_id = sfselect(sd_id, index)
    istat = sfginfo(sds_id, xname, rank, dim_sizes, num_type, attributes)
    write(*,*) 'io_read2d_I16D: name = ', name
    write(*,*) 'io_read2d_I16D: rank = ', rank
    write(*,*) 'io_read2d_I16D: dims = ', dim_sizes(1:rank)
    write(*,*) 'io_read2d_I16D: attr = ', attributes

    attr_index = sffattr(sds_id, 'idate')
    istat = sfgainfo(sds_id, attr_index, attr_name, data_type,  n_values)
    istat = sfrattr(sds_id, attr_index, idate_file)
    if( sum(abs(idate_file-idate)) == 0           .and. &
                              rank == 2           .and. &
                        +name(1:4) == +xname(1:4) .and. &
                      dim_sizes(1) == im          .and. &
                      dim_sizes(2) == jm                  ) then
       allocate(hdfr(im,jm))
       istat =  sfrdata(sds_id, start,stride,(/im,jm/),hdfr)
       if ( istat == SUCCEED ) then
          if( okdebug ) print*, 'io_read2d_I16D: Successfully retrieved '// &
               name//' from file'
       else
          print*, 'io_read2d_I16D: Failed to read '//name//' from file'
          return
       end if
       data = hdfr
       istat = sfendacc(sds_id)
       deallocate(hdfr)
       ifail = 0
       index = index+1
    end if

  end subroutine io_read2d_I16D



  subroutine io_write3d_32d(sd_id,im,labim,jm,labjm,lm,lablm,data,name,idate)
    implicit none

    ! in/out:
    integer,intent(in)                  :: im,jm,lm
    integer,intent(in)                  :: sd_id
    character(len=*),intent(in)         :: name,labim,labjm,lablm
    real,dimension(im,jm,lm),intent(in) :: data
    integer,dimension(6),intent(in)     :: idate

    ! local
    integer              :: sfcreate, sfscompress, sfwdata, sfendacc
    integer              :: sfdimid, sfsdmname, sfsnatt
    integer              :: rank = 3    ,istat
    integer,dimension(3) :: start = (/ 0,0,0 /), stride= (/1,1,1/)
    integer              :: sds_id, dimid0, dimid1, dimid2

    sds_id = sfcreate(sd_id,name, DFNT_FLOAT32, rank, (/im,jm,lm/) )
    istat =  sfsnatt(sds_id,'idate',  DFNT_INT32, 6, idate)
    dimid2 = sfdimid(sds_id, 0)
    istat = sfsdmname(dimid2,labim)
    dimid1 = sfdimid(sds_id, 1)
    istat = sfsdmname(dimid1,labjm)
    dimid0 = sfdimid(sds_id, 2)
    istat = sfsdmname(dimid0, lablm)
    !istat = sfscompress(sds_id,comp_type,comp_prm)
    istat  = sfwdata( sds_id, (/0,0,0/), (/1,1,1/), (/im,jm,lm/), &
         real(data,kind=4) )   
    istat = sfendacc(sds_id)

  end subroutine io_write3d_32d



  subroutine io_write3d_32dr(sd_id,im,labim,jm,labjm,lm,lablm, &
       data,name,idate,region)
    implicit none

    ! in/out:
    integer,intent(in)                  :: im,jm,lm
    integer,intent(in)                  :: sd_id
    character(len=*),intent(in)         :: name,labim,labjm,lablm
    real,dimension(im,jm,lm),intent(in) :: data
    integer,dimension(6),intent(in)     :: idate
    integer,intent(in)                  :: region

    ! local
    integer            :: sfcreate, sfscompress, sfwdata, sfendacc
    integer            :: sfdimid, sfsdmname, sfsnatt
    integer            :: rank = 3    ,istat
    integer,dimension(3) :: start = (/ 0,0,0 /), stride= (/1,1,1/)
    integer            :: sds_id, dimid0, dimid1, dimid2

    sds_id = sfcreate(sd_id,name, DFNT_FLOAT32, rank, (/im,jm,lm/) )
    istat =  sfsnatt(sds_id,'idate',  DFNT_INT32, 6, idate)
    istat =  sfsnatt(sds_id,'region',  DFNT_INT32, 1, region)
    dimid2 = sfdimid(sds_id, 0)
    istat = sfsdmname(dimid2,labim)
    dimid1 = sfdimid(sds_id, 1)
    istat = sfsdmname(dimid1,labjm)
    dimid0 = sfdimid(sds_id, 2)
    istat = sfsdmname(dimid0, lablm)
    !istat = sfscompress(sds_id,comp_type,comp_prm)
    istat  = sfwdata( sds_id, (/0,0,0/), (/1,1,1/), (/im,jm,lm/), &
         real(data,kind=4) )   
    istat = sfendacc(sds_id)

  end subroutine io_write3d_32dr



  subroutine io_write2d_32d(sd_id,im,labim,jm,labjm,data,name,idate)
    implicit none

    ! in/out:
    integer,intent(in)               :: im,jm
    integer,intent(in)               :: sd_id
    character(len=*),intent(in)      :: name,labim,labjm
    real,dimension(im,jm),intent(in) :: data
    integer,dimension(6),intent(in)  :: idate

    ! local
    integer            :: sfcreate, sfscompress, sfwdata, sfendacc
    integer            :: sfdimid, sfsdmname, sfsnatt
    integer            :: rank = 2
    integer            :: istat
    integer,dimension(2) :: start = (/ 0,0 /), stride= (/1,1/)
    integer            :: sds_id, dimid1, dimid2

    sds_id = sfcreate(sd_id,name, DFNT_FLOAT32, rank, (/im,jm/) )
    istat =  sfsnatt(sds_id,'idate',  DFNT_INT32, 6, idate)
    dimid2 = sfdimid(sds_id, 0)
    istat = sfsdmname(dimid2,labim)
    dimid1 = sfdimid(sds_id, 1)
    istat = sfsdmname(dimid1,labjm)
    !istat = sfscompress(sds_id,comp_type,comp_prm)
    istat  = sfwdata(sds_id, start, stride, (/im,jm/), real(data,kind=4))   
    istat = sfendacc(sds_id)

  end subroutine io_write2d_32d



  subroutine io_write2d_32dr(sd_id,im,labim,jm,labjm,data,name,idate,region)
    implicit none

    ! in/out:
    integer,intent(in)               :: im,jm
    integer,intent(in)               :: sd_id
    character(len=*),intent(in)      :: name,labim,labjm
    integer,dimension(6),intent(in)  :: idate
    real,dimension(im,jm),intent(in) :: data
    integer,intent(in)               :: region

    ! local
    integer            :: sfcreate, sfscompress, sfwdata, sfendacc
    integer            :: sfdimid, sfsdmname, sfsnatt
    integer            :: rank = 2    ,istat
    integer,dimension(2) :: start = (/ 0,0 /), stride= (/1,1/)
    integer            :: sds_id, dimid1, dimid2

    sds_id = sfcreate(sd_id,name, DFNT_FLOAT32, rank, (/im,jm/) )
    istat =  sfsnatt(sds_id,'idate',  DFNT_INT32, 6, idate)
    istat =  sfsnatt(sds_id,'region',  DFNT_INT32, 1, region)
    dimid2 = sfdimid(sds_id, 0)
    istat = sfsdmname(dimid2,labim)
    dimid1 = sfdimid(sds_id, 1)
    istat = sfsdmname(dimid1,labjm)
    !istat = sfscompress(sds_id,comp_type,comp_prm)
    istat  = sfwdata(sds_id, start, stride ,(/im,jm/) , real(data,kind=4))   
    istat = sfendacc(sds_id)

  end subroutine io_write2d_32dr



  subroutine io_read3d_32d(sd_id,im,jm,lm,data,name,index,ifail,idate)
    implicit none

    ! in/out:
    integer, intent(in)                  :: im,jm,lm
    integer, intent(in)                  :: sd_id
    real,dimension(im,jm,lm),intent(out) :: data
    character(len=*),intent(in)          :: name
    integer,intent(inout)                :: index
    integer,intent(out)                  :: ifail
    integer,dimension(6),intent(in)      :: idate

    ! local
    integer,dimension(6)           :: idate_file

    integer, parameter             :: MAX_VAR_DIMS = 32
    character(len=64)              :: xname,attr_name

    integer(kind=4)                :: rank
    integer                        :: sds_id
    integer                        :: istat, attributes, num_type
    integer                        :: n_values, data_type
    integer                        :: sffinfo, sfselect, sfginfo, sfrattr
    integer                        :: sfendacc, sfend, sfrnatt, sfrcatt
    integer                        :: sffattr, sfrdata, sfn2index
    integer                        :: sfgainfo,attr_index,lname
    integer                        :: num_ds,num_at,ids
    integer(kind=4),dimension(MAX_VAR_DIMS) :: dim_sizes

    real(kind=4),dimension(:,:,:),allocatable :: hdfr

    ifail = 1

    sds_id = sfselect(sd_id, index)    !first try the suggested index.....
    istat = sfginfo(sds_id, xname, rank, dim_sizes, num_type, attributes)
    lname = len_trim(name)
    write(*,*) 'io_read3d_32d: name = ', name,xname(1:lname)
    write(*,*) 'io_read3d_32d: rank = ', rank
    write(*,*) 'io_read3d_32d: dims = ', dim_sizes(1:rank)
    write(*,*) 'io_read3d_32d: attr = ', attributes

    attr_index = sffattr(sds_id, 'idate')
    istat = sfgainfo(sds_id, attr_index, attr_name, data_type,  n_values)
    istat = sfrattr(sds_id, attr_index, idate_file)

    write(*,*) 'io_read3d_32d: idate = ', idate,idate_file

    if(  rank == 3        .and. &
         +xname(1:lname) == +name(1:lname) .and. &
         sum(abs(idate_file-idate)) == 0   .and. &
         dim_sizes(1) == im       .and. &
         dim_sizes(2) == jm       .and. &
         dim_sizes(3) == lm      ) then

       allocate(hdfr(im,jm,lm))
       istat =    sfrdata(sds_id, (/0,0,0/),(/1,1,1/),(/im,jm,lm/),hdfr)
       if ( istat == SUCCEED ) then
          if ( okdebug ) print*, 'io_read3d_32d: Successfully retrieved '// &
               name//' from file'
       else
          print*, 'io_read3d_32d: Failed to read '//name//' from file'
          return
       end if
       data = hdfr
       istat = sfendacc(sds_id)
       deallocate(hdfr)
       ifail = 0
       index = index+1

    else

       print *, 'io_read3d_32d: Try to find '//name//' with date ',idate
       istat = sfendacc(sds_id)   !close 'wrong' ds
       istat = sffinfo(sd_id,num_ds,num_at)
       do ids = 0,num_ds-1
          sds_id = sfselect(sd_id, ids)
          istat = sfginfo(sds_id, xname, rank, dim_sizes, num_type, attributes)
          if ( +xname(1:lname) == +name(1:lname) .and.   &
               rank == 3        .and. &
               dim_sizes(1) == im       .and. &
               dim_sizes(2) == jm       .and. &
               dim_sizes(3) == lm      ) then
             attr_index = sffattr(sds_id, 'idate')
             istat = sfgainfo(sds_id, attr_index, attr_name, data_type, n_values)
             istat = sfrattr(sds_id, attr_index, idate_file)
             if ( sum(abs(idate-idate_file)) == 0 ) then
                allocate(hdfr(im,jm,lm))
                istat = sfrdata(sds_id, (/0,0,0/),(/1,1,1/),(/im,jm,lm/),hdfr)
                if ( istat == SUCCEED ) then
                   print*, 'io_read3d_32d: Successfully retrieved '//name// &
                        ' from file'
                else
                   print*, 'io_read3d_32d: Failed to read '//name//' from file'
                   return
                end if
                data = hdfr
                istat = sfendacc(sds_id)
                deallocate(hdfr)
                ifail = 0
                index = ids + 1
                return
             end if  !date fit
             istat = sfendacc(sds_id)
          end if   !name fit
       end do  !ids loop
       print*, 'io_read3d_32d: Failed to read '//name//' from file'
    end if

  end subroutine io_read3d_32d



  subroutine io_read2d_32g(sd_id,im,jm,data,name,ifail,index,idate)
    implicit none

    ! in/out:
    integer, intent(in)               :: im,jm
    integer, intent(in)               :: sd_id
    character(len=*),intent(in)       :: name
    real,dimension(im,jm),intent(out) :: data
    integer,intent(out)               :: ifail
    integer,dimension(6), optional    :: idate
    integer,optional                  :: index

    ! local
    integer,dimension(6) :: idate_file

    integer, parameter   :: MAX_VAR_DIMS = 32
    character(len=64)    :: xname,attr_name

    integer(kind=4)      :: rank
    integer              :: sds_id
    integer              :: istat, attributes, num_type, n_values, data_type
    integer              :: sffinfo, sfselect, sfginfo, sfrattr
    integer              :: sfendacc, sfend, sfrnatt, sfrcatt
    integer              :: sffattr, sfrdata, sfn2index
    integer              :: sfgainfo,attr_index,lname
    integer              :: num_ds,num_at,ids
    integer(kind=4),dimension(MAX_VAR_DIMS):: dim_sizes
    integer              :: idebug = 0

    ifail = 1

    lname = len_trim(name)
    istat = sffinfo(sd_id,num_ds,num_at)

    if ( present(index) ) then
       sds_id = sfselect(sd_id, index)   
       istat = sfginfo(sds_id, xname, rank, dim_sizes, num_type, attributes)
       if ( idebug >= 100 ) then
          write(*,*) 'io_read2d_32g: name = ', name,xname(1:lname)
          write(*,*) 'io_read2d_32g: rank = ', rank
          write(*,*) 'io_read2d_32g: dims = ', dim_sizes(1:rank)
          write(*,*) 'io_read2d_32g: attr = ', attributes
       end if
       ! check rank and name....
       names: if( +xname(1:lname) == +name(1:lname) .and.   &
            rank == 2        .and. &
            dim_sizes(1) == im       .and. &
            dim_sizes(2) == jm      ) then

          if ( present(idate) ) then
             attr_index = sffattr(sds_id, 'idate')
             istat =  &
                  sfgainfo(sds_id, attr_index, attr_name, data_type,  n_values)
             istat = sfrattr(sds_id, attr_index, idate_file)
             if ( idebug >= 100 ) then
                write(*,*) 'io_read2d_32g: idate = ', idate,idate_file
             end if
             if ( sum(abs(idate-idate_file)) == 0 ) then
                call read_it   !everything OK....proceed
                return         !and return
             end if
          else
             call read_it   !everything OK....proceed
             return         !and return
          end if
       end if names
       istat = sfendacc(sds_id)   !close 'wrong' ds
    end if

    dsloop: do ids=0,num_ds-1
       sds_id = sfselect(sd_id, ids)   
       istat = sfginfo(sds_id, xname, rank, dim_sizes, num_type, attributes)
       if ( idebug >= 100 ) then
          write(*,*) 'io_read2d_32g: name = ', name,xname(1:lname)
          write(*,*) 'io_read2d_32g: rank = ', rank
          write(*,*) 'io_read2d_32g: dims = ', dim_sizes(1:rank)
          write(*,*) 'io_read2d_32g: attr = ', attributes
       end if
       ! check rank and name....
       if ( +xname(1:lname) /= +name(1:lname) .or.   &
            rank /= 2        .or. &
            dim_sizes(1) /= im       .or. &
            dim_sizes(2) /= jm      ) then

          istat = sfendacc(sds_id)   !close 'wrong' ds
          cycle dsloop
       end if

       if ( present(idate) ) then
          attr_index = sffattr(sds_id, 'idate')
          istat = sfgainfo(sds_id, attr_index, attr_name, data_type,  n_values)
          istat = sfrattr(sds_id, attr_index, idate_file)
          if ( idebug >= 100 ) then
             write(*,*) 'io_read2d_32g: idate = ', idate,idate_file
          end if
          if( sum(abs(idate-idate_file)) /= 0 ) then 
             istat = sfendacc(sds_id)   !close 'wrong' ds
             cycle dsloop
          end if
       end if
       call read_it     !everything OK
       ! set index to next position if present...
       if ( present(index) ) index=ids+1
       return
    end do dsloop

    print *, 'io_read2d_32g: Could not find requested data set in hdf file:'
    print *, 'io_read2d_32g:   name   : ', trim(name)
    print *, 'io_read2d_32g:   shape  : ', im, jm
    if (present(idate)) print *, 'io_read2d_32g:   idate  : ', idate
    if (present(index)) print *, 'io_read2d_32g:   index  : ', index

  contains 

    subroutine read_it
      integer :: sfrdata
      real(kind=4),dimension(:,:),allocatable :: hdfr
      allocate(hdfr(im,jm))
      istat = sfrdata(sds_id, (/0,0/),(/1,1/),(/im,jm/),hdfr)
      if ( istat == SUCCEED ) then
         if ( idebug >= 100 ) &
              print*, 'io_read2d_32g: Successfully retrieved '//name//' from file'
      else
         print*, 'io_read2d_32g: Failed to read '//name//' from file'
         return
      end if
      data = hdfr
      istat = sfendacc(sds_id)
      deallocate(hdfr)
      ifail = 0
    end subroutine read_it

  end subroutine io_read2d_32g



  subroutine io_read3d_32g(sd_id,im,jm,lm,data,name,ifail,index,idate)
    implicit none

    ! in/out:
    integer, intent(in)            :: im,jm,lm
    integer, intent(in)            :: sd_id
    character(len=*),intent(in)    :: name
    real,dimension(im,jm,lm),intent(out) :: data
    integer,intent(out)            :: ifail
    integer,optional               :: index
    integer,dimension(6), optional :: idate

    ! local
    integer,dimension(6):: idate_file

    integer, parameter  :: MAX_VAR_DIMS = 32
    character(len=64)   :: xname,attr_name

    integer(kind=4)     :: rank
    integer             :: sds_id
    integer             :: istat, attributes, num_type, n_values, data_type
    integer             :: sffinfo, sfselect, sfginfo, sfrattr
    integer             :: sfendacc, sfend, sfrnatt, sfrcatt
    integer             :: sffattr, sfrdata, sfn2index
    integer             :: sfgainfo,attr_index,lname
    integer             :: num_ds,num_at,ids
    integer(kind=4),dimension(MAX_VAR_DIMS):: dim_sizes
    integer             :: idebug = 0

    ifail = 1

    lname = len_trim(name)
    istat = sffinfo(sd_id,num_ds,num_at)

    if ( present(index) ) then
       sds_id = sfselect(sd_id, index)   
       istat = sfginfo(sds_id, xname, rank, dim_sizes, num_type, attributes)
       if ( idebug >= 100 ) then
          write(*,*) 'io_read3d_32g: name = ', name,xname(1:lname)
          write(*,*) 'io_read3d_32g: rank = ', rank
          write(*,*) 'io_read3d_32g: dims = ', dim_sizes(1:rank)
          write(*,*) 'io_read3d_32g: attr = ', attributes
       end if
       !check rank and name....
       names: if ( +xname(1:lname) == +name(1:lname) .and.   &
            rank == 3        .and. &
            dim_sizes(1) == im       .and. &
            dim_sizes(2) == jm       .and. &
            dim_sizes(3) == lm      ) then

          if ( present(idate) ) then
             attr_index = sffattr(sds_id, 'idate')
             istat = sfgainfo(sds_id, attr_index, attr_name, data_type, n_values)
             istat = sfrattr(sds_id, attr_index, idate_file)
             if ( idebug >= 100 ) then
                write(*,*) 'io_read3d_32g: idate = ', idate,idate_file
             end if
             if ( sum(abs(idate-idate_file)) == 0 ) then
                call read_it   !everything OK....proceed
                return         !and return
             end if
          else
             call read_it   !everything OK....proceed
             return         !and return
          end if
       end if names
       istat = sfendacc(sds_id)   !close 'wrong' ds
    end if

    dsloop: do ids=0,num_ds-1
       sds_id = sfselect(sd_id, ids)   
       istat = sfginfo(sds_id, xname, rank, dim_sizes, num_type, attributes)
       if ( idebug >= 100 ) then
          write(*,*) 'io_read3d_32g: name = ', name,xname(1:lname)
          write(*,*) 'io_read3d_32g: rank = ', rank
          write(*,*) 'io_read3d_32g: dims = ', dim_sizes(1:rank)
          write(*,*) 'io_read3d_32g: attr = ', attributes
       end if
       !check rank and name....
       if( +xname(1:lname) /= +name(1:lname) .or.   &
            rank /= 3        .or. &
            dim_sizes(1) /= im       .or. &
            dim_sizes(2) /= jm       .or. &
            dim_sizes(3) /= lm      ) then

          istat = sfendacc(sds_id)   !close 'wrong' ds
          cycle dsloop
       end if

       if(present (idate)) then
          attr_index = sffattr(sds_id, 'idate')
          istat = sfgainfo(sds_id, attr_index, attr_name, data_type,  n_values)
          istat = sfrattr(sds_id, attr_index, idate_file)
          if ( idebug >= 100 ) then
             write(*,*) 'io_read3d_32g: idate = ', idate,idate_file
          end if
          if(sum(abs(idate-idate_file)) /=  0) then 
             istat = sfendacc(sds_id)   !close 'wrong' ds
             cycle dsloop
          end if
       end if
       call read_it     !everything OK
       if(present(index)) index=ids+1   !set index to next position if present
       return
    end do dsloop

    print *, 'io_read2d_32g: Could not find requested data set in hdf file:'
    print *, 'io_read3d_32g:   name   : ', trim(name)
    print *, 'io_read3d_32g:   shape  : ', im, jm, lm
    if (present(idate)) print *, 'io_read3d_32g:   idate  : ', idate
    if (present(index)) print *, 'io_read3d_32g:   index  : ', index

  contains 

    subroutine read_it
      integer :: sfrdata
      real(kind=4),dimension(:,:,:),allocatable :: hdfr
      allocate(hdfr(im,jm,lm))
      istat =    sfrdata(sds_id, (/0,0,0/),(/1,1,1/),(/im,jm,lm/),hdfr)
      if (istat == SUCCEED) then
         if ( idebug >= 100 ) &
              print*, 'io_read3d_32g: Successfully retrieved '// &
              name//' from file'
      else
         print*, 'io_read3d_32g: Failed to read '//name//' from file'
         return
      end if
      data = hdfr
      istat = sfendacc(sds_id)
      deallocate(hdfr)
      ifail = 0
    end subroutine read_it

  end subroutine io_read3d_32g



  subroutine io_read4d_32g(sd_id,im,jm,lm,nt,data,name,ifail,index,idate)
    implicit none
    ! in/out:
    integer, intent(in)                     :: im,jm,lm,nt
    integer, intent(in)                     :: sd_id
    character(len=*),intent(in)             :: name
    real,dimension(im,jm,lm,nt),intent(out) :: data
    integer,intent(out)                     :: ifail
    integer,dimension(6), optional          :: idate
    ! local
    integer,dimension(6):: idate_file

    integer, parameter  :: MAX_VAR_DIMS = 32
    character(len=64)   :: xname,attr_name

    integer,optional    :: index
    integer(kind=4)     :: rank
    integer             :: sds_id
    integer             :: istat, attributes, num_type, n_values, data_type
    integer             :: sfselect, sfginfo, sfrattr
    integer             :: sfendacc, sfend, sfrnatt, sfrcatt
    integer             :: sffattr, sfrdata, sfn2index
    integer             :: sfgainfo,attr_index,lname
    integer             :: sffinfo,num_ds,num_at,ids
    integer(kind=4),dimension(MAX_VAR_DIMS):: dim_sizes
    integer             :: idebug = 0

    ifail = 1

    lname = len_trim(name)
    istat = sffinfo(sd_id,num_ds,num_at)

    if ( idebug >= 100 ) then
       print *, 'io_read4d_32g: # ds & att',sd_id,num_ds,num_at
    end if

    if ( present(index) ) then
       sds_id = sfselect(sd_id, index)   
       istat = sfginfo(sds_id, xname, rank, dim_sizes, num_type, attributes)
       if ( idebug >= 100 ) then
          write(*,*) 'io_read4d_32g: name = ', name,xname(1:lname)
          write(*,*) 'io_read4d_32g: rank = ', rank
          write(*,*) 'io_read4d_32g: dims = ', dim_sizes(1:rank)
          write(*,*) 'io_read4d_32g: attr = ', attributes
       end if
       !check rank and name....
       names: if( +xname(1:lname) == +name(1:lname) .and.   &
            rank == 4        .and. &
            dim_sizes(1) == im       .and. &
            dim_sizes(2) == jm       .and. &
            dim_sizes(4) == nt       .and. &
            dim_sizes(3) == lm      ) then

          if ( present(idate) ) then
             attr_index = sffattr(sds_id, 'idate')
             istat = sfgainfo(sds_id, attr_index, attr_name, data_type, n_values)
             istat = sfrattr(sds_id, attr_index, idate_file)
             if (idebug >= 100 ) then
                write(*,*) 'io_read4d_32g: idate = ', idate,idate_file
             end if
             if( sum(abs(idate-idate_file)) == 0 ) then
                call read_it
                return         !and return
             end if
          else
             call read_it
             return         !and return
          end if
       end if names
       istat = sfendacc(sds_id)   !close 'wrong' ds
    end if

    dsloop: do ids=0,num_ds-1
       sds_id = sfselect(sd_id, ids)   
       istat = sfginfo(sds_id, xname, rank, dim_sizes, num_type, attributes)
       if (idebug >= 100 ) then
          write(*,*) 'io_read4d_32g: name = ', name,xname(1:lname)
          write(*,*) 'io_read4d_32g: rank = ', rank
          write(*,*) 'io_read4d_32g: dims = ', dim_sizes(1:rank)
          write(*,*) 'io_read4d_32g: attr = ', attributes
       end if
       !check rank and name....
       if ( +xname(1:lname) /= +name(1:lname) .or.   &
            rank /= 4        .or. &
            dim_sizes(1) /= im       .or. &
            dim_sizes(2) /= jm       .or. &
            dim_sizes(4) /= nt       .or. &
            dim_sizes(3) /= lm      ) then

          istat = sfendacc(sds_id)   !close 'wrong' ds
          cycle dsloop
       end if

       if ( present(idate) ) then
          attr_index = sffattr(sds_id, 'idate')
          istat = sfgainfo(sds_id, attr_index, attr_name, data_type,  n_values)
          istat = sfrattr(sds_id, attr_index, idate_file)
          if ( idebug >= 100 ) then
             write(*,*) 'io_read4d_32g: idate = ', idate,idate_file
          end if
          if ( sum(abs(idate-idate_file)) /= 0 ) then 
             istat = sfendacc(sds_id)   !close 'wrong' ds
             cycle dsloop
          end if
       end if
       call read_it
       return

    end do dsloop

    print*, 'io_read4d_32g: Could not find '//name//' from file'
    print*, 'io_read4d_32g: With dimensions',im,jm,lm,nt
    if ( present(idate) ) print *,'io_read4d_32g: With date..',idate

  contains 

    subroutine read_it
      integer :: sfrdata
      real(kind=4),dimension(:,:,:,:),allocatable :: hdfr
      allocate(hdfr(im,jm,lm,nt))
      istat = sfrdata(sds_id, (/0,0,0,0/),(/1,1,1,1/),(/im,jm,lm,nt/),hdfr)
      if ( istat == SUCCEED ) then
         if ( idebug >= 100 ) &
              print*, 'io_read4d_32g: Successfully retrieved '// &
              name//' from file'
      else
         print*, 'io_read4d_32g: Failed to read '//name//' from file'
         return
      end if
      data = hdfr
      istat = sfendacc(sds_id)
      deallocate(hdfr)
      ifail = 0
    end subroutine read_it

  end subroutine io_read4d_32g



  subroutine io_read4d_32d(sd_id,im,jm,lm,ntrace,data,name,index,ifail,idate)
    implicit none

    ! in/out:
    integer, intent(in)                         :: im,jm,lm,ntrace
    integer, intent(in)                         :: sd_id
    real,dimension(im,jm,lm,ntrace),intent(out) :: data
    character(len=*), intent(in)                :: name
    integer,intent(inout)                       :: index
    integer,intent(out)                         :: ifail
    integer,dimension(6),intent(in)             :: idate

    ! local
    integer,dimension(6) :: idate_file

    integer, parameter   :: MAX_VAR_DIMS = 32
    character(len=64)    :: xname,attr_name

    integer(kind=4)      :: rank
    integer              :: sds_id
    integer              :: istat, attributes, num_type, n_values, data_type
    integer              :: sffinfo, sfselect, sfginfo, sfrattr
    integer              :: sfendacc, sfend, sfrnatt, sfrcatt
    integer              :: sffattr, sfrdata, sfn2index
    integer              :: sfgainfo,attr_index,lname
    integer              :: num_ds,num_at,ids
    integer(kind=4),dimension(MAX_VAR_DIMS) :: dim_sizes

    real(kind=4),dimension(:,:,:,:),allocatable :: hdfr

    ifail = 1

    sds_id = sfselect(sd_id, index)
    istat = sfginfo(sds_id, xname, rank, dim_sizes, num_type, attributes)

    lname = len_trim(name)
    write(*,*) 'io_read4d_32d: name = ', name
    write(*,*) 'io_read4d_32d: rank = ', rank
    write(*,*) 'io_read4d_32d: dims = ', dim_sizes(1:rank)
    write(*,*) 'io_read4d_32d: attr = ', attributes

    attr_index = sffattr(sds_id, 'idate')
    istat = sfgainfo(sds_id, attr_index, attr_name, data_type,  n_values)
    istat = sfrattr(sds_id, attr_index, idate_file)
    !if(  sum(abs(idate_file-idate)) == 0        .and. &
    print*,'io_read4d_32d: ****',im,jm,lm,ntrace,dim_sizes(1:4)
    if(  &
         rank == 4        .and. &
         dim_sizes(1) == im       .and. &
         dim_sizes(2) == jm       .and. &
         dim_sizes(3) == lm       .and. &
         dim_sizes(4) == ntrace     ) then
       allocate(hdfr(im,jm,lm,ntrace))
       istat = sfrdata(sds_id, (/0,0,0,0/),(/1,1,1,1/),(/im,jm,lm,ntrace/),hdfr)
       if ( istat == SUCCEED ) then
          print*, 'io_read4d_32d: Successfully retrieved '//name//' from file'
       else
          print*, 'io_read4d_32d: Failed to read '//name//' from file'
          return
       end if

       data = hdfr
       istat = sfendacc(sds_id)
       deallocate(hdfr)
       ifail = 0
       index = index+1
    end if

  end subroutine io_read4d_32d



  subroutine io_read3d_32dr(sd_id,im,jm,lm,data,name,index,ifail,idate,region)
    implicit none

    ! in/out:
    integer, intent(in)                  :: im,jm,lm
    integer, intent(in)                  :: sd_id
    real,dimension(im,jm,lm),intent(out) :: data
    character(len=*),intent(in)          :: name
    integer,intent(inout)                :: index
    integer,intent(out)                  :: ifail
    integer,dimension(6),intent(in)      :: idate
    integer,intent(in)                   :: region

    ! local
    integer,dimension(6) :: idate_file
    integer              :: region_file
    integer, parameter   :: MAX_VAR_DIMS = 32
    character(len=64)    :: xname,attr_name

    integer              :: rank, sds_id
    integer              :: istat, attributes, num_type, n_values, data_type
    integer              :: sffinfo, sfselect, sfginfo, sfrattr
    integer              :: sfendacc, sfend, sfrnatt, sfrcatt
    integer              :: sffattr, sfrdata, sfn2index
    integer              :: sfgainfo,attr_index
    integer,dimension(MAX_VAR_DIMS) :: dim_sizes

    real(kind=4),dimension(:,:,:),allocatable :: hdfr

    ifail = 1

    sds_id = sfselect(sd_id, index)
    istat = sfginfo(sds_id, xname, rank, dim_sizes, num_type, attributes)

    !write(*,*) 'name = ', name
    !write(*,*) 'rank = ', rank
    !write(*,*) 'dims = ', dim_sizes(1:rank)
    !write(*,*) 'attr = ', attributes

    attr_index = sffattr(sds_id, 'region')
    istat = sfgainfo(sds_id, attr_index , attr_name, data_type,  n_values)
    !print*,'attr_name',attr_name
    !print*,'n_values',n_values
    istat = sfrattr(sds_id, attr_index, region_file)
    attr_index = sffattr(sds_id, 'idate')
    istat = sfgainfo(sds_id, attr_index, attr_name, data_type,  n_values)
    !print*,'attr_name',attr_name
    !print*,'n_values',n_values
    istat = sfrattr(sds_id, attr_index, idate_file)
    !print *,idate_file,idate,region_file,region
    if ( sum(abs(idate_file-idate)) == 0     .and. &
         region_file == region   .and. &
         rank == 3        .and. &
         dim_sizes(1) == im       .and. &
         dim_sizes(2) == jm       .and. &
         dim_sizes(3) == lm      ) then
       allocate(hdfr(im,jm,lm))
       istat = sfrdata(sds_id, (/0,0,0/),(/1,1,1/),(/im,jm,lm/),hdfr)
       if ( istat == SUCCEED ) then
          print*, 'io_read3d_32dr: Successfully retrieved '//name//' from file'
       else
          print*, 'io_read3d_32dr: Failed to read '//name//' from file'
          return
       end if

       data = hdfr
       istat = sfendacc(sds_id)
       deallocate(hdfr)
       ifail = 0
       index = index+1
    end if

  end subroutine io_read3d_32dr



  subroutine io_read2d_32d(sd_id,im,jm,data,name,index,ifail,idate)
    implicit none

    ! in/out:
    integer, intent(in)               :: im,jm
    integer, intent(in)               :: sd_id
    real,dimension(im,jm),intent(out) :: data
    character(len=*), intent(in)      :: name
    integer,intent(inout)             :: index
    integer,intent(out)               :: ifail
    integer,dimension(6),intent(in)   :: idate

    ! local
    integer,dimension(6) :: idate_file

    integer, parameter   :: MAX_VAR_DIMS = 32
    character(len=64)    :: xname,attr_name

    integer              :: rank, sds_id
    integer              :: istat, attributes, num_type, n_values, data_type
    integer              :: sffinfo, sfselect, sfginfo, sfrattr
    integer              :: sfendacc, sfend, sfrnatt, sfrcatt
    integer              :: sffattr, sfrdata, sfn2index
    integer              :: sfgainfo,attr_index

    integer,dimension(MAX_VAR_DIMS) :: dim_sizes

    real(kind=4),dimension(:,:),allocatable :: hdfr

    ifail = 1

    sds_id = sfselect(sd_id, index)
    istat = sfginfo(sds_id, xname, rank, dim_sizes, num_type, attributes)

    write(*,*) 'io_read2d_32d: name = ', name
    write(*,*) 'io_read2d_32d: rank = ', rank
    write(*,*) 'io_read2d_32d: dims = ', dim_sizes(1:rank)
    write(*,*) 'io_read2d_32d: attr = ', attributes

    attr_index = sffattr(sds_id, 'idate')
    istat = sfgainfo(sds_id, attr_index, attr_name, data_type,  n_values)
    istat = sfrattr(sds_id, attr_index, idate_file)
    if( &
         rank == 2        .and. &
         dim_sizes(1) == im       .and. &
         dim_sizes(2) == jm      ) then
       allocate(hdfr(im,jm))
       istat = sfrdata(sds_id, (/0,0/),(/1,1/),(/im,jm/),hdfr)
       if ( istat == SUCCEED ) then
          print*, 'io_read2d_32d: Successfully retrieved '//name//' from file'
       else
          print*, 'io_read2d_32d: Failed to read '//name//' from file'
          return
       end if

       data = hdfr
       istat = sfendacc(sds_id)
       deallocate(hdfr)
       ifail = 0
       index = index+1
    end if

  end subroutine io_read2d_32d



  subroutine io_read2d_32dr(sd_id,im,jm,data,name,index,ifail,idate,region)
    implicit none

    ! in/out:
    integer,intent(in)                :: im,jm
    integer,intent(in)                :: sd_id
    real,dimension(im,jm),intent(out) :: data
    character(len=*),intent(in)       :: name
    integer,intent(inout)             :: index
    integer,intent(out)               :: ifail
    integer,dimension(6),intent(in)   :: idate
    integer,intent(in)                :: region

    ! local
    integer,dimension(6) :: idate_file
    integer              :: region_file

    integer, parameter   :: MAX_VAR_DIMS = 32
    character(len=64)    :: xname,attr_name

    integer              :: rank, sds_id
    integer              :: istat, attributes, num_type, n_values, data_type
    integer              :: sffinfo, sfselect, sfginfo, sfrattr
    integer              :: sfendacc, sfend, sfrnatt, sfrcatt
    integer              :: sffattr, sfrdata, sfn2index
    integer              :: sfgainfo,attr_index

    integer,dimension(MAX_VAR_DIMS) :: dim_sizes

    real(kind=4),dimension(:,:),allocatable :: hdfr

    ifail = 1

    sds_id = sfselect(sd_id, index)
    istat = sfginfo(sds_id, xname, rank, dim_sizes, num_type, attributes)

    !write(*,*) 'name = ', name
    !write(*,*) 'rank = ', rank
    !write(*,*) 'dims = ', dim_sizes(1:rank)
    !write(*,*) 'attr = ', attributes

    attr_index = sffattr(sds_id, 'region')
    istat = sfgainfo(sds_id, attr_index , attr_name, data_type,  n_values)
    !print*,'attr_name',attr_name
    !print*,'n_values',n_values
    istat = sfrattr(sds_id, attr_index, region_file)
    attr_index = sffattr(sds_id, 'idate')
    istat = sfgainfo(sds_id, attr_index, attr_name, data_type,  n_values)
    !print*,'attr_name',attr_name
    !print*,'n_values',n_values
    istat = sfrattr(sds_id, attr_index, idate_file)
    !print *,idate_file,idate,region_file,region
    if(  sum(abs(idate_file-idate)) == 0        .and. &
         region_file == region   .and. &
         rank == 2        .and. &
         dim_sizes(1) == im       .and. &
         dim_sizes(2) == jm      ) then
       allocate(hdfr(im,jm))
       istat =    sfrdata(sds_id, (/0,0/),(/1,1/),(/im,jm/),hdfr)
       if ( istat == SUCCEED ) then
          print*, 'io_read2d_32dr: Successfully retrieved '//name//' from file'
       else
          print*, 'io_read2d_32dr: Failed to read '//name//' from file'
          return
       end if

       data = hdfr
       istat = sfendacc(sds_id)
       deallocate(hdfr)
       ifail = 0
       index = index+1
    end if

  end subroutine io_read2d_32dr

  
  
  function upper_case (old) result (new)
    !
    ! returns the upper-case version of the input string
    ! 
    implicit none
    character(len = *), intent(in)  ::old
    character(len = 64)             ::new
    integer :: i
    new = old
    do i=1,len_trim(old)
       if( lge (old(i:i), 'a') .and. &
            lle (old(i:i), 'z')    )   &
            new(i:i) = achar (iachar(old(i:i)) - 32)
    end do
  end function upper_case


end module io_hdf
