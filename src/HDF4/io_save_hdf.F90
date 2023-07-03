!###############################################################################
!
! contains routines to read emissions and
! to read and write the main model state from/to file
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

module io_save_hdf

  use GO, only : gol, goPr, goErr

  implicit none

  private

  public :: savehdf
  public :: readhdf
  public :: readhdfmmr
  public :: readhdfmmix
  public :: readtm3hdf

  character(len=*), parameter :: mname = 'io_save_hdf'

contains

  subroutine readtm3hdf( filename, rank, im_emis, jm_emis, lm_emis, &
                         month, field, field_name, idt )
    use dims
    use toolbox, only : escape_tm
    use io_hdf,  only : DFACC_READ, FAIL
    use io_hdf,  only : io_read2d_32g, io_read3d_32g
    implicit none

    !WP! this subroutine reads HDF emissions files on
    !    1x1 degree resolution and returns the specified field.
    !WP! input for this subroutine are:
    !WP!
    !WP! filename --> name of hdf file to be opened or read
    !MK! rank    --> rank of the dataset (1,2,3)
    !WP! im_emis --> x-DIMENSION of requested field
    !WP! jm_emis --> y-DIMENSION of requested field
    !WP! lm_emis --> z-DIMENSION of requested field
    !WP! month  --> the month to read (0 for january or simply the first field)
    !WP! field --> the field (im_emis,jm_emis,lm_emis) to put the values in
    !WP! field_name --> the name of the field to match the one in the HDF file
    !WP!
    !WP! The routine checks the dimensions and reports the succes/failure.

    ! in/out
    character(len=*),intent(in)  :: filename
    character(len=*),intent(in)  :: field_name
    integer,intent(in)           :: rank
    integer,intent(in)           :: month
    integer,intent(in)           :: im_emis, jm_emis, lm_emis
    real,dimension(im_emis,jm_emis,lm_emis),intent(out)  :: field
    integer,dimension(6),optional :: idt

    ! local
    real,dimension(:),allocatable       :: field1d   ! first dim will be 1
    real,dimension(:,:),allocatable     :: field2d
    real,dimension(:,:,:),allocatable   :: field3d

    integer :: io, sfstart, ifail, sds_id, ioo
    integer :: n_datasets, n_file_attrs
    integer :: istat, attributes, num_type
    integer :: sffinfo, sfselect, sfginfo
    integer :: sfend, sfrnatt, sfrcatt, sffattr
    integer :: idx,l
    integer :: region,ntot

    ! start

    idx=month    !FD hdf files indices run from 0 ...n-1
    io = SFSTART(filename, DFACC_READ)  !WP! assign io-number to HDF file

    if ( io /= FAIL ) then
       !print *,'readtm3hdf: ',filename,'... opened for READ access.',io
    else
       !WP!print *,'readtm3hdf: File name:',filename
       call escape_tm('readtm3hdf: Unable to open HDF file '//&
            filename//' in READTM3HDF ')
    end if

    select case(rank)
    case(1)    ! latitudional field

       allocate(field1d(jm_emis))
       if(present(idt)) call io_read2d_32g(io,im_emis,jm_emis,field1d,field_name, ifail,idate=idt)
       if(.not.present(idt)) call io_read2d_32g(io,im_emis,jm_emis,field1d,field_name, ifail,index=idx)
       if ( ifail /= 0 ) &
          call escape_tm('readtm3hdf: Error reading HDF file 1D')
       field(1,:,1)=field1d
       deallocate(field1d)

    case(2)    ! 2D surface field

       allocate(field2d(im_emis,jm_emis))
       if(present(idt)) call io_read2d_32g(io,im_emis,jm_emis,field2d,field_name, ifail,idate=idt)
       if(.not.present(idt)) call io_read2d_32g(io,im_emis,jm_emis,field2d,field_name, ifail,index=idx)
       if ( ifail /= 0 ) &
            call escape_tm('readtm3hdf: Error reading HDF file 2D')
       field(:,:,1)=field2d
       deallocate(field2d)

    case (3)   ! 3D field

       allocate(field3d(im_emis,jm_emis,lm_emis))
       if(present(idt)) call io_read3d_32g(io,im_emis,jm_emis,lm_emis,field3d,field_name, ifail,idate=idt)
       if(.not.present(idt)) call io_read3d_32g(io,im_emis,jm_emis,lm_emis,field3d,field_name, ifail,index=idx)
       if ( ifail /= 0 )  &
            call escape_tm('readtm3hdf: Error reading HDF file 3D')
       field(:,:,:)=field3d
       deallocate(field3d)

    case (23)   ! 2D field lat-pres

       allocate(field2d(jm_emis,lm_emis))
       if(present(idt)) call io_read2d_32g(io,jm_emis,lm_emis,field2d,field_name, ifail,idate=idt)
       if(.not.present(idt)) call io_read2d_32g(io,jm_emis,lm_emis,field2d,field_name, ifail,index=idx)
       if ( ifail /= 0 )  then   !might be due to (1,jm,lm) file!
          allocate(field3d(1,jm_emis,lm_emis))
          if(present(idt)) call io_read3d_32g(io,1,jm_emis,lm_emis,field3d,field_name, ifail,idate=idt)
          if(.not.present(idt)) call io_read3d_32g(io,1,jm_emis,lm_emis,field3d,field_name, ifail,index=idx)
          if ( ifail /= 0 ) &
               call escape_tm('readtm3hdf: Error reading HDF file 23')
          field(:,:,:)=field3d
          deallocate(field3d)
       else
          field(1,:,:)=field2d
          deallocate(field2d)
       end if

    case default
       print *,'readtm3hdf: Please check the rank of the required field,', &
               ' must be 1,2,3, or 23'
    end select

    istat = SFEND(io)  !close hdf file

    !! info ...
    !if ( idx == 0 ) print *,'readtm3hdf: ',filename, &
    !     ' total sum after reading the full field',sum(field)/1e9 ,' (x1e9) '
    !if ( idx > 0 ) print *,'readtm3hdf: ',filename,' ds#= ',idx+1, &
    !     ' total month sum after reading the full field', &
    !     sum(field)/1e9 ,' (x1e9) '

  end subroutine readtm3hdf



  subroutine readhdf(region,file_name)
    !
    ! data are expected in kg
    ! (typical save.hdf files from previous run...)
    ! only long-lived tracer are read, since the fields
    ! for chemistry are not yet allocated
    !
    use io_hdf
    use dims,        only : region_name, im, jm, lm, nregions, &
                            datadir, parent, adv_scheme
    use global_data, only : mass_dat
    use MeteoData  , only : m_dat
    use chem_param,  only : ntrace, ntracet
    use toolbox,     only : escape_tm
#ifdef MPI
    use mpi_comm,    only : scatter_after_read_t, stopmpi
    use mpi_const,   only : root_t, myid, ntracet_ar, ntracetloc
    use mpi_const,   only : com_trac, ierr, my_real, mpi_integer
#endif

    implicit none

    ! in/out
    integer, intent(in)             :: region
    character(len=*), intent(in)    :: file_name

    ! local
    real,dimension(:,:,:,:),pointer :: rm,rxm,rym,rzm
#ifdef secmom
    real,dimension(:,:,:,:),pointer :: rxxm,rxym,rxzm,ryym,ryzm,rzzm
#endif
    real,dimension(:,:,:),pointer   :: m

    integer                         :: io, sfstart, ifail, sds_id
    integer                         :: n_datasets, n_file_attrs
    integer                         :: istat, attributes, num_type
    integer                         :: sffinfo, sfselect, sfginfo
    integer                         :: sfend, sfrnatt, sfrcatt, sffattr
    integer,dimension(6)            :: idate_save
    integer,dimension(nregions)     :: im_file,jm_file,lm_file
    integer                         :: ntrace_file, i,j,l,n
    character(len=80)               :: msg_file
    integer                         :: ind, my_parent,imr,jmr,lmr
    real,dimension(:,:,:,:), allocatable :: field, field2
    real,dimension(:,:,:), allocatable   :: msave
    integer                         :: from_file  ! signals reading from file

    ! start

    m => m_dat(region)%data
    rm => mass_dat(region)%rm_t
    rxm => mass_dat(region)%rxm_t
    rym => mass_dat(region)%rym_t
    rzm => mass_dat(region)%rzm_t
#ifdef secmom
    rxxm => mass_dat(region)%rxxm_t
    rxym => mass_dat(region)%rxym_t
    rxzm => mass_dat(region)%rxzm_t
    ryym => mass_dat(region)%ryym_t
    ryzm => mass_dat(region)%ryzm_t
    rzzm => mass_dat(region)%rzzm_t
#endif

    imr = im(region) ; jmr = jm(region) ; lmr = lm(region)

    allocate(field(-1:imr+2,-1:jmr+2,lmr,ntracet))
    field = 0.0

    allocate(msave(imr,jmr,lmr))

#ifdef MPI
    if ( myid == root_t ) then
#endif
       from_file = 1
       allocate(field2(imr,jmr,lmr,ntracet))
       ind = len_trim(datadir)
       !io = sfstart(datadir(1:ind)//file_name//region_name(region), DFACC_READ)
       io = sfstart( file_name, DFACC_READ )
       if ( io == FAIL ) then
          my_parent = parent(region)
          if ( my_parent == 0 ) then
             !call escape_tm('readhdf: Could not open file and no parent '// &
             !     'available :'//datadir(1:ind)//file_name//region_name(region))
             call escape_tm('readhdf: Could not open file and no parent '// &
                  'available :'//file_name )
          else
             from_file = 0
             print *, 'readhdf: Trying to initialise '//region_name(region)// &
                  ' from parent....'
          end if
       end if
       if ( from_file == 1 ) then
          print*,' '
          print*,'readhdf: ',file_name,'... opened for READ access.'

          call io_read3d_32g(io,imr,jmr,lmr,msave,'m',ifail)
          if ( ifail /= 0 ) then
             call escape_tm('readhdf: Failed to read m from saveold.hdf in readhdf')
          end if

          call io_read4d_32g(io,imr,jmr,lmr,ntracet,field2,'rm',ifail)
          if( ifail /= 0 .and. ntracet == 1 ) then
             call io_read3d_32g(io,imr,jmr,lmr,field2,'rm',ifail)
          end if

          if ( ifail /= 0 ) &
               call escape_tm('readhdf: Failed to read in saveold.hdf in readhdf')
          print*, 'readhdf: Read rm, ifail = ', ifail

          field(1:imr,1:jmr,1:lmr,:)=field2 ! only transported tracers!
       end if  !from_file
#ifdef MPI
    end if  !root_t
#endif

#ifdef MPI
    ! broadcast from_file !
    call mpi_bcast(from_file, 1, mpi_integer , root_t, com_trac, ierr)
#endif

    if ( from_file == 1 ) then
#ifdef MPI
       call mpi_bcast(msave, imr*jmr*lmr, my_real, root_t, com_trac, ierr)
       call scatter_after_read_t(field, imr,jmr,lmr,2,2,0,rm,root_t)
#else
       rm = field
#endif
    else
       call init_child(region)   ! get info from parent....
    end if

    field = 0.0

    if ( from_file == 1 ) then
       if ( adv_scheme == 'slope'.or.adv_scheme.eq.'2nd_m' ) then

#ifdef MPI
          if ( myid == root_t ) then
#endif
             call io_read4d_32g(io,imr,jmr,lmr,ntracet,field2,'rxm',ifail)
             if ( ifail /= 0 .and. ntracet == 1 ) then
                call io_read3d_32g(io,imr,jmr,lmr,field2,'rxm',ifail)
             end if
             if ( ifail /= 0 ) then
                call escape_tm( &
                  'readhdf: Failed to read in rxm in saveold.hdf in readhdf')
             end if
             print*, 'readhdf: Read rxm, ifail = ', ifail

             field(1:imr,1:jmr,1:lmr,:)=field2 ! only transported tracers!
#ifdef MPI
          end if  !root_t
#endif

#ifdef MPI
          call scatter_after_read_t(field, imr,jmr,lmr,2,2,0,rxm,root_t)
#else
          rxm = field
#endif
          field = 0.0


#ifdef MPI
          if(myid == root_t) then
#endif
             call io_read4d_32g(io,imr,jmr,lmr,ntracet,field2,'rym',ifail)
             if ( ifail /= 0 .and. ntracet == 1 ) then
                call io_read3d_32g(io,imr,jmr,lmr,field2,'rym',ifail)
             end if
             if ( ifail /= 0 ) then
                call escape_tm( &
                     'readhdf: Failed to read in rym in saveold.hdf in readhdf')
             end if
             print*, 'readhdf: Read rym, ifail = ', ifail
             field(1:imr,1:jmr,1:lmr,:)=field2 ! only transported tracers!
#ifdef MPI
          end if  !root_t
#endif

#ifdef MPI
          call scatter_after_read_t(field, imr,jmr,lmr,2,2,0,rym,root_t)
#else
          rym = field
#endif
          field = 0.0

#ifdef MPI
          if ( myid == root_t ) then
#endif
             call io_read4d_32g(io,imr,jmr,lmr,ntracet,field2,'rzm',ifail)
             if ( ifail /= 0 .and. ntracet == 1 ) then
                call io_read3d_32g(io,imr,jmr,lmr,field2,'rzm',ifail)
             end if
             if ( ifail /= 0 ) then
                call escape_tm( &
                     'readhdf: Failed to read in rzm in saveold.hdf in readhdf')
             end if

             field(1:imr,1:jmr,1:lmr,:)=field2 ! only transported tracers!
             deallocate(field2)
#ifdef MPI
          end if  !root_t
#endif

#ifdef MPI
          call scatter_after_read_t(field, imr,jmr,lmr,2,2,0,rzm,root_t)
#else
          rzm = field
#endif

       end if  !slope


       ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       ! >>> second moments
       ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#ifdef secmom
       if ( adv_scheme == '2nd_m' ) then
         allocate(field2(imr,jmr,lmr,ntracet))

#ifdef MPI
   if(myid == root_t) then
#endif
      call io_read4d_32g(io,imr,jmr,lmr,ntracet,field2,'rxxm',ifail)
             if ( ifail /= 0 .and. ntracet == 1 ) then
         call io_read3d_32g(io,imr,jmr,lmr,field2,'rxxm',ifail)
      endif
             if ( ifail /= 0 ) then
                call escape_tm( &
                     'readhdf: Failed to read in rxxm in saveold.hdf in readhdf')
             end if
      print*, ' Read rxxm, ifail = ', ifail

      field(1:imr,1:jmr,1:lmr,:)=field2 ! only transported tracers!
#ifdef MPI
   endif  !root_t
#endif

#ifdef MPI
   call scatter_after_read_t(field, imr,jmr,lmr,2,2,0,rxxm,root_t)
#else
   rxxm = field
#endif
   field = 0.0


#ifdef MPI
   if(myid == root_t) then
#endif
      call io_read4d_32g(io,imr,jmr,lmr,ntracet,field2,'rxym',ifail)
             if ( ifail /= 0 .and. ntracet == 1 ) then
         call io_read3d_32g(io,imr,jmr,lmr,field2,'rxym',ifail)
      endif
             if ( ifail /= 0 ) then
                call escape_tm( &
                     'readhdf: Failed to read in rxym in saveold.hdf in readhdf')
             end if
      print*, ' Read rxym, ifail = ', ifail
      field(1:imr,1:jmr,1:lmr,:)=field2 ! only transported tracers!
#ifdef MPI
   endif  !root_t
#endif

#ifdef MPI
   call scatter_after_read_t(field, imr,jmr,lmr,2,2,0,rxym,root_t)
#else
   rxym = field
#endif
   field = 0.0

#ifdef MPI
   if(myid == root_t) then
#endif
      call io_read4d_32g(io,imr,jmr,lmr,ntracet,field2,'rxzm',ifail)
             if ( ifail /= 0 .and. ntracet == 1 ) then
         call io_read3d_32g(io,imr,jmr,lmr,field2,'rxzm',ifail)
      endif
             if ( ifail /= 0 ) then
                call escape_tm( &
                     'readhdf: Failed to read in rxzm in saveold.hdf in readhdf')
             end if
             print*, ' Read rxzm, ifail = ', ifail

      field(1:imr,1:jmr,1:lmr,:)=field2 ! only transported tracers!
#ifdef MPI
   endif  !root_t
#endif

#ifdef MPI
   call scatter_after_read_t(field, imr,jmr,lmr,2,2,0,rxzm,root_t)
#else
   rxzm = field
#endif

#ifdef MPI
   if(myid == root_t) then
#endif
      call io_read4d_32g(io,imr,jmr,lmr,ntracet,field2,'ryym',ifail)
             if ( ifail /= 0 .and. ntracet == 1 ) then
         call io_read3d_32g(io,imr,jmr,lmr,field2,'ryym',ifail)
      endif
             if ( ifail /= 0 ) then
                call escape_tm( &
                     'readhdf: Failed to read in ryym in saveold.hdf in readhdf')
             end if
      print*, ' Read ryym, ifail = ', ifail

      field(1:imr,1:jmr,1:lmr,:)=field2 ! only transported tracers!
#ifdef MPI
   endif  !root_t
#endif

#ifdef MPI
   call scatter_after_read_t(field, imr,jmr,lmr,2,2,0,ryym,root_t)
#else
   ryym = field
#endif
   field = 0.0


#ifdef MPI
   if(myid == root_t) then
#endif
      call io_read4d_32g(io,imr,jmr,lmr,ntracet,field2,'ryzm',ifail)
             if ( ifail /= 0 .and. ntracet == 1 ) then
         call io_read3d_32g(io,imr,jmr,lmr,field2,'ryzm',ifail)
      endif
             if ( ifail /= 0 ) then
                call escape_tm( &
                     'readhdf: Failed to read in ryzm in saveold.hdf in readhdf')
             end if
      print*, ' Read ryzm, ifail = ', ifail
      field(1:imr,1:jmr,1:lmr,:)=field2 ! only transported tracers!
#ifdef MPI
   endif  !root_t
#endif

#ifdef MPI
   call scatter_after_read_t(field, imr,jmr,lmr,2,2,0,ryzm,root_t)
#else
   ryzm = field
#endif
   field = 0.0

#ifdef MPI
   if(myid == root_t) then
#endif
      call io_read4d_32g(io,imr,jmr,lmr,ntracet,field2,'rzzm',ifail)
             if ( ifail /= 0 .and. ntracet == 1 ) then
         call io_read3d_32g(io,imr,jmr,lmr,field2,'rzzm',ifail)
      endif
             if ( ifail /= 0 ) then
                call escape_tm( &
                     'readhdf: Failed to read in rzzm in saveold.hdf in readhdf')
             end if
          print*, ' Read rzzm, ifail = ', ifail

      field(1:imr,1:jmr,1:lmr,:)=field2 ! only transported tracers!
      deallocate(field2)
#ifdef MPI
   endif  !root_t
#endif

#ifdef MPI
   call scatter_after_read_t(field, imr,jmr,lmr,2,2,0,rzzm,root_t)
#else
   rzzm = field
#endif

       end if  !2nd-m

#endif
       ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
       ! <<< end second moments
       ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#ifdef MPI
       if ( myid == root_t ) then
#endif
          istat = sfend(io)
          if ( istat /= FAIL ) then
             print*,'readhdf: ... file closed'
             print*,' '
          else
             call escape_tm('readhdf: ERROR in restart from HDF file')
          end if
#ifdef MPI
       end if
#endif

       ! scale in input by a possible factor.
       ! rm/msave ---> mixing ratio on file. mr*m --> something back in kilo's
#ifdef MPI
       if ( myid == root_t ) then
#endif
          print *, 'readhdf: Maximum ratio of the air masses:', &
               maxval(abs(msave/m(1:imr,1:jmr,1:lmr)))
          print *, 'readhdf: Minimum ratio of the air masses:', &
               minval(abs(msave/m(1:imr,1:jmr,1:lmr)))
#ifdef MPI
       end if
#endif

#ifdef MPI
       do n=1,ntracetloc
#else
       do n=1,ntracet
#endif
          do l = 1,lmr
             do j = 1,jmr
                do i = 1,imr
                   rm(i,j,l,n)  = rm (i,j,l,n)*m(i,j,l)/msave(i,j,l)
                   rxm(i,j,l,n) = rxm(i,j,l,n)*m(i,j,l)/msave(i,j,l)
                   rym(i,j,l,n) = rym(i,j,l,n)*m(i,j,l)/msave(i,j,l)
                   rzm(i,j,l,n) = rzm(i,j,l,n)*m(i,j,l)/msave(i,j,l)
#ifdef secmom
                   rxxm(i,j,l,n) = rxxm(i,j,l,n)*m(i,j,l)/msave(i,j,l)
                   rxym(i,j,l,n) = rxym(i,j,l,n)*m(i,j,l)/msave(i,j,l)
                   rxzm(i,j,l,n) = rxzm(i,j,l,n)*m(i,j,l)/msave(i,j,l)
                   ryym(i,j,l,n) = ryym(i,j,l,n)*m(i,j,l)/msave(i,j,l)
                   ryzm(i,j,l,n) = ryzm(i,j,l,n)*m(i,j,l)/msave(i,j,l)
                   rzzm(i,j,l,n) = rzzm(i,j,l,n)*m(i,j,l)/msave(i,j,l)
#endif
                end do
             end do
          end do
       end do

    end if   ! from_file == 1

    deallocate(msave)
    deallocate(field)
    nullify(m)
    nullify(rm)
    nullify(rxm)
    nullify(rym)
    nullify(rzm)
#ifdef secmom
    nullify(rxxm)
    nullify(rxym)
    nullify(rxzm)
    nullify(ryym)
    nullify(ryzm)
    nullify(rzzm)
#endif

  end subroutine readhdf



  subroutine init_child(region)
    !
    ! Parent is used to initilise rm, rxm,rym,rzm of region.
    ! Assumed parallel over tracers...
    !
    use dims,        only : lm, ibeg, iend, jbeg, jend, xref, yref, parent
    use global_data, only : mass_dat
    use MeteoData  , only : m_dat
    use chem_param,  only : ntrace, ntracet
    use toolbox,     only : escape_tm
    use ParTools,    only : myid, ntracetloc, root_t

    implicit none

    ! in/out
    integer, intent(in) :: region

    ! local
    real,dimension(:,:,:,:),pointer :: rm,rxm,rym,rzm   ! child
    real,dimension(:,:,:),pointer   :: m
    real,dimension(:,:,:,:),pointer :: rmp   !parent
    real,dimension(:,:,:),pointer   :: mp

    real           :: mm_parent
    integer        :: my_parent
    integer        :: n,l,jp,ip,ic,jc, xref_, yref_, i,j

    ! start

    m => m_dat(region)%data
    rm => mass_dat(region)%rm_t
    rxm => mass_dat(region)%rxm_t
    rym => mass_dat(region)%rym_t
    rzm => mass_dat(region)%rzm_t

    my_parent = parent(region)

    if ( my_parent == 0 ) then
       call escape_tm( 'init_child: Region 1 cannot be initialised')
    end if

    mp => m_dat(my_parent)%data
    rmp => mass_dat(my_parent)%rm_t

    xref_ = xref(region)/xref(my_parent)
    yref_ = yref(region)/yref(my_parent)
    do n=1,ntracetloc
       do l=1,lm(region)
          do jp=jbeg(region), jend(region)
             jc = (jp-jbeg(region))*yref_
             do ip=ibeg(region), iend(region)
                mm_parent = rmp(ip,jp,l,n)/mp(ip,jp,l)
                ic = (ip-ibeg(region))*xref_
                do j=1,yref_
                   do i=1,xref_
                      rm (ic+i,jc+j,l,n) = mm_parent*m(ic+i,jc+j,l)
                      rxm(ic+i,jc+j,l,n) = 0.0
                      rym(ic+i,jc+j,l,n) = 0.0
                      rzm(ic+i,jc+j,l,n) = 0.0
                   end do
                end do
             end do
          end do
       end do
    end do
    if ( myid == root_t ) print *, 'readhdf: Initialisation successful '

    nullify(m)
    nullify(rm)
    nullify(rxm)
    nullify(rym)
    nullify(rzm)
    nullify(mp)
    nullify(rmp)

  end subroutine init_child

  subroutine init_child_edge(region)
    !
    ! ONLY the edges!
    ! Parent is used to initilise rm, rxm,rym,rzm of region.
    ! Assumed parallel over tracers...
    !
    use dims,        only : lm, ibeg, iend, jbeg, jend, xref, yref, parent, im, jm
    use global_data, only : mass_dat, region_dat
    use MeteoData  , only : m_dat
    use chem_param,  only : ntrace, ntracet
    use toolbox,     only : escape_tm
    use ParTools,    only : myid, ntracetloc, root_t

    implicit none

    ! in/out
    integer, intent(in) :: region

    ! local
    real,dimension(:,:,:,:),pointer :: rm,rxm,rym,rzm   ! child
    real,dimension(:,:,:),pointer   :: m
    real,dimension(:,:,:,:),pointer :: rmp   !parent
    real,dimension(:,:,:),pointer   :: mp
    integer,dimension(:,:),pointer  :: edge

    real           :: mm_parent
    integer        :: my_parent
    integer        :: n,l,jp,ip,ic,jc, xref_, yref_, i,j

    ! start

    my_parent = parent(region)
    if ( my_parent == 0 ) return

    m => m_dat(region)%data
    rm => mass_dat(region)%rm_t
    rxm => mass_dat(region)%rxm_t
    rym => mass_dat(region)%rym_t
    rzm => mass_dat(region)%rzm_t
    edge => region_dat(region)%edge

    mp => m_dat(my_parent)%data
    rmp => mass_dat(my_parent)%rm_t

    xref_ = xref(region)/xref(my_parent)
    yref_ = yref(region)/yref(my_parent)
    do n=1,ntracetloc
       do l=1,lm(region)
          do j=1,jm(region)
             jp = jbeg(region) + (j-1)/yref_
             do i=1,im(region)
                if (edge(i,j) /= -1 ) cycle
                ip = ibeg(region) + (i-1)/xref_
                mm_parent = rmp(ip,jp,l,n)/mp(ip,jp,l)
                rm (i,j,l,n) = mm_parent*m(i,j,l)
                rxm(i,j,l,n) = 0.0
                rym(i,j,l,n) = 0.0
                rzm(i,j,l,n) = 0.0
             end do
          end do
       end do
    end do
    if ( myid == root_t ) print *, 'readhdf: Edges Initialisation successful '

    nullify(m)
    nullify(rm)
    nullify(rxm)
    nullify(rym)
    nullify(rzm)
    nullify(mp)
    nullify(rmp)
    nullify(edge)

  end subroutine init_child_edge


  subroutine readhdfmmr(region,file_name)
    !
    ! data are expected in mixing ratio...
    !
    use dims,        only : im, jm, lm, nregions, datadir, region_name
    use io_hdf
    use global_data, only : mass_dat
    use MeteoData  , only : m_dat
    use chem_param,  only : fscale, ntracet
    use toolbox,     only : escape_tm
#ifdef MPI
    use mpi_comm,    only : scatter_after_read_t, stopmpi
    use mpi_const,   only : root, myid, ntracet_ar, ntracetloc
#endif

    implicit none

    ! in/out
    integer, intent(in)             :: region
    character(len=*), intent(in)    :: file_name

    ! local
    real,dimension(:,:,:,:),pointer :: rm,rxm,rym,rzm
#ifdef secmom
    real,dimension(:,:,:,:),pointer :: rxxm,rxym,rxzm,ryym,ryzm,rzzm
#endif
    real,dimension(:,:,:),pointer   :: m
    integer                         :: io, sfstart, ifail, sds_id
    integer                         :: n_datasets, n_file_attrs
    integer                         :: istat, attributes, num_type
    integer                         :: sffinfo, sfselect, sfginfo
    integer                         :: sfend, sfrnatt, sfrcatt, sffattr
    integer,dimension(6)            :: idate_save
    integer,dimension(nregions)     :: im_file,jm_file,lm_file
    integer                         :: ntrace_file
    character(len=80)               :: msg_file
    integer                         :: ind,ind1,n,i,j,l,offsetn
    real,dimension(:,:,:,:), allocatable :: field,field2

    ! start

    m => m_dat(region)%data
    rm => mass_dat(region)%rm_t
    rxm => mass_dat(region)%rxm_t
    rym => mass_dat(region)%rym_t
    rzm => mass_dat(region)%rzm_t
#ifdef secmom
    rxxm => mass_dat(region)%rxxm_t
    rxym => mass_dat(region)%rxym_t
    rxzm => mass_dat(region)%rxzm_t
    ryym => mass_dat(region)%ryym_t
    ryzm => mass_dat(region)%ryzm_t
    rzzm => mass_dat(region)%rzzm_t
#endif
    allocate(field(-1:im(region)+2,-1:jm(region)+2,1:lm(region),1:ntracet))
    field=0.0
#ifdef MPI
    if ( myid == root ) then
#endif
       ! to hold all ntrace values of rm
       allocate(field2(im(region),jm(region),lm(region),ntracet))
       io = sfstart( file_name, DFACC_READ )

       if ( io == -1 ) call escape_tm('readhdfmmr: Could not open file:'//file_name)
       ind = 1  ! select 1st dataset....(0 = m)

       istat = sffinfo(io, n_datasets, n_file_attrs)
       print*, 'readhdfmmr: Number of datasets and attributes',  &
            n_datasets, n_file_attrs
       call io_read4d_32g(io,im(region),jm(region),lm(region),ntracet, &
            field2,'rm',ifail,index=ind)
       if ( ifail /= 0 .and. ntracet == 1 ) then
          call io_read3d_32g(io,im(region),jm(region),lm(region), &
               field2,'rm',ifail,index=ind)
       end if
       if ( ifail /= 0 ) then
!PB          print *, 'readhdfmmr: '//datadir(1:ind)//file_name(1:ind1-1)// &
!PB            '_'//region_name(region)//'.hdf'
!PBi
          print *, 'readhdfmmr: '//trim(datadir)//trim(file_name)// &
            '_'//region_name(region)//'.hdf'
!PBf
          call escape_tm('readhdfmmr: Failed to read fields')
       end if
       print*, 'readhdfmmr: Read rm, ifail = ', ifail

       istat = sfend(io)
       if (istat /= FAIL) then
          print*,'readhdfmmr: ... file closed'
          print*,' '
       else
          call escape_tm('readhdfmmr: ERROR in restart from HDF file')
       end if

       ! only transported tracers!
       field(1:im(region),1:jm(region),1:lm(region),:)=field2
       deallocate(field2) ! no longer needed
#ifdef MPI
    end if ! root
    call scatter_after_read_t(field,im(region),jm(region),lm(region),&
         2,2,0,rm,root)
#else
    rm = field
#endif
#ifdef MPI
    offsetn = sum(ntracet_ar(0:myid-1))
    do n=1,ntracetloc
#else
    offsetn = 0
    do n=1,ntracet
#endif
       do l = 1,lm(region)
          do j = 1,jm(region)
             do i = 1,im(region)
                rm(i,j,l,n) = rm(i,j,l,n)*m(i,j,l)/fscale(offsetn+n)
                rxm(i,j,l,n) = 0.0
                rym(i,j,l,n) = 0.0
                rzm(i,j,l,n) = 0.0
#ifdef secmom
                rxxm(i,j,l,n) = 0.0
                rxym(i,j,l,n) = 0.0
                rxzm(i,j,l,n) = 0.0
                ryym(i,j,l,n) = 0.0
                ryzm(i,j,l,n) = 0.0
                rzzm(i,j,l,n) = 0.0
#endif
             end do
          end do
       end do
    end do

    deallocate(field)

    nullify(m)
    nullify(rm)
    nullify(rxm)
    nullify(rym)
    nullify(rzm)
#ifdef secmom
    nullify(rxxm)
    nullify(rxym)
    nullify(rxzm)
    nullify(ryym)
    nullify(ryzm)
    nullify(rzzm)
#endif
    !WP! checked rm with the values in savehdf after scatter
    !    and gather on 4 processors=>OK

  end subroutine readhdfmmr

  subroutine readhdfmmix(region,file_name)
    !
    ! data are expected in mixing ratio from an mmix file
    !
    use dims,        only : im, jm, lm, nregions, datadir, region_name, parent
    use io_hdf
    use global_data, only : mass_dat
    use MeteoData  , only : m_dat
    use chem_param,  only : fscale, ntracet, names
    use toolbox,     only : escape_tm
#ifdef MPI
    use mpi_comm,    only : scatter_after_read_t, stopmpi
    use mpi_const,   only : ntracet_ar, ntracetloc
    use mpi_const,   only : com_trac, mpi_integer, root_t, ierr
#endif
    use ParTools,    only : root, myid, root_t

    implicit none

    ! in/out
    integer, intent(in)             :: region
    character(len=*), intent(in)    :: file_name

    ! local
    real,dimension(:,:,:,:),pointer :: rm,rxm,rym,rzm
    real,dimension(:,:,:),pointer   :: m
    integer                         :: io, sfstart, ifail, sds_id
    integer                         :: n_datasets, n_file_attrs
    integer                         :: istat, attributes, num_type
    integer                         :: sffinfo, sfselect, sfginfo
    integer                         :: sfend, sfrnatt, sfrcatt, sffattr
    integer,dimension(6)            :: idate_save
    integer,dimension(nregions)     :: im_file,jm_file,lm_file
    integer                         :: ntrace_file
    character(len=80)               :: msg_file
    integer                         :: ind,ind1,n,i,j,l,offsetn, from_file, my_parent
    real,dimension(:,:,:,:), allocatable :: field
    real,dimension(:,:,:),   allocatable :: field3d

    ! start

    m => m_dat(region)%data
    rm => mass_dat(region)%rm_t
    rxm => mass_dat(region)%rxm_t
    rym => mass_dat(region)%rym_t
    rzm => mass_dat(region)%rzm_t
    rm = 0.0
    rxm = 0.0
    rym = 0.0
    rzm = 0.0
    allocate(field(-1:im(region)+2,-1:jm(region)+2,1:lm(region),1:ntracet))
    field = 0.0    !09/2004 cmk: to initialise the edges!
    if ( myid == root_t ) then
       from_file = 1
       allocate(field3d(im(region),jm(region),lm(region)))
       io = sfstart(file_name, DFACC_READ)
       if ( io == FAIL ) then
          my_parent = parent(region)
          if ( my_parent == 0 ) then
             call escape_tm('readhdfmmix : Could not open file and no parent '// &
                  'available :'//file_name )
          else
             from_file = 0
             !print *, 'readhdf: Trying to initialise '//region_name(region)// &
             !     ' from parent....'
          endif
       endif
       if ( from_file == 1 ) then   ! read from file:
          !print*,' '
          !print*,'readhdfmmix: ',file_name,'... opened for READ access.'

       istat = sffinfo(io, n_datasets, n_file_attrs)

       !print*, 'readhdfmmix: Number of datasets and attributes',  &
       !     n_datasets, n_file_attrs

       do n = 1, ntracet
          call io_read3d_32(io, im(region), jm(region), lm(region), field3d, names(n), ifail)
          if (ifail == 0) then
             field(1:im(region), 1:jm(region), 1:lm(region), n) = field3d
          else
             !print *, 'Readhdfmmix:', names(n), '  not found in dataset --> 0.0 '
             field(1:im(region), 1:jm(region), 1:lm(region), n) = 0.0
          endif
       enddo

       istat = sfend(io)
       if (istat /= FAIL) then
          !print*,'readhdfmmix: ... file closed'
          !print*,' '
       else
          call escape_tm('readhdfmmix: ERROR in restart from HDF file')
       end if

       deallocate(field3d) ! no longer needed
       endif  ! from file
    endif ! root_t

#ifdef MPI
    ! broadcast from_file !
    call mpi_bcast(from_file, 1, mpi_integer , root_t, com_trac, ierr)
#endif
    if ( from_file == 1 ) then
#ifdef MPI
      call scatter_after_read_t(field,im(region),jm(region),lm(region),&
           2,2,0,rm,root)
#else
      rm = field
#endif
#ifdef MPI
      offsetn = sum(ntracet_ar(0:myid-1))
      do n=1,ntracetloc
#else
      offsetn = 0
      do n=1,ntracet
#endif
         do l = 1,lm(region)
            do j = 1,jm(region)
               do i = 1,im(region)
                  rm(i,j,l,n) = rm(i,j,l,n)*m(i,j,l)/fscale(offsetn+n)  ! from mmix to kg tracer
                  rxm(i,j,l,n) = 0.0
                  rym(i,j,l,n) = 0.0
                  rzm(i,j,l,n) = 0.0
               end do
            end do
         end do
      end do

      deallocate(field)
      call init_child_edge(region)   ! get edge info from parent....
    else   ! initialise from parent:
       call init_child(region)   ! get info from parent....
    end if

    nullify(m)
    nullify(rm)
    nullify(rxm)
    nullify(rym)
    nullify(rzm)
    !WP! checked rm with the values in savehdf after scatter
    !    and gather on 4 processors=>OK

  end subroutine readhdfmmix




  subroutine savehdf(msg,file_name)
    !------------------------------------------------------------------
    ! save all essential model parameters and fields on unit kdisk
    !
    !                                       v 8.4
    !  modified for HDF output
    !------------------------------------------------------------------
    use dims
    use global_data, only : mass_dat
    use MeteoData  , only : m_dat
    use chem_param
    use io_hdf,      only : io_write
    use io_hdf,      only : DFNT_INT32, DFNT_CHAR, DFNT_FLOAT64, DFACC_CREATE
    use datetime,    only : tstamp
    use toolbox,     only : escape_tm
#ifdef MPI
    use mpi_const
    use mpi_comm,only : barrier, gather_tracer_t
#endif
    implicit none

    ! in/out
    character(len=*),intent(in) :: file_name
    character(len=*),intent(in) :: msg

    ! local
    real,dimension(:,:,:,:),pointer :: rm,rxm,rym,rzm
#ifdef secmom
    real,dimension(:,:,:,:),pointer :: rxxm,rxym,rxzm,ryym,ryzm,rzzm
#endif
    real,dimension(:,:,:),  pointer :: m
    real,dimension(:,:,:,:), allocatable :: fieldglob4d

    integer           :: istat, sfsnatt, sfscatt, io, sfstart
    integer           :: sfcreate, sfscompress, sds_id, sfwdata
    integer           :: sfendacc, sfend, dimid3, dimid0, dimid1
    integer           :: dimid2, sfdimid, sfsdmname
    integer           :: region, ind, imr,jmr, lmr

    character(len=64) :: msg_file

    !
    ! first close and reopen mix file (makes the data available)..
    ! do region=1,nregions
    !    if (unit_mix(region) > 0 ) then
    !       istat = sfend(unit_mix(region))
    !       unit_mix(region) = sfstart('mix'//region_name(region)// &
    !           '.hdf',DFACC_WRITE)
    !    end if
    ! end do
    !

    if ( .not. splitsave ) call escape_tm('savehdf: splitsave should be true ')

    do region = 1,nregions

       !! debug ...
       !write (gol,'("io_save/savehdf: create ",a)') trim(file_name)//'_'//region_name(region)//'.hdf'; call goPr

#ifdef MPI
       if ( myid == root_t ) then
#endif
          ind = len_trim(file_name)
          io = sfstart(trim(file_name)//'_'//region_name(region)//'.hdf',DFACC_CREATE)
          !print *, 'savehdf: io unit',io, ' on file ', &
          !     file_name(1:ind)//region_name(region)

          istat = sfsnatt(io,'itau',  DFNT_INT32, 1, itau)
          istat = sfsnatt(io,'nregions',  DFNT_INT32, 1, nregions)
          istat = sfscatt(io,'region_name',  DFNT_CHAR, 6, region_name(region))
          istat = sfsnatt(io,'im',    DFNT_INT32, 1, im(region))
          istat = sfsnatt(io,'jm',    DFNT_INT32, 1, jm(region))
          istat = sfsnatt(io,'lm',    DFNT_INT32, 1, lm(region))
          istat = sfsnatt(io,'dx',    DFNT_FLOAT64, 1, dx/xref(region))
          istat = sfsnatt(io,'dy',    DFNT_FLOAT64, 1, dy/yref(region))
          istat = sfsnatt(io,'dz',    DFNT_FLOAT64, 1, dz/zref(region))
          istat = sfsnatt(io,'xbeg',  DFNT_INT32, 1, xbeg(region))
          istat = sfsnatt(io,'xend',  DFNT_INT32, 1, xend(region))
          istat = sfsnatt(io,'ybeg',  DFNT_INT32, 1, ybeg(region))
          istat = sfsnatt(io,'yend',  DFNT_INT32, 1, yend(region))
          istat = sfsnatt(io,'zbeg',  DFNT_INT32, 1, zbeg(region))
          istat = sfsnatt(io,'zend',  DFNT_INT32, 1, zend(region))
          if ( region /= 1 ) then
             istat = sfsnatt(io,'ibeg',  DFNT_INT32, 1, ibeg(region))
             istat = sfsnatt(io,'iend',  DFNT_INT32, 1, iend(region))
             istat = sfsnatt(io,'jbeg',  DFNT_INT32, 1, jbeg(region))
             istat = sfsnatt(io,'jend',  DFNT_INT32, 1, jend(region))
             istat = sfsnatt(io,'lbeg',  DFNT_INT32, 1, lbeg(region))
             istat = sfsnatt(io,'lend',  DFNT_INT32, 1, lend(region))
          end if
          istat = sfsnatt(io,'xref',  DFNT_INT32, 1, xref(region))
          istat = sfsnatt(io,'yref',  DFNT_INT32, 1, yref(region))
          istat = sfsnatt(io,'zref',  DFNT_INT32, 1, zref(region))
          istat = sfsnatt(io,'tref',  DFNT_INT32, 1, tref(region))
          istat = sfsnatt(io,'ntracet',DFNT_INT32, 1, ntracet)
          istat = sfsnatt(io,'ntrace ',DFNT_INT32, 1, ntrace)
          istat = sfsnatt(io,'nstd',DFNT_INT32, 1, nstd)
          istat = sfsnatt(io,'idate' ,DFNT_INT32, 6, idate)
          istat = sfsnatt(io,'istart',  DFNT_INT32, 1, istart)
          istat = sfsnatt(io,'ndiag', DFNT_INT32, 1, ndiag)
          istat = sfsnatt(io,'nwrite',DFNT_INT32, 1, nwrite)
          istat = sfsnatt(io,'ninst', DFNT_INT32, 1, ninst)
          istat = sfsnatt(io,'ncheck',DFNT_INT32, 1, ncheck)
          istat = sfsnatt(io,'itaui',    DFNT_INT32, 1, itaui)
          istat = sfsnatt(io,'itaue',    DFNT_INT32, 1, itaue)
          istat = sfsnatt(io,'itaut',    DFNT_INT32, 1, itaut)
          istat = sfsnatt(io,'itau0',    DFNT_INT32, 1, itau0)
          istat = sfsnatt(io,'idatei' ,  DFNT_INT32, 6, idatei)
          istat = sfsnatt(io,'idatee' ,  DFNT_INT32, 6, idatee)
          istat = sfsnatt(io,'idatet' ,  DFNT_INT32, 6, idatet)
          istat = sfsnatt(io,'idate0' ,  DFNT_INT32, 6, idate0)
          istat = sfsnatt(io,'icalendo' ,DFNT_INT32, 1, icalendo)
          istat = sfsnatt(io,'iyear0' ,  DFNT_INT32, 1, iyear0)
          istat = sfsnatt(io,'julday0' , DFNT_INT32, 1, julday0)
          istat = sfsnatt(io,'ndiagp1' , DFNT_INT32, 1, ndiagp1)
          istat = sfsnatt(io,'ndiagp2' , DFNT_INT32, 1, ndiagp2)
          istat = sfsnatt(io,'nstep'   , DFNT_INT32, 1, nstep)
          istat = sfsnatt(io,'cpu0'   ,  DFNT_FLOAT64, 1, cpu0)
          istat = sfsnatt(io,'cpu1'   ,  DFNT_FLOAT64, 1, cpu1)
          istat = sfsnatt(io,'ra'     ,  DFNT_FLOAT64, ntracet, ra)
          istat = sfsnatt(io,'fscale' ,  DFNT_FLOAT64, ntrace, fscale)
          istat = sfscatt(io,'names'  ,  DFNT_CHAR, ntrace*8, names)
          istat = sfsnatt(io,'areag'  ,  DFNT_FLOAT64, 1, areag)
          istat = sfsnatt(io,'czeta'  ,  DFNT_FLOAT64, 1, czeta)
          istat = sfsnatt(io,'czetak'  , DFNT_FLOAT64, 1, czetak)
          istat = sfscatt(io,'xlabel'  , DFNT_CHAR, 160, xlabel)
          istat = sfsnatt(io,'istd'    , DFNT_INT32, nstd, istd)
          istat = sfsnatt(io,'newyr'   , DFNT_INT32, 1, newyr)
          istat = sfsnatt(io,'newmonth', DFNT_INT32, 1, newmonth)
          istat = sfsnatt(io,'newday'  , DFNT_INT32, 1, newday)
          istat = sfsnatt(io,'newsrun' , DFNT_INT32, 1, newsrun)
          istat = sfsnatt(io,'cdebug'  , DFNT_INT32, 1, cdebug)
          istat = sfsnatt(io,'limits'  , DFNT_INT32, 1, limits)
          istat = sfsnatt(io,'at'  , DFNT_FLOAT64,lm(1)+1, at)
          istat = sfsnatt(io,'bt'  , DFNT_FLOAT64,lm(1)+1, bt)
          istat = sfscatt(io,'adv_scheme'  , DFNT_CHAR, 5, adv_scheme)
          istat = sfsnatt(io,'nsplitsteps'  , DFNT_INT32, 1, nsplitsteps)
          istat = sfscatt(io,'splitorder', DFNT_CHAR,  nsplitsteps, splitorder)
          msg_file = msg
          istat = sfscatt(io,'msg'  , DFNT_CHAR, 64, msg_file)
#ifdef MPI
       end if  ! root all PEs from here:
#endif
#ifdef MPI
       call barrier
#endif
       m => m_dat(region)%data
       rm => mass_dat(region)%rm_t
       rxm => mass_dat(region)%rxm_t
       rym => mass_dat(region)%rym_t
       rzm => mass_dat(region)%rzm_t
#ifdef secmom
       rxxm => mass_dat(region)%rxxm_t
       rxym => mass_dat(region)%rxym_t
       rxzm => mass_dat(region)%rxzm_t
       ryym => mass_dat(region)%ryym_t
       ryzm => mass_dat(region)%ryzm_t
       rzzm => mass_dat(region)%rzzm_t
#endif

#ifdef MPI
       which_par=previous_par(region)
       if ( which_par /= 'tracer' ) &
            call escape_tm('savehdf: Data should be parallel over tracers')
#endif

       imr = im(region) ; jmr = jm(region) ; lmr = lm(region)
       ! allocate global array on all PEs
       allocate(fieldglob4d(-1:imr+2,-1:jmr+2,lmr, ntracet))

#ifdef MPI
       ! false signals: do not broadcast
       call gather_tracer_t(fieldglob4d,imr,jmr,lmr,2,2,0,ntracet,rm,.false.)
#else
       fieldglob4d = rm
#endif
#ifdef MPI
       if ( myid == root_t ) then
#endif
          call io_write(io,imr,'LON'//region_name(region), &
               jmr,'LAT'//region_name(region),lmr,'HYBRID', &
               m(1:imr,1:jmr,1:lmr),'m')
          call io_write(io,imr,'LON'//region_name(region), &
               jmr,'LAT'//region_name(region),lmr, &
               'HYBRID',ntracet,'NTRACET', &
               fieldglob4d(1:imr,1:jmr,1:lmr,1:ntracet),'rm')

          !! debug ...
          !write (gol,*) 'sss1', sum(fieldglob4d**2), fieldglob4d(12,34,1,1); call goPr
#ifdef MPI
       end if
#endif
#ifdef MPI
       call barrier
       ! false signals: do not broadcast
       call gather_tracer_t(fieldglob4d,imr,jmr,lmr,2,2,0,ntracet,rxm,.false.)
#else
       fieldglob4d = rxm
#endif
#ifdef MPI
       if ( myid == root_t ) then
#endif
          call io_write(io,imr,'LON'//region_name(region), &
               jmr,'LAT'//region_name(region),lmr, &
               'HYBRID',ntracet,'NTRACET', &
               fieldglob4d(1:imr,1:jmr,1:lmr,1:ntracet),'rxm')
#ifdef MPI
       end if
#endif

#ifdef MPI
       call barrier
       ! false signals: do not broadcast
       call gather_tracer_t(fieldglob4d,imr,jmr,lmr,2,2,0,ntracet,rym,.false.)
#else
       fieldglob4d = rym
#endif

#ifdef MPI
       if ( myid == root_t ) then
#endif
          call io_write(io,imr,'LON'//region_name(region), &
               jmr,'LAT'//region_name(region),lmr, &
               'HYBRID',ntracet,'NTRACET', &
               fieldglob4d(1:imr,1:jmr,1:lmr,1:ntracet),'rym')
#ifdef MPI
       end if
       call barrier
       ! false signals: do not broadcast
       call gather_tracer_t(fieldglob4d,imr,jmr,lmr,2,2,0,ntracet,rzm ,.false.)
#else
       fieldglob4d = rzm
#endif

#ifdef MPI
       if ( myid == root_t ) then
#endif
          call io_write(io,imr,'LON'//region_name(region), &
               jmr,'LAT'//region_name(region),lmr, &
               'HYBRID',ntracet,'NTRACET', &
               fieldglob4d(1:imr,1:jmr,1:lmr,1:ntracet),'rzm')
#ifdef MPI
       end if
#endif

#ifdef secmom
! second moments
if (adv_scheme=='2nd_m') then
#ifdef MPI
   call barrier
   call  gather_tracer_t(fieldglob4d,imr,jmr, lmr, 2,2,0,ntracet, rxxm ,.false.)   ! false signals: do not broadcast
#else
   fieldglob4d = rxxm
#endif
#ifdef MPI
   if(myid == root_t) then
#endif
      call io_write(io,imr,'LON'//region_name(region),jmr,'LAT'//region_name(region),lmr, &
                        'HYBRID',ntracet,'NTRACET',fieldglob4d(1:imr,1:jmr,1:lmr,1:ntracet),'rxxm')
#ifdef MPI
   endif
#endif

#ifdef MPI
   call barrier
   call  gather_tracer_t(fieldglob4d,imr,jmr, lmr, 2,2,0,ntracet, rxym ,.false.)   ! false signals: do not broadcast
#else
   fieldglob4d = rxym
#endif

#ifdef MPI
   if(myid == root_t) then
#endif
      call io_write(io,imr,'LON'//region_name(region),jmr,'LAT'//region_name(region),lmr, &
                        'HYBRID',ntracet,'NTRACET',fieldglob4d(1:imr,1:jmr,1:lmr,1:ntracet),'rxym')
#ifdef MPI
   endif
   call barrier
   call  gather_tracer_t(fieldglob4d,imr,jmr, lmr, 2,2,0,ntracet, rxzm ,.false.)   ! false signals: do not broadcast
#else
   fieldglob4d = rxzm
#endif

#ifdef MPI
   if(myid == root_t) then
#endif
      call io_write(io,imr,'LON'//region_name(region),jmr,'LAT'//region_name(region),lmr, &
                        'HYBRID',ntracet,'NTRACET',fieldglob4d(1:imr,1:jmr,1:lmr,1:ntracet),'rxzm')
#ifdef MPI
   endif
#endif

#ifdef MPI
   call barrier
   call  gather_tracer_t(fieldglob4d,imr,jmr, lmr, 2,2,0,ntracet, ryym ,.false.)   ! false signals: do not broadcast
#else
   fieldglob4d = ryym
#endif
#ifdef MPI
   if(myid == root_t) then
#endif
      call io_write(io,imr,'LON'//region_name(region),jmr,'LAT'//region_name(region),lmr, &
                        'HYBRID',ntracet,'NTRACET',fieldglob4d(1:imr,1:jmr,1:lmr,1:ntracet),'ryym')
#ifdef MPI
   endif
#endif

#ifdef MPI
   call barrier
   call  gather_tracer_t(fieldglob4d,imr,jmr, lmr, 2,2,0,ntracet, ryzm ,.false.)   ! false signals: do not broadcast
#else
   fieldglob4d = ryzm
#endif

#ifdef MPI
   if(myid == root_t) then
#endif
      call io_write(io,imr,'LON'//region_name(region),jmr,'LAT'//region_name(region),lmr, &
                        'HYBRID',ntracet,'NTRACET',fieldglob4d(1:imr,1:jmr,1:lmr,1:ntracet),'ryzm')
#ifdef MPI
   endif
   call barrier
   call  gather_tracer_t(fieldglob4d,imr,jmr, lmr, 2,2,0,ntracet, rzzm ,.false.)   ! false signals: do not broadcast
#else
   fieldglob4d = rzzm
#endif

#ifdef MPI
   if(myid == root_t) then
#endif
      call io_write(io,imr,'LON'//region_name(region),jmr,'LAT'//region_name(region),lmr, &
                        'HYBRID',ntracet,'NTRACET',fieldglob4d(1:imr,1:jmr,1:lmr,1:ntracet),'rzzm')
#ifdef MPI
   endif
#endif


endif ! second moments
#endif



#ifdef MPI
       if ( myid == root_t ) then
         istat = sfend( io )
         !print *,'savehdf: sfend returns', istat
       end if
       call barrier
#else
       istat = sfend( io )
       !print *,'savehdf: sfend returns', istat
#endif
       deallocate(fieldglob4d)
       nullify(m)
       nullify(rm)
       nullify(rxm)
       nullify(rym)
       nullify(rzm)
#ifdef secmom
       nullify(rxxm)
       nullify(rxym)
       nullify(rxzm)
       nullify(ryym)
       nullify(ryzm)
       nullify(rzzm)
#endif

    end do  !regions...

    call tstamp(kmain,itau,msg)

    if (cdebug) then
       call tstamp(kdebug,itau,'savehdf')
    end if

  end subroutine savehdf



end module io_save_hdf
