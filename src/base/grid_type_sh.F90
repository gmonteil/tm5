! NAME
!   grid_type_sh  -  manipulate spherical harmonic field
!
! may 2002, Arjo Segers
!

module grid_type_sh


  implicit none
  
  ! --- in/out ----------------------------
  
  private
  
  public  ::  TshGridInfo
  public  ::  TshGrid

  public  ::  Init, Done, Check
  
  public  ::  SpN

  public  ::  Set
  public  ::  Truncate
  
  public  ::  sh_Pnm  !, sh_Pnm_fast
  public  ::  FourrierCoeff

  public  ::  Eval
  public  ::  EvalFast
  public  ::  Eval_Lons

  public  ::  Nabla
  
  public  ::  vod2uv
 
  ! --- const ------------------------------
  
  character(len=*), parameter  ::  mname = 'grid_type_sh'
  
  ! --- types -------------------------
  
  ! Coefficients of Associated Legendre Function
  ! * Triangular truncation: m=0,..,M=T  n=m,..,N=T
  ! * Coefficients p(m,n), stored in order:
  !     (m,n) = (0,0) (0,1) .. (0,N)  (1,1) .. (1,N) ..
  !   The storage order is similar to the order used in the EMOS library.
  !   Only coefficients for m>=0 are stored, since
  !   for real valued ALF the p(-m,n)=conj(p(m,n) .
  
  type TshGridInfo
    integer             ::  T     ! triangular truncation
    integer             ::  np    ! number of data points 
  end type TshGridInfo
  
  type TshGrid
    integer             ::  T
    complex, pointer    ::  c(:)
  end type TshGrid
  
  
  ! --- interfaces ----------------------
  
  interface Init
    module procedure shi_Init
    module procedure shi_Init_copy
    module procedure shgrid_Init
    module procedure shgrid_Init_T
    module procedure shgrid_Init_1
    module procedure shgrid_Init_1_T
  end interface
  
  interface Done
    module procedure shi_Done
    module procedure shgrid_Done
    module procedure shgrid_Done_1
  end interface
  
  interface Check
    module procedure shgrid_Check
    module procedure shgrid_Check_T_1
    module procedure shgrid_Check_T
  end interface
  
  interface Set
    module procedure shi_Set_T
    module procedure shi_Set_shi
    module procedure Set_T
    module procedure Set_c
    module procedure Set_sh
    module procedure Set_T_sh
  end interface
  
  interface Truncate
    module procedure shgrid_Truncate
  end interface
  
  interface FourrierCoeff
    module procedure shgrid_FourrierCoeff_pnm
    module procedure    shi_FourrierCoeff_pnm
    module procedure shgrid_FourrierCoeff_pnm2
  end interface

  interface Eval
    module procedure shgrid_Eval_Xm_lon
    module procedure shgrid_Eval_Pnm_lon
    module procedure shgrid_Eval_lat_lon
  end interface
  
  interface EvalFast
    module procedure shgrid_Eval_Pnm_lon_fast
  end interface
  
  interface Eval_Lons
    module procedure shi_Eval_Lons_Ulat
    module procedure shi_Eval_Lons_UPnm_1
    module procedure shgrid_Eval_Lons_Ulat
    module procedure shgrid_Eval_Lons_UPnm_1
    module procedure shgrid_Eval_Lons_Xm_1
  end interface
  
  interface Nabla
    module procedure shgrid_Nabla
  end interface
  
  interface vod2uv
    module procedure sh_vod2uv
    module procedure shi_vod2uv
  end interface
  
  

contains


  
  ! ====================================================================
  
  
  subroutine shi_Init( shi, T, status )
  
    ! --- in/out ------------------------------------
    
    type(TshGridInfo), intent(out)       ::  shi
    integer, intent(in)                  ::  T
    integer, intent(out)                 ::  status
    
    ! --- begin ----------------------------------

    ! store triangular trunction:
    shi%T = T
    
    ! number of complex data points:
    shi%np = SpN( shi%T )
    
    ! ok
    status = 0

  end subroutine shi_Init


  ! ***
  
  
  subroutine shi_Init_copy( shi, shi2, status )
  
    ! --- in/out ------------------------------------
    
    type(TshGridInfo), intent(out)       ::  shi
    type(TshGridInfo), intent(in)        ::  shi2
    integer, intent(out)                 ::  status
    
    ! --- begin ----------------------------------

    call Init( shi, shi2%T, status )

  end subroutine shi_Init_copy


  ! ***
  

  subroutine shi_Done( shi )
  
    ! --- in/out ------------------------------------
    
    type(TshGridInfo), intent(inout)       ::  shi
    
    ! --- begin ----------------------------------

    ! nothing to be done ...

  end subroutine shi_Done


  ! ***
  
  
  subroutine shi_Set_T( shi, T )
  
    ! --- in/out ------------------------------------
    
    type(TshGridInfo), intent(out)       ::  shi
    integer, intent(in)                  ::  T
    
    ! --- begin ----------------------------------

    shi%T = T

  end subroutine shi_Set_T


  ! ***
  
  
  subroutine shi_Set_shi( shi, shi2 )
  
    ! --- in/out ------------------------------------
    
    type(TshGridInfo), intent(out)       ::  shi
    type(TshGridInfo), intent(in)        ::  shi2
    
    ! --- begin ----------------------------------

    shi%T = shi2%T

  end subroutine shi_Set_shi


  ! ====================================================================
  
  
  subroutine shgrid_Init( SH )
  
    ! --- in/out ------------------------------------
    
    type(TshGrid), intent(out)       ::  SH
    
    ! --- begin ----------------------------------

!    ! check arguments ...    
!    if (associated(SH%c)) then
!      print *, 'shgrid_Init : coefficient array already associated ...'
!      stop
!    end if
    
    SH%T = -1
    nullify( SH%c )

  end subroutine shgrid_Init


  ! ***
  

  subroutine shgrid_Init_T( SH, T )
  
    ! --- in/out ------------------------------------
    
    type(TshGrid), intent(out)       ::  SH
    integer, intent(in)              ::  T
    
    ! --- begin ----------------------------------

!    ! check arguments ...    
!    if (associated(SH%c)) then
!      print *, 'shgrid_Init : coefficient array already associated ...'
!      stop
!    end if

    ! set Triangular Truncation
    if ( T < 0 ) then
      print *, 'shgrid_Init : tried to set triangular truncation ', T
      stop
    end if
    SH%T = T

    ! allocate ...
    allocate( SH%c(SpN(T)) )

    ! set to zero:
    SH%c = (0.0,0.0)
      
  end subroutine shgrid_Init_T


  ! ***
  

  subroutine shgrid_Init_1( SH )
  
    ! --- in/out ------------------------------------
    
    type(TshGrid), intent(out)       ::  SH(:)
    
    ! --- local ----------------------------------
    
    integer          ::  l
    
    ! --- begin ----------------------------------

    do l = 1, size(SH)
    
!      ! check arguments ...    
!      if ( associated(SH(l)%c) ) then
!        print *, 'shgrid_Init : coefficient array already associated ...'
!        stop
!      end if

      SH(l)%T = -1
      nullify( SH(l)%c )

    end do
      
  end subroutine shgrid_Init_1


  ! ***
  

  subroutine shgrid_Init_1_T( SH, T )
  
    ! --- in/out ------------------------------------
    
    type(TshGrid), intent(out)       ::  SH(:)
    integer, intent(in)              ::  T
    
    ! --- local ----------------------------------
    
    integer          ::  l
    
    ! --- begin ----------------------------------

    do l = 1, size(SH)
    
!      ! check arguments ...    
!      if (associated(SH(l)%c)) then
!        print *, 'shgrid_Init : coefficient array already associated ...'
!        stop
!      end if

      ! set Triangular Truncation
      if ( T < 0 ) then
        print *, 'shgrid_Init : tried to set triangular truncation ', T
        stop
      end if
      SH(l)%T = T

      ! allocate ...
      allocate( SH(l)%c(SpN(T)) )

      ! set to zero:
      SH(l)%c = (0.0,0.0)

    end do
      
  end subroutine shgrid_Init_1_T

 
  ! ***
  

  subroutine shgrid_Done( SH )
  
    ! --- in/out ------------------------------------
    
    type(TshGrid), intent(out)    ::  SH
    
    ! --- begin ----------------------------------
    
    ! deallocate if necessary ...
    if (associated(SH%c)) deallocate(SH%c)
    
  end subroutine shgrid_Done
  
  
  ! ***

  
  subroutine shgrid_Done_1( SH )
  
    ! --- in/out ------------------------------------
    
    type(TshGrid), intent(out)    ::  SH(:)
    
    ! --- local ----------------------------------
    
    integer          ::  l
    
    ! --- begin ----------------------------------

    ! deallocate if necessary ...
    do l = 1, size(SH)
      if ( associated(SH(l)%c) ) deallocate( SH(l)%c )
    end do
    
  end subroutine shgrid_Done_1
  
  
  ! ====================================================================


  ! call Check( SH )    : check wether T and c are consistent
  ! call Check( SH, T ) : check wether SH is truncated at T
  
  subroutine shgrid_Check( SH )
  
    ! --- in/out ----------------------------
    
    type(TshGrid), intent(in)   ::  SH
    
    ! --- begin ----------------------------
    

    if (SH%T<0) then

      if (associated(SH%c)) then
        print *, 'shgrid_Check : error : coefficient array associated', &
                 ' while T<0'
        stop
      end if

    else

      if (.not.associated(SH%c)) then
        print *, 'shgrid_Check : error : coefficient array not associated'
        stop
      end if

      if (size(SH%c) < SpN(SH%T) ) then
        print *, 'shgrid_Check : error : coefficient array has size ', size(SH%c), &
           'while triangular truncation ', SH%T, &
           'requires spectral number ', SpN(SH%T)
        stop
      end if

    end if

  end subroutine shgrid_Check


  ! ===


  subroutine shgrid_Check_T( SH, T )
  
    ! --- in/out ----------------------------
    
    type(TshGrid), intent(in)   ::  SH
    integer, intent(in)         ::  T
    
    ! --- begin ----------------------------
    
    if ( T < 0 ) then
      print *, 'shgrid_Check : tried to check for truncation ', T
      stop
    end if

    if ( SH%T /= T ) then
      print *, 'shgrid_Check : truncated at ', SH%T, ' instead of ', T
      stop
    end if

    if (.not.associated(SH%c)) then
      print *, 'shgrid_Check : error : coefficient array not associated'
      stop
    end if

    if (size(SH%c) < SpN(SH%T) ) then
      print *, 'shgrid_Check : error : coefficient array has size ', size(SH%c), &
         'while triangular truncation ', SH%T, &
         'requires spectral number ', SpN(SH%T)
      stop
    end if

  end subroutine shgrid_Check_T


  ! ===


  subroutine shgrid_Check_T_1( SH, T )
  
    ! --- in/out ----------------------------
    
    type(TshGrid), intent(in)   ::  SH(:)
    integer, intent(in)         ::  T
    
    ! --- local ---------------------------
    
    integer     ::  l
    
    ! --- begin ----------------------------
    
    do l = 1, size(SH)
      call Check( SH(l), T )
    end do

  end subroutine shgrid_Check_T_1


  ! ====================================================================

  ! Spectral number :  SN = (T+1)*(T+2)/2
  ! (half triangle  m=0,..,T , n=m,..,T  ;
  !  we assume real restult thus only m>=0 required)
  
  integer function SpN( TT )
  
    ! --- in/out -----------------
    
    integer, intent(in)   ::  TT   ! triangular truncation
    
    ! --- begin ------------------
    
    SpN = (TT+1)*(TT+2)/2
    
  end function SpN


  ! ====================================================================

  
  subroutine Set_T( SH, TT )
  
    ! --- in/out -------------------------
    
    type(TshGrid), intent(inout)    ::  SH
    integer, intent(in)             ::  TT
    
    ! --- local ---------------------------
    
    integer             ::  SN
    
    ! --- in/out -------------------------

    call check( SH )
    
    ! -- set triangular truncation

    ! check ...
    if ( TT < 0 ) then
      print *, 'Set : tried to set triangular truncation ', TT
      stop
    end if

    ! set:
    SH%T = TT
    
    ! -- set spectral coefficients
    
    ! number of coeff to be set:
    SN = SpN(SH%T)

    ! (de)allocate coefficient array if necessary ...
    if ( associated(SH%c) ) then
      if ( size(SH%c)/=SN ) deallocate(SH%c)
    end if
    if (.not. associated(SH%c) ) allocate(SH%c(SN))

    ! fill with zeros    
    SH%c = (0.0,0.0)
      
  end subroutine Set_T
    

  ! ===      


  subroutine Set_c( SH, TT, c )
  
    ! --- in/out -------------------------
    
    type(TshGrid), intent(inout)   ::  SH
    integer, intent(in)            ::  TT
    complex, intent(in)            ::  c(:)
    
    ! --- local ---------------------------
    
    integer             ::  SN
    
    ! --- in/out -------------------------
    
    ! init with zeros:
    call Set( SH, TT )
    
    ! number of coeff to be set:
    SN = SpN(SH%T)

    ! fill with c
    if ( size(c) /= SN ) then
      print *, 'Set_c : size of input array is ', size(c), &
        ' while triangular truncation ', TT, ' requires ', SN, &
        ' spectral coefficients'
      stop
    end if
    SH%c = c
    
  end subroutine Set_c
  
  
  ! *
  
  
  subroutine Set_T_sh( SH, TT, sh2 )
  
    ! --- in/out -------------------------
    
    type(TshGrid), intent(inout)   ::  SH
    integer, intent(in)            ::  TT
    type(TshGrid), intent(in)      ::  sh2
    
    ! --- in/out -------------------------
    
    ! init with full truncation:
    call Set( SH, sh2%T, sh2%c )
    
    ! truncate to requested number:
    call Truncate( SH, TT )
    
  end subroutine Set_T_sh
  
  
  ! *
  
  
  subroutine Set_sh( SH, sh2 )
  
    ! --- in/out -------------------------
    
    type(TshGrid), intent(inout)   ::  SH
    type(TshGrid), intent(in)      ::  sh2
    
    ! --- in/out -------------------------
    
    ! init with full truncation:
    call Set( SH, sh2%T, sh2%c )
    
  end subroutine Set_sh
  
  
  ! ====================================================================
  
  ! truncates spherical harmonic coeffiecients to lower
  ! triangular truncation.
  
  subroutine shgrid_Truncate( SH, T_new )
  
    ! --- in/out ------------------------------

    type(TshGrid), intent(inout)    ::  SH
    integer, intent(in)          ::  T_new
    
    ! --- local ----------------------------
    
    integer               ::  T, m, n
    integer               ::  k, k_new
    complex, allocatable  ::  c_new(:)
    
    ! ---- begin ----------------------------
    
    ! check arguments ...
    call Check( SH )
    if (T_new<0) then
      print *, 'shgrid_Truncate : error : tried to truncate to T=', T_new
      stop
    end if
    if (T_new > SH%T) then
      print *, 'shgrid_Truncate : error : tried to truncate to T=', T_new, &
        ' while input spherical harmonics are truncated at ', SH%T
      stop
    end if
    
    ! shorthands:
    T = SH%T

    ! truncate only if necessary ...
    if ( T_new < T ) then

      ! allocate temporary array for truncated coefficients:
      allocate( c_new(SpN(T_new)) )

      ! loop over original triangle for m <= T_new;
      k = 0
      k_new = 0
      do m = 0, T_new
        do n = m, T
          k = k + 1
          if (n <= T_new) then
            k_new = k_new + 1
            c_new(k_new) = SH%c(k)
          end if
        end do
      end do
      
      ! fill SH with truncated coefficients:
      call Set( SH, T_new, c_new )
      
      ! done
      deallocate( c_new )
      
    end if
    
  end subroutine shgrid_Truncate
    

  ! =======================================================
  
  !
  ! Evaluate associate Legendre functions at given latitude.
  ! Corner (T,T) is set to 0.0   .
  !
  ! Based on recurent formula:
  !
  !
  !   mu P(mu;m,n-1) = eps(m,n) P(mu;m,n) + eps(m,n-1) P(mu;m,n-2)
  !
  !                                          n^2 - m^2
  !      mu = sin(lat)  ,   eps(m,n) = sqrt( --------- )
  !                                          4n^2 - 1
  !
  !   P(m,n) = ( mu P(m,n-1) - eps(m,n-1) P(m,n-2) ) / eps(m,n)
  !
  !
  !   P(0,0) = 1
  !
  !
  !   P(0,1) = mu / eps(0,1) = mu sqrt(3)
  !
  !                                   
  !              2m+1   1/2  sqrt(1-mu^2)^m 
  !   P(m,m) = ( ----- )     -------------- (2m)!
  !              (2m)!            2^m m!
  !
  !                          sqrt(1-mu^2)^m 
  !          = sqrt((2m+1)!) --------------        
  !                             2^m m!
  !
  !                                    sqrt(1-mu^2)sqrt(1-mu^2)^(m-1)
  !          = sqrt((2m+1)2m(2(m-1))!) ------------------------------ 
  !                                         2 2^(m-1) m (m-1)!
  !
  !                                      cos(lat)
  !          = P(m-1,m-1) sqrt((2m+1)2m) --------
  !                                        2m
  !
  !          = P(m-1,m-1) sqrt((2m+1)) cos(lat) / sqrt(2m)
  !
  !   P(m,m+1) = mu P(m,m) / eps(m,m+1)
  !
  !

  subroutine sh_Pnm( Pnm, T, lat, status )

    ! --- in/out ----------------------------------

    real, intent(out)     ::  Pnm(:)
    integer, intent(in)   ::  T
    real, intent(in)      ::  lat      ! rad
    integer, intent(out)  ::  status

    ! --- const --------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sh_Pnm'
    
    ! --- local ------------------------------------

    integer   ::  m, n
    integer   ::  k
    real      ::  mu, rmu
    real      ::  fmm, fmmp
    real      ::  eps, eps1
    
    ! --- begin ---------------------------------------

    ! Use EMOS library:    
    !Pnm = 0.0
    !call jspleg1( Pnm, lat, T-1 )
    !return

    if ( size(Pnm) /= SpN(T) ) then
      write (*,'("ERROR - wrong size of output array:")')
      write (*,'("ERROR -   size(Pnm)   : ",i6)') size(Pnm)
      write (*,'("ERROR -   expected    : ",i6)') SpN(T)
      write (*,'("ERROR -   truncation  : ",i6)') T
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if
    
    mu   = sin(lat)
    rmu   = sqrt(1.0-mu*mu)    ! cos(lat)

    do m = 0, T-1

      if ( m == 0 ) then
      
        fmmp    = sqrt(3.0)
        k = 1; Pnm(k) = 1.0          ! (0,0)
        k = 2; Pnm(k) = fmmp*mu      ! (0,1)

      else

        fmm  = fmmp * rmu / sqrt( 2.0*m )
        fmmp = fmm * sqrt( 2*m + 3.0 )

        ! n = m and n = m+1

        k = k+1; Pnm(k) = fmm
        k = k+1; Pnm(k) = fmmp * mu

      endif

      ! leave if truncation is reached:
      if ( m == T-1 ) exit

      eps1  = 1.0 / sqrt( 2*m + 3.0 )

      do n = m+2, T

        eps = sqrt((n*n*1.0-m*m*1.0)/(4.0*n*n-1.0))
        k = k+1; Pnm(k) = ( mu*Pnm(k-1) - eps1*Pnm(k-2) ) / eps
        eps1 = eps

      end do

    end do
    
    ! set corner:
    k = k + 1; Pnm(k) = 0.0
    
    ! ok
    status = 0
    
  end subroutine sh_Pnm


  ! ***
  
  
  ! same, but without if statements in do loop
  ! (faster when compiled with pgf90)

  subroutine sh_Pnm_fast( Pnm, T, lat, status )

    ! --- in/out ----------------------------------

    real, intent(out)     ::  Pnm(:)
    integer, intent(in)   ::  T
    real, intent(in)      ::  lat      ! rad
    integer, intent(out)  ::  status

    ! --- const --------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sh_Pnm_fast'
    
    ! --- local ------------------------------------

    integer   ::  m, n
    integer   ::  k
    real      ::  mu, rmu
    real      ::  fmm, fmmp
    real      ::  eps, eps1
    
    ! --- begin ---------------------------------------

    ! Use EMOS library:    
    !Pnm = 0.0
    !call jspleg1( Pnm, lat, T-1 )
    !return

    if ( size(Pnm) /= SpN(T) ) then
      write (*,'("ERROR - wrong size of output array:")')
      write (*,'("ERROR -   size(Pnm)   : ",i6)') size(Pnm)
      write (*,'("ERROR -   expected    : ",i6)') SpN(T)
      write (*,'("ERROR -   truncation  : ",i6)') T
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if
    
    mu   = sin(lat)
    rmu   = sqrt(1.0-mu*mu)    ! cos(lat)

    ! -- first m
    
    m = 0

    fmmp    = sqrt(3.0)
    k = 1; Pnm(k) = 1.0          ! (0,0)
    k = 2; Pnm(k) = fmmp*mu      ! (0,1)

    eps1  = 1.0 / sqrt( 2*m + 3.0 )

    do n = m+2, T

      eps = sqrt((n*n*1.0-m*m*1.0)/(4.0*n*n-1.0))
      k = k+1; Pnm(k) = ( mu*Pnm(k-1) - eps1*Pnm(k-2) ) / eps
      eps1 = eps

    end do

    ! -- 0 < m < T-1
    
    do m = 1, T-2

      fmm  = fmmp * rmu / sqrt( 2.0*m )
      fmmp = fmm * sqrt( 2*m + 3.0 )

      ! n = m and n = m+1

      k = k+1; Pnm(k) = fmm
      k = k+1; Pnm(k) = fmmp * mu

      eps1  = 1.0 / sqrt( 2*m + 3.0 )

      do n = m+2, T

        eps = sqrt((n*n*1.0-m*m*1.0)/(4.0*n*n-1.0))
        k = k+1; Pnm(k) = ( mu*Pnm(k-1) - eps1*Pnm(k-2) ) / eps
        eps1 = eps

      end do

    end do
    
    ! -- m = T-1
    
    m = T-1

    fmm  = fmmp * rmu / sqrt( 2.0*m )
    fmmp = fmm * sqrt( 2*m + 3.0 )

    ! n = m and n = m+1

    k = k+1; Pnm(k) = fmm
    k = k+1; Pnm(k) = fmmp * mu

    eps1  = 1.0 / sqrt( 2*m + 3.0 )

    do n = m+2, T

      eps = sqrt((n*n*1.0-m*m*1.0)/(4.0*n*n-1.0))
      k = k+1; Pnm(k) = ( mu*Pnm(k-1) - eps1*Pnm(k-2) ) / eps
      eps1 = eps

    end do

    ! set corner:
    k = k + 1; Pnm(k) = 0.0
    
    
  end subroutine sh_Pnm_fast



  ! ====================================================================

  !                     M     N                           i m lon
  !  Eval(p,lon,lat) = sum   sum  p(m,n) P(m,n,sin(lat)) e
  !                    m=-M  n=m
  !
  !               M
  !   = c(0)  +  sum 2[Re{c(m)}cos(m lon) - Im{c(m)}sin(m lon)]  , 
  !              m=1
  !
  !              N                         
  !      c(m) = sum  p(m,n) P(m,n,sin(lat))
  !             n=m
  !
  
  real function shgrid_Eval_Xm_lon( Xm, T, lon )
  
    ! --- in/out -------------------------
    
    integer, intent(in)            ::  T
    complex, intent(in)            ::  Xm(0:T)
    real, intent(in)               ::  lon
    
    ! --- local --------------------------
    
    integer              ::  m

    ! --- begin --------------------------
    
    ! summation over triangle
    shgrid_eval_Xm_lon = Xm(0)
    do m = 1, T
      shgrid_eval_Xm_lon = shgrid_eval_Xm_lon + 2.0*( real(Xm(m))*cos(m*lon) - aimag(Xm(m))*sin(m*lon) )
    end do
    
  end function shgrid_Eval_Xm_lon
  
  
  ! ***
  
  
  real function shgrid_Eval_Pnm_lon( SH, Pnm, lon )
  
    ! --- in/out -------------------------
    
    type(TshGrid), intent(in)      ::  SH
    real, intent(in)               ::  Pnm(:)
    real, intent(in)               ::  lon
    
    ! --- local --------------------------
    
    integer              ::  T, SN
    integer              ::  m, n, k

    complex              ::  cm

    ! --- begin --------------------------
    
    ! shorthands:    
    T = SH%T
    SN = SpN(T)
    
    ! check arguments:
    call Check( SH )
    if ( size(Pnm) /= SN ) then
      print *, 'shgrid_FourrierCoeff_pnm : size(Pnm)=',size(Pnm), &
          'while',SN,'expected for T',T
      stop
    end if

    ! summation over triangle
    k = 0
    do m = 0, T
      cm = (0.0,0.0)
      do n = m, T 
        k = k + 1
        cm = cm  +  SH%c(k) * Pnm(k)
      end do
      if (m==0) then
        shgrid_eval_Pnm_lon = real(cm)
      else
        shgrid_eval_Pnm_lon = shgrid_eval_Pnm_lon + 2.0*( real(cm)*cos(m*lon) - aimag(cm)*sin(m*lon) )
      end if
    end do

  end function shgrid_Eval_Pnm_lon
  
  
  ! ***
  
  
  real function shgrid_Eval_Pnm_lon_fast( SH, Pnm, lon )
  
    ! --- in/out -------------------------
    
    type(TshGrid), intent(in)      ::  SH
    real, intent(in)               ::  Pnm(:)
    real, intent(in)               ::  lon
    
    ! --- local --------------------------
    
    integer              ::  T, SN
    integer              ::  m, n, k

    complex              ::  cm
    complex              ::  cp(size(Pnm))
    real                 ::  res

    ! --- begin --------------------------
    
    ! shorthands:    
    T = SH%T
    SN = SpN(T)
    
    ! check arguments:
    call Check( SH )
    if ( size(Pnm) /= SN ) then
      print *, 'shgrid_FourrierCoeff_pnm : size(Pnm)=',size(Pnm), &
          'while',SN,'expected for T',T
      stop
    end if

    cp = SH%c * Pnm

    ! summation over triangle
    res = sum( cp(1:T+1) )
    k = T+2
    do m = 1, T
      cm = sum( cp(k:k+T-m) )
      k = k + T-m+1
      res = res + 2.0*( real(cm)*cos(m*lon) - aimag(cm)*sin(m*lon) )
    end do
    
    shgrid_eval_Pnm_lon_fast = res
    
  end function shgrid_Eval_Pnm_lon_fast
  
  
  ! ***
  
  
  real function shgrid_Eval_lat_lon( SH, lat, lon, status )
  
    ! --- in/out -------------------------
    
    type(TshGrid), intent(in)   ::  SH
    real, intent(in)            ::  lon, lat
    integer, intent(out)        ::  status
    
    ! --- const --------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/shgrid_Eval_lat_lon'
    
    ! --- local --------------------------
    
    integer              ::  T
    real, allocatable    ::  Pnm(:)

    ! --- begin --------------------------
    
    ! check arguments:
    call Check( SH )

    ! shorthands:    
    T = SH%T
    
    ! allocate array for associated Legendre functions:
    allocate( Pnm(SpN(T)) )

    ! evaluate Legendre functions at given latitude:
    call sh_Pnm( Pnm, T, lat, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! expand summation:
    shgrid_Eval_lat_lon = Eval( SH, Pnm, lon )

    ! done
    deallocate( Pnm )
    
    ! ok
    status = 0
    
  end function shgrid_Eval_lat_lon
  
  
  
 
  ! ==============================================
  
  ! compute Fourrier coefficients:
  !
  !               T
  !   X(m;mu) =  sum  U(m,n) P(m,n;mu)  =  Xe(m;mu) + Xo(m;mu)
  !             n=|m|
  !
  !               T                            T
  !        =     sum  U(m,n) P(m,n;mu)   +    sum  U(m,n) P(m,n;mu)
  !             n=|m|                        n=|m|
  !           (n-m) even                   (n-m) odd
  !
  ! and return coefficients for both +mu and -mu using that
  !
  !   X(m,1) = X(m;+mu)  = Xe(m;mu) + Xo(m;mu)          'north'
  !
  !   X(m,2) = X(m;-mu) = Xe(m;mu) - Xo(m;mu)          'south'
  !
  ! (valid for real symmetric Pmn(mu) )
  
!  subroutine shgrid_FourrierCoeff( Xnm, Pnm, T, X )
!  
!    ! --- in/out -----------------------------------
!    
!    type(TshGrid), intent(in)   ::  Xnm, Pnm
!    integer, intent(in)      ::  T
!    complex, intent(out)     ::  X(0:T,2)
!    
!    ! --- local ---------------------------------
!    
!    complex        ::  Xe, Xo
!    integer        ::  m, n, k
!    
!    ! --- begin ----------------------------------
!    
!    ! check input ...
!    call Check( Xnm, T )
!    call Check( Pnm, T )
!  
!    ! loop over all coeff:
!    m = 0
!    n = 0
!    Xe = (0.0,0.0)
!    Xo = (0.0,0.0)
!    do k = 1, SpN(T)
!      if ( mod(n+m,2)==0 ) then
!        Xe = Xe + Xnm%c(k)*Pnm%c(k)
!      else
!        Xo = Xo + Xnm%c(k)*Pnm%c(k)
!      end if
!      ! next n
!      n = n + 1
!      ! end of current m ?
!      if ( n > T ) then
!        ! save fourrier coeff for current m :
!        X(m,1) = Xe + Xo
!        X(m,2) = Xe - Xo
!        ! next m
!        m = m + 1
!        n = m
!        Xe = (0.0,0.0)
!        Xo = (0.0,0.0)
!      end if
!    end do
!    
!    ! adhoc: ensure that X(0) is real:
!    X(0,:) = real(X(0,:)) * (1.0,0.0)
!    
!  end subroutine shgrid_FourrierCoeff
  
  
  
  
  subroutine shgrid_FourrierCoeff_pnm( Xnm, Pnm, T, X )
  
    ! --- in/out -----------------------------------
    
    type(TshGrid), intent(in)   ::  Xnm
    real, intent(in)            ::  Pnm(:)
    integer, intent(in)         ::  T
    complex, intent(out)        ::  X(0:T)
    
    ! --- local ---------------------------------
    
    complex        ::  Xe, Xo
    integer        ::  m, n, k
    
    ! --- begin ----------------------------------
    
    ! check input ...
    call Check( Xnm, T )
    if ( size(Pnm) /= SpN(T) ) then
      print *, 'shgrid_FourrierCoeff_pnm : size(Pnm)=',size(Pnm), &
          'while',SpN(T),'expected for T',T
      stop
    end if
  
    ! loop over all coeff:
    k = 0
    do m = 0, T
      X(m) = (0.0,0.0)
      do n = m, T 
        k = k + 1
        X(m) = X(m) + Xnm%c(k) * Pnm(k)
      end do
    end do
    
    ! adhoc: ensure that X(0) is real:
    X(0) = real(X(0)) * (1.0,0.0)
    
  end subroutine shgrid_FourrierCoeff_pnm
  
  
  ! ***
  
  
  pure subroutine shi_FourrierCoeff_pnm( shi_in, Xnm, shi, Pnm, X )
  
    ! --- in/out -----------------------------------
    
    type(TshGridInfo), intent(in)       ::  shi_in
    complex, intent(in)                 ::  Xnm(1:shi_in%np)
    type(TshGridInfo), intent(in)       ::  shi
    real, intent(in)                    ::  Pnm(1:shi%np)
    complex, intent(out)                ::  X(0:shi%T)
    
    ! --- const --------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/shi_FourrierCoeff_pnm'
    
    ! --- local ---------------------------------
    
    integer        ::  m, n, kx, kp
    
    ! --- begin ----------------------------------

#ifdef check_all    
    ! check input ...
    if ( shi_in%T < shi%T ) then
      write (*,'("ERROR - input truncation should not be less than output:")')
      write (*,'("ERROR -   shi_in%T    : ",i6)') shi_in%T
      write (*,'("ERROR -   shi%T       : ",i6)') shi%T
      write (*,'("ERROR in ",a)') rname; stop
    end if
    if ( size(Xnm) /= shi_in%np ) then
      write (*,'("ERROR - mismatch between size of Xnm and specified truncation:")')
      write (*,'("ERROR -   shi_in%T    : ",i6)') shi_in%T
      write (*,'("ERROR -   shi_in%np   : ",i6)') shi_in%np
      write (*,'("ERROR -   Xnm         : ",i6)') size(Xnm)
      write (*,'("ERROR in ",a)') rname; stop
    end if
    if ( size(Pnm) /= shi%np ) then
      write (*,'("ERROR - mismatch between size of Pnm and specified truncation:")')
      write (*,'("ERROR -   shi%T       : ",i6)') shi%T
      write (*,'("ERROR -   shi%np      : ",i6)') shi%np
      write (*,'("ERROR -   Pnm         : ",i6)') size(Pnm)
      write (*,'("ERROR in ",a)') rname; stop
    end if
#endif
  
    ! loop over all coeff:
    kx = 0
    kp = 0
    do m = 0, shi%T
      X(m) = (0.0,0.0)
      do n = m, shi%T 
        kx = kx + 1
        kp = kp + 1
        X(m) = X(m) + Xnm(kx) * Pnm(kp)
      end do
      kx = kx + shi_in%T - shi%T
    end do
    
    ! adhoc: ensure that X(0) is real:
    X(0) = real(X(0)) * (1.0,0.0)
    
  end subroutine shi_FourrierCoeff_pnm
  
  
  ! ***
  
  
  subroutine shgrid_FourrierCoeff_pnm2( Xnm, Pnm, T, X )
  
    ! --- in/out -----------------------------------
    
    type(TshGrid), intent(in)   ::  Xnm
    real, intent(in)            ::  Pnm(:)
    integer, intent(in)         ::  T
    complex, intent(out)        ::  X(0:T,2)
    
    ! --- local ---------------------------------
    
    complex        ::  Xe, Xo
    integer        ::  m, n, k
    
    ! --- begin ----------------------------------
    
    ! check input ...
    call Check( Xnm, T )
    if ( size(Pnm) /= SpN(T) ) then
      print *, 'shgrid_FourrierCoeff_pnm : size(Pnm)=',size(Pnm), &
          'while',SpN(T),'expected for T',T
      stop
    end if
  
    ! loop over all coeff:
    m = 0
    n = 0
    Xe = (0.0,0.0)
    Xo = (0.0,0.0)
    do k = 1, SpN(T)
      if ( mod(n+m,2)==0 ) then
        Xe = Xe + Xnm%c(k)*Pnm(k)
      else
        Xo = Xo + Xnm%c(k)*Pnm(k)
      end if
      ! next n
      n = n + 1
      ! end of current m ?
      if ( n > T ) then
        ! save fourrier coeff for current m :
        X(m,1) = Xe + Xo
        X(m,2) = Xe - Xo
        ! next m
        m = m + 1
        n = m
        Xe = (0.0,0.0)
        Xo = (0.0,0.0)
      end if
    end do
    
    ! adhoc: ensure that X(0) is real:
    X(0,:) = real(X(0,:)) * (1.0,0.0)
    
  end subroutine shgrid_FourrierCoeff_pnm2
  
  
  ! ====
  
  !
  ! Evaluate x(lat) and x(-lat) given:
  !   Xnm      : complex Associated Legendre Func. coefficients
  !   Pnm      : Associated Legendre Functions evaluated at mu=sin(lat)
  !   K        : number of lon values on a complete latitude circle;
  !              distance between to lon values is thus 2pi/K;
  !   lon0     : western boundary of lon grid (rad)
  !   Kout     : actual number of lon points from lon0 to east .
  !
  ! Evaluation at oposite latitdues at the same time saves some
  ! computation time.
  !
  !  x(t0+k2pi/K)                             ,  k=0,1,..,K0-1
  !
  !         T         i m [t0 + k 2pi/K]
  !     =  sum  X(m) e                        ,  k=0,1,..,K0-1
  !        m=-T                  
  !
  !              L such that KL >= (2T+1)
  !              k=(j+1)/L  ,  j=kL-1   ,   J=KL
  !            
  !         T         i m t0  i m (j+1)2pi/KL
  !     =  sum  X(m) e       e                    ,  j=-1,L-1,2L-1,...,(K0-1)L-1
  !        m=-T                        
  !
  !         T         im[t0+2pi/J]   i m j 2pi/J
  !     =  sum  [X(m)e            ] e             ,  j=-1,L-1,2L-1,...,(K0-1)L-1
  !        m=-T                         
  !
  !   _            
  !   X(m) = X(m) exp(i phi)    ,   phi=m[t0+2pi/J]
  !        =  (a+ib)(cos(phi)+isin(phi))
  !        = [ a cos(phi) - b sin(phi) ]  +  i [ a sin(phi) + b cos(phi) ]
  !
  ! FFT99 uses that X(-m) = conjg(X(m)), and thus requires X(m) for m>=0 only;
  ! higher fourrier coefficients are zero:
  !    X(0)  X(1)  ...  X(T)  0  0  0  ...      
  !   |---    1+floor(J/2) elements     ---|
  !
  !         o   o   x   x   x   0   0 ... 0    0   x   x
  !  m =            0   1   2                     J-2 J-1
  !
  ! Number of K longitudes on complete latitude circle,
  ! while Kout evaluations are required.
  ! Note that Kout might be both smaller and larger than K.
  ! Routine FFT99 returns J+2 evaluations:
  !    x(j=-1) x(j=0)  ...  x(j=J-1=-1) x(j=J=0)
  ! thus if Kout exceeds K+2 then the result should be taken 
  ! cyclic from the output of FFT99.
  !
  
  ! Evaluate given spherical harmonic coefficients and latitude

  subroutine shi_Eval_Lons_Ulat( shi, Xnm, lat, KK, lon_start, Kout, llgrid, status )

    ! --- in/out -------------------------
    
    type(TshGridInfo), intent(in)   ::  shi
    complex, intent(in)             ::  Xnm(1:shi%np)
    real, intent(in)                ::  lat
    integer, intent(in)             ::  KK
    real, intent(in)                ::  lon_start
    integer, intent(in)             ::  Kout
    real,intent(out)                ::  llgrid(Kout)
    integer, intent(out)            ::  status
    
    ! --- const --------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/shi_Eval_Lons_Ulat'
    
    ! --- local --------------------------

    complex   ::  X(0:shi%T)
    real      ::  Pnm(1:shi%np)
    
    ! --- begin --------------------------
    
    ! evaluate Legendre functions at given latitude:
    call sh_Pnm( Pnm, shi%T, lat, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! compute Fourrier coefficients
    ! (latitude is implicit defined in Pnm)
    call FourrierCoeff( shi, Xnm, shi, Pnm, X )

    ! evaluate at longitudes given Fourrier coefficients:
    call Eval_Lons( llgrid, X, shi%T, KK, lon_start, Kout )
    
    ! ok
    status = 0
    
  end subroutine shi_Eval_Lons_Ulat
  
  subroutine shgrid_Eval_Lons_Ulat( llgrid, Xnm, lat, KK, lon_start, Kout, status )

    ! --- in/out -------------------------
    
    type(TshGrid), intent(in)  ::  Xnm
    real, intent(in)           ::  lat
    integer, intent(in)        ::  KK
    real, intent(in)           ::  lon_start
    integer, intent(in)        ::  Kout
    real,intent(out)           ::  llgrid(Kout)
    integer, intent(out)       ::  status
    
    ! --- const --------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/shgrid_Eval_Lons_Ulat'
    
    ! --- local --------------------------

    integer                ::  T
    complex, allocatable   ::  X(:,:)
    real, allocatable      ::  Pnm(:)
    
    ! --- begin --------------------------
    
    ! timing:
    !   spleg1   65 %
    !   fourrier      25 %
    !   eval          10 %
    
    ! extract triangular truncation    
    T = Xnm%T
    
    ! allocate array for associated Legendre functions:
    allocate( Pnm(SpN(T)) )

    ! evaluate Legendre functions at given latitude:
    call sh_Pnm( Pnm, T, lat, status )
    if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if

    ! compute Fourrier coefficients
    ! (latitude is implicit defined in Pnm)
    allocate( X(0:T,2) )
    call FourrierCoeff( Xnm, Pnm, T, X )

    ! evaluate at longitudes given Fourrier coefficients:
    call Eval_Lons( llgrid, X(:,1), T, KK, lon_start, Kout )
    
    ! done
    deallocate( X )
    deallocate( Pnm )
    
    ! ok
    status = 0
    
  end subroutine shgrid_Eval_Lons_Ulat
  

  ! ***     
  

  ! Evaluate given spherical harmonic coefficients and Legendre functions

  subroutine shgrid_Eval_Lons_UPnm_1( llgrid, Xnm, Pnm, KK, lon_start, Kout, status )

    ! --- in/out -------------------------
    
    type(TshGrid), intent(in)  ::  Xnm
    real, intent(in)           ::  Pnm(:)
    integer, intent(in)        ::  KK
    real, intent(in)           ::  lon_start
    integer, intent(in)        ::  Kout
    real,intent(out)           ::  llgrid(Kout)
    integer, intent(out)       ::  status
    
    ! --- const --------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/shgrid_Eval_Lons_UPnm_1'
    
    ! --- local --------------------------

    integer                ::  T
    complex, allocatable   ::  X(:)
    
    ! --- begin --------------------------
    
    T = Xnm%T
    
    ! compute Fourrier coefficients at oposite latitudes
    ! (latitude is implicit defined in Pnm)
    allocate( X(0:T) )
    call FourrierCoeff( Xnm, Pnm, T, X )

    ! evaluate given Fourrier coefficients:    
    call Eval_Lons( llgrid, X, T, KK, lon_start, Kout )
    
    ! done
    deallocate( X )
    
    ! ok
    status = 0
    
  end subroutine shgrid_Eval_Lons_UPnm_1

  
  ! ***     
  

  ! Evaluate given spherical harmonic coefficients and Legendre functions
  !   o TXnm, Xnm        : input truncation and sh coeff
  !   o T, Pnm         : target truncation and evaluated Legende functions
  !   o KK             : number of longitudes on full circle
  !   o lon_start      : first longitude
  !   o Kout           : number of output longitudes
  !   o llgrid         : output array

  subroutine shi_Eval_Lons_UPnm_1( shi_in, Xnm, shi, Pnm, KK, lon_start, Kout, llgrid, status )

    ! --- in/out -------------------------
    
    type(TShGridInfo), intent(in)     ::  shi_in
    complex, intent(in)               ::  Xnm(shi_in%np)
    type(TShGridInfo), intent(in)     ::  shi
    real, intent(in)                  ::  Pnm(shi%np)
    integer, intent(in)               ::  KK
    real, intent(in)                  ::  lon_start
    integer, intent(in)               ::  Kout
    real,intent(out)                  ::  llgrid(Kout)
    integer, intent(out)              ::  status
    
    ! --- const --------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/sh_Eval_Lons_UPnm_1'
    
    ! --- local --------------------------

    complex       ::  X(0:shi%T)
    
    ! --- begin --------------------------

    ! compute Fourrier coefficients at oposite latitudes
    ! (latitude is implicit defined in Pnm)
    call FourrierCoeff( shi_in, Xnm, shi, Pnm, X )

    ! evaluate given Fourrier coefficients:    
    call Eval_Lons( llgrid, X, shi%T, KK, lon_start, Kout )
    
    ! ok
    status = 0
    
  end subroutine shi_Eval_Lons_UPnm_1

  

  ! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  ! Shift to first lon:
  !
  !  x(t0+k2pi/K)                             ,  k=0,1,..,KK-1
  !
  !         T         i m [t0 + k 2pi/K]
  !     =  sum  X(m) e                        ,  k=0,1,..,K0-1
  !        m=-T                  
  !
  !         T          i m t0   k 2pi/K
  !     =  sum  [X(m) e      ] e              ,  k=0,1,..,K0-1
  !        m=-T                  
  !
  ! PERFORMANCE
  !   This routine is about 7% slower than with emos fft ...
  !
  !  real(4), external ::  etime
  !  real ::  tarr(2), t0
  !  t0=etime(tarr)
  !  print *, 'time:',etime(tarr)-t0
    
  
  subroutine shgrid_Eval_Lons_Xm_1( llgrid, Xm, T, KK, lon_start, Kout )

    use grid_tools, only : pi
    use grid_singleton, only : fft
  
    ! --- in/out -------------------------
    
    integer, intent(in)        ::  T
    complex, intent(in)        ::  Xm(0:T)
    integer, intent(in)        ::  KK
    real, intent(in)           ::  lon_start
    integer, intent(in)        ::  Kout
    real,intent(out)           ::  llgrid(Kout)
    
    ! --- local --------------------------

    integer                ::  i, j, k
    integer                ::  m, n
    
    complex                ::  X(0:T)

    integer                ::  LL, JJ
    complex, allocatable   ::  C(:), F(:)
    real                   ::  fac

    ! --- begin --------------------------
    
    ! shift over west most longitude:
    !   X(m) :=  X(m) exp( i m lon0 )
    do m = 0, T
      X(m) = Xm(m) * exp( (0.0,1.0)*m*lon_start )
    end do

    ! choose L such that JJ=KK*LL >= 2T+1
    LL = 1
    do
      if ( KK*LL >= 2*T+1 ) exit
      LL = LL * 2
    end do
    
    ! number of longitudes to be evaluated:
    JJ = KK*LL

    ! fft arrays:
    allocate( C(JJ) )
    allocate( F(JJ) )
    
    ! fill coeff array:
    C = (0.0,0.0)
    C(1) = X(0)
    do m = 1, T
      C(1+m)    = X(m)
      C(JJ+1-m) = conjg(X(m))
    end do

    ! Apply fast fourrier transform.
    !
    !            1    J-1       i m lambda(j)
    !  F(j) = ------- sum C(m) e
    !         sqrt(J) m=0
    !
    !  where
    !
    !    lambda(j) = j 2pi/J = 0, dlon, 2*dlon, ...
    !                                           ____      ____
    !    C = [ X(0), X(1), .., X(T), 0, ..., 0, X(T), .., X(1) ]
    !
    ! Coeff in C are complex conj, thus imag part of F is zero.returned array is real.
    !
    F = fft( C, inv=.true. )
    fac = sqrt(JJ*1.0)

    ! Extract result.
    ! Since 'only' KK lons are evaluated (covering complete circle),
    ! select values cyclic for k>KK.
    ! Select first each group of LL elements
    do k = 0, Kout-1
      j = mod(k,KK)*LL + 1    ! 1, LL+1, 2*LL+1, ...
      llgrid(k+1) = real(F(j))*fac
    end do
    
    ! done
    deallocate( C )
    deallocate( F )
    
  end subroutine shgrid_Eval_Lons_Xm_1
  
  

  ! ==============================================
  
  
  ! compute associated legendre coeff for
  !
  !      d X           2  d X                         dX   dX
  !  ( -------- , (1-mu ) ---- )   =   r cos(theta) ( -- , -- )
  !    d lambda           d mu                        dx   dy
  !
  ! given associated legendre coeff of X(lambda,mu)
  !
  
  subroutine shgrid_Nabla( Xnm, NablaX )

    ! --- in/out -----------------------------
    
    type(TshGrid), intent(in)        ::  Xnm
    type(TshGrid), intent(inout)     ::  NablaX(2)
    
    ! --- local ------------------------
    
    integer   ::  m, n, T
    integer   ::  k
    
    ! --- begin ----------------------------

    ! extract triangular truncation
    T = Xnm%T
    
    ! allocate output arrays
    call Set( NablaX(1), T )
    call Set( NablaX(2), T )
    
    ! loop over all coeff    
    k = 0
    do m = 0, T
      do n = m, T
        k = k + 1

        ! d(ln ps)                                i m lon
        ! -------- = sum  {i m X(m,n)} P(mu;m,n) e
        !  d lon     m,n

        NablaX(1)%c(k) = (0.0,1.0) * m * Xnm%c(k)

        !      2  d(ln ps)       ~                 i m lon
        ! (1-mu ) -------- = sum X(m,n) P(mu;m,n) e
        !           d mu     m,n
        !
        !  ~      = (n+2) eps(m,n+1) X(m,n+1)                           , m<T, m=n
        !  X(m,n) = (n+2) eps(m,n+1) X(m,n+1) - (n-1) eps(m,n) X(m,n-1) , m<T, m<n<T
        !         =                           - (n-1) eps(m,n) X(m,n-1) , m<T,   n=T
        !         =                          0.0                        , m=T, n=T
        
        if ( m < T ) then
          if ( n == m ) then
            NablaX(2)%c(k) = (n+2) * epsilon(m,n+1) * Xnm%c(k+1)
          else if ( n > m .and. n < T ) then
            NablaX(2)%c(k) = (n+2) * epsilon(m,n+1) * Xnm%c(k+1) - &
                             (n-1) * epsilon(m,n)   * Xnm%c(k-1)
          else ! n == T
            NablaX(2)%c(k) = (n-1) * epsilon(m,n)   * Xnm%c(k-1)
          end if
        else
          NablaX(2)%c(k) = 0.0
        end if

      end do
    end do
    
  contains

    !             n*n - m*m  1/2
    ! epsilon = ( --------- )
    !             4*n*n - 1

    real function epsilon( m, n )

      ! --- in/out --------------------------------

      integer, intent(in) ::  m, n

      ! --- begin --------------------------------

      epsilon = sqrt((n*n-m*m)*1.0/(4*n*n-1.0))

    end function epsilon
    
  end subroutine shgrid_Nabla


  ! ==================================================================
  
  ! 
  ! Compute  (U,V) = (u,v) cos(lat)  from vorticity VO and divergence D.
  !
  ! Use relations:
  !
  !  U(m,n) = ( + delta(m,n) VO(m,n-1) + i sigma(m,n) D(m,n)  - delta(m,n+1) VO(m,n+1) ) R
  !
  !  V(m,n) = ( - delta(m,n) D(m,n-1)  + i sigma(m,n) VO(m,n) + delta(m,n+1) D(m,n+1)  ) R
  !
  ! where
  !
  !   delta(m,n) = - eps(m,n)/n     ,   eps(m,n) = sqrt( (n^2-m^2)/(4n^2-1) )
  !
  !   sigma(m,n) = - m / (n(n+1))
  !
  ! NOTE: In EMOS library, vor. and div. are assumed to be truncated at T-1;
  !   for consitency, this is used here too ...
  !   
  

  subroutine sh_vod2uv( VO, DI, U, V )
  
    use Binas, only : R => ae

    ! --- in/out -----------------------------

    type(TshGrid), intent(in)    ::  VO, DI
    type(TshGrid), intent(out)   ::  U, V

    ! --- local ------------------------------

    integer    ::  m, n, T
    integer    ::  k
    complex    ::  z

    ! --- begin ----------------------------

    ! Use EMOS library:
    !call jvod2uv( VO%c, DI%c, VO%T, U%c, V%c, VO%T )

    ! check arguments:    
    call Check( VO )
    T = VO%T
    call Check( DI, T )
    call Check( U, T )
    call Check( V, T )

    ! complex i for multiplication:
    z = (0.0, 1.0)

    ! index of coeff in array; initialise:
    k    = 0
    
    ! loop over all m:
    do m = 0, T-1

      ! (m,m)
      n = m
      k = k + 1
      if ( m == 0 ) then
        U%c(k) = (                      - delta(m,n+1)*VO%c(k+1) )*R
        V%c(k) = (                        delta(m,n+1)*DI%c(k+1) )*R
      else
        U%c(k) = ( z*sigma(m,n)*DI%c(k) - delta(m,n+1)*VO%c(k+1) )*R
        V%c(k) = ( z*sigma(m,n)*VO%c(k) + delta(m,n+1)*DI%c(k+1) )*R
      end if

      ! internal area:
      if ( m <= T-2 ) then
        do n = m+1, T-2
          k  = k + 1
          U%c(k) = (   delta(m,n)*VO%c(k-1) + z*sigma(m,n)*DI%c(k) - delta(m,n+1)*VO%c(k+1) )*R
          V%c(k) = ( - delta(m,n)*DI%c(k-1) + z*sigma(m,n)*VO%c(k) + delta(m,n+1)*DI%c(k+1) )*R
        end do

        ! (m,T-1)
        n = T-1
        k  = k + 1
        U%c(k) = (   delta(m,n)*VO%c(k-1) + z*sigma(m,n)*DI%c(k))*R
        V%c(k) = ( - delta(m,n)*DI%c(k-1) + z*sigma(m,n)*VO%c(k))*R
      end if

      ! (m,T)
      n = T
      k  = k + 1
      U%c(k) = (   delta(m,n)*VO%c(k-1) )*R
      V%c(k) = ( - delta(m,n)*DI%c(k-1) )*R

    end do

    ! (T,T)
    k  = k + 1
    U%c(k) = 0.0
    V%c(k) = 0.0

  contains

    real function delta( m, n )
      integer, intent(in)  :: m, n
      delta = - sqrt( (n*n-m*m)*1.0/(4.0*n*n-1.0) ) / (n*1.0)
    end function delta

    real function sigma( m, n ) 
      integer, intent(in)  :: m, n
      sigma = - m * 1.0/(n*(n+1.0))
    end function sigma

  end subroutine sh_vod2uv
  

  ! ***
  
  
  pure subroutine shi_vod2uv( shi_in, VO, DI, shi, U, V )
  
    use Binas, only : R => ae

    ! --- in/out -----------------------------

    type(TShGridInfo), intent(in)     ::  shi_in
    complex, intent(in)               ::  VO(shi_in%np), DI(shi_in%np)
    type(TShGridInfo), intent(in)     ::  shi
    complex, intent(out)              ::  U(shi%np), V(shi%np)

    ! --- local ------------------------------

    integer    ::  m, n
    integer    ::  kin, k
    complex    ::  z

    ! --- begin ----------------------------

    ! Use EMOS library:
    !call jvod2uv( VO%c, DI%c, VO%T, U%c, V%c, VO%T )

    ! check arguments:    
    !call Check( VO )
    !T = VO%T
    !call Check( DI, T )
    !call Check( U, T )
    !call Check( V, T )

    ! complex i for multiplication:
    z = (0.0, 1.0)

    ! index of coeff in array; initialise:
    kin = 0
    k   = 0
    
    ! loop over all m:
    do m = 0, shi%T-1

      ! (m,m)
      n = m
      kin = kin + 1
      k   = k   + 1
      if ( m == 0 ) then
        U(k) = (                      - delta(m,n+1)*VO(kin+1) )*R
        V(k) = (                        delta(m,n+1)*DI(kin+1) )*R
      else
        U(k) = ( z*sigma(m,n)*DI(kin) - delta(m,n+1)*VO(kin+1) )*R
        V(k) = ( z*sigma(m,n)*VO(kin) + delta(m,n+1)*DI(kin+1) )*R
      end if

      ! internal area:
      if ( m <= shi%T-2 ) then
        do n = m+1, shi%T-2
          kin = kin + 1
          k   = k   + 1
          U(k) = (   delta(m,n)*VO(kin-1) + z*sigma(m,n)*DI(kin) - delta(m,n+1)*VO(kin+1) )*R
          V(k) = ( - delta(m,n)*DI(kin-1) + z*sigma(m,n)*VO(kin) + delta(m,n+1)*DI(kin+1) )*R
        end do

        ! (m,T-1)
        n = shi%T-1
        kin = kin + 1
        k   = k   + 1
        U(k) = (   delta(m,n)*VO(kin-1) + z*sigma(m,n)*DI(kin))*R
        V(k) = ( - delta(m,n)*DI(kin-1) + z*sigma(m,n)*VO(kin))*R
      end if

      ! (m,T)
      n = shi%T
      kin = kin + 1
      k   = k   + 1
      U(k) = (   delta(m,n)*VO(kin-1) )*R
      V(k) = ( - delta(m,n)*DI(kin-1) )*R

      kin = kin + shi_in%T - shi%T
      
    end do

    ! (T,T)
    k  = k + 1
    U(k) = (0.0,0.0)
    V(k) = (0.0,0.0)

  contains

    pure real function delta( m, n )
      integer, intent(in)  :: m, n
      delta = - sqrt( (n*n-m*m)*1.0/(4.0*n*n-1.0) ) / (n*1.0)
    end function delta

    pure real function sigma( m, n ) 
      integer, intent(in)  :: m, n
      sigma = - m * 1.0/(n*(n+1.0))
    end function sigma

  end subroutine shi_vod2uv
  
  
  
  
end module grid_type_sh
  
