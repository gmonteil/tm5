module phys_cloudcover

  implicit none
  
  ! --- in/out ---------------------------
  
  private
  
  public  ::  cf_overhead
  
  
contains


  ! ==============================================================
  
  ! ---------------------------------------------------
  !
  ! Calculate total cloud fraction overhead the base of each layer
  ! based on random/maximum overlap assumptions
  ! Based on code provided by Rob van Dorland
  !   Peter van Velthoven - 22 November 2002
  ! 
  ! Optional arguments
  !   Arjo Segers - 25 november 2002
  !
  ! ----------------------------
  !
  ! input: 
  !   nlev  : number of vertical levels
  !   yclfr : cloud fraction (cc) per cell (0-1)
  !
  ! output:
  !   wccro: overhead cloud fraction
  !
  ! optional arguments:
  !   scheme='ecmwf'  : 'ecmwf' -> iovln=1
  !                     'other' -> iovln=0
  !   eps=1.0e-4      : cltres
  !
  ! parameters:
  !   iovln : switch
  !     1 = ecmwf (maximum random overlap assumption) scheme
  !     0 = another scheme 
  !   cltres : threshold (minimum) cloud fraction used
  !            for numerical stability (division by zero
  !            and to eliminate small unrealistic cloud fractions
  !
  ! Notes: 
  ! - Index=1 of arrays (yclfr) corresponds to model top
  ! - The clouds are supposed to be distributed homogeneously
  !   in the vertical in each layer.
  !
  ! ----------------------------------------------------

  subroutine cf_overhead( nlev, yclfr, wccro, scheme, eps )

    ! --- in/out ------------------------------

    integer, intent(in)   ::  nlev
    real, intent(in)      ::  yclfr(nlev)
    real, intent(out)     ::  wccro(nlev)
    
    character(len=*), intent(in), optional  ::  scheme
    real, intent(in), optional              ::  eps

    ! --- local ------------------------------

    real      ::   clfr0, clfr1, clfr2, ctver
    real      ::   zclear, zcloud
    integer   ::   jk

    ! --- settings -----------------------------

    !integer   ::  iovln = 0
    integer    ::  iovln = 1      ! ecmwf; maximum random overlap

    real       ::  cltres = 1.0e-4

    ! --- begin -----------------------------
    
    if ( present(scheme) ) then
      select case ( scheme )
        case ( 'ecmwf' )
          iovln = 1
        case ( 'other' )
          iovln = 0
        case default
          print *, 'Unsupported scheme "'//scheme//'".'
          stop 'FATAL BUG IN cf_overhead'
      end select
    end if
    
    if ( present(eps) ) cltres = eps
    

    select case ( iovln )

      !-----------------------------------------
      ! scheme 0: maximum overlap unless there's a 
      ! clear sky layer in between?
      !-----------------------------------------

      case ( 0 )

        clfr0 = 0.0
        clfr2 = 0.0
        ctver = 1.0
        do jk = 1, nlev
          clfr1 = yclfr(jk)
          if ( clfr1 < cltres ) then
            !----------------
            ! random overlap
            !----------------
            ctver = ctver * ( 1.0 - clfr2 )
            clfr2 = 0.0
          else
            if ( clfr0 < cltres ) then
              clfr2 = clfr1
            else
              !----------------
              ! maximum overlap
              !----------------
              clfr2 = max( clfr1,clfr2 )
            end if
          end if
          clfr0 = clfr1
          wccro(jk) = 1.0 - ctver * ( 1.0 - clfr2 )
        end do
        !ctver=ctver*(1.-clfr2)
        !wccro=1.-ctver

      !-----------------------------------------
      !     ecmwf scheme
      !-----------------------------------------
      
      case ( 1 )

        zclear = 1.0
        zcloud = 0.0
        do jk = 1, nlev
          zclear = zclear*(1.0-max(yclfr(jk),zcloud))/(1.0-min(zcloud,1.0-cltres))
          zcloud = yclfr(jk)
          wccro(jk) = 1.0 - zclear
        end do
        
      !-----------------------------------------
      ! error ...
      !-----------------------------------------
      
      case default
      
        print *, 'unknown switch',IOVLN
        stop 'FATAL BUG IN cf_overhead'
        
    end select

  end subroutine cf_overhead


  ! ***
  

  ! ---------------------------------------------------
  ! Calculate cloud fraction in joint layers based on maximum overlap
  ! which holds only for neighbouring layers !
  ! ----------------------------
  ! Peter van Velthoven - 22 November 2002
  ! ---------------------------------------------------

  subroutine join_cc_layers( YCLFR1, YCLFR2, YCLFRJOINED )

    ! --- in/out ----------------------------
    
    real, intent(in)    ::  YCLFR1, YCLFR2
    real, intent(out)   ::  YCLFRJOINED

    ! --- begin ----------------------------
    
    YCLFRJOINED  = max( YCLFR1, YCLFR2 )

  end subroutine join_cc_layers


end module phys_cloudcover
