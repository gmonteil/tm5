!
! numerical tools : random numbers
!

module Num_Random

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  gasdev
  
  
contains


  ! -------------------------------------------
  ! Returns a normally distributed deviate with zero mean and unit variance
  ! From Numerical Recipes in Fortran, 2nd edition, p.280
  !
  ! Copied from Var4D code.
  ! -------------------------------------------

  real function gasdev()

    integer, save  ::  iset = 0
    real, save     ::  gset
    real           ::  fac,rsq,v1,v2,rnd1,rnd2

    ! generate new values ?
    if ( iset == 0 ) then
      do 
        ! draw two uniform distributed random numbers:
        call random_number(rnd1)  ! [0,1]
        call random_number(rnd2)  ! [0,1]
        ! convert range:
        v1 = 2.0 * rnd1 - 1.0   ! [-1,1]
        v2 = 2.0 * rnd2 - 1.0   ! [-1,1]
        ! squared radius in 2D plane:
        rsq = v1**2 + v2**2   ! [0,2]
        ! ok if in (0,1)
        if ( (rsq > 0.0) .and. (rsq < 1.0) ) exit
      end do
      ! magic ...
      fac = sqrt( -2.0*log(rsq)/rsq )
      ! two Gaussian distributed random numbers
      gset   = v1*fac   ! store for later use
      gasdev = v2*fac   ! return value
      ! reset flag:
      iset = 1
    else
      ! now return the other one:
      gasdev = gset
      ! reset flag:
      iset = 0
    end if
    
  end function gasdev



end module Num_Random

