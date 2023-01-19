module interpol_grid

    implicit none

    public :: redis_partialcol_altitude

contains

    subroutine redis_partialcol_altitude(old_levels,new_levels,old_grid,old_prof,new_grid,new_prof,matr)
        !
        ! old_levels = Number of levels in the TM5 grid (25 or 34)
        ! old_grid is the altitude grid of TM5, given on old_levels+1 altitudes
        ! old_prof is the profile of partial column CO in #molecules cm-2 on old_levels layers
        !
        ! new_levels = Number of levels of the IASI altitude grid (<=19 levels)    
        ! new_grid  = IASI altitude grid on old_levels+1 altitudes
        ! new_prof  = new profile of partial column CO in #molecules cm-2 on new_levels layers
        ! (such that sum(new_prof)==sum(old_prof))
        !
        ! matr is the matrix of dimensions (new_levels x old_levels) such that new_prof = matr * old_prof
        !
        implicit none
        
        ! in/out
        integer, intent(in) :: old_levels, new_levels
        double precision, intent(in)    :: old_grid(old_levels+1), old_prof(old_levels), new_grid(new_levels+1)
        double precision, intent(out)   :: new_prof(new_levels), matr(new_levels, old_levels)    
        
        ! local
        integer :: i, j, last
        double precision    :: frac

        ! debug
        !write(*,*) 'old_grid = ', old_grid
        !write(*,*) 'old_prof = ', old_prof
        !write(*,*) 'new_grid = ', new_grid
        
        ! start
        new_prof = 0.
        matr     = 0.
        last     = 1
        
        last = 2
        do i = 1, new_levels
            
            ! Find boxes needed to fill the new grid            
            do j = last, old_levels+1
                if ( old_grid(j) <= new_grid(i+1) .and. old_grid(j) >= new_grid(i) ) then
                    frac = (old_grid(j)-new_grid(i))/(old_grid(j)-old_grid(j-1))
                    if ( frac > 1 ) then
                        frac = 1.
                    end if
                    matr(i,j-1) = frac
                else
                    if ( old_grid(j) >= new_grid(i) ) then
                        frac = (new_grid(i+1)-max(old_grid(j-1),new_grid(i)))/(old_grid(j)-old_grid(j-1))
                        matr(i,j-1) = frac
                        last = j
                        exit
                    end if
                end if
            end do
            !print*, ''
            
        end do
        
        new_prof = matmul(matr,old_prof)
        
        !print*,'Mass conservation? ',sum(new_prof),sum(old_prof),sum(new_prof)-sum(old_prof)
            
    end subroutine redis_partialcol_altitude


end module interpol_grid
