module tmflex

    use dims,        only : im, jm, isr, ier, jsr, jer, nregions
    use go,          only : readrc
    use global_data, only : rcf

    implicit none

    public :: initTMFlex, resetBackgroundMix, doneTMFlex
    private

    type d2struct
        logical, dimension(:,:), allocatable :: mask
    end type d2struct

    type(d2struct), dimension(nregions) :: isInDOI ! is the point in or out the domain of interest

    logical :: tmflex_apply

contains
    subroutine initTMFlex 
        integer :: status
        call readrc(rcf, 'tmflex.compute.backgrounds', tmflex_apply, status, default=.false.)
        if (tmflex_apply) call tmflex_init 
    end subroutine initTMFlex

    subroutine tmflex_init 
        use dims, only : xref, yref, dx, dy, xbeg, ybeg
        integer :: region
        integer :: ilon, ilat, status
        real    :: lon0, lon1, lat0, lat1 ! domain boundaries
        real    :: latb, late, lonb, lone ! grid-cell boundaries
        ! Define the domain outside which the mixing ratios must be reset every time-step
        call readrc(rcf, 'tmflex.lat0', lat0, status)
        call readrc(rcf, 'tmflex.lat1', lat1, status)
        call readrc(rcf, 'tmflex.lon0', lon0, status)
        call readrc(rcf, 'tmflex.lon1', lon1, status)
        do region = 1, nregions
            allocate(isInDOI(region)%mask(im(region),jm(region)))
            isInDOI(region)%mask = .False.
            do ilat = jsr(region), jer(region)
                latb = ybeg(region) + (ilat-1)*dy/yref(region)
                late = latb+dy/yref(region)
                do ilon = isr(region), ier(region)
                    lonb = xbeg(region) + (ilon-1)*dx/xref(region)
                    lone = lonb + dx/xref(region)
                    if (lonb >= lon0 .and. lone <= lon1 .and. latb >= lat0 .and. late <= lat1) then
                        ! isInDOI is set to True if the grid-cell is entirely contained in the footprint domain
                        ! (it is easier to cut the footprint than the TM5 grid-cell!)
                        ! It is however best to plan ahead and make sure that the domains are conformable
                        isInDOI(region)%mask(ilon, ilat) = .true.
                    endif
                enddo
            enddo
        enddo
    end subroutine tmflex_init 

    subroutine resetBackgroundMix(region, status)
        use global_data, only : mass_dat, region_dat
        ! reset the mixing ratios ouside the domain of interest every time-step
        ! This is done as it would be done in a chemistry module: a loss rate of 1/timestep is used
        integer, intent(in)  :: region
        integer, intent(out) :: status
        integer     :: ilon, ilat
        status = 1
        if (tmflex_apply) then
            do ilat = jsr(region), jer(region)
                do ilon = isr(region), ier(region)
#ifdef with_zoom
                    if (region_dat(region)%zoomed(ilon, ilat) /= region) cycle
#endif
                    if (.not. isInDOI(region)%mask(ilon, ilat)) then
                        mass_dat(region)%rm_t(ilon, ilat, :, 2) = 0.
#ifdef slopes
                        mass_dat(region)%rxm_t(ilon,ilat, :, 2) = 0.
                        mass_dat(region)%rym_t(ilon,ilat, :, 2) = 0.
                        mass_dat(region)%rzm_t(ilon,ilat, :, 2) = 0.
#endif
                    endif
                enddo
            enddo
        endif
        status = 0
    end subroutine resetBackgroundMix

    subroutine doneTMFlex
        integer :: region
        if (tmflex_apply) then
            do region = 1, nregions
                deallocate(isInDOI(region)%mask)
            enddo
        endif
    end subroutine doneTMFlex
end module tmflex
