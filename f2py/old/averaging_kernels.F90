module oco2

implicit none

private

public :: applyaveragingkernel

integer, parameter  :: mismatch_params = 4

contains

pure function inv_searchSorted(in_array, values) result (index_array)

    ! searchSorted = return the indices in array in_array such that
    ! in_array(idx(i)) >= value(i) > in_array(idx(i)+1). If value(i) is outside
    ! in_array then either lbound(in_array)-1 or ubound(in_array) is returned.

    double precision, intent(in)        :: in_array(:), values(:)
    integer, dimension(size(values))    :: index_array
    integer                             :: ub, lb, i, j, idx

    lb = lbound(in_array,1)
    ub = ubound(in_array,1)
    do j=lbound(values,1),ubound(values,1)
        idx = lb-1
        do i=lb,ub
            if (values(j) .ge. in_array(i)) exit
            idx = i
        end do
        index_array(j+1-lbound(values,1)) = idx
    end do

end function inv_searchSorted

pure function diff(in_array)

    double precision, intent(in), dimension(:)  :: in_array
    double precision                            :: diff(size(in_array)-1)
    integer                                     :: i, lb

    lb = lbound(in_array,1)
    do i=1,size(diff)
        diff(i) = in_array(i+lb)-in_array(i+lb-1)
    end do

end function diff

pure function inv_diff(in_array)

    double precision, intent(in), dimension(:)      :: in_array
    double precision, dimension(size(in_array)-1)   :: inv_diff
    integer                                         :: i, lb

    lb = lbound(in_array,1)
    do i=1,size(inv_diff)
        inv_diff(i) = in_array(i+lb-1)-in_array(i+lb)
    end do

end function inv_diff

pure subroutine inv_rebin(in_x_grid, in_y_vals, out_x_grid, out_y_vals, distrib_matrix, input_OK)

    ! out_x_grid is the output grid except the end-points, since
    ! the end-points of out_x_grid are the same as those of in_x_grid
    ! in_x_grid and out_x_grid are in descending order, and distrib_matrix is such that
    !     out_y_vals = (distrib_matrix) (in_y_vals)
    ! input/output
    double precision, intent(in)    :: in_x_grid(:), in_y_vals(:), out_x_grid(:)
    double precision, intent(out)   :: out_y_vals(size(out_x_grid)+1), distrib_matrix(size(out_x_grid)+1, size(in_x_grid)-1)
    logical, intent(out)            :: input_OK
    ! local
    integer                         :: in_levels, out_levels, i, j
    integer, allocatable            :: breakpoints(:), ext_breakpoints(:)
    double precision, allocatable   :: cumul(:), in_x(:), out_x(:)
    double precision                :: mass_error, lb, ub

    in_levels = size(in_x_grid)
    out_levels = size(out_x_grid)

    ! first some error checking
    input_OK = .true.
    ! check whether x_grid's are sorted in descending order and unique
    do i=lbound(in_x_grid,1),ubound(in_x_grid,1)-1
        if (in_x_grid(i) <= in_x_grid(i+1)) then
            input_OK = .false.
        end if
    end do
    do i=lbound(out_x_grid,1),ubound(out_x_grid,1)-1
        if (out_x_grid(i) <= out_x_grid(i+1)) then
            input_OK = .false.
        end if
    end do
    ! check whether the output grid is contained inside the input grid
    if (out_x_grid(lbound(out_x_grid,1)) > in_x_grid(lbound(in_x_grid,1))) then
        input_OK = .false.
    end if
    if (out_x_grid(ubound(out_x_grid,1)) < in_x_grid(ubound(in_x_grid,1))) then
        input_OK = .false.
    end if
    if ( .not. input_OK ) return
    ! now that we know the input is OK, proceed with rebinning :-)

    ! since this is a pure subroutine, use only allocatable arrays
    allocate(breakpoints(out_levels), ext_breakpoints(out_levels+2), cumul(out_levels+2), in_x(in_levels), out_x(out_levels+2))

    ! copy over data, only really needed for out_x
    in_x(:) = in_x_grid(:)
    out_x(2:size(out_x)-1) = out_x_grid(:)
    out_x(1) = in_x(1)
    out_x(size(out_x)) = in_x(size(in_x))

    breakpoints = inv_searchSorted(in_x, out_x_grid)
    ext_breakpoints(2:out_levels+1) = breakpoints(:)
    ext_breakpoints(1) = 1
    ext_breakpoints(out_levels+2) = in_levels-1

    cumul(1) = 0.0
    do i=1,size(breakpoints)
        cumul(i+1) = sum(in_y_vals(:breakpoints(i)-1) * inv_diff(in_x(:breakpoints(i))))
        cumul(i+1) = cumul(i+1) + in_y_vals(breakpoints(i)) * (in_x(breakpoints(i)) - out_x(i+1))
    end do
    cumul(size(cumul)) = sum(in_y_vals * inv_diff(in_x))
    out_y_vals = diff(cumul)/inv_diff(out_x)

    ! now calculate the distribution matrix
    distrib_matrix = 0.0 ! set all out_levels+1 x in_levels-1 elements to 0
    ! out_y_vals(i) has contributions from in_y_vals(ext_breakpoints(i)) to in_y_vals(ext_breakpoints(i+1))
    ! in_y_vals(j) lies between in_x(j) and in_x(j+1)
    ! hence, out_y_vals(i), which lines between out_x(i) and out_x(i+1), will need to examine levels in_x(ext_breakpoints(i)) and in_x(ext_breakpoints(i+1)+1)
    do i=1,out_levels+1
        do j=1,in_levels-1
            ! output level i is between out_x(i) and out_x(i+1)
            ! input level j is between in_x(j) and in_x(j+1)
            ! do something only if there is some overlap
            if ((in_x(j+1) .ge. out_x(i)) .or. (out_x(i+1) .ge. in_x(j))) then
                continue
            else
                ub = min(out_x(i),in_x(j))
                lb = max(in_x(j+1),out_x(i+1))
                distrib_matrix(i,j) = (ub-lb)/(out_x(i)-out_x(i+1))
            end if
        end do
    end do

    deallocate(breakpoints, ext_breakpoints, cumul, in_x, out_x)

end subroutine inv_rebin

subroutine applyaveragingkernel(n_obs, n_lev_mod, n_lev_oco, model_profiles, model_column_err, model_psurf, &
    oco_totalcol, oco_totalcol_err, AT, BT, oco_pres_levels, avg_ker, prior_profile, add_model_err, &
    mismatches, J_obs, departures, valid_deps)

    ! I/O (input)
    integer, intent(in)             :: n_obs, n_lev_mod, n_lev_oco
    double precision, intent(in)    :: model_profiles(n_obs, n_lev_mod), model_column_err(n_obs), model_psurf(n_obs)
    double precision, intent(in)    :: oco_totalcol(n_obs), oco_totalcol_err(n_obs), AT(n_lev_mod+1), BT(n_lev_mod+1)
    double precision, intent(in)    :: oco_pres_levels(n_obs, n_lev_oco+1), avg_ker(n_obs, n_lev_oco)
    double precision, intent(in)    :: prior_profile(n_obs, n_lev_oco)
    integer, intent(in)             :: add_model_err
    ! I/O (output)
    double precision, intent(out)   :: J_obs, departures(n_obs, n_lev_mod)
    double precision, intent(out)   :: mismatches(n_obs, mismatch_params)
    integer(1), intent(out)         :: valid_deps(n_obs)

    ! local variables
    integer                         :: i
    double precision                :: m, c, c1, c2, model_tc, tot_dtc, dep
    logical                         :: ip_ok
    ! Thread-specific variables
    double precision, allocatable   :: departures_thread(:,:), mismatches_thread(:,:), tm5_pres_levels(:)
    double precision, allocatable   :: squeezed_model_profile(:), redistrib_matrix(:,:), pres_weights(:), W(:), H(:)
    double precision                :: J_thread
    logical, allocatable            :: valid_thread(:), use_thisthread(:)

    mismatches = 0.0
    J_obs = 0.0
    departures = 0.0
    ! assume all observations are invalid, and later correct for the valid ones
    valid_deps = 0

    !$omp parallel private(departures_thread, mismatches_thread, J_thread, valid_thread, use_thisthread) &
    !$omp private(i, tm5_pres_levels, m, c, c1, c2, squeezed_model_profile, redistrib_matrix) &
    !$omp private(ip_ok, pres_weights, W, H, model_tc, tot_dtc, dep)
    allocate(departures_thread(n_obs, n_lev_mod), mismatches_thread(n_obs, mismatch_params))
    allocate(valid_thread(n_obs), use_thisthread(n_obs))
    allocate(tm5_pres_levels(n_lev_mod+1),  pres_weights(n_lev_oco))
    allocate(squeezed_model_profile(n_lev_oco), redistrib_matrix(n_lev_oco, n_lev_mod), W(n_lev_oco), H(n_lev_mod))

    departures_thread = 0.0
    mismatches_thread = 0.0
    valid_thread = .false.
    use_thisthread = .false.
    J_thread = 0.0

    !$omp do schedule(guided)
    do i = 1, n_obs
        use_thisthread(i) = .true. ! computation for this obs done on this thread

        tm5_pres_levels = AT + BT * model_psurf(i)
        ! squeeze the model profile into the instrumental grid
        m = (oco_pres_levels(i,1) - oco_pres_levels(i,n_lev_oco+1))/(tm5_pres_levels(1) - tm5_pres_levels(n_lev_mod+1))
        c1 = oco_pres_levels(i,1) - m*tm5_pres_levels(1)
        c2 = oco_pres_levels(i,n_lev_oco+1) - m*tm5_pres_levels(n_lev_mod+1)
        c = 0.5*(c1+c2)
        tm5_pres_levels = m * tm5_pres_levels + c

        ! Now re-bin the model profile onto the retrieval levels
        call inv_rebin(tm5_pres_levels, model_profiles(i,:), oco_pres_levels(i,2:n_lev_oco), &
            squeezed_model_profile, redistrib_matrix, ip_ok)

        ! if ip_ok is .false., there was a problem with the matrices, so do not use this observation, i.e.,
        ! keep valid_thread for this obs to .false.
        if (ip_ok) then
            ! For OCO-2, Chris provides normalized averaging kernels. Therefore, given a (squeezed) modeled profile, the modeled
            ! column average X^mod is
            ! X^mod = \sum_i h_i a_i x^mod_i + \sum_i h_i (1-a_i) x^pri_i
            ! where h_i = \del_i P/Psurf is the pressure weighting function, and a_i is the normalized averaging kernel
            ! It's expedient to define a vector W such that W_i = a_i h_i
            pres_weights = inv_diff(oco_pres_levels(i,:))/oco_pres_levels(i,1)
            W = pres_weights * avg_ker(i,:)
            model_tc = dot_product(W, squeezed_model_profile) + dot_product(pres_weights-W, prior_profile(i,:))
            if (add_model_err == 0) then
                tot_dtc = oco_totalcol_err(i)**2
            else
                tot_dtc = model_column_err(i)**2 + oco_totalcol_err(i)**2
            end if

            mismatches_thread(i,:) = (/ model_tc, oco_totalcol(i), tot_dtc, dble(i-1) /)
            dep = (model_tc - oco_totalcol(i))/tot_dtc
            J_thread = J_thread + 0.5 * dep * (model_tc - oco_totalcol(i))
            H = matmul(W, redistrib_matrix)
            departures_thread(i,:) = H * dep
            valid_thread(i) = .true.
        else
            mismatches_thread(i,:) = (/ 0.0d0, 0.0d0, 1.0d0, dble(i-1) /)
            departures_thread(i,:) = 0.0
        end if
    end do
    !$omp end do

    ! put all the thread results together now, executing one thread at a time
    !$omp critical
    do i = 1, n_obs
        if (use_thisthread(i)) then
            mismatches(i,:) = mismatches_thread(i,:)
            departures(i,:) = departures_thread(i,:)
            if (valid_thread(i)) valid_deps(i) = 1
        end if
    end do
    J_obs = J_obs + J_thread
    !$omp end critical
    deallocate(departures_thread, mismatches_thread, valid_thread, use_thisthread)
    deallocate(tm5_pres_levels,  pres_weights, squeezed_model_profile, redistrib_matrix, W, H)
    !$omp end parallel

end subroutine applyaveragingkernel

end module oco2

module carbonsat

implicit none

private

public :: takecolumn_createobs, applyaveragingkernel, apply_ak_perturb

interface polyval
    module procedure polyval_arr
    module procedure polyval_scal
end interface polyval

real                :: measured_totalcol = 0.0
real                :: measured_totalcol_error = 0.0
integer, parameter  :: mismatch_params = 5

contains

pure function polyval_arr(poly_coeffs, x_array) result(res_array)
    ! If poly_coeffs = [p[1], p[2], ... p[N]] and an element of res_array is x, then
    ! return p[1]*x**(N-1) + p[2]*x**(N-2) + ... + p[N]
    real(8), intent(in) :: poly_coeffs(:), x_array(:)
    real(8)             :: res_array(size(x_array))
    integer             :: i, N ! N is degree of polynomial

    N = size(poly_coeffs)
    res_array = poly_coeffs(1)
    do i = 2, N
        res_array = x_array * res_array + poly_coeffs(i)
    end do

end function polyval_arr

pure function polyval_scal(poly_coeffs, x) result(res)
    ! Same as polyval_arr, but for a single scalar x
    real(8), intent(in) :: poly_coeffs(:), x
    real(8)             :: res
    integer             :: i, N ! N is degree of polynomial

    N = size(poly_coeffs)
    res = poly_coeffs(1)
    do i = 2, N
        res = x * res + poly_coeffs(i)
    end do

end function polyval_scal

subroutine apply_ak_perturb(n_obs, n_lev, n_poly, n_total, &
    model_profiles, err_model_profiles, model_psurf, meas_tc, meas_dtc, perturb_meas, indices, AT, BT, ak_poly, prior_mix, &
    mismatches, J_obs, departures)
    ! For Chevallier's Monte Carlo method of estimating posterior errors, we need to perturb the measurements by errors consistent with
    ! MDM matrix. Since the MDM matrix is diagonal, this can be done easily, by taking a vector of normally distributed variables,
    ! multiplying them by the total error, and adding them to the measured total column. This vector os perturbations has to be the same
    ! throughout an inversion, but that is something that will be ensured by python.
    ! --------------- I/O ---------------
    integer, intent(in)     :: n_obs, n_lev, n_poly, n_total
    real(8), intent(in)     :: model_profiles(n_obs, n_lev), err_model_profiles(n_obs, n_lev), model_psurf(n_obs), prior_mix
    real(8), intent(in)     :: AT(n_lev+1), BT(n_lev+1), ak_poly(n_total, n_poly), meas_tc(n_obs), meas_dtc(n_obs), perturb_meas(n_obs)
    integer, intent(in)     :: indices(n_obs)
    ! model_profiles are the profiles of TM5 at sounding locations, output by user_output_satellite, and err_model_profiles are their errors
    ! model_psurf is the TM5 surface pressure at a sounding, prior_mix is the (uniform) prior mixing ratio for the total column retrieval
    ! AT and BT are coefficients to reconstruct the TM5 pressure levels from model_psurf, while ak_poly is the polynomial representation of the averaging kernels
    ! meas_tc and meas_dtc are the measured total column CO2 and their errors, while perturb_meas is a vector of normally distributed errors to add to meas_tc before calculating the departures
    ! indices is the ordering of the samples, which can get messed up during multi-threaded calculation
    real(8), intent(out)    :: J_obs, departures(n_obs, n_lev)
    real(8), intent(out)    :: mismatches(n_obs,mismatch_params)
    ! J_obs is the total cost function from all the observations, departures contains the departure profile (adjoint forcing) per sounding,
    ! and mismatches is a vector of [index, modelled column, measured column, total error, ordering] per sounding
    ! --------------- local variables ---------------
    integer                 :: i, idx, n_eval
    real(8), allocatable    :: departures_thread(:,:)
    real(8), allocatable    :: mismatches_thread(:,:)
    logical, allocatable    :: use_thisthread(:)
    real(8)                 :: J_thread, model_pres_levels(n_lev+1), mid_pres_levels(n_lev), col_ak(n_lev), prior_profile(n_lev)
    real(8)                 :: prior_tc, deps, pres_weights(n_lev), ret_profile(n_lev), H(n_lev)
    real(8)                 :: model_tc, tot_dtc, mod_dtc

    mismatches = 0.0
    J_obs = 0.0
    departures = 0.0

    prior_profile = prior_mix
    prior_tc = prior_mix

    !$omp parallel private(i, mismatches_thread, departures_thread, use_thisthread, J_thread) &
    !$omp private(idx, model_pres_levels, mid_pres_levels, col_ak, model_tc, mod_dtc, tot_dtc) &
    !$omp private(deps, n_eval, ret_profile, H, pres_weights)
    J_thread = 0.0
    allocate(mismatches_thread(n_obs, mismatch_params))
    allocate(departures_thread(n_obs, n_lev))
    allocate(use_thisthread(n_obs))
    mismatches_thread = 0.0
    departures_thread = 0.0
    use_thisthread = .false.
    !$omp do schedule(dynamic, 1000)
    do i = 1, n_obs
        idx = indices(i)
        model_pres_levels = AT + BT*model_psurf(i)
        pres_weights = (model_pres_levels(1:n_lev)-model_pres_levels(2:n_lev+1))/model_psurf(i)
        mid_pres_levels = 0.5*(model_pres_levels(1:n_lev) + model_pres_levels(2:n_lev+1))/model_psurf(i)
        col_ak = polyval(ak_poly(idx+1,:), mid_pres_levels) ! the indices are PYTHON indices, so their lowest value is zero, for example
        ret_profile = col_ak*(model_profiles(i,:)-prior_profile) + prior_profile
        model_tc = sum(ret_profile*pres_weights)
        H = col_ak * pres_weights
        mod_dtc = dot_product(err_model_profiles(i,:)**2, H**2)
        ! Remember, meas_dtc is the standard deviation, not the variance, so we need to square it
        tot_dtc = mod_dtc + meas_dtc(i)**2
        ! perturb_meas(i) has to be multiplied by sqrt(tot_dtc) and added to meas_tc(i)
        mismatches_thread(i,:) = (/ dble(idx), model_tc, dble(meas_tc(i)) + perturb_meas(i)*sqrt(tot_dtc), tot_dtc, dble(i-1) /)
        deps = (model_tc - meas_tc(i) - perturb_meas(i)*sqrt(tot_dtc))/tot_dtc
        J_thread = J_thread + 0.5*deps*(model_tc-meas_tc(i)-perturb_meas(i)*sqrt(tot_dtc))
        departures_thread(i,:) = H * deps
        use_thisthread(i) = .true.
    end do
    !$omp end do
    ! put all the thread results together now
    !$omp critical
    n_eval = 0
    do i = 1, n_obs
        if (use_thisthread(i)) then
            mismatches(i,:) = mismatches_thread(i,:)
            departures(i,:) = departures_thread(i,:)
            n_eval = n_eval + 1
        end if
    end do
    J_obs = J_obs + J_thread
!    write(*,'("Thread ", i2, " evaluated ", i8, " soundings")') omp_get_thread_num()+1, n_eval
    !$omp end critical
    deallocate(mismatches_thread, departures_thread, use_thisthread)
    !$omp end parallel

end subroutine apply_ak_perturb


subroutine applyaveragingkernel(n_obs, n_lev, n_poly, n_total, model_profiles, err_model_profiles, model_psurf, meas_tc, meas_dtc, indices, AT, BT, ak_poly, prior_mix, mismatches, J_obs, departures)

    use omp_lib
    ! I/O
    integer, intent(in)     :: n_obs, n_lev, n_poly, n_total
    real(8), intent(in)     :: model_profiles(n_obs, n_lev), err_model_profiles(n_obs, n_lev), model_psurf(n_obs), prior_mix
    real(8), intent(in)     :: AT(n_lev+1), BT(n_lev+1), ak_poly(n_total, n_poly), meas_tc(n_obs), meas_dtc(n_obs)
    integer, intent(in)     :: indices(n_obs)
    real(8), intent(out)    :: J_obs, departures(n_obs, n_lev)
    real(8), intent(out)    :: mismatches(n_obs,mismatch_params)
    ! local variables
    integer                 :: i, idx, n_eval
    real(8), allocatable    :: departures_thread(:,:)
    real(8), allocatable    :: mismatches_thread(:,:)
    logical, allocatable    :: use_thisthread(:)
    real(8)                 :: J_thread, model_pres_levels(n_lev+1), mid_pres_levels(n_lev), col_ak(n_lev), prior_profile(n_lev)
    real(8)                 :: prior_tc, deps, pres_weights(n_lev), ret_profile(n_lev), H(n_lev)
    real(8)                 :: model_tc, tot_dtc, mod_dtc

    mismatches = 0.0
    J_obs = 0.0
    departures = 0.0

    prior_profile = prior_mix
    prior_tc = prior_mix

    !$omp parallel private(i, mismatches_thread, departures_thread, use_thisthread, J_thread) &
    !$omp private(idx, model_pres_levels, mid_pres_levels, col_ak, model_tc, mod_dtc, tot_dtc) &
    !$omp private(deps, n_eval, ret_profile, H, pres_weights)
    J_thread = 0.0
    allocate(mismatches_thread(n_obs, mismatch_params))
    allocate(departures_thread(n_obs, n_lev))
    allocate(use_thisthread(n_obs))
    mismatches_thread = 0.0
    departures_thread = 0.0
    use_thisthread = .false.
    !$omp do schedule(dynamic, 1000)
    do i = 1, n_obs
        idx = indices(i)
        model_pres_levels = AT + BT*model_psurf(i)
        pres_weights = (model_pres_levels(1:n_lev)-model_pres_levels(2:n_lev+1))/model_psurf(i)
        mid_pres_levels = 0.5*(model_pres_levels(1:n_lev) + model_pres_levels(2:n_lev+1))/model_psurf(i)
        col_ak = polyval(ak_poly(idx+1,:), mid_pres_levels) ! the indices are PYTHON indices, so their lowest value is zero, for example
        ret_profile = col_ak*(model_profiles(i,:)-prior_profile) + prior_profile
        model_tc = sum(ret_profile*pres_weights)
        H = col_ak * pres_weights
        mod_dtc = dot_product(err_model_profiles(i,:)**2, H**2)
        ! Remember, meas_dtc is the standard deviation, not the variance, so we need to square it
        tot_dtc = mod_dtc + meas_dtc(i)**2
        mismatches_thread(i,:) = (/ dble(idx), model_tc, dble(meas_tc(i)), tot_dtc, dble(i-1) /)
        deps = (model_tc - meas_tc(i))/tot_dtc
        J_thread = J_thread + 0.5*deps*(model_tc-meas_tc(i))
        departures_thread(i,:) = H * deps
        use_thisthread(i) = .true.
    end do
    !$omp end do
    ! put all the thread results together now
    !$omp critical
    n_eval = 0
    do i = 1, n_obs
        if (use_thisthread(i)) then
            mismatches(i,:) = mismatches_thread(i,:)
            departures(i,:) = departures_thread(i,:)
            n_eval = n_eval + 1
        end if
    end do
    J_obs = J_obs + J_thread
!    write(*,'("Thread ", i2, " evaluated ", i8, " soundings")') omp_get_thread_num()+1, n_eval
    !$omp end critical
    deallocate(mismatches_thread, departures_thread, use_thisthread)
    !$omp end parallel

end subroutine applyaveragingkernel

subroutine takecolumn_createobs(n_obs, n_lev, n_poly, n_total, model_profiles, err_model_profiles, model_psurf, indices, AT, BT, ak_poly, prior_mix, mismatches, J_obs, departures)

    use omp_lib
    ! I/O
    integer, intent(in)     :: n_obs, n_lev, n_poly, n_total
    real(8), intent(in)     :: model_profiles(n_obs, n_lev), err_model_profiles(n_obs, n_lev), model_psurf(n_obs), AT(n_lev+1), BT(n_lev+1), ak_poly(n_total, n_poly), prior_mix
    integer, intent(in)     :: indices(n_obs)
    real(8), intent(out)    :: J_obs, departures(n_obs, n_lev)
    real(4), intent(out)    :: mismatches(n_obs,mismatch_params)
    ! local variables
    integer                 :: i, idx, n_eval
    real(8), allocatable    :: departures_thread(:,:)
    real(4), allocatable    :: mismatches_thread(:,:)
    logical, allocatable    :: use_thisthread(:)
    real(8)                 :: J_thread, model_pres_levels(n_lev+1), mid_pres_levels(n_lev), col_ak(n_lev), prior_profile(n_lev)
    real(8)                 :: prior_tc, deps, pres_weights(n_lev), ret_profile(n_lev), H(n_lev)
    real(4)                 :: model_tc, tot_dtc, mod_dtc

    mismatches = 0.0
    J_obs = 0.0
    departures = 0.0

    prior_profile = prior_mix
    prior_tc = prior_mix

    !$omp parallel private(i, mismatches_thread, departures_thread, use_thisthread, J_thread) &
    !$omp private(idx, model_pres_levels, mid_pres_levels, col_ak, model_tc, mod_dtc, tot_dtc) &
    !$omp private(deps, n_eval, ret_profile, H, pres_weights)
    J_thread = 0.0
    allocate(mismatches_thread(n_obs, mismatch_params))
    allocate(departures_thread(n_obs, n_lev))
    allocate(use_thisthread(n_obs))
    mismatches_thread = 0.0
    departures_thread = 0.0
    use_thisthread = .false.
    !$omp do schedule(dynamic, 1000)
    do i = 1, n_obs
        idx = indices(i)
        model_pres_levels = AT + BT*model_psurf(i)
        pres_weights = (model_pres_levels(1:n_lev)-model_pres_levels(2:n_lev+1))/model_psurf(i)
        mid_pres_levels = 0.5*(model_pres_levels(1:n_lev) + model_pres_levels(2:n_lev+1))/model_psurf(i)
        col_ak = polyval(ak_poly(idx+1,:), mid_pres_levels) ! the indices are PYTHON indices, so their lowest value is zero, for example
        ret_profile = col_ak*(model_profiles(i,:)-prior_profile) + prior_profile
        model_tc = sum(ret_profile*pres_weights)
        H = col_ak * pres_weights
        mod_dtc = dot_product(err_model_profiles(i,:)**2, H**2)
        tot_dtc = mod_dtc+measured_totalcol_error**2
        mismatches_thread(i,:) = (/ real(idx), model_tc, measured_totalcol, tot_dtc, real(i-1) /)
        deps = (model_tc - measured_totalcol)/tot_dtc
        J_thread = J_thread + 0.5*deps*(model_tc-measured_totalcol)
        departures_thread(i,:) = H * deps
        use_thisthread(i) = .true.
    end do
    !$omp end do
    ! put all the thread results together now
    !$omp critical
    n_eval = 0
    do i = 1, n_obs
        if (use_thisthread(i)) then
            mismatches(i,:) = mismatches_thread(i,:)
            departures(i,:) = departures_thread(i,:)
            n_eval = n_eval + 1
        end if
    end do
    J_obs = J_obs + J_thread
!    write(*,'("Thread ", i2, " evaluated ", i8, " soundings")') omp_get_thread_num()+1, n_eval
    !$omp end critical
    deallocate(mismatches_thread, departures_thread, use_thisthread)
    !$omp end parallel

end subroutine takecolumn_createobs

end module carbonsat

module mopitt

implicit none

private

public :: applyaveragingkernel

double precision, parameter :: m_air = 28.94E-3 ! Kg/mole for dry air
double precision, parameter :: Na = 6.02214E23  ! Avogadro's number
double precision, parameter :: grav = 9.81      ! gravitational acceleration
integer, parameter          :: mismatch_params = 3
double precision, parameter :: log_thresh = 0.1

interface ulog
    module procedure ulog_0d
    module procedure ulog_1d
end interface

interface dulog
    module procedure dulog_0d
    module procedure dulog_1d
end interface dulog

contains

pure function ulog_0d(x)

    implicit none
    double precision, intent(in)    :: x
    double precision                :: ulog_0d

    if (x .lt. log_thresh) then
        ulog_0d = tanh(x/log_thresh - 1.0) + log(log_thresh)
    else
        ulog_0d = log(x)
    end if

end function ulog_0d

pure function ulog_1d(x)

    implicit none
    double precision, intent(in)    :: x(:)
    double precision                :: ulog_1d(size(x))
    integer                         :: N, i

    N = size(x)
    do i = 1, N
        if (x(i) .lt. log_thresh) then
            ulog_1d(i) = tanh(x(i)/log_thresh - 1.0) + log(log_thresh)
        else
            ulog_1d(i) = log(x(i))
        end if
    end do

end function ulog_1d

pure function dulog_0d(x)

    implicit none
    double precision, intent(in)    :: x
    double precision                :: dulog_0d

    if (x .lt. log_thresh) then
        dulog_0d = 1.0/(log_thresh * cosh(x/log_thresh - 1.0)**2)
    else
        dulog_0d = 1.0/x
    end if

end function dulog_0d

pure function dulog_1d(x)

    implicit none
    double precision, intent(in)    :: x(:)
    double precision                :: dulog_1d(size(x))
    integer                         :: N, i

    N = size(x)
    do i = 1, N
        if (x(i) .lt. log_thresh) then
            dulog_1d(i) = 1.0/(log_thresh * cosh(x(i)/log_thresh - 1.0)**2)
        else
            dulog_1d(i) = 1.0/x(i)
        end if
    end do

end function dulog_1d

pure subroutine rebin_pressure(model_grid, model_prof, sat_grid, sat_prof, trans_mat)

    implicit none
    ! Rebin a gas profile model_prof defined on model_grid, onto a new grid sat_grid. The grids are given as pressure
    ! levels (descending in pressure, so increasing in altitude), and the profiles are mixing ratios. The matrix
    ! trans_mat should pre-multiply model_prof to give sat_prof. In general the mass of the total column need not be
    ! conserved, since there is no requirement that the end points of the two grids should match.
    !
    ! I/O
    real(8), intent(in)     :: model_grid(:), model_prof(:), sat_grid(:)
    real(8), intent(out)    :: sat_prof(size(sat_grid)-1), trans_mat(size(sat_grid)-1, size(model_grid)-1)
    ! Local
    integer                 :: i, j, model_levels, sat_levels
    real(8)                 :: frac, up, dn, dp_sat
    real(8), allocatable    :: local_model_grid(:)

    model_levels = size(model_grid)-1 ! assume that size(model_prof) = model_levels, without checking
    sat_levels = size(sat_grid)-1
    allocate(local_model_grid(model_levels+1))

    local_model_grid = model_grid
    ! extend the model grid, if needed
    if (local_model_grid(1) < sat_grid(1))                          local_model_grid(1) = sat_grid(1)
    if (local_model_grid(model_levels+1) > sat_grid(sat_levels+1))  local_model_grid(model_levels+1) = sat_grid(sat_levels+1)

    trans_mat = 0.0
    do i = 1, sat_levels
        dp_sat = sat_grid(i) - sat_grid(i+1) ! pressure thickness of satellite layer
        do j = 1, model_levels
            ! The model level is contained between pressures local_model_grid(j) and local_model_grid(j+1), whereas the satellite
            ! level is contained between pressures sat_grid(i) and sat_grid(i+1). If there is no overlap, frac = 0.
            if (local_model_grid(j) .le. sat_grid(i+1) .or. local_model_grid(j+1) .ge. sat_grid(i)) then
                frac = 0.0
            else
                up = max(local_model_grid(j+1), sat_grid(i+1)) ! lower of the two upper boundaries
                dn = min(local_model_grid(j), sat_grid(i)) ! upper of the two lower boundaries
                frac = (dn-up)/dp_sat
            end if
            trans_mat(i,j) = frac
        end do
    end do

    deallocate(local_model_grid)

    sat_prof = matmul(trans_mat, model_prof)

end subroutine rebin_pressure

subroutine applyaveragingkernel(n_obs, n_lay_mod, n_lay_sat, model_profile, model_plevs, meas_col, sigma_meas_col, sat_plevs, &
    valid_layers, prior_profiles, avg_ker, J_obs, departures, mismatches, valid_deps)

    implicit none
    ! To apply the MOPITT averaging kernel, we need
    !
    ! (1) The modeled profiles at a bunch of (x,y,t) points             --> model_profile
    ! (2) The modeled pressure levels at a bunch of (x,y,t) points      --> model_plevs
    ! (3) The observed total column at a bunch of (x,y,t) points        --> meas_col
    ! (4) The observed total column error at a bunch of (x,y,t) points  --> sigma_meas_col
    ! (4) The retrieval pressure levels at a bunch of (x,y,t) points    --> sat_plevs
    ! (5) For each retrieval, the number of "valid levels"              --> valid_layers
    ! (6) The retrieval prior profiles at a bunch of (x,y,t) points     --> prior_profiles
    ! (7) The retrieval averaging kernels at a bunch of (x,y,t) points  --> avg_ker
    !
    ! For IASI, variables from the input file were passed en toto, i.e., instead of subselecting for a region. For MOPITT,
    ! we will be a bit more clever, and subselect valid_layers, prior_profiles, sat_plevs and avg_ker for the region.
    !
    ! I/O :: input
    integer, intent(in)         :: n_obs, n_lay_mod, n_lay_sat
    real(8), intent(in)         :: model_profile(n_obs, n_lay_mod), model_plevs(n_obs, n_lay_mod+1)
    real(8), intent(in)         :: meas_col(n_obs), sigma_meas_col(n_obs)
    real(8), intent(in), target :: sat_plevs(n_obs, n_lay_sat+1)
    integer, intent(in)         :: valid_layers(n_obs)
    real(8), intent(in), target :: prior_profiles(n_obs, n_lay_sat), avg_ker(n_obs, n_lay_sat, n_lay_sat)
    ! I/O :: output
    real(8), intent(out)        :: J_obs, departures(n_obs, n_lay_mod), mismatches(n_obs, mismatch_params)
    integer(1), intent(out)     :: valid_deps(n_obs)
    ! Local variables
    integer                     :: i, l
    real(8), target             :: transfer_matrix(n_lay_sat, n_lay_mod), rebinned_model(n_lay_sat)
    real(8), pointer            :: trans_mat(:,:), prior(:), ak(:,:), satp(:), rbm(:)
    real(8)                     :: q(n_lay_sat), XCO, del_p(n_lay_sat), u(n_lay_sat), v(n_lay_sat), prefac

    mismatches = 0.0
    J_obs = 0.0
    departures = 0.0
    ! assume all observations are invalid, and later correct for the valid ones
    valid_deps = 0

    !$omp parallel do private(i, l, prior, ak, satp, trans_mat, transfer_matrix, rebinned_model, rbm) &
    !$omp private(XCO, del_p, q, u, v, prefac) reduction(+:departures, J_obs, mismatches, valid_deps) &
    !$omp schedule(dynamic, 10)
    do i = 1, n_obs
        ! choose the sub-portions of the satellite column which are valid
        l = valid_layers(i)
        prior => prior_profiles(i, 1:l)
        ak => avg_ker(i, 1:l, 1:l)
        satp => sat_plevs(i, 1:l+1)
        ! trans_mat should then be l x n_lay_mod
        trans_mat => transfer_matrix(1:l, :)
        rbm => rebinned_model(1:l)

        ! rebin the modeled profile on to the satellite grid
        call rebin_pressure(model_plevs(i,:), model_profile(i,:), satp, rbm, trans_mat)
        ! rbm is the rebinned model profile, trans_mat the transformation matrix
        !q(1:l) = exp(matmul(ak, log(rbm) - log(prior)))
        q(1:l) = exp(matmul(ak, ulog(rbm) - log(prior)))
        del_p(1:l) = satp(1:l) - satp(2:l+1)
        ! we will reuse \del p * x^pri * q, so store it
        v(1:l) = del_p(1:l) * prior * q(1:l)
        XCO = sum(v(1:l))/satp(1) ! modeled total column

        J_obs = J_obs + 0.5 * ((XCO - meas_col(i))/sigma_meas_col(i))**2 ! addition to the cost function

        ! Now calculate the adjoint forcings
        !u(1:l) = matmul(transpose(ak), v(1:l))/rbm
        u(1:l) = matmul(transpose(ak), v(1:l)) * dulog(rbm)
        prefac = (XCO - meas_col(i))/(sigma_meas_col(i)**2 * satp(1))
        departures(i,:) = prefac * matmul(transpose(trans_mat), u(1:l))

        ! Each row of mismatches contains (j-1), modeled total column, measured total column, (MDM)**2, (i-1)
        mismatches(i,:) = (/ XCO, meas_col(i), sigma_meas_col(i)**2 /) ! haven't yet figured out how i and j are different
        valid_deps(i) = 1

    end do
    !$omp end parallel do

end subroutine applyaveragingkernel

end module mopitt

module iasi

implicit none

private

public :: applyaveragingkernel

integer, parameter  :: mismatch_params = 3

contains

pure subroutine rebin_partialcols(old_grid, old_prof, new_grid, new_prof, trans_mat)
    ! old_grid is the altitude grid of TM5, given on old_levels+1 altitudes
    ! old_prof is the profile of partial column CO in #molecules cm-2 on old_levels layers
    ! new_grid  = IASI altitude grid on old_levels+1 altitudes
    ! new_prof  = new profile of partial column CO in #molecules cm-2 on new_levels layers
    ! Note that in general the masses don't have to be conserved, since there is no
    ! restriction that the end-points of the two grids should match.
    implicit none

    ! in/out
    real(8), intent(in)  :: old_grid(:), new_grid(:), old_prof(:)
    real(8), intent(out) :: new_prof(size(new_grid)-1), trans_mat(size(new_grid)-1, size(old_grid)-1)
    ! local
    integer             :: i, j, old_levels, new_levels
    double precision    :: frac, up, dn

    old_levels = size(old_grid)-1 ! assume as well that size(old_prof) = old_levels
    new_levels = size(new_grid)-1

    trans_mat = 0.0
    do i = 1, new_levels
        do j = 1, old_levels
            ! frac = fraction of the j-th old layer that is inside the i-th new layer
            ! the j-th old layer is between old_grid(j) and old_grid(j+1)
            ! the i-th new layer is between new_grid(i) and new_grid(i+1)
            if ((old_grid(j+1) .le. new_grid(i)) .or. (old_grid(j) .ge. new_grid(i+1))) then
                frac = 0.0
            else
                ! so there is some overlap...
                up = min(old_grid(j+1), new_grid(i+1))
                dn = max(old_grid(j), new_grid(i))
                frac = (up-dn)/(old_grid(j+1)-old_grid(j))
            end if
            trans_mat(i,j) = frac
        end do
    end do

    new_prof = matmul(trans_mat, old_prof)

end subroutine rebin_partialcols

pure subroutine redis_partialcol_altitude(old_grid,old_prof,new_grid,new_prof,matr)
    !
    ! Note :: This subroutine is not used any more
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
    real(8), intent(in)  :: old_grid(:), new_grid(:), old_prof(:)
    real(8), intent(out) :: new_prof(size(new_grid)-1), matr(size(new_grid)-1, size(old_grid)-1)

    ! local
    integer             :: i, j, last, old_levels, new_levels
    double precision    :: frac

    old_levels = size(old_grid)-1 ! assume as well that size(old_prof) = old_levels
    new_levels = size(new_grid)-1

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
    end do

    new_prof = matmul(matr,old_prof)

end subroutine redis_partialcol_altitude

subroutine applyaveragingkernel(n_obs, n_lev, n_lay, model_profile, model_psurf, meas_tc, meas_dtc, gph_grid, AT, BT, prior_profiles, avg_ker, mismatches, J_obs, departures, valid_deps)

    use omp_lib
    ! I/O
    integer, intent(in)     :: n_obs, n_lev, n_lay
    real(8), intent(in)     :: model_profile(n_obs,n_lev), model_psurf(n_obs), gph_grid(n_obs,n_lev+1)
    real(8), intent(in)     :: meas_tc(n_obs), meas_dtc(n_obs), AT(n_lev+1), BT(n_lev+1)
    real(8), intent(in)     :: prior_profiles(n_obs,n_lay), avg_ker(n_obs,n_lay)
    !integer, intent(in)     :: indices(n_obs)
    real(8), intent(out)    :: J_obs, departures(n_obs,n_lev)
    real(8), intent(out)    :: mismatches(n_obs,mismatch_params)
    ! We should avoid I/O with logical variables, because on the python side there is no guarantee about
    ! the value a True/False would be interpreted as. For example, with xlf, a T was 1 and an F was 0,
    ! but with ifort a T is -1 and an F is 0. So locally we can have a logical variable but the I/O must
    ! be with an integer array, say 8-bit, which will be enough.
    integer(1), intent(out) :: valid_deps(n_obs)
    ! local variables
    integer                 :: i, iasi_l_idx, model_l_idx
    real(8)                 :: IASI_alt_grid(n_lay), dp(n_lev), tm5_pres(n_lev+1)
    real(8)                 :: Mair = 28.94E-3, N_a = 6.02214E23, grav = 9.81, mult_const, model_tc_mol, var_tot, dep, model_alt_grid(n_lay+1)
    real(8), target         :: model_IASI_grid(n_lay), transfer_matrix(n_lay,n_lev)
    real(8), pointer        :: model_on_IASI(:), trans_mat(:,:)
    real(8)                 :: ppb_to_mol(n_lev), model_prof_mol(n_lev), H(n_lev), select_layers(n_lev)
    logical                 :: use_it
    ! Thread-specific variables
    real(8), allocatable    :: departures_thread(:,:)
    real(8), allocatable    :: mismatches_thread(:,:)
    real(8)                 :: J_thread
    logical, allocatable    :: valid_thread(:), use_thisthread(:)
    integer                 :: n_eval

    mismatches = 0.0
    J_obs = 0.0
    departures = 0.0
    ! assume all observations are invalid, and later correct for the valid ones
    valid_deps = 0

    ! What are the grids we need?
    ! 1. The IASI altitude grid, [0., 1000., ... 18000.], which are the lower boundaries of the 19 layers
    ! 2. The modeled geopotential heights at the pressure boundaries, 26 for a 25-layer model
    ! 3. The model altitude grid, [0., 1000., ... 18000., 200000.], which are the IASI altitudes plus the topmost geopotential height
    do i = 1, n_lay
        IASI_alt_grid(i) = (i-1)*1000.0
    end do

    mult_const = 1.0E-13*N_a/grav/Mair

    !$omp parallel private(i, departures_thread, mismatches_thread, J_thread, valid_thread, model_alt_grid) &
    !$omp private(tm5_pres, dp, ppb_to_mol, iasi_l_idx, use_it, model_l_idx, model_prof_mol, model_IASI_grid) &
    !$omp private(transfer_matrix, model_on_IASI, trans_mat, model_tc_mol, var_tot, dep, select_layers, H) &
    !$omp private(n_eval, use_thisthread)
    J_thread = 0.0
    allocate(departures_thread(n_obs,n_lev))
    allocate(mismatches_thread(n_obs,mismatch_params))
    allocate(valid_thread(n_obs))
    allocate(use_thisthread(n_obs))
    departures_thread = 0.0
    mismatches_thread = 0.0
    valid_thread = .false.
    use_thisthread = .false.
    !$omp do schedule(dynamic, 50)
    do i = 1, n_obs
        !j = indices(i)+1
        model_alt_grid(1:n_lay) = IASI_alt_grid
        model_alt_grid(n_lay+1) = gph_grid(i, n_lev+1) ! [0., 1000., ... 18000., 200000.]
        tm5_pres = AT + BT*model_psurf(i)
        dp = tm5_pres(1:n_lev) - tm5_pres(2:n_lev+1)
        ppb_to_mol = mult_const * dp ! ppb_to_mol is a conversion factor, to go from ppb to molecules per square cm
        ! Which are the valid IASI layers? The ones with positive priors
        iasi_l_idx = 1 ! lower index from which to consider IASI prior and ak
        do while (prior_profiles(i,iasi_l_idx) < -1.0)
            iasi_l_idx = iasi_l_idx + 1
        end do
        ! if an entire valid IASI layer is below TM5 ground, skip that observation
        use_it = .true.
        if (IASI_alt_grid(iasi_l_idx+1) <= gph_grid(i,1)) use_it = .false.
        if (use_it) then
            ! Which TM5 layers to use? For example, if the valid IASI lower levels are [1000., 2000., ...] and the TM5 GPH
            ! values are [-5., 52., 250., 600., 1125. 1700., ...] then we want to use the model layers starting from 600m
            ! and neglect those below. In that case, model_l_idx = 4, corresponding to iasi_l_idx = 2.
            model_l_idx = 1
            do while (gph_grid(i,model_l_idx+1) <= IASI_alt_grid(iasi_l_idx))
                model_l_idx = model_l_idx + 1
            end do
            ! Now we need to take the model_profile(model_l_idx:n_lev) and re-distribute it within model_alt_grid(iasi_l_idx:n_lay+1)
            ! Before that, we need to convert the model_profile from ppb to molecules/cm^2
            model_prof_mol = ppb_to_mol * model_profile(i,1:n_lev)
            model_IASI_grid = 0.0
            transfer_matrix = 0.0
            model_on_IASI => model_IASI_grid(iasi_l_idx:n_lay)
            trans_mat => transfer_matrix(iasi_l_idx:n_lay,model_l_idx:n_lev)
            call rebin_partialcols(gph_grid(i,model_l_idx:n_lev+1), model_prof_mol(model_l_idx:n_lev), model_alt_grid(iasi_l_idx:n_lay+1), model_on_IASI, trans_mat)
            ! Now model_on_IASI is the modeled CO profile on the (iasi_l_idx:n_lay) layers, in molecules/cm^2
            ! The relevant part of the averaging kernel is avg_ker(i,iasi_l_idx:n_lay), and that of the prior profile is prior_profiles(i,iasi_l_idx:n_lay)
            ! Use those to construct the modeled total column
            model_tc_mol = dot_product(avg_ker(i,iasi_l_idx:n_lay), model_on_IASI-prior_profiles(i,iasi_l_idx:n_lay)) + sum(prior_profiles(i,iasi_l_idx:n_lay))
            ! Assume that meas_dtc has already been inflated
            var_tot = meas_dtc(i)**2
            ! The corresponding measured IASI total column is meas_tc(i)
            dep = (model_tc_mol-meas_tc(i))/var_tot
            !J_obs = J_obs + 0.5*(model_tc_mol-meas_tc(i))*dep
            J_thread = J_thread + 0.5*(model_tc_mol-meas_tc(i))*dep
            ! What is \partial model XCO/\partial model profile in ppb?
            ! it is (avg_ker * trans_mat * ppb_to_mol * select_layers), where select_layers is a vector with ones for model layers used, zeros for layers not used
            select_layers = 0.0
            select_layers(model_l_idx:n_lev) = 1.0
            H = matmul(avg_ker(i,:), transfer_matrix)
            H = H * ppb_to_mol * select_layers
            !mismatches(i,:) = (/ real(j-1), real(model_tc_mol), real(meas_tc(i)), real(var_tot), real(i-1) /)
            !mismatches_thread(i,:) = (/ dble(j-1), model_tc_mol, meas_tc(i), var_tot, dble(i-1) /)
            mismatches_thread(i,:) = (/ model_tc_mol, meas_tc(i), dble(i-1) /)
            !departures(i,:) = H * dep
            departures_thread(i,:) = H * dep
        else
            !mismatches_thread(i,:) = (/ dble(j-1), 0.0d0, 0.0d0, 1.0d0, dble(i-1) /)
            mismatches_thread(i,:) = (/ 0.0d0, 0.0d0, dble(i-1) /)
            departures_thread(i,:) = 0.0
        end if
        valid_thread(i) = use_it
        use_thisthread(i) = .true.
    end do
    !$omp end do
    ! put all the thread results together now
    !$omp critical
    n_eval = 0
    do i = 1, n_obs
        if (use_thisthread(i)) then
            mismatches(i,:) = mismatches_thread(i,:)
            departures(i,:) = departures_thread(i,:)
            if (valid_thread(i)) valid_deps(i) = 1
            n_eval = n_eval + 1
        end if
    end do
    J_obs = J_obs + J_thread
    !$omp end critical
    deallocate(mismatches_thread, departures_thread, valid_thread, use_thisthread)
    !$omp end parallel

end subroutine applyaveragingkernel

end module iasi
