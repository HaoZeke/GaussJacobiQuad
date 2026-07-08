!> Timed multi-image CAF driver for all quadrature methods.
!>
!> Modes:
!>   single  — one rule; method=rec_caf uses node partition (speedup for rec).
!>   batch   — nbatch independent rules (α fixed, β varies); partitions rules
!>             across images for rec / gw / algo665 (speedup for all).
!>
!> Usage:
!>   gjp_bench_caf single <method> <npts> <alpha> <beta> [ntrials] [nwarm]
!>   gjp_bench_caf batch  <method> <npts> <nbatch> <alpha> <beta0> [ntrials] [nwarm]
!>
!> method: rec_caf | rec | gw | algo665
program gjp_bench_caf
use GaussJacobiQuad, only: gauss_jacobi, gauss_jacobi_batch_caf
use gjp_types, only: dp
implicit none

integer :: npts, ntrials, nwarm, it, me, nimg, nbatch, ib, narg
integer :: count0, count1, count_rate, count_max
real(dp) :: alpha, beta, beta0, elapsed, sum_t, mean_t
real(dp), allocatable :: x(:), w(:), x_ref(:), w_ref(:)
real(dp), allocatable :: xs(:, :), ws(:, :), xs_ref(:, :), ws_ref(:, :)
real(dp), allocatable :: alphas(:), betas(:)
real(dp) :: max_dx, max_dw
character(len=128) :: arg, mode
character(len=:), allocatable :: method

me = this_image()
nimg = num_images()
ntrials = 5
nwarm = 2
narg = command_argument_count()

if (narg < 1) then
    if (me == 1) call print_usage()
    error stop "need mode"
end if
call get_command_argument(1, mode)

select case (trim(mode))
case ("single")
    if (narg < 5) then
        if (me == 1) call print_usage()
        error stop "single: need method npts alpha beta"
    end if
    call get_command_argument(2, arg)
    method = trim(arg)
    call get_command_argument(3, arg)
    read (arg, *) npts
    call get_command_argument(4, arg)
    read (arg, *) alpha
    call get_command_argument(5, arg)
    read (arg, *) beta
    if (narg >= 6) then
        call get_command_argument(6, arg)
        read (arg, *) ntrials
    end if
    if (narg >= 7) then
        call get_command_argument(7, arg)
        read (arg, *) nwarm
    end if
    call run_single()
case ("batch")
    if (narg < 6) then
        if (me == 1) call print_usage()
        error stop "batch: need method npts nbatch alpha beta0"
    end if
    call get_command_argument(2, arg)
    method = trim(arg)
    call get_command_argument(3, arg)
    read (arg, *) npts
    call get_command_argument(4, arg)
    read (arg, *) nbatch
    call get_command_argument(5, arg)
    read (arg, *) alpha
    call get_command_argument(6, arg)
    read (arg, *) beta0
    if (narg >= 7) then
        call get_command_argument(7, arg)
        read (arg, *) ntrials
    end if
    if (narg >= 8) then
        call get_command_argument(8, arg)
        read (arg, *) nwarm
    end if
    call run_batch()
case default
    if (me == 1) call print_usage()
    error stop "mode must be single or batch"
end select

contains

subroutine print_usage()
    print *, "Usage:"
    print *, "  gjp_bench_caf single <method> <npts> <alpha> <beta> [ntrials] [nwarm]"
    print *, "  gjp_bench_caf batch  <method> <npts> <nbatch> <alpha> <beta0> [ntrials] [nwarm]"
    print *, "method: rec_caf|rec|gw|algo665"
    print *, "single: rec_caf node-partition; batch: multi-rule CAF for all methods"
end subroutine print_usage

subroutine run_single()
    character(len=:), allocatable :: method_ref
    allocate (x(npts), w(npts), x_ref(npts), w_ref(npts))

    ! Serial reference method (deferred-length for gauss_jacobi dummy)
    select case (trim(method))
    case ("rec_caf")
        method_ref = "rec"
    case ("gw")
        method_ref = "gw"
    case ("algo665")
        method_ref = "algo665"
    case ("rec")
        method_ref = "rec"
    case default
        error stop "single: unsupported method"
    end select

    call gauss_jacobi(npts, alpha, beta, x_ref, w_ref, method_ref)
    call gauss_jacobi(npts, alpha, beta, x, w, method)
    max_dx = maxval(abs(x - x_ref))
    max_dw = maxval(abs(w - w_ref))
    if (me == 1) then
        print '(A,A,A,I0,A,I0,A,ES12.5,A,ES12.5)', &
            "CORRECT mode=single method=", trim(method), " nimg=", nimg, &
            " npts=", npts, " max|dx|=", max_dx, " max|dw|=", max_dw
    end if
    if (max_dx > 1.0e-12_dp .or. max_dw > 1.0e-12_dp) then
        if (me == 1) print *, "FAIL single correctness"
        error stop "single correctness failed"
    end if

    do it = 1, nwarm
        call gauss_jacobi(npts, alpha, beta, x, w, method)
    end do
    sync all

    sum_t = 0.0_dp
    do it = 1, ntrials
        sync all
        call system_clock(count0, count_rate, count_max)
        call gauss_jacobi(npts, alpha, beta, x, w, method)
        call system_clock(count1, count_rate, count_max)
        sync all
        elapsed = real(count1 - count0, dp) / real(count_rate, dp)
        sum_t = sum_t + elapsed
        if (me == 1) then
            print '(A,A,1X,A,I0,1X,A,I0,1X,A,I0,1X,A,ES16.8)', &
                "TRIAL mode=single method=", trim(method), "nimg=", nimg, &
                "npts=", npts, "i=", it, "sec=", elapsed
        end if
    end do
    mean_t = sum_t / real(ntrials, dp)
    if (me == 1) then
        print '(A,A,1X,A,I0,1X,A,I0,1X,A,ES16.8)', &
            "MEAN mode=single method=", trim(method), "nimg=", nimg, &
            "npts=", npts, "sec=", mean_t
    end if
    deallocate (x, w, x_ref, w_ref)
end subroutine run_single

subroutine run_batch()
    allocate (alphas(nbatch), betas(nbatch))
    allocate (xs(npts, nbatch), ws(npts, nbatch))
    allocate (xs_ref(npts, nbatch), ws_ref(npts, nbatch))

    do ib = 1, nbatch
        alphas(ib) = alpha
        ! keep beta > -1; small positive stride
        betas(ib) = beta0 + 0.01_dp * real(ib - 1, dp)
    end do

    ! Serial reference: every rule on each image (only image 1 must match)
    do ib = 1, nbatch
        call gauss_jacobi(npts, alphas(ib), betas(ib), xs_ref(:, ib), ws_ref(:, ib), method)
        ! method is already deferred-length allocatable
    end do

    call gauss_jacobi_batch_caf(npts, nbatch, alphas, betas, xs, ws, method)
    max_dx = maxval(abs(xs - xs_ref))
    max_dw = maxval(abs(ws - ws_ref))
    if (me == 1) then
        print '(A,A,A,I0,A,I0,A,I0,A,ES12.5,A,ES12.5)', &
            "CORRECT mode=batch method=", trim(method), " nimg=", nimg, &
            " npts=", npts, " nbatch=", nbatch, &
            " max|dx|=", max_dx, " max|dw|=", max_dw
    end if
    if (max_dx > 1.0e-12_dp .or. max_dw > 1.0e-12_dp) then
        if (me == 1) print *, "FAIL batch correctness"
        error stop "batch correctness failed"
    end if

    do it = 1, nwarm
        call gauss_jacobi_batch_caf(npts, nbatch, alphas, betas, xs, ws, method)
    end do
    sync all

    sum_t = 0.0_dp
    do it = 1, ntrials
        sync all
        call system_clock(count0, count_rate, count_max)
        call gauss_jacobi_batch_caf(npts, nbatch, alphas, betas, xs, ws, method)
        call system_clock(count1, count_rate, count_max)
        sync all
        elapsed = real(count1 - count0, dp) / real(count_rate, dp)
        sum_t = sum_t + elapsed
        if (me == 1) then
            print '(A,A,1X,A,I0,1X,A,I0,1X,A,I0,1X,A,I0,1X,A,ES16.8)', &
                "TRIAL mode=batch method=", trim(method), "nimg=", nimg, &
                "npts=", npts, "nbatch=", nbatch, "i=", it, "sec=", elapsed
        end if
    end do
    mean_t = sum_t / real(ntrials, dp)
    if (me == 1) then
        print '(A,A,1X,A,I0,1X,A,I0,1X,A,I0,1X,A,ES16.8)', &
            "MEAN mode=batch method=", trim(method), "nimg=", nimg, &
            "npts=", npts, "nbatch=", nbatch, "sec=", mean_t
    end if
    deallocate (alphas, betas, xs, ws, xs_ref, ws_ref)
end subroutine run_batch

end program gjp_bench_caf
