!> Timed driver for multi-image rec_caf speedup measurement.
!> Times only the shipped gauss_jacobi(..., "rec_caf") call via system_clock.
!> Usage: gjp_bench_caf <npts> <alpha> <beta> [ntrials] [nwarm]
!> Prints machine-readable lines; image 1 only for summaries.
program gjp_bench_caf
use GaussJacobiQuad, only: gauss_jacobi
use gjp_types, only: dp
implicit none

integer :: npts, ntrials, nwarm, it, me, nimg
integer :: count0, count1, count_rate, count_max
real(dp) :: alpha, beta, elapsed, sum_t, mean_t
real(dp), allocatable :: x(:), w(:), x_ref(:), w_ref(:)
real(dp) :: max_dx, max_dw
character(len=128) :: arg
character(len=:), allocatable :: method, method_ser

me = this_image()
nimg = num_images()

ntrials = 5
nwarm = 2
if (command_argument_count() < 3) then
    if (me == 1) then
        print *, "Usage: gjp_bench_caf <npts> <alpha> <beta> [ntrials] [nwarm]"
    end if
    error stop "need at least 3 args"
end if

call get_command_argument(1, arg)
read (arg, *) npts
call get_command_argument(2, arg)
read (arg, *) alpha
call get_command_argument(3, arg)
read (arg, *) beta
if (command_argument_count() >= 4) then
    call get_command_argument(4, arg)
    read (arg, *) ntrials
end if
if (command_argument_count() >= 5) then
    call get_command_argument(5, arg)
    read (arg, *) nwarm
end if

allocate (x(npts), w(npts), x_ref(npts), w_ref(npts))
method = "rec_caf"
method_ser = "rec"  ! deferred-length for gauss_jacobi method dummy

! Correctness: multi-image rec_caf vs serial rec (on each image after compute)
call gauss_jacobi(npts, alpha, beta, x_ref, w_ref, method_ser)
call gauss_jacobi(npts, alpha, beta, x, w, method)
max_dx = maxval(abs(x - x_ref))
max_dw = maxval(abs(w - w_ref))
if (me == 1) then
    print '(A,I0,A,I0,A,I0,A,ES12.5,A,ES12.5)', &
        "CORRECT nimg=", nimg, " npts=", npts, " trials_planned=", ntrials, &
        " max|dx|=", max_dx, " max|dw|=", max_dw
end if
if (max_dx > 1.0e-12_dp .or. max_dw > 1.0e-12_dp) then
    if (me == 1) print *, "FAIL correctness"
    error stop "rec_caf disagree with rec"
end if

! Warmup (not timed)
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
        print '(A,I0,1X,A,I0,1X,A,I0,1X,A,ES16.8)', &
            "TRIAL nimg=", nimg, "npts=", npts, "i=", it, "sec=", elapsed
    end if
end do

mean_t = sum_t / real(ntrials, dp)
if (me == 1) then
    print '(A,I0,1X,A,I0,1X,A,ES16.8)', &
        "MEAN nimg=", nimg, "npts=", npts, "sec=", mean_t
end if

deallocate (x, w, x_ref, w_ref)
end program gjp_bench_caf
