!> Compare serial rec against Coarray-partitioned rec_caf on fixed cases.
!> Drives the shipped gauss_jacobi dispatch (not a reimplementation).
program test_caf_rec
use GaussJacobiQuad, only: gauss_jacobi
use gjp_rec, only: caf_root_owner, caf_owns_index
use gjp_types, only: dp
implicit none

integer, parameter :: ncase = 2
integer :: npts_list(ncase)
real(dp) :: alpha_list(ncase), beta_list(ncase)
real(dp), allocatable :: x_ser(:), w_ser(:), x_caf(:), w_caf(:)
real(dp) :: max_dx, max_dw, tol
integer :: ic, i, n, me, nimg, owner, n_owned
character(len=:), allocatable :: method_ser, method_caf
logical :: ok

me = this_image()
nimg = num_images()
ok = .true.
tol = 1.0e-12_dp

! Ownership helper sanity (partition covers 1..n exactly once)
n = 17
n_owned = 0
do i = 1, n
    owner = caf_root_owner(i, nimg)
    if (owner < 1 .or. owner > nimg) then
        if (me == 1) print*,"FAIL: caf_root_owner out of range for i=", i
        ok = .false.
    end if
    if (caf_owns_index(i, me, nimg)) n_owned = n_owned + 1
end do
! every index assigned; this image owns floor/ceil share
if (n_owned < n / nimg .or. n_owned > (n + nimg - 1) / nimg) then
    if (me == 1) print*,"FAIL: ownership count", n_owned, "for nimg=", nimg
    ok = .false.
end if

npts_list = [32, 48]
alpha_list = [0.5_dp, 0.0_dp]
beta_list = [0.5_dp, 1.5_dp]
method_ser = "rec"
method_caf = "rec_caf"

do ic = 1, ncase
    n = npts_list(ic)
    allocate (x_ser(n), w_ser(n), x_caf(n), w_caf(n))
    call gauss_jacobi(n, alpha_list(ic), beta_list(ic), x_ser, w_ser, method_ser)
    call gauss_jacobi(n, alpha_list(ic), beta_list(ic), x_caf, w_caf, method_caf)

    max_dx = maxval(abs(x_ser - x_caf))
    max_dw = maxval(abs(w_ser - w_caf))

    if (me == 1) then
        print '(A,I0,A,F6.2,A,F6.2,A,ES12.5,A,ES12.5)', &
            "case n=", n, " a=", alpha_list(ic), " b=", beta_list(ic), &
            " max|dx|=", max_dx, " max|dw|=", max_dw
        print '(A,I0,A,I0)', "  images: me=", me, " nimg=", nimg
    end if

    if (max_dx > tol .or. max_dw > tol) then
        if (me == 1) print*,"FAIL: rec vs rec_caf disagree beyond tol=", tol
        ok = .false.
    end if
    deallocate (x_ser, w_ser, x_caf, w_caf)
end do

sync all

if (.not. ok) then
    error stop "test_caf_rec failed"
end if
if (me == 1) print*,"test_caf_rec PASSED (rec vs rec_caf, nimg=", nimg, ")"
end program test_caf_rec
