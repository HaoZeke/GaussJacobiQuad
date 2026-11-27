!> fpm builds src/gjp_version.f90. fpm.toml stores major.minor.patch.
!> The module string may keep a pre-release suffix.
!> Finds fpm.toml by walking up from the executable, so this test must run
!> in-tree. An out-of-tree fpm --build-dir cannot see fpm.toml.
program test_gjp_version_fpm
use gjp_version, only: gjp_version_major, gjp_version_minor, gjp_version_patch, &
                       gjp_version_string, gjp_version_compact
implicit none
character(len=1024) :: exe, dir, path
character(len=256) :: line, quoted
character(len=32) :: triple
integer :: ios, u, i, j, level
logical :: found

call get_command_argument(0, exe, status=ios)
if (ios /= 0) error stop "argv0"
call parent_of(exe, dir)
found = .false.
do level = 1, 8
    path = trim(dir)//"/fpm.toml"
    inquire (file=trim(path), exist=found)
    if (found) exit
    call parent_of(dir, dir)
end do
if (.not. found) then
    write (*, "(A)") "exe="//trim(exe)
    error stop "fpm.toml not found"
end if
open (newunit=u, file=trim(path), status="old", action="read", iostat=ios)
if (ios /= 0) error stop "open fpm.toml"
found = .false.
do
    read (u, "(A)", iostat=ios) line
    if (ios /= 0) exit
    if (len_trim(line) >= 8) then
        if (line(1:8) == "version ") then
            found = .true.
            exit
        end if
    end if
end do
close (u)
if (.not. found) error stop "version key"
i = index(line, '"')
if (i < 1) error stop "version quote"
j = index(line(i + 1:), '"')
if (j < 1) error stop "version end quote"
quoted = line(i + 1:i + j - 1)
write (triple, '(I0,".",I0,".",I0)') gjp_version_major, gjp_version_minor, &
    gjp_version_patch
if (trim(quoted) /= trim(triple)) error stop "fpm.toml is not the release triple"
if (gjp_version_major < 0 .or. gjp_version_major > 99 .or. &
    gjp_version_minor < 0 .or. gjp_version_minor > 99 .or. &
    gjp_version_patch < 0 .or. gjp_version_patch > 99) error stop "bounds"
if (gjp_version_compact /= &
    10000 * gjp_version_major + 100 * gjp_version_minor + gjp_version_patch) &
    error stop "compact"
if (index(gjp_version_string, trim(triple)) /= 1 .and. .not. &
    (gjp_version_patch == 0 .and. &
     index(gjp_version_string//".0", trim(triple)) == 1)) &
    error stop "string triple"
print "(A)", "fpm version agrees with gjp_version"

contains

subroutine parent_of(path_in, out)
    character(len=*), intent(in) :: path_in
    character(len=*), intent(out) :: out
    integer :: n
    out = trim(path_in)
    n = len_trim(out)
    if (n < 1) then
        out = "."
        return
    end if
    do while (n > 1 .and. (out(n:n) == "/" .or. out(n:n) == char(92)))
        n = n - 1
    end do
    do while (n > 1 .and. out(n:n) /= "/" .and. out(n:n) /= char(92))
        n = n - 1
    end do
    if (n <= 1) then
        if (out(1:1) == "/") then
            out = "/"
        else
            out = "."
        end if
    else
        out = out(1:n - 1)
    end if
end subroutine parent_of

end program test_gjp_version_fpm
