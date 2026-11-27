!> fpm builds src/gjp_version.f90. The constants must match fpm.toml.
program test_gjp_version_fpm
use gjp_version, only: gjp_version_major, gjp_version_minor, gjp_version_patch, &
    gjp_version_string, gjp_version_compact
implicit none
character(len=1024) :: exe, dir, path
character(len=256) :: line, quoted
integer :: ios, u, i, j, level, major, minor, patch
logical :: found

call get_command_argument(0, exe, status=ios)
if (ios /= 0) error stop "argv0"
call parent_of(exe, dir)
found = .false.
do level = 1, 8
  path = trim(dir) // "/fpm.toml"
  inquire(file=trim(path), exist=found)
  if (found) exit
  call parent_of(dir, dir)
end do
if (.not. found) then
  write(*, "(A)") "exe=" // trim(exe)
  error stop "fpm.toml not found"
end if
open(newunit=u, file=trim(path), status="old", action="read", iostat=ios)
if (ios /= 0) error stop "open fpm.toml"
found = .false.
do
  read(u, "(A)", iostat=ios) line
  if (ios /= 0) exit
  if (len_trim(line) >= 8) then
    if (line(1:8) == "version ") then
      found = .true.
      exit
    end if
  end if
end do
close(u)
if (.not. found) error stop "version key"
i = index(line, '"')
if (i < 1) error stop "version quote"
j = index(line(i + 1:), '"')
if (j < 1) error stop "version end quote"
quoted = line(i + 1:i + j - 1)
if (trim(quoted) /= gjp_version_string) error stop "string /= fpm.toml"
call split3(trim(quoted), major, minor, patch, ios)
if (ios /= 0) error stop "parse fpm.toml version"
if (major /= gjp_version_major .or. minor /= gjp_version_minor .or. &
    patch /= gjp_version_patch) error stop "components /= fpm.toml"
if (gjp_version_major < 0 .or. gjp_version_major > 99 .or. &
    gjp_version_minor < 0 .or. gjp_version_minor > 99 .or. &
    gjp_version_patch < 0 .or. gjp_version_patch > 99) error stop "bounds"
if (gjp_version_compact /= &
    10000 * gjp_version_major + 100 * gjp_version_minor + gjp_version_patch) &
    error stop "compact"
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

subroutine split3(text, major, minor, patch, ios)
  character(len=*), intent(in) :: text
  integer, intent(out) :: major, minor, patch, ios
  character(len=len(text)) :: core
  integer :: cut, used
  core = text
  cut = scan(core, "+")
  if (cut > 0) core = core(1:cut - 1)
  cut = scan(core, "-")
  if (cut > 0) core = core(1:cut - 1)
  major = 0
  minor = 0
  patch = 0
  call take_int(core, major, used, ios)
  if (ios /= 0) return
  if (used > len_trim(core)) return
  if (core(used:used) /= ".") return
  core = core(used + 1:)
  call take_int(core, minor, used, ios)
  if (ios /= 0) return
  if (used > len_trim(core)) return
  if (core(used:used) /= ".") return
  core = core(used + 1:)
  call take_int(core, patch, used, ios)
end subroutine split3

subroutine take_int(text, value, used, ios)
  character(len=*), intent(in) :: text
  integer, intent(out) :: value, used, ios
  integer :: k
  if (len_trim(text) < 1 .or. text(1:1) < "0" .or. text(1:1) > "9") then
    ios = 1
    value = 0
    used = 1
    return
  end if
  value = 0
  k = 1
  do while (k <= len_trim(text) .and. text(k:k) >= "0" .and. text(k:k) <= "9")
    value = value * 10 + ichar(text(k:k)) - ichar("0")
    k = k + 1
  end do
  used = k
  ios = 0
end subroutine take_int

end program test_gjp_version_fpm
