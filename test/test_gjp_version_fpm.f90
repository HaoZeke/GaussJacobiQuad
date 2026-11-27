!> fpm builds src/gjp_version.f90. fpm.toml stores major.minor.patch.
!> The module string may keep a pre-release suffix.
program test_gjp_version_fpm
use gjp_version, only: gjp_version_major, gjp_version_minor, gjp_version_patch, &
                       gjp_version_string, gjp_version_compact
implicit none
character(len=1024) :: exe, dir, path
character(len=256) :: line, quoted
character(len=32) :: triple
integer :: ios, u, i, j, level, major, minor, patch
logical :: found

call expect_split("0.2.6", 0, 0, 2, 6)
call expect_split("1.2", 0, 1, 2, 0)
call expect_split("1.2.0rc1", 0, 1, 2, 0)
call expect_split("1.2.0-rc.1", 0, 1, 2, 0)
call expect_split("1.2.0.dev1", 0, 1, 2, 0)
call expect_split("1.2.3.post1", 0, 1, 2, 3)
call expect_split("1.2.0.rc1", 0, 1, 2, 0)
call expect_split("1.2.0rc.1", 0, 1, 2, 0)
call expect_split("1.2.3+gabcdef", 0, 1, 2, 3)
call expect_split("1.2.3.4", 1, 0, 0, 0)
call expect_split("1.2.3&x", 1, 0, 0, 0)
call expect_split("1.100.0", 1, 0, 0, 0)
call expect_split("1.2.100", 1, 0, 0, 0)
call expect_split("100.0.0", 1, 0, 0, 0)
call expect_split("2147483648", 1, 0, 0, 0)
call expect_split("1.2.2147483648", 1, 0, 0, 0)
call expect_split("4294967296", 1, 0, 0, 0)
call expect_split("999999999999", 1, 0, 0, 0)

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
call split3(trim(quoted), major, minor, patch, ios)
if (ios /= 0) error stop "parse fpm.toml version"
if (major /= gjp_version_major .or. minor /= gjp_version_minor .or. &
    patch /= gjp_version_patch) error stop "components /= fpm.toml"
write (triple, '(I0,".",I0,".",I0)') gjp_version_major, gjp_version_minor, &
    gjp_version_patch
if (trim(quoted) /= trim(triple)) error stop "fpm.toml is not the release triple"
if (gjp_version_major < 0 .or. gjp_version_major > 99 .or. &
    gjp_version_minor < 0 .or. gjp_version_minor > 99 .or. &
    gjp_version_patch < 0 .or. gjp_version_patch > 99) error stop "bounds"
if (gjp_version_compact /= &
    10000 * gjp_version_major + 100 * gjp_version_minor + gjp_version_patch) &
    error stop "compact"
print "(A)", "fpm version agrees with gjp_version"

contains

subroutine expect_split(text, want_ios, want_major, want_minor, want_patch)
    character(len=*), intent(in) :: text
    integer, intent(in) :: want_ios, want_major, want_minor, want_patch
    integer :: got_major, got_minor, got_patch, got_ios
    call split3(text, got_major, got_minor, got_patch, got_ios)
    if (got_ios /= want_ios) error stop "split ios"
    if (want_ios /= 0) return
    if (got_major /= want_major .or. got_minor /= want_minor .or. &
        got_patch /= want_patch) error stop "split values"
end subroutine expect_split

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
    integer :: cut, used, npart, value
    core = text
    cut = scan(core, "+")
    if (cut > 0) then
        if (cut == 1 .or. .not. meta_ok(core(cut + 1:))) then
            ios = 1
            return
        end if
        core = core(1:cut - 1)
    end if
    cut = scan(core, "-")
    if (cut > 0) then
        if (cut == 1 .or. .not. meta_ok(core(cut + 1:))) then
            ios = 1
            return
        end if
        core = core(1:cut - 1)
    end if
    major = 0
    minor = 0
    patch = 0
    npart = 0
    do
        call take_int(core, value, used, ios)
        if (ios /= 0) return
        if (value < 0 .or. value > 99) then
            ios = 1
            return
        end if
        npart = npart + 1
        if (npart > 3) then
            ios = 1
            return
        end if
        if (npart == 1) major = value
        if (npart == 2) minor = value
        if (npart == 3) patch = value
        if (used > len_trim(core)) then
            ios = 0
            return
        end if
        if (suffix_ok(core(used:))) then
            ios = 0
            return
        end if
        if (core(used:used) /= ".") then
            ios = 1
            return
        end if
        core = core(used + 1:)
    end do
end subroutine split3

subroutine take_int(text, value, used, ios)
    character(len=*), intent(in) :: text
    integer, intent(out) :: value, used, ios
    integer :: k, digit
    if (len_trim(text) < 1 .or. text(1:1) < "0" .or. text(1:1) > "9") then
        ios = 1
        value = 0
        used = 1
        return
    end if
    value = 0
    k = 1
    do while (k <= len_trim(text) .and. text(k:k) >= "0" .and. text(k:k) <= "9")
        digit = ichar(text(k:k)) - ichar("0")
        if (value > (huge(value) - digit) / 10) then
            ios = 1
            value = 0
            used = k
            return
        end if
        value = value * 10 + digit
        k = k + 1
    end do
    used = k
    ios = 0
end subroutine take_int

function meta_ok(text) result(ok)
    character(len=*), intent(in) :: text
    logical :: ok
    integer :: i, n
    character(len=1) :: ch
    n = len_trim(text)
    ok = n >= 1
    if (.not. ok) return
    do i = 1, n
        ch = text(i:i)
        if (.not. ((ch >= "0" .and. ch <= "9") .or. &
                   (ch >= "A" .and. ch <= "Z") .or. &
                   (ch >= "a" .and. ch <= "z") .or. &
                   ch == "." .or. ch == "-")) then
            ok = .false.
            return
        end if
    end do
end function meta_ok

function suffix_ok(tail) result(ok)
    character(len=*), intent(in) :: tail
    logical :: ok
    integer :: i, n, wlen
    character(len=8) :: word
    character(len=1) :: ch
    n = len_trim(tail)
    i = 1
    ok = .false.
    if (n < 1) return
    do
        if (i > n) return
        ch = tail(i:i)
        if (ch == "." .or. ch == "_") then
            i = i + 1
            if (i > n) then
                ok = .false.
                return
            end if
        end if
        wlen = 0
        word = ""
        do while (i <= n)
            ch = tail(i:i)
            if (.not. ((ch >= "A" .and. ch <= "Z") .or. (ch >= "a" .and. ch <= "z"))) exit
            wlen = wlen + 1
            if (wlen > len(word)) then
                ok = .false.
                return
            end if
            word(wlen:wlen) = lower_c(ch)
            i = i + 1
        end do
        if (wlen < 1 .or. .not. known_tag(word(1:wlen))) then
            ok = .false.
            return
        end if
        if (i <= n .and. (tail(i:i) == "." .or. tail(i:i) == "_")) i = i + 1
        do while (i <= n .and. tail(i:i) >= "0" .and. tail(i:i) <= "9")
            i = i + 1
        end do
        ok = .true.
        if (i > n) return
    end do
end function suffix_ok

function known_tag(word) result(ok)
    character(len=*), intent(in) :: word
    logical :: ok
    select case (word)
    case ("alpha", "beta", "preview", "pre", "post", "rev", "rc", "dev")
        ok = .true.
    case ("a", "b", "c", "r")
        ok = .true.
    case default
        ok = .false.
    end select
end function known_tag

function lower_c(ch) result(out)
    character(len=1), intent(in) :: ch
    character(len=1) :: out
    if (ch >= "A" .and. ch <= "Z") then
        out = char(ichar(ch) - ichar("A") + ichar("a"))
    else
        out = ch
    end if
end function lower_c

end program test_gjp_version_fpm
