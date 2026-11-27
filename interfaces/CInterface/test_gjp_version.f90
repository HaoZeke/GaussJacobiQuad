program test_gjp_version
use gjp_version, only: gjp_version_major, gjp_version_minor, gjp_version_patch, &
    gjp_version_string, gjp_version_compact, gjp_get_version
implicit none
character(len=64) :: buf
integer :: major, minor, patch, got_major, got_minor, got_patch, compact
character(len=:), allocatable :: got_string

if (command_argument_count() /= 4) error stop "expected major minor patch string"
call get_command_argument(1, buf)
read(buf, *) major
call get_command_argument(2, buf)
read(buf, *) minor
call get_command_argument(3, buf)
read(buf, *) patch
call get_command_argument(4, buf)

if (gjp_version_major /= major) error stop "major"
if (gjp_version_minor /= minor) error stop "minor"
if (gjp_version_patch /= patch) error stop "patch"
if (gjp_version_string /= trim(buf)) error stop "string"
compact = 10000 * major + 100 * minor + patch
if (gjp_version_compact /= compact) error stop "compact"

call gjp_get_version(major=got_major, minor=got_minor, patch=got_patch, string=got_string)
if (got_major /= major .or. got_minor /= minor .or. got_patch /= patch) error stop "getter"
if (got_string /= gjp_version_string) error stop "getter string"
print '(A,A,A,I0)', 'gjp version ', gjp_version_string, ' compact ', gjp_version_compact
end program test_gjp_version
