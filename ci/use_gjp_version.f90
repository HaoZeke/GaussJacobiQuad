program use_gjp_version
use gjp_version, only: gjp_version_compact, gjp_version_string, gjp_get_version
implicit none
integer :: major, minor, patch
character(len=:), allocatable :: got
if (len(gjp_version_string) < 1) error stop "empty"
if (gjp_version_compact < 0) error stop "compact"
call gjp_get_version(major, minor, patch, got)
if (got /= gjp_version_string) error stop "string"
if (10000 * major + 100 * minor + patch /= gjp_version_compact) error stop "parts"
print*,gjp_version_string, gjp_version_compact
end program use_gjp_version
