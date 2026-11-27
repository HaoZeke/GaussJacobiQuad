program use_gjp_version
use gjp_version, only: gjp_version_compact, gjp_version_string
implicit none
if (len(gjp_version_string) < 1) error stop "empty"
if (gjp_version_compact < 0) error stop "compact"
print*,gjp_version_string, gjp_version_compact
end program use_gjp_version
