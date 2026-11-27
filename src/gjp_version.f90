!> Version constants for the library.
!> Checked in so fpm, which has no configure_file, builds the same values
!> meson checks against the project version. scripts/set_version.sh updates
!> the three integers and the string together with fpm.toml and meson.build.
!>
!> Compact form matches tblite: 10000 * major + 100 * minor + patch.
!> Each component must be in 0..99. A larger value collides
!> (1.100.0 == 2.0.0, 1.2.100 == 1.3.0) and is rejected here at compile time
!> and by meson at configure time.
module gjp_version
implicit none
private
public :: gjp_version_major, gjp_version_minor, gjp_version_patch
public :: gjp_version_string, gjp_version_compact
public :: gjp_get_version

integer, parameter :: gjp_version_major = 0
integer, parameter :: gjp_version_minor = 2
integer, parameter :: gjp_version_patch = 6
character(len=*), parameter :: gjp_version_string = "0.2.6"
integer, parameter :: gjp_version_in_range = 1 / &
                      merge(1, 0, gjp_version_major >= 0 .and. gjp_version_major <= 99 .and. &
                            gjp_version_minor >= 0 .and. gjp_version_minor <= 99 .and. &
                            gjp_version_patch >= 0 .and. gjp_version_patch <= 99)
integer, parameter :: gjp_version_compact = &
                      (10000 * gjp_version_major + 100 * gjp_version_minor + gjp_version_patch) * &
                      gjp_version_in_range

contains

!> Optional arguments return the version pieces.
subroutine gjp_get_version(major, minor, patch, string)
    integer, intent(out), optional :: major
    integer, intent(out), optional :: minor
    integer, intent(out), optional :: patch
    character(len=:), allocatable, intent(out), optional :: string
    if (present(major)) major = gjp_version_major
    if (present(minor)) minor = gjp_version_minor
    if (present(patch)) patch = gjp_version_patch
    if (present(string)) string = gjp_version_string
end subroutine gjp_get_version

end module gjp_version
