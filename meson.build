project('gaussjacobiquad', ['c', 'fortran'],
       version: '0.1.0',
       meson_version: '>=0.60.0',
  default_options : [
    'warning_level=everything',
    'buildtype=release',
  ])
# IMPORTANT!! warning_level=3 passes -fimplicit-none
# Many of the older Fortran codes need implicit typing

_deps = []
_deps += dependency('lapack')

incdir = get_option('includedir')
gjp_mod_install_dir = join_paths(incdir, 'gaussjacobiquad')

# ---------------------- Core library (Fortran)

gjplib = library('GaussJacobiQuad',
              sources: [
                'src/gjp_constants.f90',
                'src/gjp_types.f90',
                'src/gjp_rec.f90',
                'src/gjp_lapack.f90',
                'src/gjp_common.f90',
                'src/gjp_imtqlx.f90',
                'src/gjp_gw.f90',
                'src/gjp_algo665.f90',
                'src/GaussJacobiQuad.f90',
              ],
              dependencies: _deps,
              install: true,
              version: meson.project_version(),
              soversion: meson.project_version().split('.')[0],
)

# ---------------------- Applications

executable('gjp_quad',
           sources: ['app/gjp_quad.f90'],
           link_with: [gjplib],
           install: true,
)

executable('gjp_quad_rec',
           sources: ['app/gjp_quad_rec.f90'],
           link_with: [gjplib],
           install: true,
)

executable('gjp_quad_gw',
           sources: ['app/gjp_quad_gw.f90'],
           link_with: [gjplib],
           install: true,
)

executable('gjp_quad_algo665',
           sources: ['app/gjp_quad_algo665.f90'],
           link_with: [gjplib],
           install: true,
)

# ---------------------- C / ISO_C_BINDING interface

gjp_cinterp = library('gjp_cinterp',
                      sources: [
                        'src/GaussJacobiQuadCCompat.f90',
                        'interfaces/CInterface/GaussJacobiQuadCInterp.c',
                      ],
                      link_with: [gjplib],
                      install: true,
                      version: meson.project_version(),
                      soversion: meson.project_version().split('.')[0],
)

install_headers('interfaces/CInterface/GaussJacobiQuadCInterp.h',
                subdir: 'gaussjacobiquad')

executable('c_cli_gjpq',
           sources: ['interfaces/CInterface/c_cli_gjpq.c'],
           link_with: [gjp_cinterp],
           install: true,
)

# ---------------------- Dependency / pkg-config export

gaussjacobiquad_dep = declare_dependency(
  link_with: [gjplib, gjp_cinterp],
  include_directories: include_directories('interfaces/CInterface'),
  dependencies: _deps,
)
meson.override_dependency('gaussjacobiquad', gaussjacobiquad_dep)

pkg = import('pkgconfig')
pkg.generate(
  name: 'gaussjacobiquad',
  description: 'Gauss-Jacobi quadrature nodes and weights (Fortran + C interface)',
  version: meson.project_version(),
  libraries: [gjplib, gjp_cinterp],
  subdirs: 'gaussjacobiquad',
  filebase: 'gaussjacobiquad',
  requires: 'lapack',
  variables: [
    'moddir=${prefix}/' + gjp_mod_install_dir,
  ],
)

# Install Fortran .mod files produced during the build (for USE statements in dependents)
meson.add_install_script(
  'sh', '-c',
  'set -eu; ' +
  'moddir="$MESON_INSTALL_DESTDIR_PREFIX/' + gjp_mod_install_dir + '"; ' +
  'mkdir -p "$moddir"; ' +
  'find "$MESON_BUILD_ROOT" -name \'*.mod\' -print0 2>/dev/null | ' +
  'while IFS= read -r -d \'\' f; do ' +
  '  base=$(basename "$f"); ' +
  '  case "$base" in ' +
  '    gjp_*.mod|gaussjacobiquad.mod|gaussjacobiquadccompat.mod) install -m 644 "$f" "$moddir/" ;; ' +
  '  esac; ' +
  'done',
)
