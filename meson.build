project('pypotlib', ['c', 'fortran'],
       version: '0.0.13',
  default_options : ['warning_level=2', 'cpp_std=c++17'])
# IMPORTANT!! warning_level=3 passes -fimplicit-none
# Many of the older Fortran codes need implicit typing

_args = [] # Extra arguments
_deps = [] # Dependencies

# ---------------------- Dependencies
_deps += dependency('lapack')

# Languages
cc = meson.get_compiler('c')

# ---------------------- Library

gjplib = library('GaussJacobiQaud',
              sources: [
                'src/gjp_constants.f90',
                'src/gjp_types.f90',
                'src/gjp_lapack.f90',
                'src/gjp_rec.f90',
                'src/gjp_gw.f90',
                'src/GaussJacobiQuad.f90',
              ],
              dependencies: _deps,
              install: true,
)

# ---------------------- Applications

gjp_quad = executable('gjp_quad',
                      sources: [
                        'app/gjp_quad.f90',
                      ],
                      link_with: [gjplib],
                      install: true,
)

gjp_quad_rec = executable('gjp_quad_rec',
                      sources: [
                        'app/gjp_quad_rec.f90',
                      ],
                      link_with: [gjplib],
                      install: true,
)

gjp_quad_gw = executable('gjp_quad_gw',
                      sources: [
                        'app/gjp_quad_gw.f90',
                      ],
                      link_with: [gjplib],
                      install: true,
)

# ---------------------- Interfaces

gjp_cinterp = library('gjp_cinterp',
                      sources: [
                        'src/GaussJacobiQuadCCompat.f90',
                        'interfaces/CInterface/GaussJacobiQuadCInterp.h',
                      ],
                      link_with: [gjplib],
                      install: true,
)

c_cli_gjpq = executable('c_cli_gjpq',
                       sources: [
                         'interfaces/CInterface/c_cli_gjpq.c',
                       ],
                       link_with: [gjp_cinterp],
                       install: true,
)
