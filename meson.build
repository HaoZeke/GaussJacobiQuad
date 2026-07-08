project('gaussjacobiquad', ['c', 'fortran'],
       version: '0.1.0',
  default_options : [
    'warning_level=everything',
    'buildtype=release',
  ])
# IMPORTANT!! warning_level=3 passes -fimplicit-none
# Many of the older Fortran codes need implicit typing

_args = [] # Extra arguments
_deps = [] # Dependencies

# Coarray Fortran for rec_caf (this_image / num_images / coarrays).
# Override with -Dfcoarray=lib when linking OpenCoarrays multi-image.
_fc = meson.get_compiler('fortran')
_caf_mode = get_option('fcoarray')
if _fc.get_id() == 'gcc'
  _args += '-fcoarray=' + _caf_mode
endif

# ---------------------- Dependencies
_deps += dependency('lapack')

# Languages
cc = meson.get_compiler('c')

# ---------------------- Library

gjplib = library('GaussJacobiQaud',
              sources: [
                'src/gjp_constants.f90',
                'src/gjp_types.f90',
                'src/gjp_rec.f90',
                'src/gjp_lapack.f90',
                'src/gjp_common.f90',
                'src/gjp_imtqlx.f90',
                'src/gjp_imtqlx_dc.f90',
                'src/gjp_gw.f90',
                'src/gjp_algo665.f90',
                'src/gjp_caf.f90',
                'src/GaussJacobiQuad.f90',
              ],
              fortran_args: _args,
              dependencies: _deps,
              install: true,
)

# ---------------------- Applications

gjp_quad = executable('gjp_quad',
                      sources: [
                        'app/gjp_quad.f90',
                      ],
                      fortran_args: _args,
                      link_with: [gjplib],
                      install: true,
)

gjp_quad_rec = executable('gjp_quad_rec',
                      sources: [
                        'app/gjp_quad_rec.f90',
                      ],
                      fortran_args: _args,
                      link_with: [gjplib],
                      install: true,
)

gjp_quad_gw = executable('gjp_quad_gw',
                      sources: [
                        'app/gjp_quad_gw.f90',
                      ],
                      fortran_args: _args,
                      link_with: [gjplib],
                      install: true,
)

gjp_quad_algo665 = executable('gjp_quad_algo665',
                      sources: [
                        'app/gjp_quad_algo665.f90',
                      ],
                      fortran_args: _args,
                      link_with: [gjplib],
                      install: true,
)

# ---------------------- Interfaces

gjp_cinterp = library('gjp_cinterp',
                      sources: [
                        'src/GaussJacobiQuadCCompat.f90',
                        'interfaces/CInterface/GaussJacobiQuadCInterp.c',
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
