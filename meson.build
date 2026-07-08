project('gaussjacobiquad', ['c', 'fortran'],
       version: '0.2.3',
  default_options : [
    'warning_level=everything',
    'buildtype=release',
  ])
# IMPORTANT!! warning_level=3 passes -fimplicit-none
# Many of the older Fortran codes need implicit typing

_args = [] # Extra arguments
_deps = [] # Dependencies
_py_pkg = get_option('python_package')

# Languages / compilers
cc = meson.get_compiler('c')
_fc = meson.get_compiler('fortran')

# Coarray Fortran for rec_caf (this_image / num_images / coarrays).
# Override with -Dfcoarray=lib when linking OpenCoarrays multi-image.
_caf_mode = get_option('fcoarray')
if _fc.get_id() == 'gcc'
  _args += '-fcoarray=' + _caf_mode
endif

# ---------------------- Dependencies
# Prefer OpenBLAS (BLAS+LAPACK). Fall back to lapack pkg-config / -l flags.
_ob = dependency('openblas', required: false)
if _ob.found()
  _deps += _ob
else
  # Honor LIBRARY_PATH / linker search path for OpenBLAS
  _oblib = cc.find_library('openblas', required: false)
  if _oblib.found()
    _deps += _oblib
  else
    _lapack = dependency('lapack', required: false)
    if _lapack.found()
      _deps += _lapack
    else
      _deps += cc.find_library('lapack', required: true)
      _deps += cc.find_library('blas', required: false)
    endif
  endif
endif

# Fortran kernel sources shared by the static/shared lib and the Python extension
_gjp_f_sources = files(
  'src/gjp_constants.f90',
  'src/gjp_types.f90',
  'src/gjp_rec.f90',
  'src/gjp_lapack.f90',
  'src/gjp_common.f90',
  'src/gjp_imtqlx.f90',
  'src/gjp_imtqlx_dc.f90',
  'src/gjp_gw.f90',
  'src/gjp_algo665.f90',
  'src/gjp_sturm.f90',
  'src/gjp_bogaert.f90',
  'src/gjp_glr.f90',
  'src/gjp_auto.f90',
  'src/gjp_caf.f90',
  'src/GaussJacobiQuad.f90',
)

# ---------------------- Native library / apps / C ABI (skip install when packaging wheel)
if not _py_pkg
  gjplib = library('GaussJacobiQaud',
                sources: _gjp_f_sources,
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

  # ---------------------- Interfaces (C/C++ ABI)

  # Install public C header
  install_headers('interfaces/CInterface/GaussJacobiQuadCInterp.h')

  gjp_cinterp = library('gjp_cinterp',
                        sources: [
                          'src/GaussJacobiQuadCCompat.f90',
                          'interfaces/CInterface/GaussJacobiQuadCInterp.c',
                        ],
                        include_directories: include_directories('interfaces/CInterface'),
                        fortran_args: _args,
                        link_with: [gjplib],
                        dependencies: _deps,
                        install: true,
  )

  # pkg-config for Rust/Python consumers
  pkg = import('pkgconfig')
  pkg.generate(
    gjp_cinterp,
    name: 'gjp_cinterp',
    description: 'GaussJacobiQuad C ABI (gauss_jacobi_rule_c)',
    version: meson.project_version(),
    filebase: 'gjp_cinterp',
    subdirs: '.',
  )

  c_cli_gjpq = executable('c_cli_gjpq',
                         sources: [
                           'interfaces/CInterface/c_cli_gjpq.c',
                         ],
                         include_directories: include_directories('interfaces/CInterface'),
                         link_with: [gjp_cinterp],
                         install: true,
  )

  test_c_rule = executable('test_c_rule',
                          sources: [
                            'interfaces/CInterface/test_c_rule.c',
                          ],
                          include_directories: include_directories('interfaces/CInterface'),
                          link_with: [gjp_cinterp],
                          install: false,
  )
  test('c_rule', test_c_rule)
endif

# ---------------------- CPython extension (meson-python / PyPI)
# Self-contained module: links Fortran kernels + C ABI into one .so (no GJP_CINTERP).
if _py_pkg
  py = import('python').find_installation(pure: false)
  py_dep = py.dependency()
  py_deps = [py_dep] + _deps
  # gfortran runtime when the C linker drives the shared object
  _gfrt = cc.find_library('gfortran', required: false)
  if _gfrt.found()
    py_deps += _gfrt
  endif
  _qm = cc.find_library('quadmath', required: false)
  if _qm.found()
    py_deps += _qm
  endif
  _m = cc.find_library('m', required: false)
  if _m.found()
    py_deps += _m
  endif

  _ext_sources = files(
    'python/src/gauss_jacobi_quad/_core.c',
    'interfaces/CInterface/GaussJacobiQuadCInterp.c',
    'src/GaussJacobiQuadCCompat.f90',
  ) + _gjp_f_sources

  # Free-threaded CPython is a different ABI: full C API + Py_mod_gil, not abi3.
  _py_is_ft = run_command(
    py, '-c',
    'import sysconfig; print(1 if sysconfig.get_config_var("Py_GIL_DISABLED") else 0)',
    check: true,
  ).stdout().strip() == '1'

  if _py_is_ft
    py.extension_module(
      '_core',
      sources: _ext_sources,
      include_directories: include_directories('interfaces/CInterface'),
      dependencies: py_deps,
      fortran_args: _args,
      c_args: ['-DGJP_FULL_API=1'],
      install: true,
      subdir: 'gauss_jacobi_quad',
      link_language: 'fortran',
    )
  else
    py.extension_module(
      '_core',
      sources: _ext_sources,
      include_directories: include_directories('interfaces/CInterface'),
      dependencies: py_deps,
      fortran_args: _args,
      # Stable ABI: one abi3 wheel for CPython >= 3.9 (see Py_LIMITED_API in _core.c)
      limited_api: '3.9',
      install: true,
      subdir: 'gauss_jacobi_quad',
      link_language: 'fortran',
    )
  endif


  py.install_sources(
    [
      'python/src/gauss_jacobi_quad/__init__.py',
      'python/src/gauss_jacobi_quad/cli.py',
    ],
    pure: false,
    subdir: 'gauss_jacobi_quad',
  )
endif
