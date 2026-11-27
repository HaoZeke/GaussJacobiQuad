project('gaussjacobiquad', ['c', 'fortran'],
       version: '0.2.6',
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
# Prefer FlexiBLAS: one wrapper over BLAS, LAPACK and LAPACKE whose
# backend is chosen at run time, which is how the distributions and the
# clusters that carry it expect a consumer to link. OpenBLAS comes next,
# then a plain lapack pkg-config file, then bare -l flags.
_blas_backend = ''
_fb = dependency('flexiblas', required: false)
if not _fb.found()
  _fb = cc.find_library('flexiblas', required: false)
endif
if _fb.found()
  _deps += _fb
  _blas_backend = 'flexiblas'
else
  _ob = dependency('openblas', required: false)
  if _ob.found()
    _deps += _ob
    _blas_backend = 'openblas'
  else
    # Honor LIBRARY_PATH / linker search path for OpenBLAS
    _oblib = cc.find_library('openblas', required: false)
    if _oblib.found()
      _deps += _oblib
      _blas_backend = 'openblas'
    else
      _lapack = dependency('lapack', required: false)
      if _lapack.found()
        _deps += _lapack
        _blas_backend = 'lapack'
      else
        _deps += cc.find_library('lapack', required: true)
        _deps += cc.find_library('blas', required: false)
        _blas_backend = 'lapack'
      endif
    endif
  endif
endif

# gjp_lapack.f90 and gjp_imtqlx_dc.f90 call LAPACKE_dsteqr and
# LAPACKE_dsyev unconditionally, so this is not optional. Probing it
# with required: false moved the failure from configure to link, and
# only where the link forbids undefined symbols: a toolchain without
# -Wl,--no-undefined produced a library with LAPACKE_* left dangling,
# which fails at load in whatever program picks it up. FlexiBLAS carries
# LAPACKE itself, so it needs no second library.
if _blas_backend != 'flexiblas'
  _lapacke = dependency('lapacke', required: false)
  if _lapacke.found()
    _deps += _lapacke
  else
    _deps += cc.find_library('lapacke', required: true)
  endif
endif

# C executables still link with cc even when link_language is fortran.
# libgfortran then needs libquadmath from the gfortran private libdir.
fs = import('fs')
_qm = run_command(_fc, '-print-file-name=libquadmath.so', check: false)
if _qm.returncode() == 0
  _qm_path = _qm.stdout().strip()
  if fs.is_absolute(_qm_path) and fs.exists(_qm_path)
    _deps += cc.find_library('quadmath', dirs: [fs.parent(_qm_path)], required: false)
  endif
endif

summary({'BLAS backend': _blas_backend}, section: 'Dependencies')

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
  # C header is generated. The Fortran module is src/gjp_version.f90 so fpm
  # (no configure_file) and meson compile the same constants.
  # Packing matches tblite: 10000*major + 100*minor + patch, each in 0..99.
  _py = find_program('python3', required: true)
  _gjp_parsed = run_command(
    _py,
    files('scripts/parse_gjp_version.py'),
    '--check-fortran',
    files('src/gjp_version.f90'),
    meson.project_version(),
    check: true,
  )
  _gjp_nums = _gjp_parsed.stdout().strip().split()
  _gjp_major = _gjp_nums[0]
  _gjp_minor = _gjp_nums[1]
  _gjp_patch = _gjp_nums[2]
  _gjp_ver_data = configuration_data()
  _gjp_ver_data.set('GJP_VERSION_MAJOR', _gjp_major)
  _gjp_ver_data.set('GJP_VERSION_MINOR', _gjp_minor)
  _gjp_ver_data.set('GJP_VERSION_PATCH', _gjp_patch)
  _gjp_ver_data.set('GJP_VERSION_STRING', meson.project_version())
  gjp_version_h = configure_file(
    input: 'interfaces/CInterface/gjp_version.h.in',
    output: 'gjp_version.h',
    configuration: _gjp_ver_data,
    install: true,
    install_dir: get_option('includedir'),
  )
  # Compiler-specific, same layout tblite uses for installed .mod files.
  _mod_sub = meson.project_name() + '/' + _fc.get_id() + '-' + _fc.version()

  gjplib = library('GaussJacobiQuad',
                sources: [_gjp_f_sources, files('src/gjp_version.f90')],
                fortran_args: _args,
                dependencies: _deps,
                install: true,
  )

  gaussjacobiquad_dep = declare_dependency(
    link_with: gjplib,
    # .mod files are written to the library private build dir. Export that
    # directory so a subproject consumer can `use GaussJacobiQuad`.
    include_directories: gjplib.private_dir_include(),
    dependencies: _deps,
  )

  test_dep_mod = executable('test_dep_mod',
    sources: files('test/test_dep_mod.f90'),
    dependencies: gaussjacobiquad_dep,
    fortran_args: _args,
    link_language: 'fortran',
  )
  test('dep_mod', test_dep_mod)

  meson.add_install_script(
    _py,
    files('scripts/install_gjp_mod.py'),
    get_option('includedir') + '/' + _mod_sub,
  )

  # ---------------------- Applications

  gjp_quad = executable('gjp_quad',
                        sources: [
                          'app/gjp_quad.f90',
                        ],
                        fortran_args: _args,
                        link_with: [gjplib],
                        dependencies: _deps,
                        link_language: 'fortran',
                        install: true,
  )

  gjp_quad_rec = executable('gjp_quad_rec',
                        sources: [
                          'app/gjp_quad_rec.f90',
                        ],
                        fortran_args: _args,
                        link_with: [gjplib],
                        dependencies: _deps,
                        link_language: 'fortran',
                        install: true,
  )

  gjp_quad_gw = executable('gjp_quad_gw',
                        sources: [
                          'app/gjp_quad_gw.f90',
                        ],
                        fortran_args: _args,
                        link_with: [gjplib],
                        dependencies: _deps,
                        link_language: 'fortran',
                        install: true,
  )

  gjp_quad_algo665 = executable('gjp_quad_algo665',
                        sources: [
                          'app/gjp_quad_algo665.f90',
                        ],
                        fortran_args: _args,
                        link_with: [gjplib],
                        dependencies: _deps,
                        link_language: 'fortran',
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
                        include_directories: [
                          include_directories('interfaces/CInterface'),
                          gjplib.private_dir_include(),
                        ],
                        fortran_args: _args,
                        link_with: [gjplib],
                        dependencies: _deps,
                        install: true,
  )

  # Header dir and the installed module dir. subdirs '.' is ${includedir}
  # (gjp_version.h). _mod_sub is the gjp_version.mod directory.
  pkg = import('pkgconfig')
  pkg.generate(
    gjp_cinterp,
    name: 'gjp_cinterp',
    description: 'GaussJacobiQuad C ABI (gauss_jacobi_rule_c)',
    version: meson.project_version(),
    filebase: 'gjp_cinterp',
    subdirs: ['.', _mod_sub],
    variables: {'moduledir': '${includedir}/' + _mod_sub},
  )

  c_cli_gjpq = executable('c_cli_gjpq',
                         sources: [
                           'interfaces/CInterface/c_cli_gjpq.c',
                         ],
                         include_directories: include_directories('interfaces/CInterface'),
                         link_with: [gjp_cinterp],
                         dependencies: _deps,
                         link_language: 'fortran',
                         install: true,
  )

  test_c_rule = executable('test_c_rule',
                          sources: [
                            'interfaces/CInterface/test_c_rule.c',
                          ],
                          include_directories: include_directories('interfaces/CInterface'),
                          link_with: [gjp_cinterp],
                          dependencies: _deps,
                          link_language: 'fortran',
                          install: false,
  )
  test('c_rule', test_c_rule)

  test_c_version = executable('test_c_version',
    sources: ['interfaces/CInterface/test_c_version.c', gjp_version_h],
    c_args: [
      '-DEXPECT_MAJOR=' + _gjp_major,
      '-DEXPECT_MINOR=' + _gjp_minor,
      '-DEXPECT_PATCH=' + _gjp_patch,
      '-DEXPECT_STRING="' + meson.project_version() + '"',
    ],
  )
  test('c_version', test_c_version)

  test_gjp_version = executable('test_gjp_version',
    sources: 'interfaces/CInterface/test_gjp_version.f90',
    fortran_args: _args,
    link_with: [gjplib],
    dependencies: _deps,
    link_language: 'fortran',
  )
  test('fortran_version', test_gjp_version, args: [
    _gjp_major, _gjp_minor, _gjp_patch, meson.project_version(),
  ])

  test('version_parse', _py, args: [files('test/test_gjp_version_parse.py')])
  test('version_agree', _py, args: [files('test/test_gjp_version_agree.py')])
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
