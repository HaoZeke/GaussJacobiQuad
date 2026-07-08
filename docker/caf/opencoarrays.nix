# OpenCoarrays runtime (libcaf_mpi + caf/cafrun) against MPICH.
# gfortran ≥14 rejects OpenMPI for CAF; always pair with mpich.
{
  lib,
  stdenv,
  fetchFromGitHub,
  cmake,
  gfortran,
  mpich,
  perl,
  bash,
  pkg-config,
  writeShellScriptBin,
}:

let
  cafrunWrapper = writeShellScriptBin "cafrun" ''
    set -euo pipefail
    export PATH="${mpich.bin}/bin:$PATH"
    N=""
    ARGS=()
    while [[ $# -gt 0 ]]; do
      case "$1" in
        -n|-np|--nimages) N="$2"; shift 2 ;;
        *) ARGS+=("$1"); shift ;;
      esac
    done
    if [[ -z "$N" ]]; then
      echo "cafrun: usage: cafrun -n <images> <program> [args...]" >&2
      exit 2
    fi
    exec ${mpich.bin}/bin/mpiexec -n "$N" "''${ARGS[@]}"
  '';
in
stdenv.mkDerivation rec {
  pname = "opencoarrays";
  # 2.10.2 lacks gfortran ≥15 CAF ABI symbols (get_from_remote, register_accessor).
  # Pin a post-2.10.2 commit that implements them.
  version = "2.10.2-74a5d0a";

  src = fetchFromGitHub {
    owner = "sourceryinstitute";
    repo = "OpenCoarrays";
    rev = "74a5d0ac3f2d6ee5f985a13478a87f309a444b07";
    hash = "sha256-7ZntIRpvDN0HTnZgBGnOhJOX6GKeu5k0jv6Esuois1M=";
  };

  nativeBuildInputs = [
    cmake
    gfortran
    perl
    pkg-config
  ];

  buildInputs = [
    mpich
    mpich.dev
    mpich.bin
  ];

  cmakeFlags = [
    "-DCAF_ENABLE_FAILED_IMAGES=OFF"
    "-DCMAKE_C_COMPILER=${mpich.dev}/bin/mpicc"
    "-DCMAKE_Fortran_COMPILER=${mpich.dev}/bin/mpifort"
    "-DMPI_C_COMPILER=${mpich.dev}/bin/mpicc"
    "-DMPI_Fortran_COMPILER=${mpich.dev}/bin/mpifort"
    "-DMPIEXEC_EXECUTABLE=${mpich.bin}/bin/mpiexec"
    "-DMPI_HOME=${mpich.dev}"
  ];

  buildPhase = ''
    runHook preBuild
    cmake --build . --target caf_mpi caf_mpi_static -j$NIX_BUILD_CORES
    runHook postBuild
  '';

  dontCheckForBrokenSymlinks = true;

  installPhase = ''
    runHook preInstall
    mkdir -p $out/lib $out/bin $out/include

    # cmake may link libs directly into $out/lib during the build
    for f in $out/lib/libcaf_mpi.so* $out/lib/libcaf_mpi.a; do
      [ -e "$f" ] || true
    done
    find . -name 'libcaf_mpi.so*' -o -name 'libcaf_mpi.a' | while read -r f; do
      cp -a "$f" "$out/lib/" 2>/dev/null || true
    done

    find .. -name 'libcaf*.h' -exec cp -a {} $out/include/ \; 2>/dev/null || true

    cp ${cafrunWrapper}/bin/cafrun $out/bin/cafrun

    cat > $out/bin/caf <<EOF
    #!${bash}/bin/bash
    set -euo pipefail
    export PATH="${gfortran}/bin:${mpich.dev}/bin:${mpich.bin}/bin:\$PATH"
    export LIBRARY_PATH="$out/lib:${mpich}/lib\''${LIBRARY_PATH:+:\$LIBRARY_PATH}"
    export LD_LIBRARY_PATH="$out/lib:${mpich}/lib\''${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}"
    exec ${gfortran}/bin/gfortran -fcoarray=lib "\$@" \
      -L$out/lib -lcaf_mpi \
      -L${mpich}/lib -lmpifort -lmpi
    EOF
    sed -i 's/^    //' $out/bin/caf
    chmod +x $out/bin/caf $out/bin/cafrun

    if [ ! -e $out/lib/libcaf_mpi.so ] && ls $out/lib/libcaf_mpi.so.* >/dev/null 2>&1; then
      ln -sf "$(basename "$(ls -1 $out/lib/libcaf_mpi.so.* | head -1)")" $out/lib/libcaf_mpi.so
    fi

    rm -f $out/include/opencoarrays.mod
    rm -rf $out/include/OpenCoarrays-* 2>/dev/null || true

    test -e $out/lib/libcaf_mpi.so -o -e $out/lib/libcaf_mpi.a
    test -x $out/bin/cafrun
    runHook postInstall
  '';

  meta = with lib; {
    description = "OpenCoarrays multi-image CAF runtime (MPI) for gfortran";
    homepage = "https://github.com/sourceryinstitute/OpenCoarrays";
    license = licenses.bsd3;
    platforms = platforms.linux;
  };
}
