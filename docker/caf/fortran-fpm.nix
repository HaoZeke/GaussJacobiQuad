# Fortran Package Manager (fortran-lang/fpm) — NOT the Ruby fpm gem in nixpkgs.
{
  lib,
  stdenv,
  fetchurl,
}:

stdenv.mkDerivation rec {
  pname = "fortran-fpm";
  version = "0.13.0";

  src = fetchurl {
    url = "https://github.com/fortran-lang/fpm/releases/download/v${version}/fpm-${version}-linux-x86_64-gcc-12";
    hash = "sha256-LlCoZJrEE2uc9rY8XRsBrfJsAu7cGPMuFiYcDOZV1G0=";
  };

  dontUnpack = true;

  installPhase = ''
    runHook preInstall
    mkdir -p $out/bin
    cp $src $out/bin/fpm
    chmod +x $out/bin/fpm
    runHook postInstall
  '';

  meta = with lib; {
    description = "Fortran Package Manager (fortran-lang)";
    homepage = "https://fpm.fortran-lang.org/";
    license = licenses.mit;
    platforms = [ "x86_64-linux" ];
    mainProgram = "fpm";
  };
}
