{
  description = "GaussJacobiQuad multi-image CAF toolchain (gfortran + MPICH + OpenCoarrays)";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
    }:
    flake-utils.lib.eachSystem [ "x86_64-linux" "aarch64-linux" ] (
      system:
      let
        pkgs = import nixpkgs {
          inherit system;
          config.allowUnfree = false;
        };

        opencoarrays = pkgs.callPackage ./opencoarrays.nix { };
        # nixpkgs `fpm` is the Ruby packaging tool — use fortran-lang fpm
        fortran-fpm = pkgs.callPackage ./fortran-fpm.nix { };

        # Toolchain contents for the CI/dev image
        cafTools = pkgs.buildEnv {
          name = "gjp-caf-tools";
          paths = with pkgs; [
            bash
            coreutils
            gnugrep
            gnused
            gawk
            findutils
            gnumake
            which
            git
            gfortran
            gcc
            binutils
            pkg-config
            mpich
            mpich.bin
            mpich.dev
            opencoarrays
            fortran-fpm
            openblas
            lapack
            python3
            python3Packages.pytest
            python3Packages.numpy
            cacert
          ];
          pathsToLink = [
            "/bin"
            "/lib"
            "/include"
            "/share"
            "/etc"
          ];
        };

        # OCI image via nixpkgs dockerTools (load with: docker load < result)
        ociImage = pkgs.dockerTools.buildLayeredImage {
          name = "ghcr.io/haozeke/gaussjacobiquad-caf";
          tag = "latest";
          created = "now";
          contents = [ cafTools ];
          config = {
            Cmd = [ "${pkgs.bash}/bin/bash" ];
            WorkingDir = "/work";
            Env = [
              "PATH=${cafTools}/bin"
              "LIBRARY_PATH=${cafTools}/lib"
              "LD_LIBRARY_PATH=${cafTools}/lib"
              "PKG_CONFIG_PATH=${cafTools}/lib/pkgconfig"
              "SSL_CERT_FILE=${pkgs.cacert}/etc/ssl/certs/ca-bundle.crt"
              # fpm / gfortran CAF defaults used by CI scripts
              "GJP_CAF_FCFLAGS=-fcoarray=lib"
              "GJP_CAF_LDFLAGS=-lcaf_mpi -lmpifort -lmpi"
            ];
            Volumes = {
              "/work" = { };
            };
          };
          fakeRootCommands = ''
            mkdir -p ./work
          '';
        };

        # Streamed image for large CI (pipe into docker load)
        ociStream = pkgs.dockerTools.streamLayeredImage {
          name = "ghcr.io/haozeke/gaussjacobiquad-caf";
          tag = "latest";
          contents = [ cafTools ];
          config = {
            Cmd = [ "${pkgs.bash}/bin/bash" ];
            WorkingDir = "/work";
            Env = [
              "PATH=${cafTools}/bin"
              "LIBRARY_PATH=${cafTools}/lib"
              "LD_LIBRARY_PATH=${cafTools}/lib"
              "PKG_CONFIG_PATH=${cafTools}/lib/pkgconfig"
              "SSL_CERT_FILE=${pkgs.cacert}/etc/ssl/certs/ca-bundle.crt"
              "GJP_CAF_FCFLAGS=-fcoarray=lib"
              "GJP_CAF_LDFLAGS=-lcaf_mpi -lmpifort -lmpi"
            ];
          };
        };

      in
      {
        packages = {
          default = opencoarrays;
          opencoarrays = opencoarrays;
          fortran-fpm = fortran-fpm;
          tools = cafTools;
          oci = ociImage;
          oci-stream = ociStream;
        };

        devShells.default = pkgs.mkShell {
          packages = [
            pkgs.gfortran
            pkgs.mpich
            pkgs.mpich.bin
            pkgs.mpich.dev
            opencoarrays
            fortran-fpm
            pkgs.openblas
            pkgs.pkg-config
            pkgs.python3
            pkgs.python3Packages.pytest
            pkgs.python3Packages.numpy
          ];
          shellHook = ''
            export LIBRARY_PATH="${opencoarrays}/lib:${pkgs.mpich}/lib:${pkgs.openblas}/lib''${LIBRARY_PATH:+:$LIBRARY_PATH}"
            export LD_LIBRARY_PATH="${opencoarrays}/lib:${pkgs.mpich}/lib:${pkgs.openblas}/lib''${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
            export PATH="${opencoarrays}/bin:${pkgs.mpich.bin}/bin:${fortran-fpm}/bin:$PATH"
            export GJP_CAF_FCFLAGS="-fcoarray=lib"
            export GJP_CAF_LDFLAGS="-L${opencoarrays}/lib -lcaf_mpi -L${pkgs.mpich}/lib -lmpifort -lmpi"
            echo "GJP CAF shell: gfortran+MPICH+OpenCoarrays+fortran-fpm (multi-image via cafrun)"
          '';
        };

        checks.opencoarrays-probe =
          pkgs.runCommand "opencoarrays-probe"
            {
              nativeBuildInputs = [
                pkgs.gfortran
                pkgs.mpich
                pkgs.mpich.bin
                opencoarrays
              ];
            }
            ''
              cat > probe.f90 <<'EOF'
              program probe
                print *, "me=", this_image(), "nimg=", num_images()
              end program
              EOF
              export PATH="${opencoarrays}/bin:${pkgs.mpich.bin}/bin:$PATH"
              export LD_LIBRARY_PATH="${opencoarrays}/lib:${pkgs.mpich}/lib"
              gfortran -fcoarray=lib probe.f90 -o probe \
                -L${opencoarrays}/lib -lcaf_mpi \
                -L${pkgs.mpich}/lib -lmpifort -lmpi
              # MPICH 5 may abort on finalize after a good multi-image run
              set +e
              cafrun -n 2 ./probe > $out 2>&1
              set -e
              cat $out
              grep -q "nimg=.*2" $out
            '';
      }
    );
}
