//! Locate libgjp_cinterp (meson product of GaussJacobiQuad).
//!
//! Search order:
//! 1. `GJP_CINTERP` — path to the shared library file
//! 2. `GJP_CINTERP_DIR` / `GJP_PREFIX` — directory containing libgjp_cinterp
//! 3. pkg-config `gjp_cinterp` if available
//! 4. Default linker search (user must pass `-L` via RUSTFLAGS)

use std::env;
use std::path::PathBuf;

fn main() {
    println!("cargo:rerun-if-env-changed=GJP_CINTERP");
    println!("cargo:rerun-if-env-changed=GJP_CINTERP_DIR");
    println!("cargo:rerun-if-env-changed=GJP_PREFIX");
    println!("cargo:rerun-if-env-changed=PKG_CONFIG_PATH");

    if let Ok(libfile) = env::var("GJP_CINTERP") {
        let p = PathBuf::from(&libfile);
        if let Some(dir) = p.parent() {
            println!("cargo:rustc-link-search=native={}", dir.display());
        }
        // libgjp_cinterp.so → link name gjp_cinterp
        println!("cargo:rustc-link-lib=dylib=gjp_cinterp");
        // Fortran runtime often needed when the .so was not linked fully
        link_fortran_runtime();
        return;
    }

    for key in ["GJP_CINTERP_DIR", "GJP_PREFIX"] {
        if let Ok(dir) = env::var(key) {
            let d = PathBuf::from(&dir);
            let libdir = if d.join("lib").is_dir() {
                d.join("lib")
            } else {
                d
            };
            println!("cargo:rustc-link-search=native={}", libdir.display());
            println!("cargo:rustc-link-lib=dylib=gjp_cinterp");
            link_fortran_runtime();
            return;
        }
    }

    // pkg-config
    if std::process::Command::new("pkg-config")
        .args(["--exists", "gjp_cinterp"])
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
    {
        let libs = std::process::Command::new("pkg-config")
            .args(["--libs", "gjp_cinterp"])
            .output()
            .expect("pkg-config --libs");
        let s = String::from_utf8_lossy(&libs.stdout);
        for tok in s.split_whitespace() {
            if let Some(dir) = tok.strip_prefix("-L") {
                println!("cargo:rustc-link-search=native={dir}");
            } else if let Some(lib) = tok.strip_prefix("-l") {
                println!("cargo:rustc-link-lib={lib}");
            }
        }
        return;
    }

    // Fallback: hope the linker finds it (CI sets LIBRARY_PATH)
    println!("cargo:rustc-link-lib=dylib=gjp_cinterp");
    link_fortran_runtime();
}

fn link_fortran_runtime() {
    // gfortran runtime; harmless if already linked into the shared lib
    println!("cargo:rustc-link-lib=dylib=gfortran");
    println!("cargo:rustc-link-lib=dylib=m");
}
