{pkgs ? import <nixpkgs> {} }:

with pkgs; mkShell.override { stdenv = stdenvNoCC; } {
  nativeBuildInputs = [ pkgconfig ninja doxygen ];
  buildInputs = [ gsl eigen];
}
