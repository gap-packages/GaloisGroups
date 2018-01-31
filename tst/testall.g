#
# GaloisGroups: Computing Galois Groups using GAP and PARI
#
# This file runs package tests. It is also referenced in the package
# metadata in PackageInfo.g.
#
LoadPackage( "GaloisGroups" );

TestDirectory(DirectoriesPackageLibrary( "GaloisGroups", "tst" ),
  rec(exitGAP := true));

FORCE_QUIT_GAP(1); # if we ever get here, there was an error
