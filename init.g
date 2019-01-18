#
# GaloisGroups: Computing Galois Groups using GAP and PARI
#
# Reading the declaration part of the package.
#
#_PATH_SO:=Filename(DirectoriesPackagePrograms("GaloisGroups"), "GaloisGroups.so");
#if _PATH_SO <> fail then
#    LoadDynamicModule(_PATH_SO);
#fi;
#Unbind(_PATH_SO);

ReadPackage( "GaloisGroups", "gap/short_cosets.gd");
ReadPackage( "GaloisGroups", "gap/transitive_lattice.gd");
ReadPackage( "GaloisGroups", "gap/relative_invariants.gd");
ReadPackage( "GaloisGroups", "gap/GaloisGroups.gd");
