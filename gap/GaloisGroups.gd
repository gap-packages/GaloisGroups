#! @Chapter GaloisGroups
#
#! @Section Main function

#! @Arguments p
#! @Returns a list
#! @Description
#!  Returns information about the Galois group of the polynomial
#!  <A>p</A>. Currently, the output format is a list of length 2
#!  where the first element is the index of the group in the
#!  transitive group table.
#! @BeginExampleSession
#! gap> LoadPackage("GaloisGroups");
#! gap> Galois(x^5-2);
#! [ 3, (3,4) ]
#! @EndExampleSession
DeclareGlobalFunction("Galois");
