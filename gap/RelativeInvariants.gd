########################################################################
##
##  GaloisGroups package
##
##  Copyright 2018-2019
##    Bill Allombert <bill.allombert@math.u-bordeaux.fr>
##    Vincent Delecroix <vincent.delecroix@math.cnrs.fr>
##    Markus Pfeiffer <markus.pfeiffer@morphism.de>
##
## Licensed under the GPL 2 or later.
##
########################################################################

#! @Chapter RelativeInvariants
#
#! @Section Relative invariant for descents H &lt; G

#! @Arguments G H
#! @Returns pairs of integers
#! @Description
#!  Return the pair (n, dim) so that n is the minimal degree
#!  so that there exists a polynomial in Q[x1, ..., xd] of
#!  degree n that is <A>H</A> invariant but not G invariant. The
#!  second returned argument dim is the dimension of
#!  the algebra of <A>G</A>-invariant polynomials on the
#!  basis of <A>H</A>-invariant ones.
#! @BeginExampleSession
#! gap> T63 := TransitiveGroup(6,3);;
#! gap> T61 := TransitiveGroup(6,1);;
#! gap> RelativeInvariantMinimalDegree(T63, T61);
#! [ 3, 3 ]
#! @EndExampleSession
#! @BeginExampleSession
#! gap> T65 := TransitiveGroup(6,5);;
#! gap> T62 := TransitiveGroup(6,2);;
#! gap> RelativeInvariantMinimalDegree(T65, T62^(2,4));
#! [ 2, 2 ]
#! @EndExampleSession
#! @BeginExampleSession
#! gap> T69 := TransitiveGroup(6,9);;
#! gap> T63 := TransitiveGroup(6,4);;
#! gap> RelativeInvariantMinimalDegree(T69, T63^(2,6,4));
#! [ 2, 1 ]
#! @EndExampleSession
#!
#! If <A>H</A> is not a subgroup of <A>G</A> then an error is raised
#! @BeginExampleSession
#! gap> RelativeInvariantMinimalDegree(T65, T62);
#! @EndExampleSession
DeclareGlobalFunction( "RelativeInvariantMinimalDegree" );

#! @Arguments d
#! @Returns the lattice of descent of transitive subgroups
#!
DeclareGlobalFunction( "GaloisDescentTable" );

#TODO: is this really needed!?
# the G-orbit shouuld be enough
#! @Arguments G H q
#! @Returns the list of monomials
DeclareGlobalFunction( "AllMonomials" );

DeclareGlobalFunction( "IndicesToPolynomial" );

DeclareGlobalFunction( "PrintGaloisDescentTable" );
