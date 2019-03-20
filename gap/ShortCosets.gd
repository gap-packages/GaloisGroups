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

#! @Chapter ShortCosets
#
#! @Section Computation with short cosets

DeclareGlobalFunction( "NaiveShortCosets" );

#! @Arguments G H t
#! @Returns list
#! @Description
#!  Returns the list of cosets of G / C_G(H) such that H^s
#!  contains t. For Galois group consideration, t is either
#!  the Frobenius action on the roots (in case p-adic roots
#!  are considered) or the complex conjugation.
#! @BeginExampleSession
#! gap> T6_5 := TransitiveGroup(6,5);;
#! gap> T6_2 := TransitiveGroup(6,2)^(2,4);;
#! gap> ShortCosets(T6_5, T6_2, (1,4)(2,5)(3,6));
#! [ (1,5,3), (1,6)(2,3)(4,5), (2,4,6) ]
#! gap> ShortCosets(T6_5, T6_2, (1,4,3,6,5,2));
#! [  ]
#! @EndExampleSession
#! Here is a very huge descent where the short coset can dramatically
#! reduce the search.
#! @BeginExampleSession
#! gap> T10_45 := TransitiveGroup(10, 45);;
#! gap> T10_35 := TransitiveGroup(10, 35)^(2,8,7)(3,5)(4,9);;
#! gap> Index(T10_45, T10_35);
#! 2520
#! gap> Size(ShortCosets(T10_45, T10_35, (1,8)(2,9)(3,7)(4,10)));
#! 24
#! gap> Size(ShortCosets(T10_45, T10_35, (1,8,10)(2,7,6)(3,5,9)));
#! 9
#! @EndExampleSession
DeclareGlobalFunction( "ShortCosets" );
