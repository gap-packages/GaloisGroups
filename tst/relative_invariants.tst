gap> PartitionsByDegree(3);
[ [ [ 2, 1 ] ], [  ], [ [ 1, 1, 1 ] ] ]
gap> PartitionsByDegree(4);
[ [ [ 3, 1 ] ], [ [ 2, 2 ] ], [ [ 2, 1, 1 ] ], [  ], [  ], [ [ 1, 1, 1, 1 ] ] 
 ]
gap> S5 := SymmetricGroup(5);;
gap> A5 := AlternatingGroup(5);;
gap> RelativeInvariantMinimalDegree(S5,A5);
[ 10, 1 ]
gap> MonomialMinimalDegree(S5,A5,PartitionsByDegree(5));
[ [ 1 ], [ 2 ], [ 3 ], [ 4 ] ]
gap> G := TransitiveGroup(11, 4);;
gap> H := TransitiveGroup(11, 2);;
gap> RelativeInvariantMinimalDegree(G, H);
[ 2, 4 ]
gap> MonomialMinimalDegree(G,H,PartitionsByDegree(11));
[ [ 1, 2 ] ]
