#
# TODO: Document
#
# here is an example
#
# gap> A13 := TransitiveGroup(13, 8);
# gap> L13 := TransitiveGroup(13, 7);
# gap> t := (1,4,6,2)(3,7,8,5,13,11,10,9);
# gap> ShortCosets(A13, L13, t);
# [ (1,2,3,10,8)(4,9,13,5)(7,12), (1,2,8)(4,7,12,5)(9,10,13,11), (1,2,3,8,6)(4,7,12,9,13,11,5), (1,2,8,6)(3,13,9,10)(4,5)(7,12) ]

# WARNING: This is slow, don't use it in production code.
#          Its mainly used for testing
InstallGlobalFunction(NaiveShortCosets,
function(G, H, t)
    return Filtered( RightTransversal(G, H), x -> t in H^x );
end);

InstallGlobalFunction(ShortCosets,
function(G, H, t)
  local CGt, R, cc, g, CHt, A;
  if not (t in G) then
    Error("t must be an element of G");
  fi;
  if not IsSubgroup(G, H) then
    Error("H must be a subgroup of G");
  fi;
  CGt := Centralizer(G, t);
  R := [];
  for cc in ConjugacyClasses(H) do
    g := RepresentativeAction(G, Representative(cc), t);
    if g <> fail then
      CHt := Centralizer(H^g, t);
      A := RightTransversal(CGt, CHt);
      Append(R, List(A, x -> g*x));
    fi;
  od;
  return R;
end);
