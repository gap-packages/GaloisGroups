# this test takes too long!
# gap> A13 := TransitiveGroup(13,8);;
# gap> L13 := TransitiveGroup(13,7);;
# gap> t := (1,4,6,2)(3,7,8,5,13,11,10,9);;
# gap> sc := Set(ShortCosets(A13, L13, t), x -> RightCoset(L13, x));;
# gap> nsc := Set(NaiveShortCosets(A13, L13, t), x -> RightCoset(L13, x));;
# gap> sc = nsc;
# true
gap> A7 := TransitiveGroup(7,6);;
gap> L7 := TransitiveGroup(7,5);;
gap> ForAll([1..20], function(x)
>  local t, sc, nsc;
>  t := Random(A7);
>  sc := Set(ShortCosets(A7, L7, t), x -> RightCoset(L7, x));;
>  nsc := Set(NaiveShortCosets(A7, L7, t), x -> RightCoset(L7, x));;
>  return sc = nsc;
>  end);
true
