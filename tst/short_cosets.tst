gap> A13 := TransitiveGroup(13,8);;
gap> L13 := TransitiveGroup(13,7);;
gap> t := (1,4,6,2)(3,7,8,5,13,11,10,9);;
gap> sc := Set(ShortCosets(A13, L13, t), x -> RightCoset(L13, x));;
gap> nsc := Set(NaiveShortCosets(A13, L13, t), x -> RightCoset(L13, x));;
gap> sc = nsc;
true
