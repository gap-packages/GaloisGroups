gap> chk := function(a)
> local Pa, b;
> Pa := GaloisTestPolynomials[a];
> Print (a,": ",Size(Pa), "\n");
> for b in [1..Size(Pa)] do
>   if Galois(Pa[b])[1] <> b then
>     Print("wrong a=", a, " b=", b, "\n");
>   fi;
> od;
> end;;
gap> chk(2);
2: 1
gap> chk(3);
3: 2
gap> chk(4);
4: 5
gap> chk(5);
5: 5
gap> chk(6);
6: 16
gap> chk(7);
7: 7
gap> chk(8);
8: 50
gap> chk(9);
9: 34
gap> chk(10);
10: 45
