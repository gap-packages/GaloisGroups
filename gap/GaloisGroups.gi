#############################################################################
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
#############################################################################

#
# GaloisGroups: Computing Galois Groups using GAP and PARI via
# Stauduhar method

########################################################
# Function calls to PARI/GP

PARIInitialise();

GPToPerm := function(l)
  return AsPermutation(Transformation(l));
end;

PermToGP := function(p, l)
  return Permuted([1..l], p^-1);
end;

# Call the pari getroots function
# return an approximation of the roots of the polynomial P
# (for now hardcoded with precision 100)
#
# Examples:
#
# gap> GetRoots(x^5-2);
# PARI([[1891, 1036*y + 1033, 3146*y + 1697, 23*y + 1720, 2133*y + 3166]~, y^2 + y + 2, 3169])
#
PARIGetRoots := function(P)
  return PARI_CALL2("getroots", PARIPolynomial(P), INT_TO_PARI_GEN(100));
end;

# Perform a random Tschirnausen transformation on P
#
# Examples: (though output is random)
#
# gap> PARITschirnhaus(x^5-2);
# x^5+10*x^4+40*x^3-40*x^2-940*x-2142
PARITschirnhaus := function(P)
  local Q;
  Q := PARI_CALL1("tschirnhaus", PARIPolynomial(P));
  return UnivariatePolynomial(Rationals, PARI_GEN_GET_DATA(Q));
end;

# Perform a random Tschirnausen transformation on P where Q are the roots
# Return [P1,Q1] where P1 is the new polynomial and Q1 is the roots of P1
# ordered using the same embedding order as for Q.
PARITschirnhausen := function(P, Q)
  local U;
  U := PARI_VEC_TO_LIST(PARI_CALL2("Tschirnhausen", PARIPolynomial(P), Q));
  return [UnivariatePolynomial(Rationals, PARI_GEN_GET_DATA(U[1])), U[2]];
end;


# Return the action fo the Frobenius on the roots as a permutation
PARIGetFrobenius := function(Q)
  return GPToPerm(PARI_GEN_GET_DATA(PARI_CALL1("getperm", Q)));
end;

# Check whether the resolvant is squarefree
PARICosets_squarefree := function(C, K, Q, P)
  return PARI_GEN_GET_DATA(PARI_CALL4("cosets_squarefree", PARI_VECVECSMALL(C), PARI_VECVECSMALL(K), Q, PARIPolynomial(P))) = 1;
end;


# Check whether the resolvant has a rational root
PARICosets4 := function(C, K, Q, P)
  local r;
  r:= PARI_GEN_GET_DATA(PARI_CALL4("cosets4", PARI_VECVECSMALL(C), PARI_VECVECSMALL(K), Q, PARIPolynomial(P)));
  if IsList(r) then
    return GPToPerm(r);
  else
    return r;
  fi;
end;

# Examples:
#
# gap> FlatMonomial([[1,2,3]]);
#[  1, 2, 3 ]
# gap> FlatMonomial([[1,2,3],[4,5,6]]);
# [ 1, 2, 3, 4, 4, 5, 5, 6, 6 ]
#
FlatMonomial := function(p)
  local i, j, k, q;
  q := [];
  for i in [1..Size(p)] do
    for j in [1..Size(p[i])] do
      for k in [1..i] do
        Add(q, p[i][j]);
      od;
    od;
  od;
  return q;
end;

#####################################################################
# Main function

# Compute the Galois group by doing descent (Stauduhar method)
#
# Input: P a polynomial
#
# Examples:
#
# gap> GaloisDescentStauduhar(x^5-2);
# [ 3, (3,4) ]
# gap> GaloisDescentStauduhar(x^5+20*x+16);
# [ 4, () ]
GaloisDescentStauduhar := function(P)
  local Ta, d, S, T, s, Q, FC, C, l, a, rho, tau, sigma, name, gen, list, b, bloc, K, frob, G, H, blocGP, av, U, Qo;
  d := Degree(P);
  if d = 1 then return [ 1, () ]; fi;
  Ta := GaloisDescentTable(Degree(P));
  av := PARIGetAvma();
  Q := PARIGetRoots(P);
  Qo := PARI_GEN_TO_OBJ(Q);
  frob := PARIGetFrobenius(Q);
  T := List(Ta, x->x[3]);
  if IsSquareInt(Discriminant(P)) then
    a := Length(T)-1; #A_n
  else
    a := Length(T);   #S_n
  fi;
  sigma := ();
  repeat
    list := T[a];
    rho := 0;
    for l in list do
      b := l[1]; tau := l[2]; bloc := List(l[3], x->x[2]);
      # current group
      G := TransitiveGroup(d,a)^sigma;
      # subgroup to be tested
      H := TransitiveGroup(d,b)^(tau^-1*sigma);

      if not IsSubgroup(G, H) then
        Error("wrong assumption for descent\n");
      fi;

      # construct the resolvant
      # NOTE: each element of bloc is a valid resolvant (ie a H-invariant but
      #       not G-invariant polynomial). For now, we only consider bloc[1].
      #Print("#I G = ", G, " ~ T(", d, ",", a, ")\n");
      #Print("#I H = ", H, " ~ T(", d, ",", b, ")\n");
      S := OnTuplesSets(bloc[1], sigma);
      K := List(RightTransversal(H,Stabilizer(H,S,OnTuplesSets)),x->FlatMonomial(OnTuplesSets(S,x)));;
      FC := List(RightTransversal(G, H), c->PermToGP(c, d));
      while not PARICosets_squarefree(FC, K, Q, P) do
      #  Print("#I Applying Tschirnhausen transform\n");
        U := PARITschirnhausen(P, Q);
        P := U[1]; Q := U[2];
      od;
      C := List(ShortCosets(G, H, frob), c->PermToGP(c, d));
      rho := PARICosets4(C, K, Q, P);
      if IsPerm(rho) then
        sigma := tau^-1*sigma*rho;
        a := b;
        break;
      fi;
    od;
  until rho = 0;
  PARISetAvma(av);
  return [a,sigma,Qo];
end;

InstallGlobalFunction(Galois, function(P)
  return GaloisDescentStauduhar(P);
end);
