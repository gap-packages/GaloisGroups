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
#
# Implementations
#
# Right now, the descent is done manually. At each step, the
# user should choose between a coset or fail. This step should
# be replaced with PARI/GP function testing integrality of
# polynomials.

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

# Perform a random Tschirnaus transformation on P
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

PARIGetFrobenius := function(Q)
  return GPToPerm(PARI_GEN_GET_DATA(PARI_CALL1("getperm", Q)));
end;

# evaluation of resolvent... what is this function actually doing?
PARICosets3 := function(C, K, Q, P)
  local r;
  r := PARI_GEN_GET_DATA(PARI_CALL4("cosets3", PARI_VECVECSMALL(C), PARI_VECVECVECSMALL(K), Q, PARIPolynomial(P)));
  if IsList(r) then
    return GPToPerm(r);
  else
    return r;
  fi;
end;

PARICosets_squarefree := function(C, K, Q, P)
  return PARI_GEN_GET_DATA(PARI_CALL4("cosets_squarefree", PARI_VECVECSMALL(C), PARI_VECVECVECSMALL(K), Q, PARIPolynomial(P))) <> 0;
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
  local Ta, d, T, s, Q, FC, C, l, a, rho, tau, sigma, name, gen, list, b, bloc, K, frob, G, H, blocGP;
  Ta := GaloisDescentTable(Degree(P));
  Q := PARIGetRoots(P);
  frob := PARIGetFrobenius(Q);
  d := Degree(P);
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
      Print("#I G = ", G, " ~ T(", d, ",", a, ")\n");
      Print("#I H = ", H, " ~ T(", d, ",", b, ")\n");
      K := Orbit(H, OnTuplesSets(bloc[1], sigma), OnTuplesSets);
      K := List(K, FlatMonomial);
      FC := List(RightTransversal(G, H), c->PermToGP(c, d));
      Print("#I resolvant ", IndicesToPolynomial(K, d), "\n");
      Print("#I cosets ", FC, "\n");
      if PARICosets_squarefree(FC, K, Q, P) then
        Print("#I Applying Tschirnhausen transform\n");
        return GaloisDescentStauduhar(PARITschirnhaus(P));
      fi;
      # short cosets of G/H (wrt the Frobenius permutation)
      C := List(ShortCosets(G, H, frob), c->PermToGP(c, d));
      Print("#I short cosets ", C, "\n");
      rho := PARICosets3(C, K, Q, P);
      if IsPerm(rho) then
        sigma := tau^-1*sigma*rho;
        a := b;
        break;
      fi;
    od;
  until rho = 0;
  return [a,sigma];
end;

InstallGlobalFunction(Galois, function(P)
  return GaloisDescentStauduhar(P);
end);

