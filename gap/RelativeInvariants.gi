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

InstallGlobalFunction(RelativeInvariantMinimalDegree,
function(G, H)
  local fG, fH, vG, vH, n;
  if MovedPoints(G) <> MovedPoints(H) then
    Error("G and H have different domains");
  fi;
  if not IsSubgroup(G, H) then
    Error("H must be a subgroup of G");
  fi;
  fG := MolienSeries(NaturalCharacter(G));
  fH := MolienSeries(NaturalCharacter(H));
  n := 1;
  vG := ValueMolienSeries(fG, n);
  vH := ValueMolienSeries(fH, n);
  while vG = vH do
    n := n + 1;
    vG := ValueMolienSeries(fG, n);
    vH := ValueMolienSeries(fH, n);
  od;
  return [n, vH - vG];
end);

# Convert a vector (typically a partition obtained from
# PartitionsByDegree) and convert it into a set
# partition. The vector [v0, v1, ..., vk] is converted
# to [[1..v1], [v1+1..(v1+v2)], ...].
#
# We assume that if there is one component vi which is
# zero, then all the following one are as well.
#
# Examples:
#
# gap> SizesToSetPartition([3,1,1]);
# [ [ 1 ], [ 2 ] ]
# gap> SizesToSetPartition([1,1,2,1]);
#
SizesToSetPartition := function(p)
  local j,k,part;
  k := 1;
  part := [];
  for j in [2..Size(p)] do
    if p[j] = 0 then
      break;
    fi;
    Add(part, List([k..(k+p[j]-1)]));
    k := k + p[j];
  od;
  return part;
end;


# Return the Partitions of d ordered by the degree of the associated
# monomial.
#
PartitionsByDegree := function(d)
  local i, p, s, L;
  L := [];
  for i in [1..(d*(d-1)/2)] do
    Add(L, []);
  od;
  for p in Partitions(d) do
    s := 0;
    for i in [2..Size(p)] do
      s := s + (i-1) * p[i];
    od;
    if s <> 0 then
      Add(L[s], p);
    fi;
  od;
  return L;
end;

# Return a monomial of minimal degree that gives
# a polynomial in Q[x1,x2,...,xd] that is invariant
# under H but not under G
#
# The output is given as a set partitions [ I1, I2, ...]
# that has to be interpreted as x_{I1}^1 x_{I2}^2 ...
# where I1 is the product of the variables indexed
# by I1.
MonomialMinimalDegree := function(G, H, P)
  local d,p,q,s,n,sG,sH,S;
  d := NrMovedPoints(G);
  s := RelativeInvariantMinimalDegree(G, H);
  n := s[1];

  # try the naive monomial
  for p in P[n] do
    q := SizesToSetPartition(p);
    sG := Size(G) / Size(Stabilizer(G, q, OnTuplesSets));
    sH := Size(H) / Size(Stabilizer(H, q, OnTuplesSets));
    if sG <> sH then
      return q;
    fi;
  od;

  # random
  S := SymmetricGroup(d);
  while true do
    for p in P[n] do
      q := SizesToSetPartition(p);
      q := OnTuplesSets(q, Random(S));
      sG := Size(G) / Size(Stabilizer(G, q, OnTuplesSets));
      sH := Size(H) / Size(Stabilizer(H, q, OnTuplesSets));
      if sG <> sH then
        return q;
      fi;
    od;
  od;
end;

InstallGlobalFunction(AllMonomials,
function(G, H, q)
  local o,k,L;
  # o := Orb(G, q, OnTuplesSets);
  # Enumerate(o);
  # return FindSuborbits(o, H);
  o := Orbit(G, q, OnTuplesSets);
  L := [];
  for k in Orbits(H, o, OnTuplesSets) do
    Add(L, [Size(k), Minimum(k)]);
  od;
  return L;
end);

InstallGlobalFunction(IndicesToPolynomial,
function(I, n)
  local R,i,j,p,x,m;
  R := PolynomialRing(Rationals, List([1..n], i->Concatenation("x", String(i))));
  x := IndeterminatesOfPolynomialRing(R);
  p := 0;
  for i in I do
    m := 1;
    for j in i do
      m := m * x[j];
    od;
    p := p + m;
  od;
  return p;
end);


# Compute data associated to the transitive subgroups in SymmetricGroup(d)
# the i-th list corresponds to the group TransitiveGroup(d,i) as follows
#   line[1] : name
#   line[2] : list of generators
#   line[3] : data corresponding to maximal subgroups
#
# the data for maximal subgroups is a list of the following form
#
#   subg[1] : index of the subgroup in the databse
#   subg[2] : a permutation such that conjugates
#   subg[3] : a list of pairs [n, t] where n is an integer and t
#             an ordered set partition

CachedGaloisDescentTable := [];

InstallGlobalFunction(GaloisDescentTable,
function(d)
  local L,t,P,i,j,G,gr,subg,subgs,k,H,cH,q;
  if not IsInt(d) or d <= 1 then
    Error("invalid d parameter");
  fi;
  if IsBound(GaloisDescentTables[d]) then
    return GaloisDescentTables[d];
  fi;
 if IsBound(CachedGaloisDescentTable[d]) then
    return CachedGaloisDescentTable[d];
  fi;

  L := [];
  t := TransitiveMaximalSubgroups(d);
  P := PartitionsByDegree(d);
  if Size(t[2]) <> NrTransitiveGroups(d) or Size(t[3]) <> NrTransitiveGroups(d) then
    Error("wrong assumption");
  fi;
  for i in [1..NrTransitiveGroups(d)] do
    G := TransitiveGroup(d,i);
    gr := [];

    # 1 group string
    Add(gr, ViewString(G));

    # 2. generators
    Add(gr, GeneratorsOfGroup(G));

    # 3. list of subgroups
    subgs := [];
    for j in [1..Size(t[2][i])] do
      subg := [];
      k := t[2][i][j];
      H := t[3][i][j][1];
      cH := t[3][i][j][2];
      if H^cH <> TransitiveGroup(d,k) then
        Error("pb");
      fi;

      # Ignore the case Sign(G) <> Sign(H)
      if SignPermGroup(G) <> SignPermGroup(H) then
        continue;
      fi;

      # 3.a conjugation
      Add(subg, k);
      Add(subg, cH);

      # 3.b monomials
      q := MonomialMinimalDegree(G, H, P);
      Add(subg, AllMonomials(G, H, q));

      Add(subgs, subg);
    od;

    Add(gr, subgs);
    Add(L, gr);
  od;

  CachedGaloisDescentTable[d] := L;
  return L;
end);

InstallGlobalFunction(PrintGaloisDescentTable,
function(d)
  local T,a,G,b,H,l,K;
  T := GaloisDescentTable(d);
  if NrTransitiveGroups(d) <> Size(T) then
    Error("PROUT");
  fi;
  for a in [1..Size(T)] do
    Print("T(", d, ",", a, " = ", T[a][1], ")\n");
    G := TransitiveGroup(d, a);
    for l in T[a][3] do
      b := l[1];
      H := TransitiveGroup(d,b)^(l[2]^-1);
      K := Orbit(H, l[3][1][2], OnTuplesSets);
      K := List(K, FlatMonomial);
      Print("  T(", d, ",", b, "): ", IndicesToPolynomial(K,d), "\n");
    od;
  od;
end);


#######################################################################
# Old code that might be recycled later on
#
#ListAppendToNoSpace := function(filename, L)
#  local x,first;
#  first := true;
#  AppendTo(filename, "[");
#  for x in L do
#    if not first then
#      AppendTo(filename, ",");
#    fi;
#    if IsString(x) and Size(x) <> 0 then
#      AppendTo(filename, "\"", x, "\"");
#    elif IsList(x) then
#      ListAppendToNoSpace(filename, x);
#    else
#      AppendTo(filename, x);
#    fi;
#    first := false;
#  od;
#  AppendTo(filename, "]");
#end;
#
## Write canonical factorization for transitive subgroups
## of Sd into files.
##
## The format is the following: the line i corresponds to the
## group Transitive(d,i)
# PrintGaloisDescentTable := function(T, filename...)
#  local line;
#  if Size(filename) = 0 then
#    filename := "*stdout*";
#  elif Size(filename) = 1 then
#    filename := filename[1];
#  else
#    Error("wrong argument lists");
#  fi;
#  PrintTo(filename);   # delete everything
#
#  if IsInt(T) then
#    T := GaloisDescentTable(T);
#  fi;
#
#  for line in T do
#    ListAppendToNoSpace(filename, line);
#    AppendTo(filename, "\n");
#  od;
# end;
#
#PermToGP := function(p,l)
#  return Permuted([1..l],p^-1);
#end;
#
#
#PrintGaloisDescentTableGP := function(T, filename...)
#  local d, line, G;
#  if Size(filename) = 0 then
#    filename := "*stdout*";
#  elif Size(filename) = 1 then
#    filename := filename[1];
#  else
#    Error("wrong argument lists");
#  fi;
#  PrintTo(filename);   # delete everything
#
#  if IsInt(T) then
#    d := T;
#    T := GaloisDescentTable(T);
#  else
#    G := Group(T[Size(T)][2]); # should be S_d
#    d := NrMovedPoints(G);
#  fi;
#
#  for line in T do
#    line[2] := List(line[2], x->PermToGP(x,d));
#    line[3] := List(line[3], x->[x[1], PermToGP(x[2],d), List(x[3], y->FlatMonomial(y[2]))]);
#    ListAppendToNoSpace(filename, line);
#    AppendTo(filename, "\n");
#  od;
#end;
