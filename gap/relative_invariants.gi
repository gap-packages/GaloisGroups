# Return the pair (n, dim) so that n is the minimal degree
# so that there exists a polynomial in Q[x1, ..., xd] of
# degree n that is H invariant but not G invariant. The
# second returned argument dim is the dimension of
# the algebra of G-invariant polynomials on the
# basis of H-invariant ones.
#

PermToGP := function(p,l)
  return Permuted([1..l],p);
end;

InstallGlobalFunction(RelativeInvariantMinimalDegree,
function(G, H)
  local fG, fH, vG, vH, n;
  if NrMovedPoints(G) <> NrMovedPoints(H) then
    Error("G and H have different domains");
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
# Examples
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

# Return a monomial of minimal degree that gives
# a polynomial in Q[x1,x2,...,xd] that is invariant
# under H but not under G
#
# The output is given as a set partitions [ I1, I2, ...]
# that has to be interpreted as x_{I1}^1 x_{I2}^2 ...
# where I1 is the product of the variables indexed
# by I1.
MonomialMinimalDegree := function(G, H, P)
  local d,p,q,s,n,dim,sG,sH,S;
  d := NrMovedPoints(G);
  s := RelativeInvariantMinimalDegree(G, H);
  n := s[1];
  dim := s[2];

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
    Add(L, Minimum(k));
  od;
  return L;
end);

# Write canonical factorization for transitive subgroups
# of Sd into files.
InstallGlobalFunction(GaloisDescentTable,
function(d,filename)
  local first,ffirst,i,j,k,g,t,q,qq,cH,P,G,H;
  PrintTo(filename);   # delete everything
  t := TransitiveMaximalSubgroups(d);
  P := PartitionsByDegree(d);
  for i in [1..NrTransitiveGroups(d)] do
    G := TransitiveGroup(d,i);
    AppendTo(filename, "[\"", ViewString(G), "\",[");
    first := true;
    for g in GeneratorsOfGroup(G) do
      if not first then
        AppendTo(filename,",");
      fi;
      AppendTo(filename, PermToGP(g,d));
      first := false;
    od;
    AppendTo(filename, "]");
    first := true;
    if Size(t[2][i]) <> 0 then
      AppendTo(filename, ",");
    fi;
    for j in [1..Size(t[2][i])] do
      k := t[2][i][j];
      # Ignore the case Sd -> Ad
      if i = NrTransitiveGroups(d) and k = i-1 then
        continue;
      fi;

      G := TransitiveGroup(d, i);
      H := t[3][i][j][1];
      cH := t[3][i][j][2];
      if H^cH <> TransitiveGroup(d,k) then
        Error("pb");
      fi;
      q := MonomialMinimalDegree(G, H, P);
      Print("q=",q,"\n");
      if not first then
        AppendTo(filename , ",");
      fi;
      AppendTo(filename, "[", k, ",", PermToGP(cH, d), ",[");
      ffirst := true;
      for qq in AllMonomials(G,H,q) do
        if not ffirst then
          AppendTo(filename, ",");
        fi;
        AppendTo(filename, FlatMonomial(qq));
        ffirst := false;
      od;
      AppendTo(filename, "]]");
      first := false;
    od;
    AppendTo(filename, "]\n");
  od;
end);
