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

PermToGP := function(p,l)
  return Permuted([1..l],p^-1);
end;

GPGetRoots := function(P)
  return PARI_CALL2("getroots",PARIPolynomial(P),INT_TO_PARI_GEN(100));
end;

GPTschirnhaus := function(P)
  local Q;
  Q := PARI_CALL1("tschirnhaus",PARIPolynomial(P));
  return UnivariatePolynomial(Rationals,PARI_GEN_GET_DATA(Q));
end;

GPGetFrobenius := function(Q)
  return GPToPerm(PARI_GEN_GET_DATA(PARI_CALL1("getperm",Q)));
end;

GPCosets3 := function(C,K,Q,P)
  local r;
  r := PARI_GEN_GET_DATA(PARI_CALL4("cosets3",PARI_VECVECSMALL(C),PARI_VECVECVECSMALL(K),Q,PARIPolynomial(P)));
  if IsList(r) then
    return GPToPerm(r);
  else
    return r;
  fi;
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

#####################################################################
# Main function

#InstallGlobalFunction( GaloisGroup2,
#function(p)
#  local d, a, b, sigma, tau, rho, frob, G, H, blocks, C, i, j, s, k, O, subgroup;
#
#  d := Degree(p);
#  if d <> 7 then
#    Error("not implemented yet");
#  fi;
#
#  # the Frobenius should be computed from PARI/GP
#  frob := (1,3,5,2,4,7,6);
#
#  if SignPerm(frob) = -1 then
#    a := NrTransitiveGroups(d);
#  else
#    a := NrTransitiveGroups(d) - 1;
#  fi;
#  sigma := ();
#
#  Print("frob=",frob,"\n");
#  while true do
#    G := TransitiveGroup(d,a)^sigma;
#    i := fail;
#    for subgroup in T[a][4] do
#      b := subgroup[1];
#      tau := subgroup[2];
#      blocks := subgroup[3];
#      H := TransitiveGroup(d,b)^(tau^-1*sigma);
#
#      Print("a = ",a," b = ", b, "\n");
#      Print("sigma=",sigma,"\n");
#      Print("tau=",tau,"\n");
#
#      C := ShortCosets(G, H, frob);
#      i := -1;
#      for i in [1..Size(C)] do
#        rho := C[i];
#        Print("  ", i, " -> ", rho, "\n");
#        for j in [1..Size(blocks)] do
#          s := blocks[j][1];
#          k := blocks[j][2];
#          O := Orbit(H^rho, OnTuplesSets(k, sigma*rho), OnTuplesSets);
#          if Size(O) <> s then
#            Error("wrong coset size");
#          fi;
#          Print("  fac", j, "=",SortedList(O), "\n");
#        od;
#      od;
#
#      # here we choose rho by asking to the user... should
#      # use PARI/GP instead by computing the potential factors
#      # of the resolvent obtained from O.
#      i := -1;
#      while not (i in [1..Size(C)] or i = fail) do
#        i := InputFromUser("choice> ");
#      od;
#
#      if i <> fail then
#        rho := C[i];
#        break;
#      fi;
#    od;
#
#    if i = -1 or i = fail then
#      break;
#    else
#      a := b;
#      sigma := tau^-1*sigma*rho;
#    fi;
#  od;
#
#  return [T[a][1], a, TransitiveGroup(d,a)^sigma];
#end );


GaloisFromTable:= function(P,Ta)
  local d, T, s, Q, C, l, a, rho, tau, sigma, name, gen, list, b, bloc, K, frob, G, H, blocGP;
  Q := GPGetRoots(P);
  frob := GPGetFrobenius(Q);
  d := Degree(P);
  T := List(Ta, x-> x[3]);
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
      G := TransitiveGroup(d,b)^(tau^-1*sigma);
      H := TransitiveGroup(d,a)^sigma;
      C := List(ShortCosets(H, G, frob), c->PermToGP(c, d));
      K := List(bloc,bl->Orbit(G, OnTuplesSets(bl, sigma), OnTuplesSets));
      K := List(K, y->SortedList(List(y, FlatMonomial)));
      rho := GPCosets3(C, K, Q, P);
      if IsPerm(rho) then
        sigma := tau^-1*sigma*rho;
        a := b;
        break;
      elif rho = 1 then
        Print("#I Applying Tschirnhausen transform\n");
        return GaloisFromTable(GPTschirnhaus(P), Ta);
      fi;
    od;
  until rho = 0;
  return [a,sigma];
end;

InstallGlobalFunction(Galois, function(P)
  return GaloisFromTable(P, GaloisDescentTable(Degree(P)));
end);

