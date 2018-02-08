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

#####################################################################
# table of groups and subgroups for degree 7
# This has been generated with the function  GaloisDescentTable
# from the GaloisGroups package.
T := [
["C(7) = 7",1,[(1,2,3,4,5,6,7)],[]],
["D(7) = 7:2",-1,[(1,2,3,4,5,6,7),(1,6)(2,5)(3,4)],[[1,(),[[ 7, [ [ 1 ], [ 2 ] ] ],[ 7, [ [ 1 ], [ 7 ] ] ]]]]],
["F_21(7) = 7:3",1,[(1,2,3,4,5,6,7),(1,2,4)(3,6,5)],[[1,(),[[ 7, [ [ 1, 2 ] ] ],[ 7, [ [ 1, 3 ] ] ],[ 7, [ [ 1, 4 ] ] ]]]]],
["F_42(7) = 7:6",-1,[(1,2,3,4,5,6,7),(1,3,2,6,4,5)],[[3,(),[[ 21, [ [ 1 ], [ 2 ] ] ],[ 21, [ [ 1 ], [ 4 ] ] ]]],[2,(),[[ 7, [ [ 1, 2 ] ] ],[ 7, [ [ 1, 3 ] ] ],[ 7, [ [ 1, 4 ] ] ]]]]],
["L(7) = L(3,2)",1,[(1,2,3,4,5,6,7),(1,2)(3,6)],[[3,(),[[ 21, [ [ 1, 2, 3 ] ] ],[ 7, [ [ 1, 2, 6 ] ] ]]]]],
["A7",1,[(1,2,3,4,5,6,7),(5,6,7)],[[5,(),[[ 28, [ [ 1, 2, 3 ] ] ],[ 7, [ [ 1, 2, 4 ] ] ]]],[5,(2,5)(3,4,6),[[ 28, [ [ 1, 2, 3 ] ] ],[ 7, [ [ 1, 2, 4 ] ] ]]]]],
["S7",-1,[(1,2,3,4,5,6,7),(1,2)],[[4,(2,4,7),[[ 14, [ [ 1, 2, 3 ] ] ],[ 21, [ [ 1, 2, 4 ] ] ]]]]],
];

#######################################################
# functions that can alternatively be found in TransGrp
if LoadPackage("TransGrp") = fail then
  TransitiveGroup := function(d, a)
    return Group(T[a][3]);
  end;

  NrTransitiveGroups := function(d)
    return Size(T);
  end;
fi;

#####################################################################
# Main function

InstallGlobalFunction( GaloisGroup2,
function(p)
  local d, a, b, sigma, tau, rho, frob, G, H, blocks, C, i, j, s, k, O, subgroup;

  d := Degree(p);
  if d <> 7 then
    Error("not implemented yet");
  fi;

  # the Frobenius should be computed from PARI/GP
  frob := (1,3,5,2,4,7,6);

  if SignPerm(frob) = -1 then
    a := NrTransitiveGroups(d);
  else
    a := NrTransitiveGroups(d) - 1;
  fi;
  sigma := ();

  Print("frob=",frob,"\n");
  while true do
    G := TransitiveGroup(d,a)^sigma;
    i := fail;
    for subgroup in T[a][4] do
      b := subgroup[1];
      tau := subgroup[2];
      blocks := subgroup[3];
      H := TransitiveGroup(d,b)^(tau^-1*sigma);

      Print("a = ",a," b = ", b, "\n");
      Print("sigma=",sigma,"\n");
      Print("tau=",tau,"\n");

      C := ShortCosets(G, H, frob);
      i := -1;
      for i in [1..Size(C)] do
        rho := C[i];
        Print("  ", i, " -> ", rho, "\n");
        for j in [1..Size(blocks)] do
          s := blocks[j][1];
          k := blocks[j][2];
          O := Orbit(H^rho, OnTuplesSets(k, sigma*rho), OnTuplesSets);
          if Size(O) <> s then
            Error("wrong coset size");
          fi;
          Print("  fac", j, "=",SortedList(O), "\n");
        od;
      od;

      # here we choose rho by asking to the user... should
      # use PARI/GP instead by computing the potential factors
      # of the resolvent obtained from O.
      i := -1;
      while not (i in [1..Size(C)] or i = fail) do
        i := InputFromUser("choice> ");
      od;

      if i <> fail then
        rho := C[i];
        break;
      fi;
    od;

    if i = -1 or i = fail then
      break;
    else
      a := b;
      sigma := tau^-1*sigma*rho;
    fi;
  od;

  return [T[a][1], a, TransitiveGroup(d,a)^sigma];
end );

