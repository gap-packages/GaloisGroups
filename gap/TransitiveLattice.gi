# Super ugly and inefficient
InstallGlobalFunction(TransitiveMaximalSubgroups,
function(deg)
    local omega
        , g
        , Sym
        , TMCR
        , edges
        , i, j, k
        , graph
        , subtrees
        , todo
        , conjugate
        , trans;
    if not IsPosInt(deg) then
        Error("deg has to be a positive integer");
    fi;

    Sym := SymmetricGroup(deg);
    trans := AllTransitiveGroups(NrMovedPoints, deg);

    omega := [1..deg];
    subtrees := [];
    conjugate := [];
    edges := [];

    todo := [Length(trans)]; # This is the symmetric group

    while not IsEmpty(todo) do
        k := Remove(todo);
        Info(InfoTransitiveMaximal, 5,
             "Looking at ", k, " (", trans[k], ")");

        if not IsBound(subtrees[k]) then
            TMCR := Filtered(MaximalSubgroupClassReps(TransitiveGroup(deg, k)), x -> IsTransitive(x, omega));
            subtrees[k] := [];
            conjugate[k] := [];
            for g in TMCR do
                i := TransitiveIdentification(g);
                Info(InfoTransitiveMaximal, 5,
                     "rep act: ", g, " ", trans[i]);
                j := RepresentativeAction(Sym, g, trans[i]);
                Add(subtrees[k], i);
                Add(conjugate[k], [g, j]);
            od;

            Append(edges, List(subtrees[k], x -> [k, x]));
            Append(todo, subtrees[k]);
        fi;
    od;

    graph := DigraphByEdges(edges);
    return [graph, subtrees, conjugate];
end);

InstallGlobalFunction(TestTransitiveMaximals,
function(deg)
    local t, i, j, g;

    t := TransitiveMaximalSubgroups(deg);

    for i in [1..Length(t[3])] do
        for j in [1..Length(t[3][i])] do
            if TransitiveGroup(deg, t[2][i][j]) <> t[3][i][j][1] ^ (t[3][i][j][2]) then
                Error("Group failed verification.");
            fi;
        od;
    od;
    return true;
end);


InstallGlobalFunction(TransitiveMaximalSubgroups2,
function(deg)
    local grps, edges, todo, i, j, n, graph, t;

    grps := AllTransitiveGroups(NrMovedPoints, deg);
    n := Length(grps);
    graph := [];
    edges := [];

    for i in [n,n-1..1] do
        graph[i] := [];
        for j in [i-1,i-2..1] do
            graph[i][j] := ContainedConjugates(SymmetricGroup(deg), grps[i], grps[j]);
            for t in graph[i][j] do
                Add(edges, [i,j]);
            od;
        od;
    od;
    return [DigraphByEdges(edges), graph];
end);

# Return the list of indices [G:H] for maximal subgroups H that strictly contains
# the Stabilizer(G, 1).
# TODO: we actually want *all* subgroups sandwiched between Stabilizer(G,1) and G
#       and not only the maximal ones... don't know how to do
DegreesMaximalStabilizer := function(G)
  local S,H,L,i,g;
  L := [];
  S := Stabilizer(G, 1);
  for H in MaximalSubgroupClassReps(G) do
    if Size(H) <= Size(S) then
      continue;
    fi;

    i := Index(G, H);
    for g in RightTransversal(G, Normalizer(G,H)) do
      if IsSubgroup(H^g, S) then
          Add(L, i);
      fi;
    od;
  od;
  return SortedList(L);
end;

# TODO: we actually want *all* subgroups and not only the maximal ones!
StabilizerData := function(d)
  local G;
  return List(AllTransitiveGroups(NrMovedPoints, d), G->DegreesMaximalStabilizer(G));
end;

