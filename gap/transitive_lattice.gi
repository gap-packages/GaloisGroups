# Super ugly and inefficient
InstallGlobalFunction(TransitiveMaximalSubgroups,
function(deg)
    local omega, g, G, TMCR, tree, edges, i, j, k, graph, subtrees, todo, conjugate, trans;
    
    omega := [1..deg];
    subtrees := [];
    conjugate := [];
    tree := [];

    G := SymmetricGroup(deg);
    trans := AllTransitiveGroups(NrMovedPoints, deg);

    # Note that it is important to specify the domain!
    tree := [G];
    edges := [];
    todo := [TransitiveIdentification(G)];
    while not IsEmpty(todo) do
        k := Remove(todo);
        Print("Looking at ", k, " (", trans[k], ")\n");
        
        if not IsBound(subtrees[k]) then
            TMCR := Filtered(MaximalSubgroupClassReps(TransitiveGroup(deg, k)), x -> IsTransitive(x, omega));
            subtrees[k] := [];
            conjugate[k] := [];
            for g in TMCR do
                i := TransitiveIdentification(g);
                Print("rep act: ", g, " ", trans[i], "\n");
                j := RepresentativeAction(SymmetricGroup(deg), g, trans[i]);
                # j := ContainingConjugates(SymmetricGroup(deg), g, trans[i]);
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

# FIXME: Reference is algorithm 4.6 from Klüners/Fieker "Computing Galois Groups"
InstallGlobalFunction(BlockLadder,
function(grp, omega, a)
    local  b, i, H;
        
end);

