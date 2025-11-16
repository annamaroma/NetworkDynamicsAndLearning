function min_cut = extra_units_capacity(nodes, G, maxIter)
% algotirhm idea:
%  - find current min cut value and save in vecIndex all cuts that achieve it
%  - for those cuts, collect the set of edges crossing the cut
%  - find edges common to all min-cuts and increment their capacity by 1
%  - save the min cut value into vector min_cut for each iteration
% stop when when maxIter is reached
 
min_cut = zeros(maxIter,1);

for iter = 1:maxIter
    [U,V,C] = compute_cuts(nodes, G);
    minC = min(C); %find minimum cut value
    min_cut(iter) = minC; %store min cut value for this iteration

    vecIndex = find(C==minC);
    min_cuts = length(vecIndex);
    fprintf('\nIteration %d: minimum cut = %d (found in %d combos)\n', iter, minC, min_cuts);

    % count links among min-cuts
    linksLists = cell(min_cuts,1);
    for j = 1:min_cuts
        comb_U = U{vecIndex(j)};
        comb_V = V{vecIndex(j)};
        links = [];
        for a = 1:numel(comb_U)
            for b = 1:numel(comb_V)
                eidx = findedge(G, comb_U(a), comb_V(b));
                if eidx > 0
                    links = [links; eidx];
                end
            end
        end
        linksLists{j} = unique(links);
    end

    % find common links among all min-cut edge sets
    if isempty(linksLists)
        fprintf('No min-cut combos found, stopping.\n');
        break;
    end
    commonEdges = linksLists{1};
    for k = 2:numel(linksLists)
        commonEdges = intersect(commonEdges, linksLists{k});
    end

    if isempty(commonEdges)
        % No single edge common to all min-cuts. Choose one edge to increment
        % Strategy: select the edge that appears most frequently across the
        % min-cut edge lists.
        allEdges = [];
        for ii = 1:numel(linksLists)
            allEdges = [allEdges; linksLists{ii}(:)];
        end
        uniqE = unique(allEdges);
        counts = zeros(size(uniqE));
        for ii = 1:numel(uniqE)
            counts(ii) = sum(allEdges==uniqE(ii));
        end
        [~, maxIdx] = max(counts);
        chosenEdge = uniqE(maxIdx);
        fprintf('No common link among all min-cuts at iteration %d. Chosen edge %d appears in %d cuts -> incrementing it.\n', iter, chosenEdge, counts(maxIdx));
        commonEdges = chosenEdge;
    end

    % increment capacity of a single chosen common link by 1
    % selection strategy: choose the first common edge (changeable)
    chosenEdge = commonEdges(1);
    G.Edges.Weight(chosenEdge) = G.Edges.Weight(chosenEdge) + 1;
    % also update original w vector if present and aligned with edges
    if exist('w','var') && numel(w) >= chosenEdge
        w(chosenEdge) = w(chosenEdge) + 1;
    end
    fprintf('  Incremented chosen edge %d, new weight = %d\n', chosenEdge, G.Edges.Weight(chosenEdge));
end


