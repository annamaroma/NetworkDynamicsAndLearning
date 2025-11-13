function [used_combinations_cell, unused_combinations_cell, capacities] = compute_cuts(V, G)

v=[];
for i=2:size(V,1)-1
    v = [v V(i)];
end

n = length(v);

used_combinations_cell = {1};
unused_combinations_cell = {[5,2,3,4]};

for k = 1:n
    C_k = nchoosek(v, k);
    
    num_combos_k = size(C_k, 1);
    
    C_k_with_one = [ones(num_combos_k, 1), C_k];

    C_k_cell = mat2cell(C_k_with_one, ones(num_combos_k, 1), k+1);
    used_combinations_cell = [used_combinations_cell; C_k_cell];

    unused_k_cell = cell(num_combos_k, 1);
    
    for i = 1:num_combos_k
        current_combo = C_k(i, :);
        
        unused_elements = setdiff(v, current_combo);
        
        unused_with_five = [5, unused_elements];
        
        unused_k_cell{i} = unused_with_five;
    end
    
    unused_combinations_cell = [unused_combinations_cell; unused_k_cell];
end

capacities = [];

for i = 1:size(used_combinations_cell,1)
    comb_U = used_combinations_cell{i};
    comb_V = unused_combinations_cell{i};
    comb_cap = 0;
    for j = 1:size(comb_U,2)
        for k = 1:size(V,1)-size(comb_U,2)
            edge_idx = findedge(G,comb_U(j), comb_V(k));
            if edge_idx > 0
                comb_cap = comb_cap + G.Edges.Weight(edge_idx);
            end
        end
    end
    capacities = [capacities; comb_cap];
end

% Display the results side-by-side in a table for clarity
disp('Combined results table:');
T = table(used_combinations_cell, unused_combinations_cell, capacities);
T.Properties.VariableNames = {'U', 'V', 'C'};
disp(T);