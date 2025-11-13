function y = create_cut_r(V, U, index)

v=[];
for i=2:size(V,1)
    v = [v;V(i)];
end

n = length(v);

all_combinations_cell = {};

for k = 1:n
    C_k = nchoosek(v, k);
    C_k_cell = mat2cell(C_k, ones(size(C_k, 1), 1), size(C_k, 2));
    all_combinations_cell = [all_combinations_cell; C_k_cell];
end

% Display the result
disp('All non-empty subsets:');
disp(all_combinations_cell);

size(all_combinations_cell,1)

