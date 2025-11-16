clear all
clc

% EXERCISE 2

s = [1,1,1,1,1,2,2,2,2,3,3,3,4,4,5,6,6,7,8,9,9,9,9,9];
t = [2,3,4,5,6,3,4,5,6,4,5,6,5,6,6,7,15,8,9,10,11,12,13,14];

names={'n1','n2','n3','n4','n5','n6','n7','n8','n9','n10','n11','n12','n13','n14','n15'};
G = graph(s,t);
G.Nodes.Name = names(:);
numberOfNodes = numnodes(G);
W = adjacency(G);

%compute transition matrix
somma_colonne = sum(W,1);
P = W./somma_colonne; 
P(:,somma_colonne==0) = 0; 


%plot(G);

%% exercise 2.a
% compute Katz centrality 
beta = 0.15;
mi = ones(numberOfNodes,1)/numberOfNodes;
z_0 = mi;
tol = 1e-9;

% calculate spectral radius 
lambda = abs(eigs(W, 1, 'largestabs'));


z_katz = katz_centrality(W, lambda, beta, mi, z_0, tol);
fprintf('Katz centrality computed:\n');
for i = 1:length(z_katz)
    fprintf('Node %d: %.6f\n', i, z_katz(i));
end

%% exercise 2.b
z_pageRank = pageRank_centrality(P, beta, mi, z_0, tol);
fprintf('PageRank centrality computed:\n');
for i = 1:length(z_pageRank)
    fprintf('Node %d: %.6f\n', i, z_pageRank(i));
end

%grafico che plotta  
%figure; 
%bar(z_pageRank); 
%ylabel('PageRank centrality value'); 
%title('PageRank Centrality Values for each Node'); 
%xticks(1:numberOfNodes); 
%xticklabels(G.Nodes.Name);

%% exercise 2.c
%barplot that shows z_katz and z_pageRank
Z = [z_katz(:) z_pageRank(:)];

%figure;
%bar(Z, 'grouped');
%ylabel('centrality reached');
%title('Comparison between Katz and PageRank Centrality Algorithms');
%legend('Katz Centrality Vector', 'PageRank Centrality Vector');
%xticks(1:numberOfNodes);

%% exercise 2.d
beta_vec = [0 1/4 1/2 3/4 1];
beta_len = length(beta_vec);
Z_out = zeros(numberOfNodes, beta_len);

for i = 1:beta_len
    beta = beta_vec(i);
    z = pageRank_centrality(P,beta,mi,z_0,tol);
    Z_out(:, i) = z;
end

Z_out;
 
figure;
bar(Z_out, 'grouped');
ylabel('centrality reached');
title('PageRank Centrality Algorithm with different beta');
legend('beta = 0', 'beta = 1/4', 'beta = 1/2', 'beta = 3/4', 'beta = 1');
xticks(1:numberOfNodes);


