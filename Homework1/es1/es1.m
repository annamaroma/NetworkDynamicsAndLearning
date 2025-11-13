clear all
clc

s = [1,1,2,3,2,4];
t = [2,3,3,4,5,5];
w = [3,3,1,3,3,2]
names={'o','a','b','c','d'};
G = digraph(s,t,w,names);

% plot(G);

V = [1;2;3;4;5];

%% exercise 1.a
[U,V,C]=compute_cuts(V, G);

%minC = min(C);
vecIndex = find(C==min(C));
for i =1:size(vecIndex,1)
    fprintf('Minimum cut value: %d\n', min(C));
    U_minCut = U(vecIndex(i))
    V_minCut = V(vecIndex(i))
end

%% exercise 1.b
%first tempt choose x = 4 to be added to the arch cd , the one that both minCuts have in common

