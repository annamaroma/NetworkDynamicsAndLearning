clear all
clc

s = [1,1,2,3,2,4];
t = [2,3,3,4,5,5];
w = [3,3,1,3,3,2];
names={'o','a','b','c','d'};
G = digraph(s,t,w,names);

% plot(G);

nodes = [1;2;3;4;5];

%% exercise 1.a
% compute cuts
[U,V,C] = compute_cuts(nodes, G);

%find and display minimum cut value and corresponding sets U and V
minC = min(C);
fprintf('Minimum cut value: %d\n', minC);

vecIndex = find(C==minC);
for i =1:size(vecIndex,1)
    U_minCut = U{vecIndex(i)};
    V_minCut = V{vecIndex(i)};
end

%% exercise 1.b

maxIter = 50;
%min_cut = extra_units_capacity(nodes, G, maxIter);

%fprintf('\nmin_cut history (per x>0):\n');
%disp(min_cut);

% bar plot to visualize final min_cut history
%bar(min_cut);
%xlabel('x value');
%ylabel('Minimum Cut Value');
%title('Minimum Cut Value over x>0');

%% exercise 1.c
% add a new link with unit capacity from node 3 to node 5
s = [1,1,2,3,2,4,1];
t = [2,3,3,4,5,5,5];
w = [3,3,1,3,3,2,1];
names={'o','a','b','c','d'};
G = digraph(s,t,w,names);

%compute cuts again
[U,V,C] = compute_cuts(nodes, G);

%find the new minimum cut value and corresponding sets U and V
minC = min(C);
fprintf('Minimum cut value: %d\n', minC);

vecIndex = find(C==minC);
for i =1:size(vecIndex,1)
    U_minCut = U{vecIndex(i)};
    V_minCut = V{vecIndex(i)};
end

maxIter = 100;
min_cut = extra_units_capacity(nodes, G, maxIter);

%fprintf('\nmin_cut history (per x>0):\n');
%disp(min_cut);

% bar plot to visualize final min_cut history
%bar(min_cut);
%xlabel('x value');
%ylabel('Minimum Cut Value');
%title('Minimum Cut Value over x>0');