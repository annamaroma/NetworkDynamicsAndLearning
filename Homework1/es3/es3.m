clear all
clc

% EXERCISE 3

% B = traffic.mat
% Le righe di B corrispondono ai nodi della rete
% Le colonne di B corrispondono agli archi della rete
% La i-esima colonna di B:
%   - contiene 1 nella riga corrispondente al nodo di partenza dell'arco e_i
%   - contiene -1 nella riga corrispondente al nodo di arrivo dell'arco e_i
% Ogni nodo rappresenta un'intersezione tra autostrade
% (e parte dell'area circostante)
B = load('traffic.mat');
B = B.traffic;

[n,m] = size(B); %n=#nodi m=#archi
s = zeros(1,m);
t = zeros(1,m);

for j = 1:m
    s(j) = find(B(:,j) == 1); %partenza
    t(j) = find(B(:,j) == -1); %arrivo
end

% ogni arco ha un maximum flow capacity C
C = load('capacities.mat');
C = C.capacities;

%ogni arco ha un tempo di percorrenza minimo che il guidatore impiega quando la strada è vuota
%cioè quando ci metto meno tempo? quando non c'è traffico
l = load('traveltime.mat');
l = l.traveltime;
w = l;

G = digraph(s,t,w);
names = {'n1','n2','n3','n4','n5','n6','n7','n8','n9','n10','n11','n12','n13','n14','n15','n16','n17'};
numberOfNodes = numnodes(G);
G.Nodes.Name = names(:);
plot(G);


% exercise 3.a
%find the shortest path between node 1 and node 17
[path, pathLength] = shortestpath(G, 'n1', 'n17', 'Method', 'positive');
fprintf('Shortest path from n1 to n17: \n');
disp(path);
fprintf('Path length: %d\n', pathLength);

% exercise 3.b
% find the maximum flow between node 1 and node 17
V = [1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17];
[U, V, C] = compute_cuts(V, G);

