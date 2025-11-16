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


%% exercise 3.a
%find the shortest path between node 1 and node 17
[path, pathLength] = shortestpath(G, 'n1', 'n17', 'Method', 'positive');
fprintf('Shortest path from n1 to n17: \n');
disp(path);
fprintf('Path length: %d\n', pathLength);

%% exercise 3.b
% find the maximum flow between node 1 and node 17
G_cap = digraph(s, t, C);
G_cap.Nodes.Name = names(:);
mf = maxflow(G_cap, 1, 17);
fprintf('Maximum flow from n1 to n17: %d\n', mf);

%% exercise 3.c
% Given the flow vector in flow.mat, compute the vector ν satisfying Bf = ν.
% In the following, we assume that the exogenous inflow is zero in all the nodes except for node 1,
% for which ν1 has the same value computed in the point (c), and node 17, for which ν17 = −ν1.
f = load('flow.mat');
f = f.flow;
nu = B * f;
disp('Vector nu satisfying Bf = nu:');
disp(nu);

%% exercise 3.d 
objective = @(f) sum(f .* l ./ (1 - f./C)); % Objective function

f_init = f; 

% Constraints
Aeq = B;
beq = nu;
% Inequality constraints
A = [];
b = [];
% Bounds: 0 <= f < C
lb = zeros(m,1);
epsilon = 1e-7;
ub = C - epsilon;

options = optimoptions('fmincon','Display','off','Algorithm','interior-point');
[f_star_fmincon, cost] = fmincon(objective, f_init, A, b, Aeq, beq, lb, ub, [], options);

fprintf('Social optimum f*:\n');
disp(f_star_fmincon);
fprintf('Total cost (total travel time):\n');
disp(cost);

%% exercise 3.e
fprintf('Computing Wardrop equilibrium with fmincon...\n');
wardrop_objective = @(f) sum(-C.*l.*log(1-f./C));

options = optimoptions('fmincon','Display','off','Algorithm','interior-point');
[f_wardrop_fmincon, wardrop_cost] = fmincon(wardrop_objective, f_init, A, b, Aeq, beq, lb, ub, [], options);

fprintf('Wardrop equilibrium f^(0):\n');
disp(f_wardrop_fmincon);

% The value of the objective function for Wardrop is not the travel time
fprintf('Value of Wardrop objective function:\n');
disp(wardrop_cost);

% To compare with social optimum, we need to compute the total travel time for the Wardrop equilibrium flow
total_travel_time_wardrop_fmincon = sum(f_wardrop_fmincon .* l ./ (1 - f_wardrop_fmincon./C));
fprintf('Total travel time for Wardrop equilibrium flow:\n');
disp(total_travel_time_wardrop_fmincon);

%% exercise 3.f
fprintf('\n--- Exercise 3.f (marginal-cost tolls) ---\n');
% Introduce tolls such that the toll on link e is omega_e = f*_e * tau'_e(f*_e)
% where tau_e(f) = l_e / (1 - f/C_e).

f_star = f_star_fmincon;

% compute derivative tau'(f) = l ./ (C .* (1 - f./C).^2)
tau_prime = l ./ (C .* (1 - f_star./C).^2);

% toll per unit flow (constant) at social optimum
omega = f_star .* tau_prime;

% Now compute the Wardrop equilibrium when users perceive cost = tau(f) + omega_e (omega is constant)
% Use the potential: sum_e 
%   Phi_e(f_e) = \int_0^{f_e} tau_e(s) ds = -C*l*log(1 - f/C)
% So objective = sum_e Phi_e(f_e) + sum_e omega_e * f_e

wardrop_with_tolls_obj = @(x) sum(-C .* l .* log(1 - x ./ C) + omega .* x);

% initial guess: previous Wardrop or f_init
f0_wardrop = f_wardrop_fmincon;
if any(f0_wardrop <= 0) || any(f0_wardrop >= ub)
    f0_wardrop = f_init;
end

options = optimoptions('fmincon','Display','off','Algorithm','interior-point');
[f_wardrop_tolls, wcost_tolls] = fmincon(wardrop_with_tolls_obj, f0_wardrop, A, b, Aeq, beq, lb, ub, [], options);

fprintf('\nWardrop equilibrium with marginal-cost tolls:\n');
disp(f_wardrop_tolls);

% Compute total travel times (without tolls) for comparison
total_travel_time_wardrop_tolls = sum(f_wardrop_tolls .* l ./ (1 - f_wardrop_tolls ./ C));
fprintf('Total travel time for Wardrop with tolls: %.8f\n', total_travel_time_wardrop_tolls);
fprintf('Total travel time for social optimum (previously computed): %.8f\n', cost);

if diff_norm < 1e-6
    fprintf('Wardrop with tolls matches the social optimum within tolerance.\n');
else
    fprintf('Wardrop with tolls differs from social optimum (see norm above).\n');
end

%% exercise 3.g
fprintf('\n--- Exercise 3.g (system cost = total additional travel time) ---\n');
% Define system cost psi_e(f_e) = f_e*(tau_e(f_e) - l_e)
objective_psi = @(x) sum( x .* ( l ./ (1 - x ./ C) - l ) );

% initial guess: use previous social optimum
f0_psi = f_star_fmincon;
if any(f0_psi <= lb) || any(f0_psi >= ub)
    f0_psi = f_init;
end

options = optimoptions('fmincon','Display','off','Algorithm','interior-point');
[f_star_psi, cost_psi] = fmincon(objective_psi, f0_psi, A, b, Aeq, beq, lb, ub, [], options);

fprintf('System optimum f* (for additional travel time objective):\n');
disp(f_star_psi);
fprintf('Total additional travel time (system cost) at f*: %.8f\n', cost_psi);

% Construct toll vector omega* = f*_e * tau'_e(f*_e) - l_e
tau_prime_psi = l ./ (C .* (1 - f_star_psi ./ C).^2);
omega_star = f_star_psi .* tau_prime_psi - l;

fprintf('Constructed tolls omega* (first 10 edges):\n');
for i=1:min(10,length(omega_star))
    fprintf('  edge %d: omega* = %.8f\n', i, omega_star(i));
end

% Compute Wardrop equilibrium under tolls omega* (minimize potential + omega'*f)
wardrop_with_omega_obj = @(x) sum(-C .* l .* log(1 - x ./ C) + omega_star .* x);

f0_w_omega = f_wardrop_fmincon;
if any(f0_w_omega <= lb) || any(f0_w_omega >= ub)
    f0_w_omega = f_init;
end

[f_wardrop_omega, wcost_omega] = fmincon(wardrop_with_omega_obj, f0_w_omega, A, b, Aeq, beq, lb, ub, [], options);

fprintf('\nWardrop equilibrium under constructed tolls:\n');
disp(f_wardrop_omega);

% Compare system cost values (additional travel time)
system_cost_wardrop_omega = objective_psi(f_wardrop_omega);
fprintf('System cost (additional travel time) at Wardrop with tolls: %.8f\n', system_cost_wardrop_omega);
fprintf('System cost at system optimum: %.8f\n', cost_psi);

if norm_diff_psi < 1e-6
    fprintf('Wardrop under constructed tolls matches the system optimum within tolerance.\n');
else
    fprintf('Wardrop under constructed tolls differs from system optimum (see norm above).\n');
end

