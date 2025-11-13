function z = pageRank_centrality(P, beta, mi, z_0, tol)

% PAGE RANK

maxIter = 10000;

z_prev = z_0;
for iter = 1:maxIter
    z_current = (1-beta)*P*z_prev + beta*mi;

    if norm(z_current - z_prev, inf) < tol
        break;
    end

    z_prev = z_current;
end

if iter == maxIter
    warning('katz_centrality:NoConvergence', 'Maximum iterations reached without convergence');
end

% normalize z
z = z_current / sum(z_current);