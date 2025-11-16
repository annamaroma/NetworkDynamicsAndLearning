function z = katz_centrality(W,lambda, beta, mi, z_0, tol)

% KATZ CENTRALITY 

maxIter = 10000;
W_transpose = W';
z_prev = z_0;
alpha = (1-beta)/lambda;
for iter = 1:maxIter
    
    z_current = alpha*W_transpose*z_prev + beta*mi;

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

end

