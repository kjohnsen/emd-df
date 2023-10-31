function [x,P,stats] = solver_EMDDF_CVX(h,x0,D,C,beta,kappa,mu)

% solver_EMDDF_beckman.m
% 
% This function solves the EMD-DF problem under standard formulation:
%
% min_{x,P,u}
%         0.5   * sum_square(h - D*x)
%       + beta  * norm(x,1) 
%       + kappa * sum( C(:)'*P(:) )
%       - mu    * u
% s.t.  sum(P,2)  <= x
%       sum(P,1)' <= x0
%       sum(P(:)) == u
%       u         <= sum(x)
%       u         <= sum(x0)
% 
% Input definitions:
% h     Observation vector, size = m x 1
% x0    Vectorized prediction signal, size = n x 1
% D     Observation matrix, size = m x n
% C     Transport cost matrix, size = n x n
% beta  Sparsity parameter
% kappa Dynamical parameter
% mu    Transport capacity promoting parameter
% 
% Output definitions:
% x     Vectorized signal estimate, size = n x 1
% stats Struct containing optimization statistics
% P     Optimal transport (EMD) flows
% 
% Notes:
% 1) Software dependency: CVX (http://cvxr.com/cvx/)
% 2) mu should be set sufficiently high, so that the sufficient transport
%    capacity is present. A rough heuristic is mu = k*max(beta,kappa),
%    where k is in the range [1,5].
% 
% Copyright 2018. John Lee.
% Georgia Institute of Technology. Sensory Information Processing Lab.

N = numel(x0);

% Solve 
tic
cvx_begin quiet
% cvx_solver mosek
cvx_solver SDPT3
cvx_precision high
    variable x(N,1) nonnegative;
    variable P(N,N) nonnegative;
    variable u nonnegative;
    minimize( 0.5 * sum_square(h(:) - D*x) ...
              + beta * norm(x,1) ...
              + kappa * C(:)'*P(:) ...
              - mu*u ...
             );
    subject to
        sum(P,2)  <= x(:);
        sum(P,1)' <= x0(:);
        sum(P(:)) == u;
        u         <= sum(x(:));
        u         <= sum(x0(:));
cvx_end

stats.obj     =   0.5   * sum_square(h - D*x) ...
                + beta  * norm(x,1) ...
                + kappa * sum(C(:)'*P(:)) ...
                - mu    * u;
stats.runtime = toc;
stats.nbriter = cvx_slvitr;

end