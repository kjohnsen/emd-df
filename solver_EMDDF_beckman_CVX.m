function [x,m_x,m_y,stats] = solver_EMDDF_beckman_CVX(h,x0,D,beta,kappa,mu)

% solver_EMDDF_beckman.m
% 
% This function solves the EMD-DF problem under a Beckman formulation:
%
% min_{x,m_x,m_y,v0,v1,u}
%         0.5   * sum_square(h - D*x)
%       + beta  * norm(x,1) 
%       + kappa * sum( norms([m_x';m_y'],2) )
%       - mu    * u
% s.t.  D_x*m_x + D_y*m_y + v1(:) - v0(:) == 0;
%       sum(v0) == u;
%       sum(v1) == u;
%       v1 <= x(:);
%       v0 <= x0(:);
%       u <= sum(x0(:));
%       u <= sum(x(:));
% 
% Input definitions:
% h     Observation vector, size = m x 1
% x0    Vectorized prediction signal, size = n x 1
% D     Observation matrix, size = m x n
% beta  Sparsity parameter
% kappa Dynamical parameter
% mu    Transport capacity promoting parameter
% 
% Output definitions:
% x     Vectorized signal estimate, size = n x 1
% stats Struct containing optimization statistics
% m_x   Horizontal transport flux
% m_y   Vertical transport flux
% 
% Notes:
% 1) Software dependency: CVX (http://cvxr.com/cvx/)
% 2) D_x and D_y are defined as linear divergence operators acting on
%    spatial flux terms. We employ zero-flux boundary conditions.
% 3) mu should be set sufficiently high, so that the sufficient transport
%    capacity is present. A rough heuristic is mu = k*max(beta,kappa),
%    where k is in the range [1,5].
% 
% Copyright 2018. John Lee.
% Georgia Institute of Technology. Sensory Information Processing Lab.

N = numel(x0);
n = sqrt(N); % assumes input is square image

% Create divergence matrices
D_x = eye(n^2) - circshift(eye(n^2),1); D_x(1,n^2) = 0; D_x(:,n:n:n^2) = zeros(n^2,n);
D_y = eye(n^2) - circshift(eye(n^2),n); D_y(:,end-n+1:end) = zeros(n^2,n);

% Solve 
tic
cvx_begin quiet
% cvx_solver mosek
cvx_solver SDPT3
cvx_precision high
    variable x(N,1) nonnegative;
    variable m_x(N,1);
    variable m_y(N,1);
    variable v0(N,1) nonnegative;
    variable v1(N,1) nonnegative;
    variable u nonnegative;
    minimize( 0.5 * sum_square(h(:) - D*x) ...
              + beta * norm(x,1) ...
              + kappa * sum(norms([m_x';m_y'],2))...
              - mu*u ...
             );
    subject to
        D_x*m_x + D_y*m_y + v1(:) - v0(:) == 0;
        sum(v0) == u;
        sum(v1) == u;
        v1 <= x(:);
        v0 <= x0(:);
        u <= sum(x0(:));
        u <= sum(x(:));
cvx_end

stats.obj     =   0.5   * sum_square(h - D*x) ...
                + beta  * norm(x,1) ...
                + kappa * sum(norms([m_x';m_y'],2)) ...
                - mu    * u;
stats.runtime = toc;
stats.nbriter = cvx_slvitr;

end