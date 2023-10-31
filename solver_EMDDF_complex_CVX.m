function [s, F, stats] = solver_EMDDF_complex_CVX(y, x0, A, R, lambda, gamma, slack_mult)
%SOLVER_EMDDF_COMPLEX_CVX Complex variant of EMD-DF implemented with CVX
%
% Solve the complex variant of EMD-DF:
%   min_{z, F, u}   0.5    * ||y - As||_2^2
%                 + lambda * ||s||_1
%                 + gamma  * sum(R(:).*F(:))
%                 - mu     * u
%   subject to    sum(F, 1)' <= x0;
%                 sum(F, 2) <= zr_p + zr_n + zi_p + zr_n;
%                 sum(F(:)) == u;
%                 sum(zr_p + zr_n + zi_p + zr_n) >= u;
%                 sum(x0) >= u;
%
% Inputs
%   y            measurement vector
%   x0           prediction
%   A            measurement matrix
%   R            flow cost matrix
%   lambda       sparsity parameter
%   gamma        dynamics parameter
%   slack_mult   slack parameter (notation from paper: mu = slack_mult * gamma)
%
%
% Outputs
%   s           solution
%   F           EMD flow matrix
%   stats       optimization statistics
%
% Requires the CVX optimization package: http://cvxr.com/cvx/
%
% Author: Nicholas Bertrand
% Georgia Institute of Technology
% Sensory Information Processing Lab

A = [A, -A, 1i*A, -1i*A];

N = numel(x0);
ind = x0 > 0;
x0_active = abs(x0(ind));
K = sum(ind);
R = R(:, ind);

h_tic = tic();
cvx_begin quiet
    % cvx_precision low
    variable xp(N,1) nonnegative;
    variable xn(N,1) nonnegative;
    variable yp(N,1) nonnegative;
    variable yn(N,1) nonnegative;
    variable F(N,K)  nonnegative;
    variable u nonnegative;
    minimize( 0.5*sum_square_abs(y - A*[xp; xn; yp; yn]) ...
              + lambda*norm([xp; xn; yp; yn], 1)...
              ... + lambda*sum(norms([xp xn yp yn], 2, 2))...
              + gamma*R(:)'*F(:) ...
              - slack_mult*gamma*u );
    subject to
              sum(F, 1)' <= x0_active;
              sum(F, 2) <= xp + xn + yp + yn;
              sum(F(:)) == u;
              sum(xp + xn + yp + yn) >= u;
              sum(x0_active) >= u;
cvx_end

s = (xp - xn) + 1i*(yp - yn);

          
stats.obj     = 0.5*sum_square_abs(y - A*[xp; xn; yp; yn]) ...
                + lambda*norm([xp; xn; yp; yn], 1)...
                + gamma*R(:)'*F(:) ...
                - slack_mult*gamma*u;
stats.runtime = toc(h_tic);
stats.nbriter = cvx_slvitr;
