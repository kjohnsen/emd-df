function x = solver_BPDN_CVX(y, A, lambda)
N = size(A, 2);
cvx_begin quiet
cvx_solver SDPT3
cvx_precision high
    variable x(N, 1) complex;
    minimize( 0.5 * sum_square_abs(y(:) - A*x) ...
              + lambda * norm(x,1) ...
             );
cvx_end