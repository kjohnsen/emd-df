% demo.m
% 
% This script demonstates CVX solvers of 
% (1) the (standard) EMD-DF and
% (2) the Beckman-EMD-DF
% which were proposed in:
%     Earth Mover's Distance as a Dynamics Regularizer for Sparse Signal Tracking
%     N. Bertrand, A. Charles, J. Lee, P. Dunn and C.J. Rozell. June 2018. Submitted.
% 
% Copyright 2018. John Lee.
% Georgia Institute of Technology. Sensory Information Processing Lab.

clearvars

%% Define problem parameters

% Parameters
n = 24; % n x n temporal frame
nbr_frames = 2;
K = 0.05; % sparsity fraction
B = 1; % Maximum distance of pixels between frames
noise_sigma = 0.05; % zero-mean gaussian noise
M = 0.15; % measurement fraction

% Generate EMD matrices
[xg,yg] = meshgrid(1:n,1:n); % create coordintes
xg = xg(:); yg = yg(:); % vectorize
Xg = xg * ones(1,n^2); Yg = yg * ones(1,n^2);
r = sqrt((Xg-Xg').^2 + (Yg-Yg').^2); % create distance matrix (Euclidean)

% EMD-DF parameters
lambda = 2e-1;
kappa = 1e-2;
mu = 3e-1;

%% Simulate toy problem

N = n^2;

% Simulate Problem
rng(1); % make repeatible
x_gt = simulate_pixels(n, nbr_frames, ceil(K*N), B);
[y, A] = take_gaussian_meas(x_gt, ceil(M*N), noise_sigma^2);

%% Run Algorithms

% EMD-DF
[x_EMDDF,~,stats_EMDDF] = solver_EMDDF_CVX(y(:,2),x_gt(:,1),A(:,:,2),r,lambda,kappa,mu);

% Beckman's formulation of EMD-DF
[x_BmEMDDF,~,~,stats_BmEMDDF] = solver_EMDDF_beckman_CVX(y(:,2),x_gt(:,1),A(:,:,2),lambda,kappa,mu);

%% Compare

compute_rMSE = @(x,x_gt) norm(x(:)-x_gt(:))^2/norm(x_gt(:))^2;
disp(['EMD-DF: relative MSE (rMSE) = ' num2str(compute_rMSE(x_EMDDF,x_gt(:,2))) ', runtime = ' num2str(stats_EMDDF.runtime)]);
disp(['Beckman EMD-DF: relative MSE (rMSE) = ' num2str(compute_rMSE(x_BmEMDDF,x_gt(:,2))) ', runtime = ' num2str(stats_BmEMDDF.runtime)]);

compute_RMSE = @(x1,x2) sqrt(sum((x1(:)-x2(:)).^2)/numel(x1));
disp(['Difference between solutions, Root MSE (RMSE) = ' num2str(compute_RMSE(x_EMDDF,x_BmEMDDF))]);

figure(1);
subplot(221); imagesc(reshape(x_gt(:,1),n,n)); axis square; title('Prediction'); colorbar;
subplot(222); imagesc(reshape(x_gt(:,2),n,n)); axis square; title('Ground Truth'); colorbar;
subplot(223); imagesc(reshape(x_EMDDF,n,n)); axis square; title(['EMD-DF, rMSE=' num2str(compute_rMSE(x_EMDDF,x_gt(:,2))) ', runtime = ' num2str(stats_EMDDF.runtime)]); colorbar;
subplot(224); imagesc(reshape(x_BmEMDDF,n,n)); axis square; title(['Beckman-EMD-DF, rMSE=' num2str(compute_rMSE(x_BmEMDDF,x_gt(:,2))) ', runtime = ' num2str(stats_BmEMDDF.runtime)]); colorbar;


%% Auxiliary Functions

function [y, G] = take_gaussian_meas(x, M, noise_var)
N = size(x, 1);
T = size(x, 2);
G = randn(M, N, T)/sqrt(M);
y = zeros(M, T);
for kk = 1:T
    y(:, kk) = G(:, :, kk)*x(:, kk) + sqrt(noise_var)*randn(M, 1);
end
end