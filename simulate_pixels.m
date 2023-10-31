function [x] = simulate_pixels(n,F,K,B)

% n = frame width
% F = number of frames
% K = number of targets
% B = EMD budget, i.e., maximum distance that pixels can move

N = n^2;

% % Special case: K == 0.5*N
% if K==0.5*N
%     x = zeros(N,F);
%     x(1:2:end,1:2:end) = 1;
%     x(2:2:end,2:2:end) = 1;
%     return;
% end
% 
% % Special case: K > 0.5*N
% % Strategy: then flip 1s and 0s
% if K >0.5*N, K_org = K; K = N - K; end

% Special case: K = N (paint everything 1)
if K == N, x = ones(N,F); return; end

for attempt = 1:10000 % max attempts
    
    % Reset
    empty_support_restart = 0;
    x = zeros(N,F);

    % Generate first frame
    x(randperm(N,K)',1) = 1;

    % Generate subsequent frames
    for f = 2:F
        targets = find(x(:,f-1)==1);
        for t = 1:K
            % Identify each target
            [tar_x,tar_y] = ind2sub([n,n],targets(t));
            % Find support for each target
            [rr cc] = meshgrid(1:n);
            S = sqrt((rr-tar_y).^2+(cc-tar_x).^2)<=B;
            % Remove existing targets from support
            S = max( S-reshape(x(:,f),n,n) , 0 );
            new_sup = find(S == 1);
            if isempty(new_sup), empty_support_restart = 1; break; end
            % Place target randomly in new support
            new_loc = new_sup(randperm(length(new_sup),1));
            x(new_loc,f) = 1;
        end
        if empty_support_restart == 1, break; end 
    end
    if empty_support_restart == 0, break; end % Terminate upon success
end

if empty_support_restart == 1, error(['Too dense! Attempts =' num2str(attempt)]); end

end