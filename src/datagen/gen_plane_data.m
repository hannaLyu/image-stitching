
function [X,Xi,Xo] = gen_plane_data(N)
    % inilers percentage
    p = 0.75;
    % noise
    sigma = 0.01;
    
    % make it pseudo-random
    rand('twister', 2222);
    randn('state', 2222);

    % build a plane passing through three points 
    X = 2*(rand(3, 3)-0.5);
    v1 = X(:, 2) - X(:, 1); 
    v1 = v1/norm(v1);
    v2 = X(:, 3) - X(:, 1);
    v2 = v2/norm(v2);
    X0 = 2*(rand(3, 1)-0.5);

    % generate a set of points spread in a cube
    Ni = round(p*N);
    No = N-Ni;

    % inliers
    lambda1 = 2*(rand(1, Ni)-0.5);
    lambda2 = 2*(rand(1, Ni)-0.5);
    Xi = repmat(X0, [1 Ni]) + ...
        repmat(lambda1, [3 1]).*repmat(v1, [1 Ni]) + ...
        repmat(lambda2, [3 1]).*repmat(v2, [1 Ni]); 

    % and add some noise
    Xi = Xi + sigma*randn(3, Ni);

    % outliers
    Xo = 2*(rand(3, No)-0.5) + repmat(X0, [1 No]);
    X = [Xi Xo];

    % scrample (just in case...)
    [dummy ind] = sort(rand(1, N));
    X = X(:, ind);
end
