
function X = gen_line_data(N)
    % inilers percentage
    p = 0.15;
    % noise
    sigma = 0.05;

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data Generation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % make it pseudo-random
    rand('twister', 2222);
    randn('state', 2222);

    % line parameters
    m = -0.8;
    q = 0.3;

    % generate a set of points 
    Ni = round(p*N);
    No = N-Ni;

    % inliers
    X1i = 3 * 2*(rand(1, Ni)-0.5);
    X2i = m*X1i+q;

    % and add some noise
    X1i = X1i + sigma*randn(1, Ni);
    X2i = X2i + sigma*randn(1, Ni);

    % outliers
    X1o = 3 * 2*(rand(1, No)-0.5);
    X2o = 3 * 2*(rand(1, No)-0.5);

    X1 = [X1i X1o];
    X2 = [X2i X2o];

    % scrample (just in case...)
    [dummy ind] = sort(rand(1, N));
    X1 = X1(:, ind);
    X2 = X2(:, ind);

    % form the input data pairs
    X = [X1; X2];
end
