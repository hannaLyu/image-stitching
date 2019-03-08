
function result = ransac_routine(x, ransac)
    iter = 1;
    n = size(x,2);
    while iter < ransac.maxiter
        % sample
        id = randperm(n, ransac.minimumset);
        % estimate
        params = ransac.estt_fun(x(:,id));
        % consensus
        errs = ransac.eval_fun(x, params);
        % verify
        inliers = errs < ransac.threshold;
        num_inlier = sum(inliers);
        % update
        if num_inlier > sum(ransac.inliers)
            ransac.inliers = inliers;
            pin = num_inlier / n;
            ransac.maxiter = round(log(1-ransac.pinlier)/log(1-pin^(ransac.minimumset)));
        end
        iter = iter + 1;
    end
    % refinement
    result.params = ransac.estt_fun(x(:,ransac.inliers));
    result.inliers = ransac.inliers;
end