function result = ransac_routine_homo(x1, x2, ransac)
    iter = 1;
    n = size(x1,2);
    while iter < ransac.maxiter
        % sample
        id = randperm(n, ransac.minimumset);
        % estimate
        params = ransac.estt_fun(x1(:,id), x2(:,id));
        if isempty(params)
            continue;
        end
        % consensus
        errs = ransac.eval_fun(params, x1, x2);
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
    result.params = ransac.estt_fun(x1(:,ransac.inliers),x2(:,ransac.inliers));
    result.inliers = ransac.inliers;
end