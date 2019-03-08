
function params = line_estimation(ph)
    if size(ph,2) == 2
        % minimum set, use cross product
        params = cross(ph(:,1),ph(:,2));
    else
        A = ph';
        [~,~,V] = svd(A);
        params = V(:,end);
    end
end