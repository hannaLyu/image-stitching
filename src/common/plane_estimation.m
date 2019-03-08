
function params = plane_estimation(ph)
    A = ph';
    [~,~,V] = svd(A);
    params = V(:,end);
end
