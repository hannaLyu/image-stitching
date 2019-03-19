function error = reprojectionError(H,p1,p2)
    error = zeros(1,size(p1,2));
    %% assume p1 and p2 are already in homogeous coordianates
    p1proj = H*p1;
    p2proj = H\p2;
    %% to inhomogeneous coordinate
    p1proj = p1proj ./ p1proj(3,:);
    p2proj = p2proj ./ p2proj(3,:);

    err1 = p1proj(1:2,:) - p2(1:2,:);   
    err2 = p2proj(1:2,:) - p1(1:2,:);

    err1_norm = sqrt(err1(1,:).^2+err1(2,:).^2);%
    err2_norm = sqrt(err2(1,:).^2+err2(2,:).^2);%

    error = err1_norm + err2_norm;
end