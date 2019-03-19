function H = HestWithNormalization(x1, x2)
    if rank(x1) ~= 3 || rank(x2) ~= 3
        H = [];
        return;
    end
    [x1h, T1] = point_normalization(x1);
    [x2h, T2] = point_normalization(x2);

    A = zeros(size(x1h,2)*2,9);
    n = size(x1h,2);
    A(1:2:n*2-1,:) = [zeros(n,3) -x1h(1,:)' -x1h(2,:)' -x1h(3,:)' x1h(1,:)'.*x2h(2,:)' x1h(2,:)'.*x2h(2,:)' x1h(3,:)'.*x2h(2,:)'];
    A(2:2:2*n,:) = [x1h(1,:)' x1h(2,:)' x1h(3,:)' zeros(n,3) -x1h(1,:)'.*x2h(1,:)' -x1h(2,:)'.*x2h(1,:)' -x1h(3,:)'.*x2h(1,:)'];
%     try
        [~,~,V] = svd(A);
%     catch
%         error('ss');
%     end
    h1 = V(:,end);
    h1 = h1./h1(end);
    mask = [1 2 3;4 5 6;7 8 9];
    Hn = h1(mask);
    H = inv(T2)*Hn*T1;
end