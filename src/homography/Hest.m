function H = Hest(x1h, x2h)
    A = zeros(size(x1h,2)*2,9);
    n = size(x1h,2);
    A(1:2:n*2-1,:) = [zeros(n,3) -x1h(1,:)' -x1h(2,:)' -x1h(3,:)' x1h(1,:)'.*x2h(2,:)' x1h(2,:)'.*x2h(2,:)' x1h(3,:)'.*x2h(2,:)'];
    A(2:2:2*n,:) = [x1h(1,:)' x1h(2,:)' x1h(3,:)' zeros(n,3) -x1h(1,:)'.*x2h(1,:)' -x1h(2,:)'.*x2h(1,:)' -x1h(3,:)'.*x2h(1,:)'];
    
    [~,~,V] = svd(A);
    
    h1 = V(:,end);
    h1 = h1./h1(end);
    mask = [1 2 3;4 5 6;7 8 9];
    H = h1(mask);
end