function [inlier1,inlier2,inlierId] = ransac(X1,X2)
number1=size(X1,2);
number2=size(X2,2);
prob=0.99;
mini=4;
iter=1;
sigma =6;
pretotal=0;
itereal=1e3;
while iter<itereal
    idx1 = randperm(number1,4);
    x1 = X1(:,idx1);
    [T1,x1h] = normalization(x1);
%     idx2 = randperm(number2,mini);
    x2 = X2(:,idx1);
    [T2,x2h] = normalization(x2);
    H1=Hest(x1h,x2h);
    H=inv(T2)*H1*T1;
    
    X1t = H*X1;
    X2t = H\X2;
    X1t = X1t./X1t(3,:);
    X2t = X2t./X2t(3,:);
    distance = sqrt(diag((X1t - X2)'*(X1t - X2)))+ sqrt(diag((X2t - X1)'*(X2t - X1)));
    distance=transpose(distance);
    inlierIdx = find(distance<sigma);%find inlier which colse to line
    inlierNum = length(inlierIdx);  %caculate inlier'number
    if inlierNum>pretotal         %find best line
        pretotal=inlierNum;
        inlierId=inlierIdx;
        pin=inlierNum/number1;
        itereal = round(log(1-prob)/log(1-pin^(mini)));
        disp(H);
    end
  
iter=iter+1;
end
inlier1= X1(:,inlierId(1,:));
inlier2=X2(:,inlierId(1,:));


   
end






