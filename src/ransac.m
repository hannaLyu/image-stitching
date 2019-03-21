function inlierNum = ransac(X1,X2)
number1=size(X1,2);
number2=size(X2,2);
X1= tohomogeneous(X1);
X2= tohomogeneous(X2);
distance=zeros(1,number1);
sigma =1;
j=1;
pretotal=0;  %points which meet requirement
iter=500;
while (j<=iter)
   idx1 = randperm(number1,4);
   x1 = X1(:,idx1); 
  [T1,x1h] = normalization(x1);
   idx2 = randperm(number2,4);
   x2 = X2(:,idx2);
   [T2,x2h] = normalization(x2);
    H1=Hest(x1h,x2h);
   H=inv(T2)*H1*T1;
   
  for i=1:iter
   a=H* X1(:,i);
   b=inv(H)*X2(:,i);
   distance(:,i)=norm(X2(:,i)-a)+norm(X1(:,i)-b);
  end
   inlierIdx = find(distance<sigma);%find inlier which colse to line
   inlierNum = length(inlierIdx);  %caculate inlier'number
   if inlierNum>pretotal         %find best line  
      inlierNum=pretotal;
      
    end       
  
   j=j+1;
end




end

