function [inliers,best]=tohomogeneous(X)
number = size(X,2); % total points
sigma =0.1;   %threshold
pretotal=0;  %points which meet requirement
iter=500;
type='plane';
if strcmp(type,'line')
for i=1:iter
   idx = randperm(number,2);
   sample = X(:,idx); 
   x = sample(:, 1);
   y = sample(:, 2);
  
    %coordinates of two points
   p1=[x(1);y(1);1];
   p2=[x(2);y(2);1];
    %line which two points decided
   line=cross(p1,p2); 
   line=transpose(line);
  distance=abs(line * [X;ones(1,size(X,2))] )./ sqrt(line(1)^2+line(2)^2);
    inlierIdx = find(distance<sigma);%find inlier which colse to line
    inlierNum = length(inlierIdx);  %caculate inlier'number
     if inlierNum>pretotal         %find best line  
         pretotal=inlierNum;
         best=line;   
         mask=distance;
    end  
end
elseif strcmp(type,'plane')
for i=1:iter
   idx = randperm(number,3);
   sample = X(:,idx); 
   x = sample(:, 1);
   y = sample(:, 2);
   z = sample(:, 3);
   p1=[x(1);y(1);z(1);1];
   p2=[x(2);y(2);z(2);1];
   p3=[x(3);y(3);z(3);1]; %caculate plane function
  A=transpose([p1 p2 p3]);
%   plane=null(A)
[~,~,V]=svd(A);
plane=V(:,end);
   plane=transpose(plane);
   distance=abs(plane*[X; ones(1,size(X,2))])./sqrt(plane(1)^2+plane(2)^2+plane(3)^2); 
   inlierIdx = find(distance<sigma);%find inlier which colse to line
   inlierNum = length(inlierIdx);  %caculate inlier'number
     if inlierNum>pretotal           
         pretotal=inlierNum;
         best=plane;
         mask=distance;
    end  
end
else
end
inliers=mask<sigma;

end

