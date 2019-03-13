function [result,linepara]=tohomogeneous(X)
number = size(X,2); % total points
sigma =0.1;   %threshold
 pretotal=0;  %points which meet requirement
iter=500;


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
         bestline=line;        
    end  
end

  
result=distance<sigma;
linepara=bestline;
end

