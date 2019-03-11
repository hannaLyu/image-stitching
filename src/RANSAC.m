
X = gen_line_data(500);
 figure;
 plot(X(1,:),X(2,:),'o');hold on; % show points
number = size(X,2); % total points
sigma = 0.3;   %threshold
 pretotal=0;  %points which meet requirement
iter=1000;


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
%caculate distance between every point and current line  
distance=abs(line(1:2)*(X - repmat(sample(:,1),1,number)));
    inlierIdx = find(distance<sigma);%find inlier which colse to line
    inlierNum = length(inlierIdx);  %caculate inlier'number
     if inlierNum>pretotal           
         pretotal=inlierNum;
         bestline=line;        
    end  
 end

   distance=distance<sigma;    
hold on;
k=1;
for i=1:length(distance)
    if distance(i)
        inliers(1,k) = X(1,i);
        k=k+1;
        plot(X(1,i),X(2,i),'+');
    end
end
