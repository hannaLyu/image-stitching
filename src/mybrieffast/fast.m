function [corner,result]= fast(img)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% img=imread('library2.jpg');
% imshow(img)
% img = im2double(img);
[m n]=size(img);
result=zeros(m,n);

t=40;   
for i=4:m-3
    for j=4:n-3
        Ic=img(i,j);    
        Inn=[img(i-3,j) img(i-3,j+1) img(i-2,j+2) img(i-1,j+3) img(i,j+3) img(i+1,j+3) img(i+2,j+2) img(i+3,j+1) img(i+3,j) img(i+3,j-1) img(i+2,j-2) img(i+1,j-3) img(i,j-3) img(i-1,j-3) img(i-2,j-2) img(i-3,j-1)];
        if abs(Inn(1)-Ic)<t && abs(Inn(9)-Ic)<t
           continue; 
        end   
        I1_5_9_13=[abs(Inn(1)-Ic)>t abs(Inn(5)-Ic)>t abs(Inn(9)-Ic)>t abs(Inn(13)-Ic)>t];
        if sum(I1_5_9_13)>=3
            ind=find(abs(Inn-Ic)>t);         
            if length(ind)>=9
               result(i,j) = sum(abs(Inn-Ic)); 
           
            end
         
            
        end
    end
end
for i=4:m-3
    for j=4:n-3
        if result(i,j)~=0
            if max(max(result(i-2:i+2,j-2:j+2)))==result(i,j)               
                [img(i-3,j), img(i-3,j+1), img(i-2,j+2), img(i-1,j+3), img(i,j+3), img(i+1,j+3), img(i+2,j+2), img(i+3,j+1), ...
                 img(i+3,j), img(i+3,j-1), img(i+2,j-2), img(i+1,j-3), img(i,j-3), img(i-1,j-3), img(i-2,j-2), img(i-3,j-1)]= ...
                 deal(255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255);
            end
        end
    end
end

hsize = 5;
[x,y]=meshgrid(-hsize:hsize,-hsize:hsize);
x = x(:);y=y(:);
id = x == 0 & y == 0;
x(id) = []; y(id) = [];
 results = zeros(m,n);
for i = 4:m-3
    for j = 4:n-3
        if result(i,j) == 0
            continue;
        end
             nn = [i+x j+y];
        valid = nn(:,1) > 0 & nn(:,1) <= m & nn(:,2)>0 & nn(:,2) <= n;
        nn = nn(valid,:);
        ind0 = sub2ind([m,n],nn(:,1),nn(:,2));
        maxval = max(result(ind0));
        if result(i,j) > maxval
            results(i,j) = 1;
        else
            result(i,j) = 0;
        end
    end
end
% [row, col] = find(results == 1);
% figure
% imshow(img),title('fast corner with bresenham circle'),
% hold on
% plot(col,row, 'r*','MarkerSize',5),
% hold off
 [mm,nn]=size(result);
a=0;
for i=1:mm
    for j=1:nn
    if result(i,j)~=0
        a=a+1;
      corner(a,:)=[i,j];
      
    end
    end

end



end

