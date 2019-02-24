function I = fast_origion
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
img=imread('library2.jpg');
imshow(img)
close all;
[m,n]=size(img);
result=zeros(m,n);

t=40;   
for i=4:m-3
    for j=4:n-3
        Ic=img(i,j);    
        Inn=[img(i-3,j) img(i-3,j+1) img(i-2,j+2) img(i-1,j+3) img(i,j+3) img(i+1,j+3) img(i+2,j+2) img(i+3,j+1) img(i+3,j) img(i+3,j-1) img(i+2,j-2) img(i+1,j-3) img(i,j-3) img(i-1,j-3) img(i-2,j-2) img(i-3,j-1)];
            ind=find(abs(Inn-Ic)>t);         
            if length(ind)>=12
               result(i,j) = length(ind);      
            end
    end
end

% [row,col] = find(result==1);
% figure
% I=imshow(img)
% hold on;

% for i=4:m-3
%     for j=4:n-3
%         if result(i,j)~=0
%             if max(max(result(i-2:i+2,j-2:j+2)))==result(i,j)  
%                 [img(i-3,j), img(i-3,j+1), img(i-2,j+2), img(i-1,j+3), img(i,j+3), img(i+1,j+3), img(i+2,j+2), img(i+3,j+1), ...
%                  img(i+3,j), img(i+3,j-1), img(i+2,j-2), img(i+1,j-3), img(i,j-3), img(i-1,j-3), img(i-2,j-2), img(i-3,j-1)]= ...
%                  deal(255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255);
%             end
%         end
%     end
% end
result=nonmaxima_suppression(result, 3);
[row,col] = find(result==1);
imshow(img);hold on;
plot(col(:),row(:),'rx','MarkerSize',5)

end

