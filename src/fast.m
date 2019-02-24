function I= fast
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

img=imread('library2.jpg');
imshow(img)

[m n]=size(img);
result=zeros(m,n);

t=40;   
% for proof of concept, it is ok to use loop in matlab. However, for and while are extremely slow in Matlab.
% It is a good ideal to use vectorization whenever possible.
for i=4:m-3
    for j=4:n-3
        Ic=img(i,j);    
        Inn=[img(i-3,j) img(i-3,j+1) img(i-2,j+2) img(i-1,j+3) ...
             img(i,j+3) img(i+1,j+3) img(i+2,j+2) img(i+3,j+1) ...
             img(i+3,j) img(i+3,j-1) img(i+2,j-2) img(i+1,j-3) ...
             img(i,j-3) img(i-1,j-3) img(i-2,j-2) img(i-3,j-1)];
        if abs(Inn(1)-Ic)<t && abs(Inn(9)-Ic)<t
           continue; 
        end   
        I1_5_9_13=[abs(Inn(1)-Ic)>t abs(Inn(5)-Ic)>t abs(Inn(9)-Ic)>t abs(Inn(13)-Ic)>t];
        if sum(I1_5_9_13)>=3
            ind=find(abs(Inn-Ic)>t);         
            if length(ind)>=12
               result(i,j) = sum(abs(Inn-Ic));      
            end
        end
    end
end
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
result=nonmaxima_suppression(result, 4);
[row,col] = find(result==1);
imshow(img);hold on;
plot(col(:),row(:),'rx','MarkerSize',5)
% figure;
% I=imshow(img);
end

