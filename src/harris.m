function [img] =harris
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
I=imread('checkerboard.jpg');
figure
imshow(I)
if size(I,3) == 3
    gray=rgb2gray(I);
else
    gray = I;
end
gray = im2double(gray);
hx= fspecial('sobel');
hy=hx';
gradx = filter2(hx,gray,'same'); 
grady=filter2(hy,gray,'same'); 
gradx = abs(gradx);
grady=abs(grady);
grad=gradx+grady;
subplot(2,2,1);imshow(gray);title('original');
subplot(2,2,2);imshow(gradx,[]);title(' horizon gradient');
subplot(2,2,3);imshow(grady,[]);title('vertical gradient');
subplot(2,2,4);imshow(grad,[]);title('sobel gradient');
% 2.
Ix=gradx.^2;
Iy=grady.^2;
Ixy=gradx.^grady;

% 3.
h=fspecial('gaussian',[3,1],1);
w=h*h';
A=imfilter(Ix,w);
B=imfilter(Iy,w);
C=imfilter(Ixy,w);

% 4.5
k=0.04;
RMax=0;
height=size(gray,1);
width=size(gray,2);
R=zeros(height,width);
for h=1:height
    for w=1:width
         M=[A(h,w) C(h,w);C(h,w) B(h,w)]; 
        R(h,w)=det(M) - k*(trace(M))^2;   
        if(R(h,w)>RMax)
            RMax=R(h,w);                  
        end
    end
end

% 6.

a=0.01;
R_corner=(R>=(a*RMax)).*R;
figure;imshow(R_corner);
hsize = 5;
[x,y]=meshgrid(-hsize:hsize,-hsize:hsize);
x = x(:);y=y(:);
id = x == 0 & y == 0;
x(id) = []; y(id) = [];
result = zeros(height,width);
for h = 1:height
    for w = 1:width
        if R_corner(h,w) == 0
            continue;
        end
        nn = [h+x w+y];
        valid = nn(:,1) > 0 & nn(:,1) <= height & nn(:,2)>0 & nn(:,2) <= width;
        nn = nn(valid,:);
        ind0 = sub2ind([height,width],nn(:,1),nn(:,2));
        maxval = max(R(ind0));
        if R(h,w) > maxval
            result(h,w) = 1;
        else
            R(h,w) = 0;
        end
    end
end

[row, col] = find(result == 1);

figure
img=imshow(I),title('my-Harris'),
hold on
plot(col,row, 'ro','MarkerSize',10),
hold off
frame = getframe(gca);
im2 = frame2im(frame);
[imind,cm] = rgb2ind(im2,256); 

end

