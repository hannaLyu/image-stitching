clc;close all;clear all;

%% the first example code to load image, display image, and do some basic image filtering

% addpath to the data folder
addpath ../data/

% load an image
im1 = imread('checkerboard.jpg');

% display an image
figure;imshow(im1)

% convert img from rgb to gray scale
if size(im1,3) == 3
    imgray = rgb2gray(im1);
else
    imgray = im1;
end
figure;imshow(imgray);

% convert image from uint8 to double, make nothing but convert 255 - 1, 0-0
imdouble = im2double(imgray);
figure;imshow(imdouble);

% img filtering
imhedge = edge(imdouble,'Sobel',[],'horizontal');
imvedge = edge(imdouble,'Sobel',[],'vertical');

% combine two images to one large image along column direction, 1-row
% direction, 2-column direction, 3-depth direction
imcombine = cat(2, imhedge,imvedge);
figure;imshow(imcombine);

% apply a Gaussian filter, just image bluring, hsize is the size of the
% gaussian kernel, 1 is the sigma (recall the 2D Gaussian function.)
g = fspecial('gaussian', 3, 1);% try play with this two values to see what will hapen
imblur = imfilter(imgray, g, 'replicate');
figure;imshow(imblur);

% save image
output_path = '../results/';
filename = strcat(output_path,'edge.jpg');
imwrite(imcombine,filename);

% harris
Ix = edge(imdouble,'Sobel',[],'horizontal');
Iy = edge(imdouble,'Sobel',[],'vertical');

Ixx = Ix.*Ix;
Iyy = Iy.*Iy;
Ixy = Ix.*Iy;

g = fspecial('gaussian', 3, 1);% try play with this two values to see what will hapen
Ixx = imfilter(Ixx,g,'replicate');
Iyy = imfilter(Iyy,g,'replicate');
Ixy = imfilter(Ixy,g,'replicate');

C = Ixx.*Iyy - Ixy.^2 - 0.04.*(Ixx+Iyy).^2;
% figure;imshow(C,[]);

threshold = 0.3*max(C(:));
Cb = imbinarize(C,threshold);
figure;imshow(Cb,[]);

% nonmaxima suppression
hsize = 5;
[X,Y]=meshgrid(-hsize:hsize,-hsize:hsize);
X = X(:);
Y = Y(:);
id = X == 0 & Y == 0;
X(id) = [];Y(id) = [];

corners = zeros(size(C,1)*size(C,2),2);
K = 0;
for i = 1:size(C,1)
    Xi = X + i;
    vXi = Xi > 0 & Xi <= size(C,1);
    for j = 1:size(C,2)
        if Cb(i,j) ~= 0
            Yi = Y + j;
            vYi = Yi > 0 & Yi <= size(C,2);
            v = vXi & vYi;
            ind0 = sub2ind(size(C),Xi(v),Yi(v));
            maxv = max(C(ind0));
            if maxv < C(i,j)
                K = K + 1;
                corners(K,:) = [i,j];
            else
                C(i,j) = 0;
            end
        end
    end
end
corners(K+1:end,:)=[];
figure;imshow(imdouble,[]);hold on;
dtured = [153/255,0,0];
plot(corners(:,2),corners(:,1),'o','MarkerSize',5,'MarkerEdgeColor',dtured);

