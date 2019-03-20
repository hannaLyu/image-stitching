clc; close all; clear all;
% assignment for homography estimation using DLT
cd(fileparts(which(mfilename)));

addpath(genpath('./'));

%% task 1
N = 4;
x1 = rand(2,N);
x1h = tohomogeneous(x1);

H = [1, 2, 3;4, 5, 6;3,8,9];
x2 = H*x1h;
x2h = x2 ./ x2(3,:);

% implement this your self
H1 = Hest(x1h, x2h);

%% task 2
% play with real image

% im1 = imread('../data/Tiles_perspective_distort.png');
% im1 = im2double(im1);
% imshow(im1);
% [x,y] = ginput(4);
% x1 = round([x y]');
% x1h = tohomogeneous(x1);
% 
% im2 = imread('../data/Tiles_perspective_undistort.png');
% im2 = im2double(im2);
% imshow(im2);
% [x,y] = ginput(4);
% x2 = round([x y]');
% x2h = tohomogeneous(x2);
% 
% % implement this your self
% H1 = Hest(x1h, x2h);
% 
% % im3 = warpping(im1,H1);
% % tform = projective2d(H1');
% % im3 = imwarp(im1,tform);
% im3 = WarpNViewMod(H1,im1,im2);
% 
% figure
% imshow(im3);

%% task 3: normalization
[x1hn, T1] = point_normalization(x1h);
[x2hn, T2] = point_normalization(x2h);
Hn = Hest(x1hn, x2hn);
H3 = inv(T2)*Hn*T1;

%% task 4: ransac+homography+normalization
N = 20;
x1 = rand(2,N);
x1h = tohomogeneous(x1);

H = [1, 2, 3;4, 5, 6;3,8,9];
x2 = H*x1h;
x2h = x2 ./ x2(3,:);

x3h = rand(2,10);
x4h = rand(2,10);
x3h = tohomogeneous(x3h);
x4h = tohomogeneous(x4h);
x5h = [x1h x3h];
x6h = [x2h x4h];

ransac.pinlier = 0.99;
ransac.estt_fun = @HestWithNormalization;%plane_estimation
ransac.eval_fun = @reprojectionError;%dist2plane
ransac.maxiter = 1e6;
ransac.threshold = 0.1;
ransac.inliers = [];
ransac.minimumset = 4;
result = ransac_routine_homo(x5h, x6h, ransac);

%% task 5: work with real images
% addpath to the data folder
addpath ../data/

addpath(genpath('./features'));
addpath(genpath('./matching'));


% load an image
im1 = imread('../data/pier/1.jpg');
im2 = imread('../data/pier/2.jpg');
im3 = imread('../data/pier/3.jpg');

[s11,s12,s13] = size(im1);
[s21,s22,s23] = size(im2);
[s31,s32,s33] = size(im3);

% display an image
% figure;imshow(im1)

% convert img from rgb to gray scale
if size(im1,3) == 3
    im1gray = rgb2gray(im1);
else
    im1gray = im1;
end

if size(im2,3) == 3
    im2gray = rgb2gray(im2);
else
    im2gray = im2;
end

if size(im3,3) == 3
    im3gray = rgb2gray(im3);
else
    im3gray = im3;
end

% convert image from uint8 to double, make nothing but convert 255 - 1, 0-0
im1d = im2double(im1gray);
im2d = im2double(im2gray);
im3d = im2double(im3gray);

% fast1 = fast_corner_detector(im1d, 500);
% fast2 = fast_corner_detector(im2d, 500);
% fast3 = fast_corner_detector(im3d, 500);

fast1 = Harris(im1d, 1000);
fast2 = Harris(im2d, 1000);
fast3 = Harris(im3d, 1000);

% brief_pattern;
load('./features/brief_pattern256.mat')
descriptors1 = extract_brief_descriptor(im1d,fast1,pattern);
descriptors2 = extract_brief_descriptor(im2d,fast2,pattern);
descriptors3 = extract_brief_descriptor(im3d,fast3,pattern);

matchingpairs = brief_matching(descriptors2, descriptors3);
id = ~isinf(matchingpairs(:,2));

ransac.pinlier = 0.99;
ransac.estt_fun = @HestWithNormalization;%plane_estimation
ransac.eval_fun = @reprojectionError;%dist2plane
ransac.maxiter = 1e3;
ransac.threshold = 6;
ransac.inliers = [];
ransac.minimumset = 4;
x1h = tohomogeneous(fast2(matchingpairs(id,1),[2,1])');
x2h = tohomogeneous(fast3(matchingpairs(id,2),[2,1])');
result = ransac_routine_homo(x1h, x2h, ransac);
x1h = x1h(:,result.inliers);
x2h = x2h(:,result.inliers);

imshow1 = cat(2, im2, im3);
figure;imshow(imshow1);hold on;

plot(fast2(:,2),fast2(:,1), 'ro','MarkerSize',5);
plot(fast3(:,2)+size(im1d,2),fast3(:,1), 'bo','MarkerSize',5);


shift = size(im1d,2);
cmap = jet(32);
k = 1;
for i = 1:size(x1h,2)
    ptdraw = [x1h(2,i), x1h(1,i);
              x2h(2,i), x2h(1,i)+shift];
    plot(ptdraw(:,2),ptdraw(:,1),'LineStyle','-','LineWidth',0.5,'Color',cmap(k,:));
    k = mod(k+1,32);if k == 0 k = 1;end
end


im12 = WarpAndBlend(result.params,im2,im3);
figure;imshow(im12);

