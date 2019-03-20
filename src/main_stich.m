clc;clear all;close all;
% load an image
im1 = imread('../data/hill/1.jpg');
im2 = imread('../data/hill/2.jpg');
im3 = imread('../data/hill/3.jpg');


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

img1 = single(im1gray);
img2 = single(im2gray);
img3 = single(im3gray);

[feature1,descriptor1] = vl_sift(img1) ;
[feature2,descriptor2] = vl_sift(img2) ;
[feature3,descriptor3] = vl_sift(img3) ;

% tic
% [matches, scores] = vl_ubcmatch(descriptor1, descriptor2,1.5);
% toc
% knn matching
tic
ds2 = single(descriptor2);
ds1 = single(descriptor1);
kdtree = vl_kdtreebuild(ds2);
matches = zeros(size(descriptor1,2),2);
k = 1;
thresh = 2;
for i = 1:size(descriptor1,2)
    [index, distance] = vl_kdtreequery(kdtree, ds2, ds1(:,i), ...
                                       'NumNeighbors', 2, 'MaxComparisons', 15);
    if distance(1)*thresh > distance(2) continue; end
    matches(k,:) = [i,index(1)];
    k = k + 1;
end
matches(k:end,:) = [];
matches = matches';
toc

ransac.pinlier = 0.99;
ransac.estt_fun = @HestWithNormalization;%plane_estimation
ransac.eval_fun = @reprojectionError;%dist2plane
ransac.maxiter = 1e3;
ransac.threshold = 12;
ransac.inliers = [];
ransac.minimumset = 4;

x1h = tohomogeneous(feature1(1:2,matches(1,:)));
x2h = tohomogeneous(feature2(1:2,matches(2,:)));
result = ransac_routine_homo(x1h, x2h, ransac);
x1h = x1h(:,result.inliers);
x2h = x2h(:,result.inliers);

imshow1 = cat(2, im1, im2);
figure;imshow(imshow1);hold on;

plot(x1h(1,:),x1h(2,:), 'ro','MarkerSize',5);
plot(x2h(1,:)+size(img1,2),x2h(2,:), 'bo','MarkerSize',5);

shift = size(img1,2);
cmap = jet(32);
k = 1;
for i = 1:size(x1h,2)
    ptdraw = [x1h(2,i), x1h(1,i);
              x2h(2,i), x2h(1,i)+shift];
    plot(ptdraw(:,2),ptdraw(:,1),'LineStyle','-','LineWidth',0.5,'Color',cmap(k,:));
    k = mod(k+1,32);if k == 0 k = 1;end
end

im12 = WarpAndBlend(result.params,im1,im2);
figure;imshow(im12);

% tic
% ds2 = single(descriptor2);
% ds1 = single(descriptor3);
% kdtree = vl_kdtreebuild(ds2);
% matches = zeros(size(descriptor3,2),2);
% k = 1;
% thresh = 2;
% for i = 1:size(descriptor3,2)
%     [index, distance] = vl_kdtreequery(kdtree, ds2, ds1(:,i), ...
%                                        'NumNeighbors', 2, 'MaxComparisons', 15);
%     if distance(1)*thresh > distance(2) continue; end
%     matches(k,:) = [i,index(1)];
%     k = k + 1;
% end
% matches(k:end,:) = [];
% matches = matches';
% toc
% 
% ransac.pinlier = 0.99;
% ransac.estt_fun = @HestWithNormalization;%plane_estimation
% ransac.eval_fun = @reprojectionError;%dist2plane
% ransac.maxiter = 1e3;
% ransac.threshold = 12;
% ransac.inliers = [];
% ransac.minimumset = 4;
% 
% x1h = tohomogeneous(feature3(1:2,matches(1,:)));
% x2h = tohomogeneous(feature2(1:2,matches(2,:)));
% result = ransac_routine_homo(x1h, x2h, ransac);
% x1h = x1h(:,result.inliers);
% x2h = x2h(:,result.inliers);
% 
% im123 = WarpAndBlend(result.params,im3,im12);
% figure;imagesc(im123);

