clear all;
close all;
clc;
VLFEATROOT='/Users/xiaohu/Documents/MATLAB/image-stitching/src/3rdparty/vlfeat-0.9.21';
run(strcat(VLFEATROOT,'/toolbox/vl_setup'));

I1=imread('pf_desk.jpg');
I2=imread('pf_stand.jpg');
[s11,s12,s13] = size(I1);
[s21,s22,s23] = size(I2);

if s11 * s12 ~= s21*s22
    if s11*s12<s21*s22
        I2 = imresize(I2,[s11,s12]);
    else
        I1 = imresize(I1,[s21,s22]);
    end
end

if size(I1,3) == 3
    im1gray = rgb2gray(I1);
else
    im1gray = I1;
end

if size(I2,3) == 3
    im2gray = rgb2gray(I2);
else
    im2gray = I2;
end

img1 = single(im1gray);
img2 = single(im2gray);

[feature1,descriptor1] = vl_sift(img1) ;
[feature2,descriptor2] = vl_sift(img2) ;
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

imshow1 = cat(2, I1, I2);
figure;imshow(imshow1);hold on;
% h1 = vl_plotframe(feature1(:,:)) ;
plot(feature1(1,:),feature1(2,:), 'ro','MarkerSize',5);
plot(feature2(1,:)+size(I1,2),feature2(2,:), 'bo','MarkerSize',5);
shift = size(I1,2);
cmap = jet(32);
k = 1;
for i = 1:size(matches,2)
    ptdraw = [feature1(2,matches(1,i)), feature1(1,matches(1,i));
              feature2(2,matches(2,i)), feature2(1,matches(2,i))+shift];
    plot(ptdraw(:,2),ptdraw(:,1),'LineStyle','-','LineWidth',0.5,'Color',cmap(k,:));
    k = mod(k+1,32);
    if k == 0 k = 1; end
end