clear all;
close all;
clc;
addpath('mybrieffast/');
addpath('precompute/');

I1=imread('chickenbroth_01.jpg');
I2=imread('chickenbroth_02.jpg');
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

[corner1,res1]=fast(im1gray);
[corner2,res2]=fast(im2gray);

img1 = im2double(im1gray);
img2 = im2double(im2gray);

% 
% corner1 = fast_corner_detector(im1gray, 500);
% corner2 = fast_corner_detector(im2gray, 500);

brief_pattern;
% descriptor1= extract_brief_descriptor(img1,corner1,pattern);
% descriptor2= extract_brief_descriptor(img2,corner2,pattern);
descriptor1 = extractdes(img1,corner1,pattern);
descriptor2 = extractdes(img2,corner2,pattern);

% matching_pairs = bruteforce(descriptor1,descriptor2);
matchingpairs = brief_matching(descriptor1, descriptor2);

imshow1 = cat(2, img1, img2);
figure;imshow(imshow1);hold on;

plot(corner1(:,2),corner1(:,1), 'ro','MarkerSize',5);
plot(corner2(:,2)+size(I1,2),corner2(:,1), 'bo','MarkerSize',5);
shift = size(I1,2);
cmap = jet(32);
k = 1;
for i = 1:size(matchingpairs,1)
    if ~isinf(matchingpairs(i,2))
        ptdraw = [corner1(matchingpairs(i,1),1), corner1(matchingpairs(i,1),2);
                  corner2(matchingpairs(i,2),1), corner2(matchingpairs(i,2),2)+shift];
        plot(ptdraw(:,2),ptdraw(:,1),'LineStyle','-','LineWidth',0.5,'Color',cmap(k,:));
        k = mod(k+1,32);
        if k == 0 k = 1;
    end
    end
end