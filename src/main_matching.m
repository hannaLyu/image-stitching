clc;close all;clear all;

%% the first example code to load image, display image, and do some basic image filtering

% addpath to the data folder
addpath ../data/

addpath(genpath('./features'));
addpath(genpath('./matching'));


% load an image
im1 = imread('pf_desk.jpg');
im2 = imread('pf_stand.jpg');
[s11,s12,s13] = size(im1);
[s21,s22,s23] = size(im2);

if s11 * s12 ~= s21*s22
    if s11*s12<s21*s22
        im2 = imresize(im2,[s11,s12]);
    else
        im1 = imresize(im1,[s21,s22]);
    end
end

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

% convert image from uint8 to double, make nothing but convert 255 - 1, 0-0
im1double = im2double(im1gray);
im2double = im2double(im2gray);

fast1 = fast_corner_detector(im1double, 500);
fast2 = fast_corner_detector(im2double, 500);
% fast1 = detectFASTFeatures(im1double);
% fast2 = detectFASTFeatures(im2double);
% fast1 = fast1.selectStrongest(500);
% fast2 = fast2.selectStrongest(500);
% fast1 = fast1.Location;
% fast2 = fast2.Location;
% fast1 = fliplr(fast1);
% fast2 = fliplr(fast2);

brief_pattern;
descriptors1 = extract_brief_descriptor(im1double,fast1,pattern);
descriptors2 = extract_brief_descriptor(im2double,fast2,pattern);

matchingpairs = brief_matching(descriptors1, descriptors2);

imshow1 = cat(2, im1, im2);
figure;imshow(imshow1);hold on;

plot(fast1(:,2),fast1(:,1), 'ro','MarkerSize',5);
plot(fast2(:,2)+size(im1double,2),fast2(:,1), 'bo','MarkerSize',5);


shift = size(im1double,2);
cmap = jet(32);
k = 1;
for i = 1:size(matchingpairs,1)
    if ~isinf(matchingpairs(i,2))
        ptdraw = [fast1(matchingpairs(i,1),1), fast1(matchingpairs(i,1),2);
                  fast2(matchingpairs(i,2),1), fast2(matchingpairs(i,2),2)+shift];
        plot(ptdraw(:,2),ptdraw(:,1),'LineStyle','-','LineWidth',0.5,'Color',cmap(k,:));
        k = mod(k+1,32);if k == 0 k = 1;end
    end
end



