clc;close all;clear all;

%% the first example code to load image, display image, and do some basic image filtering

% addpath to the data folder
addpath ../data/

% load an image
im1 = imread('right.jpg');

% display an image
figure;imshow(im1)

% convert img from rgb to gray scale
imgray = rgb2gray(im1);
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
g = fspecial('gaussian', 50, 50);% try play with this two values to see what will hapen
imblur = imfilter(imgray, g, 'replicate');
figure;imshow(imblur)
