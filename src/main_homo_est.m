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



