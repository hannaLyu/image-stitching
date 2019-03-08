clc;close all;clear all;

addpath(genpath('./'));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% either load or call gen_line_data()
% load ../data/LineData
% X = gen_line_data(500);

% for plane
[X,Xi,Xo] = gen_plane_data(200);

ph = tohomogeneous(X);

ransac.pinlier = 0.99;
ransac.estt_fun = @plane_estimation;%plane_estimation
ransac.eval_fun = @dist2plane;%dist2plane
ransac.maxiter = 1e6;
ransac.threshold = 0.1;
ransac.inliers = [];
ransac.minimumset = 2;
result = ransac_routine(ph, ransac);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Results Visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_plane_case(X,result,Xi,Xo);







