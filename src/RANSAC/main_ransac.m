clc;close all;clear all;
addpath(genpath('./'));

% addpath(genpath('./'));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% either load or call gen_line_data()
% load ../data/LineData
% X = gen_line_data(500);


%for line
% X = gen_line_data(500);
% figure;
% plot(X(1,:),X(2,:),'o');hold on; % show points
%for plane
[X,Xi,Xo] = gen_plane_data(500);
figure;
plot3(X(1,:),X(2,:),X(3,:),'o');hold on; % show points
 
 
[inliers,linepara]=tohomogeneous(X)
result.params = linepara;
result.inliers = inliers;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Results Visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %plot_line_case(X,result);
 plot_plane_case(X,result,Xi,Xo);







