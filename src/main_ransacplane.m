clc;close all;clear all;

% addpath(genpath('./'));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% either load or call gen_line_data()
% load ../data/LineData
% X = gen_line_data(500);


%for line
[X,Xi,Xo] = gen_plane_data(500);
 figure;
figure;plot3(X(1,:),X(2,:),X(3,:),'o');hold on; % show points
 
 [result]=tohomogeneousplane(X);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Results Visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 plot_plane_case(X,result,Xi,Xo);







