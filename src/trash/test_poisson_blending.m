clc;
close all;
clear all;
% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end

imgpath = '../../data/pier/';  

% Load images.
buildingScene = imageDatastore(imgpath);

im1 = readimage(buildingScene, 1);
im2 = readimage(buildingScene, 2);

% mask = uint8(ones(size(im1,1),size(im1,2)));
% mask(:,1:round(size(im1,2)*0.5)) = 0;
% imblend = laplacianBlend(im1, im2, mask);
% imblend = TwoBandBlend(imblend);

figure;imshow(im1);
% figure;imshow(imblend);
imshow(im1);hold on;
% [x,y] = ginput();x = x';y = y';

load('test.mat');

% 
x = [x x(1)];
y = [y y(1)];
plot(x,y,'r-');

figure;imshow(im2);
% [anchorx,anchory] = ginput(1);
hold on;plot(x-x(1)+anchorx,y-y(1)+anchory,'r-');


% 
bw = poly2mask(x,y,size(im1,1),size(im1,2));
figure
imshow(bw);hold on;
boundary = findBoundary(bw);
% [row, col] = find(boundary~=0);plot(col,row,'ro');
% [row, col] = find(xor(bw,boundary~=0));plot(col,row,'b.');
boundary = boundary ~= 0;
region = bw;%xor(bw,boundary);

for i = 1:size(im1,3)   
    source = im2double(im1(:,:,i));
    target = im2double(im2(:,:,i));
    imblend = poissonBlend(source,target,region,boundary);
end

% poisson blending
function boundary = findBoundary(bw)
    boundary = zeros(size(bw));
    for i = 2:size(bw,1)-1
        for j = 2:size(bw,2)-1
            block = bw(i-1:i+1,j-1:j+1);
            if block(5) == 1 && sum(block([1,2,3,4,6,7,8,9])==0)~=0
                boundary(i,j) = 1;
            end
        end
    end
end

function imblend = poissonBlend(source,target,mask,boundary)
    m = size(source,1);n = size(source,2);
    N = m*n;
    A = spdiags([-4.*ones(N,1) ones(N,1) ones(N,1) ones(N,1) ones(N,1)],[0,1,-1,m,-m],N,N);
    ii = m:m:N-m;
    jj = ii+1;
    A((ii-1)*N+jj) = 0;
    ii = m:m:N-m;
    A((ii.*N+ii)) = 0;
    % gradient of source
    ds = A * source(:);
    % gradient of target
    dt = A * target(:);
    
    %
%     b1 = dt(~mask);
%     internal = xor(mask,boundary);
%     b1 = dt(internal);
%     b2 = dt(boundary);
% 
%     [row,col] = find(boundary);
%     nnl = [row, col-1];nn1ind = sub2ind(size(source),nnl(:,1),nnl(:,2));
%     nnr = [row, col+1];nnrind = sub2ind(size(source),nnr(:,1),nnr(:,2));
%     nnd = [row+1,col];nndind = sub2ind(size(source),nnd(:,1),nnd(:,2));
%     nnu = [row-1,col];nnuind = sub2ind(size(source),nnu(:,1),nnu(:,2));
%     
%     neightbors = source(nn1ind) + source(nnrind) + + source(nndind) + + source(nnuind);
% 
%     b2 = b2 + neightbors;
%     b = [b1;b2];
%     
%     ind1 = 1:N;
%     ind1 = ind1(mask(:));
%     Ab = A(ind1,ind1);

    b = dt;
    [row,col] = find(~boundary);
    nnl = [row, col-1];
    nn1ind = sub2ind(size(source),nnl(:,1),nnl(:,2));
    nnr = [row, col+1];
    nnrind = sub2ind(size(source),nnr(:,1),nnr(:,2));
    nnd = [row+1,col];
    nndind = sub2ind(size(source),nnd(:,1),nnd(:,2));
    nnu = [row-1,col];
    nnuind = sub2ind(size(source),nnu(:,1),nnu(:,2));
end