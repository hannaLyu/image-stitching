clc;
close all;
clear all;
% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end

imgpath = '../../data/poission_blending/';  

% Load images.
buildingScene = imageDatastore(imgpath);

im1 = readimage(buildingScene, 9);
im2 = readimage(buildingScene, 19);

% im1=imresize(im1,0.8);
% im2=imresize(im2,0.8);


% mask = uint8(ones(size(im1,1),size(im1,2)));
% mask(:,1:round(size(im1,2)*0.5)) = 0;
% imblend = laplacianBlend(im1, im2, mask);
% imblend = TwoBandBlend(imblend);

figure;imshow(im1);
% figure;imshow(imblend);
imshow(im1);hold on;
[x,y] = ginput();x = x';y = y';

% load('test.mat');
x(x<2) = 2;
x(x>size(im1,2)-1) = size(im1,2)-1;
y(y<2) = 2;
y(y>size(im1,1)-1) = size(im1,1)-1;
% 
x = [x x(1)];
y = [y y(1)];
plot(x,y,'r-');

figure;imshow(im2);
[anchorx,anchory] = ginput(1);
hold on;plot(x-x(1)+anchorx,y-y(1)+anchory,'r-');


% 
bwsrc = poly2mask(x,y,size(im1,1),size(im1,2));
figure
imshow(bwsrc);hold on;
boundarysrc = findBoundary(bwsrc);
% [row, col] = find(boundary~=0);plot(col,row,'ro');
% [row, col] = find(xor(bw,boundary~=0));plot(col,row,'b.');
boundarysrc = boundarysrc ~= 0;
% region = bw;%xor(bw,boundary);

bwdst = zeros(size(im2,1),size(im2,2));
[row,col] = find(bwsrc);
bwdst((col+round(anchorx-x(1))-1).*size(im2,1)+(row+round(anchory-y(1)))) = 1;
bwdst = bwdst == 1;

figure
imshow(bwdst);hold on;
boundarydst = findBoundary(bwdst);
boundarydst = boundarydst ~= 0;

A = precomputeA(bwdst, boundarydst);
im1 = im2double(im1);
im2 = im2double(im2);
source = (rgb2gray(im1));
target = (rgb2gray(im2));
for i = 1:size(im1,3)   
    source = (im1(:,:,i));
    target = (im2(:,:,i));
    g = precomputeFlow(source,target,bwsrc,boundarysrc, bwdst, boundarydst, A);
    imblend(:,:,i) = poissonBlend(source,target,bwsrc,boundarysrc, bwdst, boundarydst, A, g);
end
imshow(imblend)

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

function A = precomputeA(bwdst, boundarydst)
    int = xor(bwdst, boundarydst);
    N = sum(int(:));

    m = size(bwdst,1);
    n = size(bwdst,2);
    [row, col] = find(int);
    map1((col-1).*m+row) = 1:N;
    % establish
    
    nnu = [row-1,col];n1 = nnu(:,1)>0;  nuint=n1;nuint(n1) = int((nnu(n1,2)-1).*m+nnu(n1,1)) == 1;
    nnd = [row+1,col];n2 = nnd(:,1)<=m; ndint=n2;ndint(n2) = int((nnd(n2,2)-1).*m+nnd(n2,1)) == 1;
    nnl = [row,col-1];n3 = nnl(:,2)>0;  nlint=n3;nlint(n3) = int((nnl(n3,2)-1).*m+nnl(n3,1)) == 1;
    nnr = [row,col+1];n4 = nnr(:,2)<=n; nrint=n4;nrint(n4) = int((nnr(n4,2)-1).*m+nnr(n4,1)) == 1;
    
    Np = double(n1) + double(n2) + double(n3) + double(n4);

%     A = sparse(N,N);
%     for i = 1:N
%         A(i,i) = Np(i);
%         if nuint(i) == 1
%             A(i,(map1((nnu(i,2)-1).*m+nnu(i,1)))) = -1;
%         end
%         if ndint(i) == 1
%             A(i,(map1((nnd(i,2)-1).*m+nnd(i,1)))) = -1;
%         end
%         if nlint(i) == 1
%             A(i,(map1((nnl(i,2)-1).*m+nnl(i,1)))) = -1;
%         end
%         if nrint(i) == 1
%             A(i,(map1((nnr(i,2)-1).*m+nnr(i,1)))) = -1;
%         end
%     end
%     A = sparse(A);

    A1 = sparse(N,N);
    id = 1:N;id1 = id - 1;
    A1(id1.*N+id) = Np;
    A2 = sparse(id(nuint), map1((nnu(nuint,2)-1).*m+nnu(nuint,1)), -1*ones(sum(nuint),1), N,N);
    A3 = sparse(id(ndint), map1((nnd(ndint,2)-1).*m+nnd(ndint,1)), -1*ones(sum(ndint),1), N,N);
    A4 = sparse(id(nlint), map1((nnl(nlint,2)-1).*m+nnl(nlint,1)), -1*ones(sum(nlint),1), N,N);
    A5 = sparse(id(nrint), map1((nnr(nrint,2)-1).*m+nnr(nrint,1)), -1*ones(sum(nrint),1), N,N);
%     A1(id1(nuint).*N+map1((nnu(nuint,2)-1).*m+nnu(nuint,1))) = -1;
%     A1(id1(ndint).*N+map1((nnd(ndint,2)-1).*m+nnd(ndint,1))) = -1;
%     A1(id1(nlint).*N+map1((nnl(nlint,2)-1).*m+nnl(nlint,1))) = -1;
%     A1(id1(nrint).*N+map1((nnr(nrint,2)-1).*m+nnr(nrint,1))) = -1;
    A = A1 + A2 + A3 + A4 + A5;

%     A1 = spdiags([Np -1*ones(N,1) -1*ones(N,1) -1*ones(N,1) -1*ones(N,1)],[0,-1,1,-m,m],N,N);
%     id = [find(~nuint);find(~ndint);find(~nlint);find(~nrint)];
%     for i = 1:length(id)
%         colid = find(A1(id(i),:));
%         colid(colid == id(i)) = [];
%         valid = nuint(colid) | ndint(colid) | nlint(colid) | nrint(colid);
%         A1(id(i),colid(~valid)) = 0;
%     end
    
%     id = find(~ndint);
%     for i = 1:length(id)
%         if id(i)+1 <= N
%             A(id(i),id(i)+1) = 0;
%         end
%     end
% 
%     id = find(~nlint);
%     for i = 1:length(id)
%         if id(i)-m > 0
%             A(id(i),id(i)-m) = 0;
%         end
%     end
%     
%     id = find(~nrint);
%     for i = 1:length(id)
%         if id(i)+m <= N
%             A(id(i),id(i)+m) = 0;
%         end
%     end
end

function g = precomputeFlow(source,target,bwsrc,boundarysrc, bwdst, boundarydst, A)
    int = xor(bwdst, boundarydst);
    
    m = size(target,1);
    n = size(target,2);
    [row, col] = find(int);
    
    nnu = [row-1,col];n1 = nnu(:,1)>0;  
    nnd = [row+1,col];n2 = nnd(:,1)<=m; 
    nnl = [row,col-1];n3 = nnl(:,2)>0;  
    nnr = [row,col+1];n4 = nnr(:,2)<=n; 

    % 
    g1 = target(int);
    g1u = zeros(length(g1),1);
    g1u(n1) = g1(n1) - target((nnu(n1,2)-1).*m+nnu(n1,1));
    
    g1d = zeros(length(g1),1);
    g1d(n2) = g1(n2) - target((nnd(n2,2)-1).*m+nnd(n2,1));
    
    g1l = zeros(length(g1),1);
    g1l(n3) = g1(n3) - target((nnl(n3,2)-1).*m+nnl(n3,1));
    
    g1r = zeros(length(g1),1);
    g1r(n4) = g1(n4) - target((nnr(n4,2)-1).*m+nnr(n4,1));
    
%     % boundary
%     nubint=n1;nubint(n1) = boundarydst((nnu(n1,2)-1).*m+nnu(n1,1)) == 1;
%     ndbint=n2;ndbint(n2) = boundarydst((nnd(n2,2)-1).*m+nnd(n2,1)) == 1;
%     nlbint=n3;nlbint(n3) = boundarydst((nnl(n3,2)-1).*m+nnl(n3,1)) == 1;
%     nrbint=n4;nrbint(n4) = boundarydst((nnr(n4,2)-1).*m+nnr(n4,1)) == 1;
    % for target
%     g1 = target(int);
% %     vpq1 = A*g1;
%     vpq1(nubint) = vpq1(nubint) - target((nnu(nubint,2)-1).*m+nnu(nubint,1));
%     vpq1(ndbint) = vpq1(ndbint) - target((nnd(ndbint,2)-1).*m+nnd(ndbint,1));
%     vpq1(nlbint) = vpq1(nlbint) - target((nnl(nlbint,2)-1).*m+nnl(nlbint,1)); 
%     vpq1(nrbint) = vpq1(nrbint) - target((nnr(nrbint,2)-1).*m+nnr(nrbint,1));

    % vpq
    int2 = xor(bwsrc, boundarysrc);
    g2 = source(int2);
    vpq2 = A*g2;
    
    [row2, col2] = find(int2);

    nnu2 = [row2-1,col2];n1 = nnu2(:,1)>0;  
    nnd2 = [row2+1,col2];n2 = nnd2(:,1)<=size(bwsrc,1); 
    nnl2 = [row2,col2-1];n3 = nnl2(:,2)>0;  
    nnr2 = [row2,col2+1];n4 = nnr2(:,2)<=size(bwsrc,2); 

    m2 = size(bwsrc,1);
    
    g2u = zeros(length(g2),1);
    g2u(n1) = g2(n1) - source((nnu2(n1,2)-1).*m2+nnu2(n1,1));
    
    g2d = zeros(length(g2),1);
    g2d(n2) = g2(n2) - source((nnd2(n2,2)-1).*m2+nnd2(n2,1));
    
    g2l = zeros(length(g2),1);
    g2l(n3) = g2(n3) - source((nnl2(n3,2)-1).*m2+nnl2(n3,1));
    
    g2r = zeros(length(g2),1);
    g2r(n4) = g2(n4) - source((nnr2(n4,2)-1).*m2+nnr2(n4,1));
    
% 
%     nubint2=n1;nubint2(n1) = boundarysrc((nnu2(n1,2)-1).*m2+nnu2(n1,1)) == 1;
%     ndbint2=n2;ndbint2(n2) = boundarysrc((nnd2(n2,2)-1).*m2+nnd2(n2,1)) == 1;
%     nlbint2=n3;nlbint2(n3) = boundarysrc((nnl2(n3,2)-1).*m2+nnl2(n3,1)) == 1;
%     nrbint2=n4;nrbint2(n4) = boundarysrc((nnr2(n4,2)-1).*m2+nnr2(n4,1)) == 1;
% 
%     vpq2(nubint2) = vpq2(nubint2) - source((nnu2(nubint2,2)-1).*m2+nnu2(nubint2,1));
%     vpq2(ndbint2) = vpq2(ndbint2) - source((nnd2(ndbint2,2)-1).*m2+nnd2(ndbint2,1));
%     vpq2(nlbint2) = vpq2(nlbint2) - source((nnl2(nlbint2,2)-1).*m2+nnl2(nlbint2,1)); 
%     vpq2(nrbint2) = vpq2(nrbint2) - source((nnr2(nrbint2,2)-1).*m2+nnr2(nrbint2,1));

%     id = abs(vpq1) > abs(vpq2);
%     g = vpq2;%vpq2*0.5+vpq1*0.5;
%     g(id) = vpq1(id);
    %
    
    %% mixed
    id  = abs(g1u) > abs(g2u);gmu = g2u;gmu(id) = g1u(id);
    id  = abs(g1d) > abs(g2d);gmd = g2d;gmd(id) = g1d(id);
    id  = abs(g1l) > abs(g2l);gml = g2l;gml(id) = g1l(id);
    id  = abs(g1r) > abs(g2r);gmr = g2r;gmr(id) = g1r(id);

    g = gmu + gmd + gml + gmr;
end

function imblend = poissonBlend(source,target,bwsrc,boundarysrc, bwdst, boundarydst, A, g)
    int = xor(bwdst, boundarydst);
    N = sum(int(:));
    
    m = size(target,1);
    n = size(target,2);
    [row, col] = find(int);
    
    nnu = [row-1,col];n1 = nnu(:,1)>0;  
    nnd = [row+1,col];n2 = nnd(:,1)<=m; 
    nnl = [row,col-1];n3 = nnl(:,2)>0;  
    nnr = [row,col+1];n4 = nnr(:,2)<=n; 

    % boundary
    nubint=n1;nubint(n1) = boundarydst((nnu(n1,2)-1).*m+nnu(n1,1)) == 1;
    ndbint=n2;ndbint(n2) = boundarydst((nnd(n2,2)-1).*m+nnd(n2,1)) == 1;
    nlbint=n3;nlbint(n3) = boundarydst((nnl(n3,2)-1).*m+nnl(n3,1)) == 1;
    nrbint=n4;nrbint(n4) = boundarydst((nnr(n4,2)-1).*m+nnr(n4,1)) == 1;
    
    b1 = zeros(N,1);
    b1(nubint) = b1(nubint) + target((nnu(nubint,2)-1).*m+nnu(nubint,1));
    b1(ndbint) = b1(ndbint) + target((nnd(ndbint,2)-1).*m+nnd(ndbint,1));
    b1(nlbint) = b1(nlbint) + target((nnl(nlbint,2)-1).*m+nnl(nlbint,1));
    b1(nrbint) = b1(nrbint) + target((nnr(nrbint,2)-1).*m+nnr(nrbint,1));

    % for target
%     g1 = target(int);
%     vpq1 = A*g1;
%    
%     vpq1(nubint) = vpq1(nubint) - target((nnu(nubint,2)-1).*m+nnu(nubint,1));
%     vpq1(ndbint) = vpq1(ndbint) - target((nnd(ndbint,2)-1).*m+nnd(ndbint,1));
%     vpq1(nlbint) = vpq1(nlbint) - target((nnl(nlbint,2)-1).*m+nnl(nlbint,1)); 
%     vpq1(nrbint) = vpq1(nrbint) - target((nnr(nrbint,2)-1).*m+nnr(nrbint,1));
% 
%     % vpq
%     int2 = xor(bwsrc, boundarysrc);
%     g = source(int2);
%     vpq2 = A*g;
%     
%     [row2, col2] = find(int2);
% 
%     nnu2 = [row2-1,col2];n1 = nnu2(:,1)>0;  
%     nnd2 = [row2+1,col2];n2 = nnd2(:,1)<=size(bwsrc,1); 
%     nnl2 = [row2,col2-1];n3 = nnl2(:,2)>0;  
%     nnr2 = [row2,col2+1];n4 = nnr2(:,2)<=size(bwsrc,2); 
% 
%     m2 = size(bwsrc,1);
% 
%     nubint2=n1;nubint2(n1) = boundarysrc((nnu2(n1,2)-1).*m2+nnu2(n1,1)) == 1;
%     ndbint2=n2;ndbint2(n2) = boundarysrc((nnd2(n2,2)-1).*m2+nnd2(n2,1)) == 1;
%     nlbint2=n3;nlbint2(n3) = boundarysrc((nnl2(n3,2)-1).*m2+nnl2(n3,1)) == 1;
%     nrbint2=n4;nrbint2(n4) = boundarysrc((nnr2(n4,2)-1).*m2+nnr2(n4,1)) == 1;
% 
%     vpq2(nubint2) = vpq2(nubint2) - source((nnu2(nubint2,2)-1).*m2+nnu2(nubint2,1));
%     vpq2(ndbint2) = vpq2(ndbint2) - source((nnd2(ndbint2,2)-1).*m2+nnd2(ndbint2,1));
%     vpq2(nlbint2) = vpq2(nlbint2) - source((nnl2(nlbint2,2)-1).*m2+nnl2(nlbint2,1)); 
%     vpq2(nrbint2) = vpq2(nrbint2) - source((nnr2(nrbint2,2)-1).*m2+nnr2(nrbint2,1));
% 
%     id = abs(vpq1) > abs(vpq2);
%     vpq = vpq2;
% %     vpq(id) = vpq1(id);
    %
    b = b1 + g;

    x = A\b;

    imblend(~int) = target(~int);
    imblend(int) = x;
    imblend = reshape(imblend,m,n);

%     imblend = (imblend - min(imblend(:))) ./ (max(imblend(:))-min(imblend(:)));
end