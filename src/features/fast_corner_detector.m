function varargout = fast_corner_detector(I, num)
%     clc;close all;clear all;

    %% the first example code to load image, display image, and do some basic image filtering
    % addpath to the data folder
%     addpath ../data/
% 
%     % load an image
%     im1 = imread('library2.jpg');
% 
%     % display an image
%     figure;imshow(im1)
% 
%     % convert img from rgb to gray scale
%     if size(im1,3) == 3
%         imgray = rgb2gray(im1);
%     else
%         imgray = im1;
%     end
%     figure;imshow(imgray);
% 
%     % convert image from uint8 to double, make nothing but convert 255 - 1, 0-0
%     imdouble = im2double(imgray);

% clc;  close all; clear all;
% %% loading the vidoe to extract the images
% video = VideoReader('rhinos.avi'); % here I put the name of the video
% vidWidth = video.Width;
% vidHeight = video.Height;
% frame = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),'colormap',[]);
% i =1;
% while hasFrame(video)
%     frame(i).cdata = readFrame(video);
%     i = i+1;
% end
% img11 = frame(10).cdata ; % choose the 1st frame to compare
% if length(size(img11)) == 3
% 	img1 =  uint8(img11(:,:,2));
% else
% 	img1 = uint8(img11);
%    
% end
% img22 = frame(12).cdata; % choose the second frame to compare
% if length(size(img22)) == 3
% 	img2 =  uint8(img22(:,:,2));
% else
% 	img2 = uint8(img22);
% end
%     imgray = rgb2gray(img11);
%     imdouble = im2double(imgray);figure;imshow(img11);
%     I = imdouble;
    
    [m,n] = size(I);
    tbl = [[-3 0]; [-3 1]; [-2 2]; [-1 3]; [0 3]; [1 3]; [2 2]; [3 1]; ...
           [3 0]; [3 -1]; [2 -2]; [1 -3]; [0 -3]; [-1 -3]; [-2 -2]; [-3 -1]];% l x 1
    l = size(tbl,1);
    
    threshold = 0.2;
    hsize = 3;
    st = max(hsize,3);
    
    [X,Y] = meshgrid(st+1:m-st,st+1:n-st);
    X = X(:)';% 1 x mn
    Y = Y(:)';% 1 x mn
    mn = length(X);
    X1 = repmat(X,l,1);% l x mn
    Y1 = repmat(Y,l,1);% l x mn
    
    tblrepx = repmat(tbl(:,1), 1, mn);
    tblrepy = repmat(tbl(:,2), 1, mn);
    
    X1 = X1 + tblrepx;
    Y1 = Y1 + tblrepy;
    
    X1 = X1(:);Y1 = Y1(:);
    
    ind0 = sub2ind([m,n],X1,Y1);
    ind1 = sub2ind([m,n],X,Y);
    tic
    score = zeros(m,n);
    for i = 1:1:mn
        il = i*l;
        val = I(ind0(il-l+1:il));
        diffval = (val - I(ind1(i)));
        s = length(find(abs(diffval) > threshold));
        if s > 9
            score(ind1(i)) = sum(abs(diffval));
        end
    end
    toc
    % nonmaxima suppression
    [x,y]=meshgrid(-hsize:hsize,-hsize:hsize);
    x = x(:);y=y(:);
    id = x == 0 & y == 0;
    x(id) = []; y(id) = [];
    ll = length(x);
    
    [row, col] = find(score ~= 0);
    mn2 = length(row);
    
    X2 = repmat(row',ll,1);% l x mn
    Y2 = repmat(col',ll,1);% l x mn
    
    repxx = repmat(x, 1, mn2);
    repyy = repmat(y, 1, mn2);
    
    X2 = X2 + repxx;
    Y2 = Y2 + repyy;
    
    X2 = X2(:);Y2 = Y2(:);
    
    ind2 = sub2ind([m,n],X2,Y2);
    ind4 = sub2ind([m,n],row,col);
    for i = 1:1:mn2
        if score(ind4(i)) == 0 continue; end
        ill = i*ll;
        val = score(ind2(ill-ll+1:ill));
        maxscore = max((val));
        if score(ind4(i)) >= maxscore
            score(ind2(ill-ll+1:ill)) = 0;
        else
            score(ind4(i)) = 0;
        end
    end
    [row, col] = find(score ~= 0);
    ind5 = sub2ind([m,n],row,col);
    [~,maxid] = sort(score(ind5),'descend');
    if num > length(row)
        num = length(row);
    end
    selid = maxid(1:num);
    corners = [row(selid) col(selid)];
    varargout{1} = corners;

%     corners = detectFASTFeatures(imdouble);
%     plot(corners.selectStrongest(500));
%     hold on;
%     plot(,row(maxid(:)), 'ro','MarkerSize',10),

end
