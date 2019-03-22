clc;clear all;close all;

% imgpath = '../data/pier/';  

% Load images.
% buildingScene = imageDatastore(imgpath);

% Load images.
buildingDir = fullfile(toolboxdir('vision'), 'visiondata', 'building');
buildingScene = imageDatastore(buildingDir);

% Display images to be stitched
montage(buildingScene.Files);

numImages = numel(buildingScene.Files);

% Initialize variable to hold image sizes.
imageSize = zeros(numImages,2);

%% feature detection and description
seq.features = cell(numImages,1);
seq.descriptors = cell(numImages,1);
seq.imgs = cell(numImages,1);
for n = 1:numImages
    im1 = readimage(buildingScene, n);
    % convert img from rgb to gray scale
    if size(im1,3) == 3
        im1gray = rgb2gray(im1);
    else
        im1gray = im1;
    end
    img1 = single(im1gray);
    [feature1,descriptor1] = vl_sift(img1);
    seq.features{n} = feature1;
    seq.descriptors{n} = descriptor1;
    seq.imgs{n} = im1;
    imageSize(n,:) = size(im1gray);
end

%% pair matching
pairs = nchoosek(1:numImages,2);
mactheScore = zeros(numImages,numImages);
ss = cell(size(pairs,1),1);
tt = cell(size(pairs,1),1);
weights = zeros(size(pairs,1),1);
for i = 1:size(pairs,1)
    descriptor1 = seq.descriptors{pairs(i,1)};
    descriptor2 = seq.descriptors{pairs(i,2)};
    [matches, scores] = vl_ubcmatch(descriptor1, descriptor2,1.6);
    mactheScore(pairs(i,1),pairs(i,2)) = size(matches,2);
    mactheScore(pairs(i,2),pairs(i,1)) = size(matches,2);

    ss{i} = num2str(pairs(i,1));
    tt{i} = num2str(pairs(i,2));
    weights(i) = 1./(size(matches,2)+1e-3);
end

%% treat as a graphy, find the minimum spanning tree
G = graph(ss,tt,weights);
p = plot(G,'EdgeLabel',G.Edges.Weight, 'MarkerSize', 5);

nl = p.NodeLabel;
p.NodeLabel = '';
xd = get(p, 'XData');
yd = get(p, 'YData');
text(xd, yd, nl, 'FontSize',20, 'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','middle')
n2 = p.EdgeLabel;
% minimum spanning tree
[T, pred] = minspantree(G);
highlight(p,T,'EdgeColor','r','LineWidth',5);

% find the short path
minCost = 1e6;
minRoot = [];
minPath = {};
for i = 1:1:size(T.Nodes.Name,1)
    s = T.Nodes.Name{i};
    totalCost = 0;
    currRoot = s;
    currPath = {};
    k = 1;
    for j = 1:1:size(T.Nodes.Name,1)
        if i == j
            continue;
        else
            t = T.Nodes.Name{j};
            [P,d] = shortestpath(G,s,t);
            totalCost = totalCost + d;
            currPath{k} = P;
            k = k+1;
        end
    end
    
    if totalCost < minCost
        minCost = totalCost;
        minRoot = currRoot;
        minPath = currPath;
    end
end

%% now, use the minimum spanning tree to do image stiching
Hs = zeros(3,3,numImages);% homography to base
flagpair = zeros(numImages,numImages);
flagsingle = zeros(1,numImages);
minRoot = str2num(minRoot);
Hs(:,:,minRoot) = eye(3);% root is base
flagsingle(minRoot) = 1;
for i = 1:1:size(minPath,2)
    for j = 1:1:size(minPath{i},2)-1
        ii = str2num(minPath{i}{j});
        jj = str2num(minPath{i}{j+1});

        if flagpair(ii,jj) == 0
            % not yet estimated
            descriptor1 = seq.descriptors{jj};
            descriptor2 = seq.descriptors{ii};
            matches = vl_ubcmatch(descriptor1, descriptor2,1.6);
            H = estimateTransformation(seq.features{jj},seq.features{ii},matches);
            flagpair(ii,jj) = 1;
            flagpair(jj,ii) = 1;            
        end

        if flagsingle(jj) == 0 && flagsingle(ii) == 1
            Hs(:,:,jj) = Hs(:,:,ii)*H;
            flagsingle(jj) = 1;
        end
    end
end

for i = 1:numImages
    tforms(i) = maketform('projective',Hs(:,:,i)');
end

% Compute the output limits  for each transform
for i = 1:numel(tforms)
    [~,xlim(i,:),ylim(i,:)]=imtransform(zeros(imageSize(i,1),imageSize(i,2)),tforms(i));
end

maxImageSize = max(imageSize);

% Find the minimum and maximum output limits
xMin = min([1; xlim(:)]);
xMax = max([maxImageSize(2); xlim(:)]);

yMin = min([1; ylim(:)]);
yMax = max([maxImageSize(1); ylim(:)]);

% Create a 2-D spatial reference object defining the size of the panorama.
xLimits = [xMin xMax];
yLimits = [yMin yMax];

% Width and height of panorama.
mask0 = imtransform(zeros(imageSize(minRoot,1),imageSize(minRoot,2)), tforms(minRoot), 'XData',xLimits,'YData',yLimits);

width  = size(mask0,2);
height = size(mask0,1);

% Initialize the "empty" panorama.
panorama = zeros([height width 3]);
weights = zeros([height width]);
% blender = vision.AlphaBlender('Operation', 'Binary mask', ...
%     'MaskSource', 'Input port');
margin = 20;
se = strel('square',30);
% Create the panorama.
for i = 1:numImages

    I = readimage(buildingScene, i);

    % Transform I into the panorama.
    warpedImage = imtransform(I, tforms(i), 'XData',xLimits,'YData',yLimits,'FillValues',zeros(size(I,3),1));

%%%%%%%%%%%%%%%%%%%%% Alpha blend %%%%%%%%%%%%%%%%%%%%%%%%%
%     % Generate a binary mask.
%     mask = imtransform(ones(size(I,1),size(I,2)), tforms(i), 'XData',xLimits,'YData',yLimits,'FillValues',0)>0;
%     % Overlay the warpedImage onto the panorama.
%     panorama = step(blender, panorama, warpedImage, mask);
%%%%%%%%%%%%%%%%%%%%% sinple blend using distance %%%%%%%%%%%%%%%%%%%%%%%%%
    mask = imtransform(ones(size(I,1),size(I,2)), tforms(i), 'XData',xLimits,'YData',yLimits,'FillValues',0)>0;
    dist1 = zeros(size(I,1),size(I,2));%dist1(round(size(I,1)*0.5),round(size(I,2)*0.5)) = 1;
    dist1(1:margin,:) = 1;dist1(end-margin:end,:) = 1;dist1(:,1:margin) = 1;dist1(:,end-margin:end) = 1;
    dist1 = bwdist(dist1,'chessboard');%maxdist = max(dist1(:));dist1 = (maxdist+1) - dist1;
    dist1 = imdilate(dist1,se);
    dist1t=imtransform(dist1,tforms(i),'XData',xLimits,'YData',yLimits,'FillValues',0) + 1e-3;
    [panorama,weights] = simpleBlend(warpedImage, panorama, weights, mask, dist1t);

%%%%%%%%%%%%%%%%%%%%% laplacian blend using distance %%%%%%%%%%%%%%%%%%%%%%%%%
%     mask = ones(size(I,1),size(I,2))*1;
%     mask = imtransform(mask, tforms(i), 'XData',xLimits,'YData',yLimits,'FillValues',0);
%     mask((rgb2gray(panorama) ~= 0)&mask) = 1;
%     panorama = laplacianBlend(warpedImage, panorama, mask);
end

figure
imshow(panorama)

function panorama = simpleBlend(im1, panorama, weights, mask, dist)
    

    if size(im1,3) == 1
        panorama = combine(im1, panorama, mask, dist);
    else
        panorama(:,:,1) = combine(im1(:,:,1), panorama(:,:,1), mask, dist);
        panorama(:,:,2) = combine(im1(:,:,2), panorama(:,:,2), mask, dist);
        panorama(:,:,3) = combine(im1(:,:,3), panorama(:,:,3), mask, dist);
    end
end

function im12 = combine(im1, im2, mask, dist)
    im12 = im2;
    mask12 = im12==0;
    mask1 = mask12 & mask;
    im12(mask1) = im1(mask1);
    mask12 = im12~=0;
    mask1 = mask12 & mask;
    weight = dist(mask1)./(max(dist(mask1)));
    im12(mask1) = uint8(weight.*double(im1(mask1))+(1-weight).*double(im12(mask1)));
end

function H = estimateTransformation(feature1,feature2,matches)
    ransac.pinlier = 0.99;
    ransac.estt_fun = @HestWithNormalization;%plane_estimation
    ransac.eval_fun = @reprojectionError;%dist2plane
    ransac.maxiter = 1e3;
    ransac.threshold = 6;
    ransac.inliers = [];
    ransac.minimumset = 4;

    x1h = tohomogeneous(feature1(1:2,matches(1,:)));
    x2h = tohomogeneous(feature2(1:2,matches(2,:)));
    result = ransac_routine_homo(x1h, x2h, ransac);
    H = result.params;
end