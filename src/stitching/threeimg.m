clc;clear;
figure;hold on;
imgpath = 'M:\Documents\stitching\stitching\hill\hill';  
buildingScene = imageDatastore(imgpath);
montage(buildingScene.Files);
hold on
numImages = numel(buildingScene.Files);


patten = brief_pattern_generator();
% feature detect
seq.features = cell(numImages,1);
seq.descriptors = cell(numImages,1);
seq.imgs = cell(numImages,1);
for n=1:numImages
im1 = readimage(buildingScene, n);
  if size(im1,3) == 3
        im1gray = rgb2gray(im1);
    else
        im1gray = im1;
    end
 img1 = double(im1gray);
corner = fast(img1);
descriptor = extractdes(img1,corner,patten);
    seq.corner{n} = transpose(corner);
    seq.descriptors{n} = descriptor;
    seq.imgs{n} = im1;
    imageSize(n,:) = size(im1gray);
end
%check answer
subplot(2,2,1);
imshow(seq.imgs{1,1});
hold on

plot(seq.corner{1,1}(2,:),seq.corner{1,1}(1,:),'o');
subplot(2,2,2);
imshow(seq.imgs{2,1});
hold on
plot(seq.corner{1,2}(2,:),seq.corner{1,2}(1,:),'o');
subplot(2,2,3);
imshow(seq.imgs{3,1});
hold on
plot(seq.corner{1,3}(2,:),seq.corner{1,3}(1,:),'o');

%match
seq.matchingpair=cell(numImages,1);
pairs = nchoosek(1:numImages,2);
for i = 1:size(pairs,1)
    descriptor1 = seq.descriptors{pairs(i,1)};
    descriptor2 = seq.descriptors{pairs(i,2)};
    matchingpair = bruteforce(descriptor1, descriptor2);
    [I,J]=find(isinf(matchingpair));
    matchingpair(I,:)=[];
    matchingpair=transpose(matchingpair);
    seq.matchingpair{i}=matchingpair;
end

%ransac
seq.matchpoints=cell(numImages,1);
seq.inlier1=cell(numImages,1);
seq.inlier2=cell(numImages,1);
seq.inlierId=cell(numImages,1);
seq.x1=cell(numImages,1);
seq.x2=cell(numImages,1);
for i=1:size(pairs,1)
    corner1=seq.corner{pairs(i,1)};
    corner2=seq.corner{pairs(i,2)};
    seq.matchpoint1{i}=corner1(:,seq.matchingpair{i}(1,:));
    seq.matchpoint2{i}=corner2(:,seq.matchingpair{i}(2,:));
end

for i=1:size(pairs,1)
    
    x1=tohomogeneous(seq.matchpoint1{i});
    x2=tohomogeneous(seq.matchpoint2{i});
    [inlier1,inlier2,inlierId]= ransac(x1,x2);
   seq.x1{i}=x1;
   seq.x2{i}=x2;
    seq.inlier1{i}=inlier1;
    seq.inlier2{i}=inlier2;
    seq.inlierId{i}=inlierId;
end

% check answer
% subplot(2,2,1);
imshow1 = cat(2, seq.imgs{1,1}, seq.imgs{2,1});
figure;imshow(imshow1);hold on;
plot(seq.inlier1{1,1}(2,:),seq.inlier1{1,1}(1,:), 'ro','MarkerSize',3);
plot(seq.inlier2{1,1}(2,:)+size(seq.imgs{1,1},2),seq.inlier2{1,1}(1,:), 'bo','MarkerSize',3);

shift = size(seq.imgs{1},2);
cmap = jet(32);
k = 1;
for i = 1:size(seq.inlierId{1,1},2)
    
    if ~isinf(seq.inlierId{1}(1,i))
        ptdraw = [seq.x1{1}(1,seq.inlierId{1}(1,i)), seq.x1{1}(2,seq.inlierId{1}(1,i));
                  seq.x2{1}(1,seq.inlierId{1}(1,i)), seq.x2{1}(2,seq.inlierId{1}(1,i))+shift];
        plot(ptdraw(:,2),ptdraw(:,1),'LineStyle','-','LineWidth',0.5,'Color',cmap(k,:));
        k = mod(k+1,32);
        if k == 0 k = 1;
        end
    end
    end

% subplot(2,2,2);
imshow1 = cat(2, seq.imgs{1,1}, seq.imgs{3,1});
figure;imshow(imshow1);hold on;
plot(seq.inlier1{2,1}(2,:),seq.inlier1{2,1}(1,:), 'ro','MarkerSize',3);
plot(seq.inlier2{2,1}(2,:)+size(seq.imgs{2,1},2),seq.inlier2{2,1}(1,:), 'bo','MarkerSize',3);

shift = size(seq.imgs{2},2);
cmap = jet(32);
k = 1;
for i = 1:size(seq.inlierId{2,1},2)
    
    if ~isinf(seq.inlierId{2}(1,i))
        ptdraw = [seq.x1{2}(1,seq.inlierId{2}(1,i)), seq.x1{2}(2,seq.inlierId{2}(1,i));
                  seq.x2{2}(1,seq.inlierId{2}(1,i)), seq.x2{2}(2,seq.inlierId{2}(1,i))+shift];
        plot(ptdraw(:,2),ptdraw(:,1),'LineStyle','-','LineWidth',0.5,'Color',cmap(k,:));
        k = mod(k+1,32);
        if k == 0 k = 1;
        end
    end
    end

% subplot(2,2,3);
imshow1 = cat(2, seq.imgs{2,1}, seq.imgs{3,1});
figure;imshow(imshow1);hold on;
plot(seq.inlier1{3,1}(2,:),seq.inlier1{3,1}(1,:), 'ro','MarkerSize',3);
plot(seq.inlier2{3,1}(2,:)+size(seq.imgs{3,1},2),seq.inlier2{3,1}(1,:), 'bo','MarkerSize',3);

shift = size(seq.imgs{3},2);
cmap = jet(32);
k = 1;
for i = 1:size(seq.inlierId{3,1},2)
    
    if ~isinf(seq.inlierId{3}(1,i))
        ptdraw = [seq.x1{3}(1,seq.inlierId{3}(1,i)), seq.x1{3}(2,seq.inlierId{3}(1,i));
                  seq.x2{3}(1,seq.inlierId{3}(1,i)), seq.x2{3}(2,seq.inlierId{3}(1,i))+shift];
        plot(ptdraw(:,2),ptdraw(:,1),'LineStyle','-','LineWidth',0.5,'Color',cmap(k,:));
        k = mod(k+1,32);
        if k == 0 k = 1;
        end
    end
end

%find H Matrix
seq.H=cell(numImages,1);


for i=1:numImages
seq.inlier1{i}(1:2,:)=flipud(seq.inlier1{i}(1:2,:));
seq.inlier2{i}(1:2,:)=flipud(seq.inlier2{i}(1:2,:));
H =Hrecacl(seq.inlier1{i},seq.inlier2{i});
seq.H{i}=H;
end

% check answer

I1 = WarpAndBlend(seq.H{1},seq.imgs{1},seq.imgs{2});
subplot(2,2,1);
imshow(I1);
subplot(2,2,2);
I2 = WarpAndBlend(seq.H{2},seq.imgs{1},seq.imgs{3});
imshow(I2);
subplot(2,2,3);
I3 = WarpAndBlend(seq.H{3},seq.imgs{2},seq.imgs{3});
imshow(I3);

% compute out put limits for each transform
seq.tform=cell(numImages,1);
seq.T=cell(numImages,1);
for i = 1:numel(seq.H)
    H=seq.H{i};
    img=seq.imgs{i};
    tform = maketform('projective',H);
    [B,xlim(i,:),ylim(i,:)] = imtransform(img,tform); 
    tforms(i)=tform;
end


% find center image
avgXLim = mean(xlim, 2);
[~, idx] = sort(avgXLim);
centerIdx = floor((numel(seq.tform)+1)/2);
centerImageIdx = idx(centerIdx);

%recompute
Tinv=inv(tforms(centerImageIdx).tdata.T);
for i = 1:numel(tforms)
img=seq.imgs{i};
    tforms(i).tdata.T = tforms(i).tdata.T * Tinv;
     [~,xlim(i,:),ylim(i,:)] = imtransform(img,tforms(i)); 
end

maxImageSize = max(imageSize);

% Find the minimum and maximum output limits
xMin = min([1; xlim(:)]);
xMax = max([maxImageSize(2); xlim(:)]);

yMin = min([1; ylim(:)]);
yMax = max([maxImageSize(1); ylim(:)]);
blender = vision.AlphaBlender('Operation', 'Binary mask', ...
    'MaskSource', 'Input port');
% Create a 2-D spatial reference object defining the size of the panorama.
xLimits = [xMin xMax];
yLimits = [yMin yMax];

% Width and height of panorama.
width  = round(xMax - xMin);
height = round(yMax - yMin);

% Initialize the "empty" panorama.

panorama = zeros([height width 3] ,'like', I);


% Create the panorama.
for i = 1:numImages

    I = readimage(buildingScene, i);

    % Transform I into the panorama.
    warpedImage = imtransform(I, tforms(i), 'XData',xLimits,'YData',yLimits,'FillValues',zeros(size(I,3),1));
    mask = imtransform(zeros(imageSize(1,1),imageSize(1,2)), tforms(1), 'XData',xLimits,'YData',yLimits);
    panorama = step(blender, panorama, warpedImage, mask);

end

figure
imshow(panorama);
