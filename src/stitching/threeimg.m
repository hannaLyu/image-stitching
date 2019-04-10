clc;clear;
figure;hold on;
 imgpath = 'M:\Documents\stitching\stitching\hill\hill';  
buildingScene = imageDatastore(imgpath);
montage(buildingScene.Files);
hold on;
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
seq.H{i}=inv(H);
end

% check answer based on img1

% Compute the output limits  for each transform
[~,xlim(1,:),ylim(1,:)]=imtransform(zeros(size(seq.imgs{2},1),size(seq.imgs{2},2)),maketform('projective',eye(3)'));
[~,xlim(2,:),ylim(2,:)]=imtransform(zeros(size(seq.imgs{2},1),size(seq.imgs{2},2)),maketform('projective',seq.H{1}'));
[~,xlim(3,:),ylim(3,:)]=imtransform(zeros(size(seq.imgs{3},1),size(seq.imgs{3},2)),maketform('projective',seq.H{2}'));
xMin = min([1; xlim(:)]);
xMax = max(xlim(:));

yMin = min([1; ylim(:)]);
yMax = max(ylim(:));

% Create a 2-D spatial reference object defining the size of the panorama.
xLimits = [xMin xMax];
yLimits = [yMin yMax];
% warpedImage=cell(3,1);
warpedImage1 = imtransform(seq.imgs{1}, maketform('projective',eye(3)), 'XData',xLimits,'YData',yLimits,'FillValues',zeros(3,1));
warpedImage2 = imtransform(seq.imgs{2}, maketform('projective',seq.H{1}'), 'XData',xLimits,'YData',yLimits,'FillValues',zeros(3,1));
warpedImage3 = imtransform(seq.imgs{3}, maketform('projective',seq.H{2}'), 'XData',xLimits,'YData',yLimits,'FillValues',zeros(3,1));

imshow(warpedImage1);
imshow(warpedImage2);
imshow(warpedImage3);
%  I1 = WarpAndBlend(seq.H{1},warpedImage1,warpedImage2);
%  imshow(I1);
%  hold on;
%  I2=WarpAndBlend(seq.H{2},I1,warpedImage3);
%  imshow(I2);

for i = 1:numImages
    H=seq.H{i};
    img=seq.imgs{i};
    tform = maketform('projective',inv(H));
%     [~,xlim(i,:),ylim(i,:)] = imtransform(img,tform); 
    tforms(i)=tform;
end



mask = imtransform(zeros(imageSize(1,1),imageSize(1,2)), tforms(1), 'XData',xLimits,'YData',yLimits);

width  = size(mask,2);
height = size(mask,1);
margin=20;
panorama = zeros([height width 3]);
weights = zeros([height width]);
% for i = 1:numImages
    I = readimage(buildingScene, 1);
    % Transform I into the panorama.
  
    dist1 = zeros(size(I,1),size(I,2));%dist1(round(size(I,1)*0.5),round(size(I,2)*0.5)) = 1;
    dist1(1:margin,:) = 1;dist1(end-margin:end,:) = 1;dist1(:,1:margin) = 1;dist1(:,end-margin:end) = 1;
    dist1 = bwdist(dist1,'euclidean');%maxdist = max(dist1(:));dist1 = (maxdist+1) - dist1;
    dist1t=imtransform(dist1,tforms(1),'XData',xLimits,'YData',yLimits,'FillValues',0) + 1e-3;
    [panorama,weights] = simpleBlend(warpedImage1, panorama, weights, dist1t);
    I = readimage(buildingScene, 2);
     dist2 = zeros(size(I,1),size(I,2));%dist1(round(size(I,1)*0.5),round(size(I,2)*0.5)) = 1;
    dist2(1:margin,:) = 1;dist2(end-margin:end,:) = 1;dist2(:,1:margin) = 1;dist2(:,end-margin:end) = 1;
    dist2 = bwdist(dist2,'euclidean');%maxdist = max(dist1(:));dist1 = (maxdist+1) - dist1;
    dist2t=imtransform(dist2,tforms(2),'XData',xLimits,'YData',yLimits,'FillValues',0) + 1e-3;
    [panorama,weights] = simpleBlend(warpedImage2, panorama, weights, dist2t);
    I = readimage(buildingScene, 3);
     dist3 = zeros(size(I,1),size(I,2));%dist1(round(size(I,1)*0.5),round(size(I,2)*0.5)) = 1;
    dist3(1:margin,:) = 1;dist3(end-margin:end,:) = 1;dist3(:,1:margin) = 1;dist3(:,end-margin:end) = 1;
    dist3 = bwdist(dist3,'euclidean');%maxdist = max(dist1(:));dist1 = (maxdist+1) - dist1;
    dist3t=imtransform(dist3,tforms(3),'XData',xLimits,'YData',yLimits,'FillValues',0)+1e-3 ;
    [panorama,weights] = simpleBlend(warpedImage3, panorama, weights, dist3t);
% end
weights = weights+1e-3;
panorama = panorama ./ weights;
figure;
imshow(panorama);






for i = 1:numImages-1
    H=seq.H{i};
    img=seq.imgs{i};
    tform = maketform('projective',inv(H));
    [~,xlim(i,:),ylim(i,:)] = imtransform(img,tform); 
    tforms(i)=tform;
end

function [panorama,weights] = simpleBlend(im1, panorama, weights, dist)
    weight = dist./(max(dist(:)));
    weights = weights + weight;
    if size(im1,3) == 1
        panorama = combine(im1, panorama, weight);
    else
        panorama(:,:,1) = combine(im1(:,:,1), panorama(:,:,1), weight);
        panorama(:,:,2) = combine(im1(:,:,2), panorama(:,:,2), weight);
        panorama(:,:,3) = combine(im1(:,:,3), panorama(:,:,3), weight);
    end
end

function im12 = combine(im1, im2, weight)
    im12 = im2;
    im12 = im12 + weight.*im2double(im1);
end

   
 
 


