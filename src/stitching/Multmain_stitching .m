
clc;clear ;
% corner detection
I1=imread('right.jpg');
I2=imread('left.jpg');
[s11,s12,s13] = size(I1);
[s21,s22,s23] = size(I2);


I1=imresize(I1,[300,400]);
I2=imresize(I2,[300,400]);


% if s11 * s12 ~= s21*s22
%     if s11*s12<s21*s22
%         I2 = imresize(I2,[s11,s12]);
%     else
%         I1 = imresize(I1,[s21,s22]);
%     end
% end

im1gray = rgb2gray(I1);
im2gray=rgb2gray(I2);
corner1=fast(I1);
figure;hold on;
corner2=fast(I2);

%extract discriptor

patten=brief_pattern_generator;
descriptor1=extractdes(im1gray,corner1,patten);
descriptor2=extractdes(im2gray,corner2,patten);
matchingpairs=bruteforce(descriptor1,descriptor2);
imshow1 = cat(2, I1, I2);
figure;imshow(imshow1);hold on;
plot(corner1(:,2),corner1(:,1), 'ro','MarkerSize',3);
plot(corner2(:,2)+size(I1,2),corner2(:,1), 'bo','MarkerSize',3);

matchingpairs = matchingpairs(~isinf(matchingpairs(:,2)),:);

shift = size(I1,2);
cmap = jet(32);
k = 1;
for i = 1:size(matchingpairs,1)
    if ~isinf(matchingpairs(i,2))
        ptdraw = [corner1(matchingpairs(i,1),1), corner1(matchingpairs(i,1),2);
                  corner2(matchingpairs(i,2),1), corner2(matchingpairs(i,2),2)+shift];
        plot(ptdraw(:,2),ptdraw(:,1),'LineStyle','-','LineWidth',0.5,'Color',cmap(k,:));
        k = mod(k+1,32);
        if k == 0 k = 1;
    end
    end
end


%ransac
matchPoints1 = transpose(corner1(matchingpairs(:,1), :));
matchPoints2 = transpose(corner2(matchingpairs(:,2), :));
x1 = tohomogeneous(matchPoints1);
x2 =tohomogeneous(matchPoints2);

[inlier1,inlier2,inlierId]= ransac(x1,x2);

imshow2 = cat(2, I1, I2);
figure;imshow(imshow2);hold on;
plot(inlier1(2,:),inlier1(1,:), 'ro','MarkerSize',3);
plot(inlier2(2,:)+size(I1,2),inlier2(1,:), 'bo','MarkerSize',3);

shift = size(I1,2);
cmap = jet(32);
k = 1;
for i = 1:size(inlierId,2)
    if ~isinf(inlierId(1,i))
        ptdraw = [x1(1,inlierId(1,i)), x1(2,inlierId(1,i));
                  x2(1,inlierId(1,i)), x2(2,inlierId(1,i))+shift];
        plot(ptdraw(:,2),ptdraw(:,1),'LineStyle','-','LineWidth',0.5,'Color',cmap(k,:));
        k = mod(k+1,32);
        if k == 0 k = 1;
    end
    end
end

%stitching
inlier1(1:2,:) = flipud(inlier1(1:2,:));
inlier2(1:2,:) = flipud(inlier2(1:2,:));
H =Hrecacl(inlier1,inlier2);
I3 = WarpAndBlend(H,I1,I2);
% figure
% imshow(I3);

I4=imread('hill\hill\3.JPG');
[s11,s12,s13] = size(I3);
[s21,s22,s23] = size(I4);


I3=imresize(I3,[300,400]);
I4=imresize(I4,[300,400]);

im3gray = rgb2gray(I3);
im4gray=rgb2gray(I4);
corner3=fast(I3);
figure;hold on;
corner4=fast(I4);
patten=brief_pattern_generator;
descriptor3=extractdes(im3gray,corner3,patten);
descriptor4=extractdes(im4gray,corner4,patten);
matchingpairs=bruteforce(descriptor3,descriptor4);
imshow3 = cat(2, I3, I4);
figure;imshow(imshow3);hold on;
plot(corner3(:,2),corner3(:,1), 'ro','MarkerSize',3);
plot(corner4(:,2)+size(I3,2),corner4(:,1), 'bo','MarkerSize',3);

matchingpairs = matchingpairs(~isinf(matchingpairs(:,2)),:);


shift = size(I3,2);
cmap = jet(32);
k = 1;
for i = 1:size(matchingpairs,1)
    if ~isinf(matchingpairs(i,2))
        ptdraw = [corner3(matchingpairs(i,1),1), corner3(matchingpairs(i,1),2);
                  corner4(matchingpairs(i,2),1), corner4(matchingpairs(i,2),2)+shift];
        plot(ptdraw(:,2),ptdraw(:,1),'LineStyle','-','LineWidth',0.5,'Color',cmap(k,:));
        k = mod(k+1,32);
        if k == 0 k = 1;
    end
    end
end

%ransac
matchPoints3 = transpose(corner3(matchingpairs(:,1), :));
matchPoints4 = transpose(corner4(matchingpairs(:,2), :));
x3 = tohomogeneous(matchPoints3);
x4 =tohomogeneous(matchPoints4);

[inlier3,inlier4,inlierId]= ransac(x3,x4);

imshow4 = cat(2, I3, I4);
figure;imshow(imshow4);hold on;
plot(inlier3(2,:),inlier3(1,:), 'ro','MarkerSize',3);
plot(inlier4(2,:)+size(I3,2),inlier4(1,:), 'bo','MarkerSize',3);

shift = size(I3,2);
cmap = jet(32);
k = 1;
for i = 1:size(inlierId,2)
    if ~isinf(inlierId(1,i))
        ptdraw = [x3(1,inlierId(1,i)), x3(2,inlierId(1,i));
                  x4(1,inlierId(1,i)), x4(2,inlierId(1,i))+shift];
        plot(ptdraw(:,2),ptdraw(:,1),'LineStyle','-','LineWidth',0.5,'Color',cmap(k,:));
        k = mod(k+1,32);
        if k == 0 k = 1;
    end
    end
end

%stitching
inlier3(1:2,:) = flipud(inlier3(1:2,:));
inlier4(1:2,:) = flipud(inlier4(1:2,:));
H =Hrecacl(inlier3,inlier4);
I5 = WarpAndBlend(H,I3,I4);
figure;
imshow(I5);
