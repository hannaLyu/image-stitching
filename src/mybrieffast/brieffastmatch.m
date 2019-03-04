img1=imread('chickenbroth_01.jpg');
img2=imread('chickenbroth_02.jpg');
[s11,s12,s13] = size(img1);
[s21,s22,s23] = size(img2);

if s11 * s12 ~= s21*s22
    if s11*s12<s21*s22
        img2 = imresize(img2,[s11,s12]);
    else
        img1 = imresize(img1,[s21,s22]);
    end
end

corner1=fast(img1);
corner2=fast(img2);
patten=brief_pattern_generator;
descriptor1=extractdes(img1,corner1,patten);
descriptor2=extractdes(img2,corner2,patten);
matchingpairs=bruteforce(descriptor1,descriptor2);
imshow1 = cat(2, img1, img2);
figure;imshow(imshow1);hold on;

plot(corner1(:,2),corner1(:,1), 'ro','MarkerSize',3);
plot(corner2(:,2)+size(img1,2),corner2(:,1), 'bo','MarkerSize',3);

shift = size(im1,2);
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