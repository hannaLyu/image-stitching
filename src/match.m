addpath('features/');
addpath('precompute/');

I1=imread('chickenbroth_01.jpg');
I2=imread('chickenbroth_02.jpg');
[s11,s12,s13] = size(I1);
[s21,s22,s23] = size(I2);

if s11 * s12 ~= s21*s22
    if s11*s12<s21*s22
        I2 = imresize(I2,[s11,s12]);
    else
        I1 = imresize(I1,[s21,s22]);
    end
end

if size(I1,3) == 3
    im1gray = rgb2gray(I1);
else
    im1gray = I1;
end

if size(I2,3) == 3
    im2gray = rgb2gray(I2);
else
    im2gray = I2;
end


 corner1=fast(im1gray);
 corner2=fast(im2gray);
% 
% corner1 = fast_corner_detector(im1gray, 500);
% corner2 = fast_corner_detector(im2gray, 500);

img1 = im2double(im1gray);
img2= im2double(im2gray);

brief_pattern;
% descriptor1= extract_brief_descriptor(img1,corner1,pattern);
% descriptor2= extract_brief_descriptor(img2,corner2,pattern);
descriptor1 = brief_descriptor(img1,corner1,pattern);
descriptor2 = brief_descriptor(img2,corner2,pattern);

% matching_pairs = bruteforce(descriptor1,descriptor2);
matchingpairs = brief_matching(descriptor1, descriptor2);

imshow1 = cat(2, img1, img2);
figure;imshow(imshow1);hold on;

plot(corner1(:,2),corner1(:,1), 'ro','MarkerSize',5);
plot(corner2(:,2)+size(I1,2),corner2(:,1), 'bo','MarkerSize',5);
