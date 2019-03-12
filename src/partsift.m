clear all
run('E:\Aѧϰ\Special course\SIFT\vlfeat-0.9.21-bin\vlfeat-0.9.21/toolbox/vl_setup');

I=vl_impattern('roofs1');
image(I);
I = single(rgb2gray(I));
[f,d] = vl_sift(I) ;
perm = randperm(size(f,2)) ;
sel = perm(1:50) ;
h1 = vl_plotframe(f(:,sel)) ;
h2 = vl_plotframe(f(:,sel)) ;
set(h1,'color','k','linewidth',3) ;
set(h2,'color','y','linewidth',2) ;
h3 = vl_plotsiftdescriptor(d(:,sel),f(:,sel)) ;
set(h3,'color','g') ;

% imga=vl_impattern('roofs1');
% imgb=vl_impattern('roofs2');


Ia =imread('pf_desk.jpg');
Ib = imread('pf_stand.jpg');
[s11,s12,s13] = size(Ia);
[s21,s22,s23] = size(Ib);

if s11 * s12 ~= s21*s22
    if s11*s12<s21*s22
        Ib = imresize(Ib,[s11,s12]);
    else
        Ia= imresize(Ia,[s21,s22]);
    end
end

if size(Ia,3) == 3
    im1gray = rgb2gray(Ia);
else
    im1gray = Ia;
end

if size(Ib,3) == 3
    im2gray = rgb2gray(Ib);
else
    im2gray = Ib;
end

img1 = single(im1gray);
img2 = single(im2gray);
[fa, da] = vl_sift(img1) ;
[fb, db] = vl_sift(img2) ;
[matches, scores] = vl_ubcmatch(da, db) ;
matches=transpose(matches);

imshow1 = cat(2, Ia, Ib);
figure;imshow(imshow1);hold on;
im1double=im2double(rgb2gray(Ia));
 plot(fa(1,1:end),fa(2,1:end), 'ro','MarkerSize',5);
 plot(fb(1,1:end)+size(im1double,2),fb(2,1:end), 'bo','MarkerSize',5);

shift = size(im1double,2);
cmap = jet(32);
k = 1;
for i = 1:size(matches,1)
    if ~isinf(matches(i,2))
        ptdraw = [fa(2,matches(i,1)), fa(1,matches(i,1));
                  fb(2,matches(i,2)), fb(1,matches(i,2))+shift];
        plot(ptdraw(:,2),ptdraw(:,1),'LineStyle','-','LineWidth',0.5,'Color',cmap(k,:));
        k = mod(k+1,32);if k == 0 k = 1;end
    end
end
