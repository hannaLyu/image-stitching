
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

imga=vl_impattern('roofs1');
imgb=vl_impattern('roofs2');


Ia = single(rgb2gray(imga));
Ib = single(rgb2gray(imgb));
[fa, da] = vl_sift(Ia) ;
[fb, db] = vl_sift(Ib) ;
[matches, scores] = vl_ubcmatch(da, db) ;
matches=transpose(matches);

imshow1 = cat(2, imga, imgb);
figure;imshow(imshow1);hold on;
imadouble=im2double(rgb2gray(imga));
 plot(fa(1,1:end),fa(2,1:end), 'ro','MarkerSize',5);
 plot(fb(1,1:end)+size(imadouble,2),fb(2,1:end), 'bo','MarkerSize',5);

shift = size(imadouble,2);
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
