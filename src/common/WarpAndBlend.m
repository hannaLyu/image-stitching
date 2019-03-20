function Out = WarpAndBlend(H,ImL,ImR)
    T = maketform('projective',H');
    % do homography transformation, find xmin, xmax, ymin, ymax
    [~,XData,YData]=imtransform(ImL,T,'FillValues',zeros(size(ImL,3),1));

    % find out output image size
    XData=[floor(XData(1)) ceil(XData(2))];
    YData=[floor(YData(1)) ceil(YData(2))];
    nX = [min(0,XData(1)) max(size(ImR,2),XData(2))];
    nY = [min(0,YData(1)) max(size(ImR,1),YData(2))];

    [ImLH]=imtransform(ImL,T,'XData',nX,'YData',nY,'FillValues',zeros(size(ImL,3),1));
    [ImRH]=imtransform(ImR,maketform('affine',eye(3)),'XData',nX,'YData',nY,'FillValues',zeros(size(ImL,3),1));

    % do the same for a mask
    Mask1=ones(size(ImL,1),size(ImL,2));
    Mask2=ones(size(ImR,1),size(ImR,2));

    Mask1=imtransform(Mask1,T,'XData',nX,'YData',nY,'FillValues',0)>0;
    Mask2=imtransform(Mask2,maketform('affine',eye(3)),'XData',nX,'YData',nY,'FillValues',0)>0;
    Maskavg = Mask1 & Mask2;
    Out = uint8(zeros(size(ImLH)));
    
%     blender = vision.AlphaBlender('Operation', 'Binary mask', ...
%         'MaskSource', 'Input port');
% 
%     Out = step(blender, ImLH, ImRH, Mask2);

    if size(ImL,3) == 1
        Out = combine(ImLH, ImRH, Mask1, Mask2, Maskavg);
    else
        Out(:,:,1) = combine(ImLH(:,:,1), ImRH(:,:,1), Mask1, Mask2, Maskavg);
        Out(:,:,2) = combine(ImLH(:,:,2), ImRH(:,:,2), Mask1, Mask2, Maskavg);
        Out(:,:,3) = combine(ImLH(:,:,3), ImRH(:,:,3), Mask1, Mask2, Maskavg);
    end
end

function im12 = combine(im1, im2, mask1, mask2, mask12)
    im12 = uint8(zeros(size(im1)));
    im12(mask1) = im1(mask1);
    im12(mask2) = im2(mask2);
    im12(mask12) = uint8(0.5*im1(mask12)+0.5*im2(mask12));
end