function corners = Harris(im, num)
    % harris
    Ix = edge(im,'Sobel',[],'horizontal');
    Iy = edge(im,'Sobel',[],'vertical');

    Ixx = Ix.*Ix;
    Iyy = Iy.*Iy;
    Ixy = Ix.*Iy;

    g = fspecial('gaussian', 3, 1);% try play with this two values to see what will hapen
    Ixx = imfilter(Ixx,g,'replicate');
    Iyy = imfilter(Iyy,g,'replicate');
    Ixy = imfilter(Ixy,g,'replicate');

    C = Ixx.*Iyy - Ixy.^2 - 0.04.*(Ixx+Iyy).^2;
    % figure;imshow(C,[]);

    threshold = 0.3*max(C(:));
%     Cb = imbinarize(C,threshold);
%     figure;imshow(Cb,[]);
    C(C(:)<threshold) = 0;

    C = nonmaxima_suppression(C, 3);
    [y,x] = find(C>0);

    [~,id] = sort(C,'descend');
    if num > length(y)
        corners = [y x];
    else
        corners = [y(id(1:num)) x(id(1:num))];
    end
end