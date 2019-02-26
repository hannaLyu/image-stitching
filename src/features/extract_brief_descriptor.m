function varargout = extract_brief_descriptor(I,corners,pattern)
    % 1st pad image
    [M,N] = size(I);
    hsize = 10;
    Ienlarge = zeros(M+2*hsize+1,N+2*hsize+1);
    Ienlarge(hsize+2:M+hsize+1, hsize+2:N+hsize+1) = I;
    [MM,NN] = size(Ienlarge);
    
    corners(:,1) = corners(:,1) + hsize + 1;
    corners(:,2) = corners(:,2) + hsize + 1;
    
    g = fspecial('gaussian', [3,1], 1);% try play with this two values to see what will hapen
    I = imfilter(I,g,'replicate');% horizontal
    I = imfilter(I,g','replicate');% vertical
    
    brief_desp = zeros(size(corners,1), size(pattern,2),'logical');
    % 2nd compare series
    for i = 1:size(corners,1)
        %
        u1 = corners(i,1) + pattern(1,:);
        v1 = corners(i,2) + pattern(2,:);
        
        u2 = corners(i,1) + pattern(3,:);
        v2 = corners(i,2) + pattern(4,:);
        try
            ind0 = sub2ind([MM,NN],u1,v1);
            ind1 = sub2ind([MM,NN],u2,v2);
        catch
            error('sss');
        end
        val0 = Ienlarge(ind0)';
        val1 = Ienlarge(ind1)';
        desp = val0 < val1;
        brief_desp(i,:) = desp;
    end
    varargout{1} = brief_desp;
end