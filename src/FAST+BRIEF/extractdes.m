function descriptor = extractdes(im,corner,patten)
 [M,N] = size(im);
  g = fspecial('gaussian', [3,1], 1);
    I = imfilter(im,g,'replicate');
    I = imfilter(im,g','replicate');
        descriptor = zeros(size(corner,1), size(patten,2),'logical');

    for i = 1:size(corner,1)
       
        u1 = corner(i,1) + patten(1,:);
        v1 = corner(i,2) + patten(2,:);
        
        u2 = corner(i,1) + patten(3,:);
        v2 = corner(i,2) + patten(4,:);
        
        % out of border, skip
        invalid = find(u1 < 1 | u1 > M | v1 < 1 | v1 > N,1);
        
        if isempty(invalid)
            ind1 = sub2ind([M,N],u1,v1);
            ind2 = sub2ind([M,N],u2,v2);
        else
            continue;
            error('sss');
        end
        desp1 = im(ind1)';
        desp2 = im(ind2)';
        desp = desp1 < desp2;
        descriptor(i,:) = desp;
    end
end

