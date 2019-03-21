function imblend = TwoBandBlend(source, target, mask)
    M = floor(log2(max(size(source))));
    for i = 1:size(source,3)
        im1p = cell(M,1);
        im2p = cell(M,1);
        im1p{1} = source(:,:,i);
        im2p{1} = target(:,:,i);
        mp = cell(M,1);
        mp{1} = mask;
        % Gaussian pyramid
        for n = 2 : M
            % downsample image
            im1p{n} = imresize(im1p{n-1}, 0.5);
            im2p{n} = imresize(im2p{n-1}, 0.5);
            % downsample blending mask
            mp{n} = imresize(mp{n-1}, 0.5, 'bilinear');
        end
 
        % Laplician pyramid
        for n = 1 : M-1
            im1p{n} = im1p{n} - imresize(im1p{n+1}, [size(im1p{n},1), size(im1p{n},2)]);
            im2p{n} = im2p{n} - imresize(im2p{n+1}, [size(im2p{n},1), size(im2p{n},2)]);   
        end   
 
        % Multi-band blending Laplician pyramid
        for n = 1 : M
            imp{n} = im1p{n} .* mp{n} + im2p{n} .* (1-mp{n});
        end
 
        % Laplician pyramid reconstruction
        im = imp{M};
        for n = M-1 : -1 : 1
            im = imp{n} + imresize(im, [size(imp{n},1) size(imp{n},2)]);
        end
    end
end

function pyr = GaussianPyr(im, layer)
    pyr = cell(layer,1);
    pyr{1} = im;
    for i = 2:layer
        img = imgaussfilt(pyr{i-1}, 3, 'Padding','replicate');
        pyr{i} = imresize(img,[round(size(pyr{i-1},1)*0.5), round(size(pyr{i-1},2)*0.5)],'bilinear');
    end
end

function lappyr = LaplacianPyr(pyr, layer)
    lappyr = cell(layer,1);
    lappyr{layer} = pyr{layer};
    for i = layer-1:-1:1
        lappyr{i} = -imresize(pyr{i+1},[size(pyr{i},1),size(pyr{i},2)],'bilinear') + pyr{i};
    end
end

