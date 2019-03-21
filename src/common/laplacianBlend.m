function imblend = laplacianBlend(source, target, mask)
    layer = 5;
    for i = 1:size(source,3)
        im1 = source(:,:,i);im1 = im2double(im1);
        im2 = target(:,:,i);im2 = im2double(im2);
        
        pyr1 = GaussianPyr(im1, layer);
        pyr2 = GaussianPyr(im2, layer);
    
        lappyr1 = LaplacianPyr(pyr1, layer);
        lappyr2 = LaplacianPyr(pyr2, layer);
        
        mask = double(mask);
        pyr3 = GaussianPyr(mask, layer);

        mergedlappyr = cell(layer);
        for j = 1:layer
            mergedlappyr{j} = pyr3{j}.*lappyr1{j} + (1-pyr3{j}).*lappyr2{j};
        end

        % collapse
        imb = mergedlappyr{layer};
        for j = layer-1:-1:1
            imb = imresize(imb,[size(mergedlappyr{j},1),size(mergedlappyr{j},2)],'bilinear') + mergedlappyr{j};
        end
        imblend(:,:,i) = imb;
    end
end
function pyr = GaussianPyr(im, layer)
%     kernel = fspecial('gaussian',15,30);
    pyr = cell(layer,1);
    img = imgaussfilt(im, 0.1, 'Padding','replicate');
    pyr{1} = img;
    for i = 2:layer
        img = imgaussfilt(pyr{i-1}, 1, 'Padding','replicate');
%         img = imfilter(pyr{i-1},kernel,'replicate');
        pyr{i} = imresize(img,[round(size(pyr{i-1},1)*0.5), round(size(pyr{i-1},2)*0.5)],'bilinear');
    end
    
    % debug
%     figure
%     for i = 1:layer
%         subplot(1,layer,i);imshow(pyr{i});
%     end
end

function lappyr = LaplacianPyr(pyr, layer)
    lappyr = cell(layer,1);
    lappyr{layer} = pyr{layer};
    for i = layer-1:-1:1
        lappyr{i} = -imresize(pyr{i+1},[size(pyr{i},1),size(pyr{i},2)],'bilinear') + pyr{i};
    end
    figure
    for i = 1:layer
        subplot(1,layer,i);imshow(lappyr{i});
    end
end

