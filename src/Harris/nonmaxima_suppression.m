function result = nonmaxima_suppression(score, hsize)
    % hsize = 5;
    [x,y]=meshgrid(-hsize:hsize,-hsize:hsize);% create a hsize x hsize index template
    [height,width]=size(score);% fetch image size
    x = x(:);y=y(:);% from a 2D matrix to a vector
    id = x == 0 & y == 0;% index of centering pixel
    x(id) = []; y(id) = [];% remove the index of the center pixel 
    result = zeros(height,width);% result
    for h = 1:height
        for w = 1:width
            if score(h,w) == 0% if no response, skip, there is a faster way: featch non zero index first, then traverse 
                continue;
            end
            nn = [h+x w+y];% index array of nn pixel
            valid = nn(:,1) > 0 & nn(:,1) <= height & nn(:,2)>0 & nn(:,2) <= width;% filter out pixels out of boarder
            nn = nn(valid,:);
            % sub to ind, image is a 2D array, so you can visit by I(i,j) or by I(j*height+i) since 
            % in Matlab, array is column-majored.
            ind0 = sub2ind([height,width],nn(:,1),nn(:,2));
            maxval = max(score(ind0));% maximum value
            if score(h,w) > maxval
                result(h,w) = 1;
            else
                result(h,w) = 0;
            end
        end
    end
end
