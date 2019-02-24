function result = nonmaxima_suppression(score, hsize)
    % hsize = 5;
    [x,y]=meshgrid(-hsize:hsize,-hsize:hsize);
    [height,width]=size(score);
    x = x(:);y=y(:);
    id = x == 0 & y == 0;
    x(id) = []; y(id) = [];
    result = zeros(height,width);
    for h = 1:height
        for w = 1:width
            if score(h,w) == 0
                continue;
            end
            nn = [h+x w+y];
            valid = nn(:,1) > 0 & nn(:,1) <= height & nn(:,2)>0 & nn(:,2) <= width;
            nn = nn(valid,:);
            ind0 = sub2ind([height,width],nn(:,1),nn(:,2));
            maxval = max(score(ind0));
            if score(h,w) > maxval
                result(h,w) = 1;
            else
                result(h,w) = 0;
            end
        end
    end
end