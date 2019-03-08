
function p = fromhomogeneous(ph)
    p = ph(1:2,:) ./ ph(3,:);
end
