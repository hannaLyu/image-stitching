
function p = fromhomogeneous(ph)
    p = ph(1:end-1,:) ./ ph(end,:);
end
