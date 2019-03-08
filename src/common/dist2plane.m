
function d = dist2plane(ph, params)
    params = params ./ norm(params(1:3));
    d = abs(ph' * params);
end