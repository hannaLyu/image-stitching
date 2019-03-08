
function d = dist2line(ph, params)
    params = params ./ norm(params(1:2));
    d = abs(ph' * params);
end