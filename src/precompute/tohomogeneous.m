function ph = tohomogeneous(p)
    ph = [p;ones(1,size(p,2))];
end