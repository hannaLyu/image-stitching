function patten = brief_pattern_generator()
    clear all;close all;
    type = 'Gaussianiid';
    hsize = 9;
    bits = 256;% 256, 512
    fsize = 2*hsize;
    if strcmp(type, 'Uniform')
        X1 = round(rand(1,bits).*(fsize)-hsize);
        Y1 = round(rand(1,bits).*(fsize)-hsize);
        X2 = round(rand(1,bits).*(fsize)-hsize);
        Y2 = round(rand(1,bits).*(fsize)-hsize);
    elseif strcmp(type, 'Gaussianiid')
        Sigma = fsize^2/25;
        sigma = fsize/5;
        X1 = round(normrnd(0,sigma,[1,bits]));
        Y1 = round(normrnd(0,sigma,[1,bits]));
        X2 = round(normrnd(0,sigma,[1,bits]));
        Y2 = round(normrnd(0,sigma,[1,bits]));
        
        X1(X1>hsize) = hsize;X1(X1<-hsize) = -hsize;
        X2(X2>hsize) = hsize;X2(X2<-hsize) = -hsize;
        Y1(Y1>hsize) = hsize;Y1(Y1<-hsize) = -hsize;
        Y2(Y2>hsize) = hsize;Y2(Y2<-hsize) = -hsize;
    else
        
    end

    
    
    [Xx,Yy] = meshgrid(-hsize:1:hsize,-hsize:1:hsize);
    Xx = Xx(:);Yy = Yy(:);
    plot(Xx, Yy, '.');hold on;
    for i = 1:bits
        plot([X1(i), X2(i)],[Y1(i), Y2(i)],'-');
    end
    patten=([X1;Y1;X2;Y2]);
end