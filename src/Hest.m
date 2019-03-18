function H1 = Hest(x1h,x2h)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



x1h1=transpose(x1h(:,1));
x1h2=transpose(x1h(:,2));
x1h3=transpose(x1h(:,3));
x1h4=transpose(x1h(:,4));

x2h1=transpose(x2h(:,1));
x2h2=transpose(x2h(:,2));
x2h3=transpose(x2h(:,3));
x2h4=transpose(x2h(:,4));

A=[zeros(1,3) -x1h1 x2h1(2)*x1h1;
    x1h1 zeros(1,3) -x2h1(1)*x1h1;
    zeros(1,3) -x1h2 x2h2(2)*x1h2;
    x1h2 zeros(1,3) -x2h2(1)*x1h2;
    zeros(1,3) -x1h3 x2h3(2)*x1h3;
    x1h3 zeros(1,3) -x2h3(1)*x1h3;
    zeros(1,3) -x1h4 x2h4(2)*x1h4;
    x1h4 zeros(1,3) -x2h4(1)*x1h4;];




 [U,S,V] = svd(A);
h=transpose(V(:,end));
H1=[h(1,1:3);h(1,4:6);h(1,7:9)];


end

