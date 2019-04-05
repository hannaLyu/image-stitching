function H =Hrecacl(inlier1,inlier2)
number=size(inlier1,2);

for i=1:number
    p1=transpose(inlier1(:,i));
    p2=transpose(inlier2(:,i));
    pre1(i,:)=[zeros(1,3) -p1 p2(2)*p1] ;
    pre2(i,:)=[p1 zeros(1,3) -p2(1)*p1];
    
end
A(1:2:2*number,:)=pre1;
A(2:2:2*number,:)=pre2;


 [U,S,V] = svd(A);
h=transpose(V(:,end));
% h=h./h(9);
H=[h(1,1:3);h(1,4:6);h(1,7:9)];
H=H./H(3,3);
end

