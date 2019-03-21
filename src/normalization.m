function [T,x_hat] = normalization(p)

pbar=mean(p,2);
pdiff=p(1:2,:);
mdist=mean(sqrt(diag(pdiff'*pdiff)));
scale=sqrt(2)/mdist;
pnormalized(1:2,:)=scale.*pdiff;
pnormalized(3,:)=1;
T=[scale 0 -pbar(1)*scale;
    0 scale -pbar(2)*scale;
    0 0 1];
x_hat=[T*p(:,1) T*p(:,2) T*p(:,3) T*p(:,4)];
end

