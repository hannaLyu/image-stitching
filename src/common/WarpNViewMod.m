function Out = WarpNViewMod(H,ImL,ImR)



T = maketform('projective',H');
%% do homography transformation, find xmin, xmax, ymin, ymax
[ImH,XData,YData]=imtransform(ImL,T,'FillValues',[0;0;0]);

%% do the same for a mask
Mask=ones(size(ImL,1),size(ImL,2));%ones(480,640);
MaskH=imtransform(Mask,T,'FillValues',0)>0;

%% find out output image size
XData=[floor(XData(1)) ceil(XData(2))];
YData=[floor(YData(1)) ceil(YData(2))];
nX=max(size(ImR,2),XData(2))-min(0,XData(1));
nY=max(size(ImR,1),YData(2))-min(0,YData(1));

OrigoWarp = [0 0];
OrigoRaw = [0 0];

if XData(1) > 0
    OrigoRaw(1) = XData(1);
else
    OrigoRaw(1) = -XData(1);
end
if YData(1) > 0
    OrigoRaw(2) = YData(1);
else
    OrigoRaw(2) = -YData(1);
end

M=ones(nY,nX)*1e-3;
M(OrigoWarp(2)+(1:size(ImH,1)),OrigoWarp(1) + (1:size(ImH,2))) = MaskH+1e-3;
M([OrigoRaw(2)+(1:size(ImR,1))],[OrigoRaw(1)+(1:size(ImR,2))])=...
    M([OrigoRaw(2)+(1:size(ImR,1))],[OrigoRaw(1)+(1:size(ImR,2))])+1;
Mask=zeros(nY,nX,3);
Mask(:,:,1)=M;
Mask(:,:,2)=M;
Mask(:,:,3)=M;

Out=zeros(nY,nX,3);
Out(OrigoWarp(2)+(1:size(ImH,1)),OrigoWarp(2) + (1:size(ImH,2)),:)=...
    double(ImH)./Mask(OrigoWarp(2)+(1:size(ImH,1)),OrigoWarp(2) + (1:size(ImH,2)),:);
Out([OrigoRaw(2)+(1:size(ImR,1))],[OrigoRaw(1)+(1:size(ImR,2))],:)=...
    Out([OrigoRaw(2)+(1:size(ImR,1))],[OrigoRaw(1)+(1:size(ImR,2))],:)+...
    double(ImR)./Mask([OrigoRaw(2)+(1:size(ImR,1))],[OrigoRaw(1)+(1:size(ImR,2))],:);

% %% if x>0, than origo x is x, otherwise
% % Origo=[-XData(1) -YData(1)];
% Origo = [XData(1) YData(1)];
% if Origo(1) < 0
%     Origo(1) = Origo(1)*-1;
% end
% if Origo(2) < 0
%     Origo(2) = Origo(2)*-1;
% end
% 
% M=ones(nY,nX)*1e-3;
% M(1:size(ImH,1),1:size(ImH,2))=MaskH+1e-3;
% M([Origo(2)+(1:size(ImR,1))],[Origo(1)+(1:size(ImR,2))])=...
%     M([Origo(2)+(1:size(ImR,1))],[Origo(1)+(1:size(ImR,2))])+1;
% Mask=zeros(nY,nX,3);
% Mask(:,:,1)=M;
% Mask(:,:,2)=M;
% Mask(:,:,3)=M;
% 
% 
% Out=zeros(nY,nX,3);
% 
% Out(1:size(ImH,1),1:size(ImH,2),:)=double(ImH)./Mask(1:size(ImH,1),1:size(ImH,2),:);
% 
% Out([Origo(2)+(1:size(ImR,1))],[Origo(1)+(1:size(ImR,2))],:)=...
% Out([Origo(2)+(1:size(ImR,1))],[Origo(1)+(1:size(ImR,2))],:)+...
% double(ImR)./Mask([Origo(2)+(1:size(ImR,1))],[Origo(1)+(1:size(ImR,2))],:);
% 
% Out=Out/255;
% % imagesc(Out)
