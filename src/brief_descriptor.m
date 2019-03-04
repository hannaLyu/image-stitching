function [descriptor] =brief_descriptor(I,corner,patten)
% I = im2double(I);  


  [m,n]=size(I);
  g = fspecial('gaussian', [3,1], 1);
  I = imfilter(I,g,'replicate');
  I = imfilter(I,g','replicate');

 descriptor = zeros(size(corner,1), size(patten,2),'logical');
 for  i = 1:size(corner,1)
         u1=corner(i,1)+patten(1,:);
         v1=corner(i,2)+patten(2,:);
         
         u2=corner(i,1)+patten(3,:);
         v2=corner(i,2)+patten(4,:);
      try
            ind1=sub2ind(size(I),u1,v1);
            ind2=sub2ind(size(I),u2,v2);
      catch
            error('sss');
      end
        desp1=I(ind1)';
        desp2=I(ind2)';
        desp=desp1<desp2;
        descriptor(i,:)=desp;      
end
 end

 



