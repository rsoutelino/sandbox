function  s=smoo(a,bad,wid,flat);  
%   function  s=smoo(a,bad,wid,flat);  
%  'bad' is the bad data flag. 
%  'wid' is the size of the square matrix.
%  'flat' controls the weights. if 'flat' is less than one,
%  the elements in the 'wid x wid' window tend to have 
%  similar weights (more smooth). If 'flat' is larger than 
%  one, more weigth is given to the central region (less smooth).
%
%  calculate smoothing function
r=floor(wid/2);
[x,y]=meshgrid(-flat:flat/r:flat,-flat:flat/r:flat);
z=(erf(sqrt(x.^2+y.^2))).^2;z=-z+max(z(:));z=z/max(z(:));
%  pre-load matrices
b=ones(size(a));s=zeros(size(a));c=zeros(size(a));
%
gud=a~=bad;
ngud=(~gud);
a=a.*gud;
b=b.*gud;
%
for i=-r:1:r,
 ii=i+r+1;
 for j=-r:1:r,
  jj=j+r+1;
  s=s+roam(a,i,j)*z(ii,jj);
  c=c+roam(b,i,j)*z(ii,jj);
 end
end
cz=c==0;c=c+cz;
s=(s./c);s=s+bad*ngud;
clear a b c cz gud ngud bad size flat r x y z i ii j jj
