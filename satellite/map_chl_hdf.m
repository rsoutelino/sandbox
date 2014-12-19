% Oceano Leste II:  01/02/2005 a 12/03/2005
% Cruzeiro Abrolhos 2: 12 a 22/03/2005
% -22S a -17S
% -41W a -35W
% chl=base^((slope*l3m_data)+intercept);
clear all;close all;
base=10;
slope=0.015;
intercept=-2;
x=linspace(-180,180,8640);
y=linspace(-90,90,4320);
g=(x>=-42)&(x<=-32);
x0=find(x==min(x(g)));
x1=find(x==max(x(g)));
lon=x(x0:x1);
g=(y>=-22)&(y<=-10);
y0=find(y==min(y(g)));
y1=find(y==max(y(g)));
lat=y(y0:y1);
%A20050332005040.L3m_8D_CHLO_4  
f=['20050332005040';'20050492005056';'20050652005072';'20050412005048'; ...
   '20050572005064';'20050732005080'];
d=[ '02 Feb 2005 to 09 Feb 2005';
    '10 Feb 2005 to 17 Feb 2005';
    '18 Feb 2005 to 25 Feb 2005';
    '26 Feb 2005 to 05 Mar 2005';
    '06 Mar 2005 to 13 Mar 2005';
    '14 Mar 2005 to 21 Mar 2005'];



k=0;
for i=1:2,
 for j=1:3,
  k=k+1;
  file=['A' f(k,:) '.L3m_8D_CHLO_4'];
  c=double(flipud(hdfread(file,'l3m_data','Index',{[1.0 1.0],[1.0 1.0],[4320.0 8640.0]})));
  c=c(y0:y1,x0:x1);
  g=c~=255;
  l10c=(slope*c)+intercept;
  figure(1)
  imagesc(lon,lat,l10c.*double(g)-999*double(~g),[-1.5 0]);
  axis('xy','equal','tight');title(d(k,:));degx;degy;colorbar('h')
%    subplot(3,2,k),imagesc(lon,lat,l10c.*double(g)-999*double(~g),[-1.5 0]);
%    axis('xy','equal','tight');title(d(k,:));degx;degy;colorbar('h')
 end
end
