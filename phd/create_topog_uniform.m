%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ETOPO 1 - based topography
% Creates sea bottom grids for ROMS runs - PhD work
% Eastern brazilian coast
% Creates real topography and semi-idealized topography removing
% seamounts, banks and other steep topographic features
% Rafael Soutelino, June, 2010
% rsoutelino@gmail.com	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc; warning off

%%% BASIC SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loading info for ETOPO 1 database: 
url = ['ftp://ftp.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/bedrock/grid_registered/netcdf/'];
etopo1file = ['/Users/rsoutelino/rsoutelino/misc/etopo1_3901x6601_atl_sul.xyz'];
filename = ['etopo1_1981x1801_atl_sul']; l1 = 6601; l2 = 3901;
newfile = 'n';
mycolorbar = ['cc = colorbar(''west''); pos = get(cc,''position''); set(cc,''position'',[pos(1) pos(2) pos(3)/3 pos(4)])'];

% declaring geographical limits
lonlim = [-54 -20]; latlim = [-32 -4]; m_proj('mercator','lon',lonlim,'lat',latlim);
load costa.mat; %  m_gshhs_f('save','costa.mat');

% topography editting parameters
shelfsize  = 20; % fixed number of grid points for the shelf (when removing the banks) => depends upon original topography
shelfbreak = 70; % assuming this depth for the shelf break (when removing the banks)
win = 51; % window size for topography smoothing
load cont_ocean.txt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% LOADING VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if newfile == 's'
   load(etopo1file);
   eval(['etopo1 = ',filename,';']);
   xb = etopo1(:,1); yb = etopo1(:,2); zb = etopo1(:,3);
   xb = reshape(xb,l1,l2); yb = reshape(yb,l1,l2); zb = reshape(zb,l1,l2); 
   xb = xb'; yb = yb'; zb = zb';
   f1 = find(yb(:,1) >= latlim(1) & yb(:,1) <= latlim(2));
   f2 = find(xb(1,:) >= lonlim(1) & xb(1,:) <= lonlim(2));
   xb = xb(f1,f2); yb = yb(f1,f2); zb = zb(f1,f2); 

else 
   load etopo1.mat
   xb = lon; yb = lat; zb = topo; [xb,yb] = meshgrid(xb,yb);
   f1 = find(yb(:,1) >= latlim(1) & yb(:,1) <= latlim(2));
   f2 = find(xb(1,:) >= lonlim(1) & xb(1,:) <= lonlim(2));
   xb = xb(f1,f2); yb = yb(f1,f2); zb = zb(f1,f2); 
end

% trying to smooth it before computations
[a,b] = size(zb); zb1 = zb;
disp('Initial smoothing....')
zb = smoo2(zb,-9999,51,1);

%%% CREATING SEMI-IDEALIZED TOPOGRAPHY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% finding a sample profile 
fex = find(yb(:,1)==-15.6);
zex = zb(fex,:); xex = xb(fex,:);
zex = zex(find(zex<=0)); xex = xex(find(zex<=0));

% fitting it into a polynome to avoid noise
zfit = zex(20:end); xfit = xex(20:end);
pol = polyfit(xfit,zfit,10); 
zfit = polyval(pol,xfit);
zex = [zex(1:20-1) zfit];
zex = zex(1:end-170); xex = xex(1:end-170);

% scanning parallels ----------------------------------------
tic
for k = 1:a;
%     figure(1); pl = plot(xb(k,:),yb(k,:),'w');

   xk = xb(k,:); zk = zb(k,:);   

   figure(2)
   plot(xk,zk,'g'); hold on; axis([lonlim -6500 1000])
   plot([lonlim],[0 0],'r--'); % grid on;
   fcont = find(zk >= 0); fcoast = fcont(end);
   plot([xk(fcoast) xk(fcoast)],[-6500 500],'r--');

   zw =  zk(fcoast:b); xw = xk(fcoast:b);    
   dif = length(zw) - length(zex);

   if dif == 0
      zk(fcoast:b) = zw;
   elseif dif < 0
      ex = abs(dif);
      zk(fcoast:b) = zex(1:[length(zex)-ex]);
   else
      zend = zex(end).*ones(1,dif);
      zaux = [zex zend];
      zk(fcoast:b) = zaux;
   end
   
   plot(xk,zk)
   LAT = num2str(yb(k,1)); title(LAT,'fontweight','bold');
   figure(2); hold off
%     figure(1); set(pl,'visible','off')
   zb(k,:) = zk;
end
close all
toc

figure; subplot(211); pcolor(xb,yb,zb1); shading flat; hold on;  caxis([-6000 10])
axis('equal'); axis([lonlim latlim])
contour(xb,yb,zb1,[0 0],'k'); contour(xb,yb,zb1,[-200 -1000],'w');
subplot(212); pcolor(xb,yb,zb); shading flat; hold on; caxis([-6000 10])
axis('equal'); axis([lonlim latlim])
contour(xb,yb,zb1,[0 0],'k'); contour(xb,yb,zb1,[-200 -1000],'w');
%  colorbar('southoutside')

STOP

%%% FANCY PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

zb2 = zb; f = find(zb1 > 50 & xb < -30); 

for k = f
   zb2(k) = zb1(k);
end

figure(1); set(gcf,'color','w')
pcolor(xb,yb,zb1); shading flat; hold on;  caxis([-6000 2450]); colorbar; colormap(cont_ocean)
axis('equal'); axis([lonlim latlim])
%  contour(xb,yb,zb1,[0 0],'k'); contour(xb,yb,zb1,[-200 -1000],'w');
xlabel('longitude [^\circ W]','fontweight','bold')
ylabel('latitude [^\circ S]','fontweight','bold')

figure(2); set(gcf,'color','w')
pcolor(xb,yb,zb2); shading flat; hold on; caxis([-6000 2450]); colorbar; colormap(cont_ocean)
axis('equal'); axis([lonlim latlim])
%  contour(xb,yb,zb1,[0 0],'k'); contour(xb,yb,zb1,[-200 -1000],'w');
xlabel('longitude [^\circ W]','fontweight','bold')
ylabel('latitude [^\circ S]','fontweight','bold')

figure(3)
set(gcf,'color','w')
surf(xb,yb,zb1); shading interp; hold on
%  contour(xb,yb,zb1,[0 0],'k'); contour(xb,yb,zb1,[-200 -1000],'w');
caxis([-6000 2450]); colormap(cont_ocean)
axis([lonlim latlim -6000 1000])
view([30 30])
daspect([1 1 3000]) 
h = light;
set(h,'Position',[10 -60 40000]);
material dull
%  text(-36,-17.7,'Banco Hot Spur','color','w','fontsize',10,'fontweight','bold')
%  text(-37.7,-16,'Banco Royal Charlote','color','w','fontsize',10,'fontweight','bold')
%  text(-39.5,-18.6,'Banco de Abrolhos','color','w','fontsize',10,'fontweight','bold')
%  text(-37.5,-21.5,'Cadeia Vitoria-Trindade','color','w','fontsize',10,'fontweight','bold')
%  text(-40,-13,100,'Salvador','color','w','fontsize',10,'fontweight','bold')
%  text(-42.5,-22.2,100,'C. Sao Tome','color','w','fontsize',10,'fontweight','bold')
%  text(-41,-20,100,'Vitoria','color','w','fontsize',10,'fontweight','bold')
%  xlabel('longitude [^\circ W]','fontweight','bold')
%  ylabel('latitude [^\circ S]','fontweight','bold')
%  zlabel('depth [m]','fontweight','bold')

figure(4)
set(gcf,'color','w')
surf(xb,yb,zb2); shading flat; hold on
%  contour(xb,yb,zb1,[0 0],'k'); contour(xb,yb,zb1,[-200 -1000],'w');
caxis([-6000 2450]); colormap(cont_ocean)
axis([lonlim latlim -6000 1000])
view([30 30])
daspect([1 1 3000]) 
h = light;
set(h,'Position',[10 -60 40000]);
material dull
%  text(-36,-17.7,'Banco Hot Spur','color','w','fontsize',10,'fontweight','bold')
%  text(-37.7,-16,'Banco Royal Charlote','color','w','fontsize',10,'fontweight','bold')
%  text(-39.5,-18.6,'Banco de Abrolhos','color','w','fontsize',10,'fontweight','bold')
%  text(-37.5,-21.5,'Cadeia Vitoria-Trindade','color','w','fontsize',10,'fontweight','bold')
%  text(-40,-13,100,'Salvador','color','w','fontsize',10,'fontweight','bold')
%  text(-42.5,-22.2,100,'C. Sao Tome','color','w','fontsize',10,'fontweight','bold')
%  text(-41,-20,100,'Vitoria','color','w','fontsize',10,'fontweight','bold')
%  xlabel('longitude [^\circ W]','fontweight','bold')
%  ylabel('latitude [^\circ S]','fontweight','bold')
%  zlabel('depth [m]','fontweight','bold')





























