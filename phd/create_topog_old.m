%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ETOPO 1 - based topography
% Creates sea bottom grids for ROMS runs - PhD work
% Eastern brazilian coast
% Creates real topography and semi-idealized topography removing
% seamounts, banks and other steep topographic features
% Rafael Soutelino, June, 2010
% rsoutelino@gmail.com	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc

%%% BASIC SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loading info for ETOPO 1 database: 
url = ['ftp://ftp.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/bedrock/grid_registered/netcdf/'];
etopo1file = ['/Users/rsoutelino/rsoutelino/misc/etopo1_1981x1801_atl_sul.xyz'];
filename = ['etopo1_1981x1801_atl_sul']; l1 = 1801; l2 = 1981;
newfile = 'n';
mycolorbar = ['cc = colorbar(''west''); pos = get(cc,''position''); set(cc,''position'',[pos(1) pos(2) pos(3)/3 pos(4)])'];

% declaring geographical limits
lonlim = [-45 -25]; latlim = [-22.8 -10]; m_proj('mercator','lon',lonlim,'lat',latlim);
load costa.mat; %  m_gshhs_f('save','costa.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loading variables
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
   f1 = find(yb(:,1) >= latlim(1) & yb(:,1) <= latlim(2));
   f2 = find(xb(1,:) >= lonlim(1) & xb(1,:) <= lonlim(2));
   xb = xb(f1,f2); yb = yb(f1,f2); zb = zb(f1,f2); 
end

%  figure; set(gcf,'color','w');
%  m_pcolor(xb,yb,zb);shading flat; hold on                    
%  m_usercoast('costa.mat','patch',[.9 .9 .9]); m_grid; eval(mycolorbar);
%  m_contour(xb,yb,zb,[-100 -200 -500 -1000 -1500 -2000 -3000 -4000],'k');

%%% REMOVING FEATURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% trying to smooth it before computations
zb = smoo2(zb,-9999,9,0.5);

% removing continent topography
f = find(zb > 10); zb(f) = 10;
zb2 = zb; zb3 = zb;

% removing semounts and VT Ridge generically --------------------------------------------
[a,b] = size(zb);
figure(1); pcolor(xb,yb,zb); shading flat; hold on; 
axis('equal'); axis([lonlim latlim])
contour(xb,yb,zb,[0 0],'k');

% scanning parallels ----------------------------------------
% general location of seamounts
int = find(yb(:,1) >= -23 & yb(:,1) <= -15);
tic
for k = int(1):int(end);
   figure(1); pl = plot(xb(k,:),yb(k,:),'w');

   x1 = xb(k,:); z1 = zb(k,:);

   figure(2)
   plot(x1,z1); hold on; axis([lonlim -6500 100])
   plot([lonlim],[0 0],'g--'); grid on;

   f = find(z1 <= -2000); zfit = z1(f); xfit = x1(f); xdeep = x1(f(1)-1:end); zdeep = z1(f(1)-1:end);
   xshallow = x1(1:f(1)-2); zshallow = z1(1:f(1)-2);
   pol = polyfit(xfit,zfit,4); 
   zfitted = polyval(pol,xfit); 
   zdeepfitted = interp1(xfit,zfitted,xdeep);
   plot(xdeep,zdeepfitted,'r');
   LAT = num2str(yb(k,1)); title(LAT,'fontweight','bold');
   z2 = [zshallow zdeepfitted];
   f = find(z1-z2 > 200); zb2(k,f) = z2(f)-200;
%     pause
   figure(2); hold off
   figure(1); set(pl,'visible','off')
   plot(xb(k,f),yb(k,f),'w.'); 
end
close all
toc

%  figure(2); pcolor(xb,yb,zb2); shading flat; hold on; 
%  axis('equal'); axis([lonlim latlim])
%  contour(xb,yb,zb,[0 0],'k');

zb3 = zb2;
figure(1); pcolor(xb,yb,zb2); shading flat; hold on; 
axis('equal'); axis([lonlim latlim])
contour(xb,yb,zb,[0 0],'k');
%  


%  % scanning meridians -------------------------------------------
% general location of seamounts
int = find(xb(1,:) >= -37.1 & xb(1,:) <= -27);
tic
for k = int(1):int(end);
   figure(1); pl = plot(xb(:,k),yb(:,k),'w');

   y1 = yb(:,k); z1 = zb(:,k);

   figure(2)
   plot(y1,z1); hold on; axis([latlim -6500 100])
   plot([latlim],[0 0],'g--'); grid on;

   f = find(z1 <= -2000); 
   if isempty(f) == 0
   zfit = z1(f); yfit = y1(f); ydeep = y1(f(1):end); zdeep = z1(f(1):end);
   if f(1) == 1
      yshallow = []; zshallow = [];;
   else
      yshallow = y1(1:f(1)-1); zshallow = z1(1:f(1)-1);
   end
   pol = polyfit(yfit,zfit,4); 
   zfitted = polyval(pol,yfit); 
   zdeepfitted = interp1(yfit,zfitted,ydeep);
   plot(ydeep,zdeepfitted,'r');
   LON = num2str(xb(1,k)); title(LON,'fontweight','bold');
   z2 = [zshallow ; zdeepfitted];
   f = find(z1-z2 > 200); zb3(f,k) = z2(f)-200;
   end
%     pause
   figure(2); hold off
   figure(1); set(pl,'visible','off')
   plot(xb(f,k),yb(f,k),'w.'); 
end
close all
toc


stop

% removing banks --------------------------------------------------------------------

































