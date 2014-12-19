%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCRIPT FOR VIZUALIZATION OF AVISO ALTIMETRY 
%      PRODUCT AND COMPUTING MEAN FIELDS
% - HORIZONTAL MAPS
% - Rafael Soutelino - rsoutelino@gmail.com
% - last update: August, 2009
%
% - Most parameters to be changed are at 'SETTINGS' section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%%% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataset = 'LESTE'

pathname = ['/Users/rsoutelino/rsoutelino/phd/data/aviso/',dataset,'/MATFILES/']
figdir = ['/Users/rsoutelino/rsoutelino/phd/data/figures/']

lonlim = [-46 -30]; latlim = [-25 -9];
m_proj('mercator','lat',latlim,'lon',lonlim);
[zb,xb,yb] = m_tbase([lonlim latlim]);
m_gshhs_l('save','costa.mat');
load costa.mat;
load /Users/rsoutelino/rsoutelino/mymatlab/lado/m_map/redblue;
% basic figures settings
PRINT = 'n'; % flag to print figures or not
DASPECT = [1 1 1];
DASPECT2 = [1 250 1];
s = 7; % scale size for vectors
dx = 1; dy = 1; % subsampling fator for vectors
ax = ['axis([min(xsec) max(xsec) -1000 0])'];

% END OF SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

YY = []; DAY = []; MON = []; U = []; V = []; Uab = []; Vab = []; Ucf = []; Vcf = [];

for k = 7:9
   matfile = ['absvel_200',num2str(k),'.mat']
   load([pathname,matfile]);
   yy = str2num(['200',num2str(k)]);
   yy = yy*ones(size(day));
   YY = [YY;yy];

   day = datestr(day,6);
   mon = day(:,1:2); mon = str2num(mon);
   day = day(:,4:5); day = str2num(day);

   DAY = [DAY;day]; 
   MON = [MON;mon];
   U = [U;u];
   V = [V;v];
   
   % making a temporal series of BC at two different locations
   fab1 = near(lat,-18,1);
   fab2 = near(lon,-36.9,1);
   fcf1 = near(lat,-23,1);
   fcf2 = near(lon,-40.3,1);

   uab = squeeze(u(:,fab2,fab1)); vab = squeeze(v(:,fab2,fab1));
   Uab = [Uab;uab]; Vab = [Vab;vab]; 
    
   ucf = squeeze(u(:,fcf2,fcf1)); vcf = squeeze(v(:,fcf2,fcf1));
   Ucf = [Ucf;ucf]; Vcf = [Vcf;vcf];
end

tt = datenum(YY,MON,DAY);

disp(' ')
disp('==> .................. COMPUTING TOTAL MEAN FIELD ....................')
disp(' ')

um = squeeze(nanmean(U(:,:,:))); um = um';
vm = squeeze(nanmean(V(:,:,:))); vm = vm';


%  lon2 = lon(1):0.1:lon(end);
%  lat2 = lat(1):0.1:lat(end);
[lon,lat] = meshgrid(lon,lat);
%  [lon2,lat2] = meshgrid(lon2,lat2);

%  um = griddata(lon,lat,um,lon2,lat2,'linear');
%  vm = griddata(lon,lat,vm,lon2,lat2,'linear');
mag = sqrt(um.^2 + vm.^2);

disp(' ')
disp('==> ............... MASKING CONTINENTAL SHELF .................')
disp(' ')

zb2 = griddata(xb,yb,zb,lon,lat);
f = find(zb2 > -1000);
um(f) = nan; vm(f) = nan; mag(f) = nan;

%  figure
%  pcolor(lon,lat,mag); shading flat; colorbar; hold on
%  colormap(redblue(32:end,:));
%  str = ['(1:dx:end,1:dy:end)'];
%  eval(['quiver(lon' str ',lat' str ',um' str ',vm' str ',3,''k'')']);
%  plot(ncst(:,1),ncst(:,2),'k');
%  plot(lon(fab1,fab2),lat(fab1,fab2),'*b','markersize',10)
%  plot(lon(fcf1,fcf2),lat(fcf1,fcf2),'*b','markersize',10)
%  contour(xb,yb,zb,[-100 -200 -500 -1000],'k');
%  caxis([4 20])
%  axis([lonlim latlim])
%  daspect([DASPECT])

figure
set(gcf,'color','w')
m_pcolor(lon,lat,mag); shading interp; colorbar; hold on
colormap(redblue(32:end,:));
str = ['(1:dx:end,1:dy:end)'];
eval(['m_quiver(lon' str ',lat' str ',um' str ',vm' str ',s,''k'',''linewidth'',1)']);
m_contourf(xb,yb,zb,[-100 -200 -500 -1000 -2000],'k');
m_usercoast('costa.mat','patch',[0 0 0])
m_plot(lon(fab1,fab2),lat(fab1,fab2),'ob','markersize',10,'markerfacecolor','y','markeredgecolor','k')
m_plot(lon(fcf1,fcf2),lat(fcf1,fcf2),'ob','markersize',10,'markerfacecolor','y','markeredgecolor','k')
%  set(h,'facecolor',[.5 .5 .5])
caxis([4 20])
%  axis([lonlim latlim])
%  daspect([DASPECT])
m_grid('box','fancy')
eval(['print -dpng ',figdir,'media_aviso_omar09.png'])

figure
set(gcf,'color','w')
subplot(211)
plot([tt(1) tt(end)],[0 0],'k'); hold on; grid on
plot(tt,Vab,'k','linewidth',1);
plot(tt,mean(Vab)*ones(size(tt)),'r','linewidth',2)
axis([tt(1) tt(end) -50 50])
%  xlabel('Janeiro de 2007 - Agosto de 2009','fontweight','bold')
ylabel('V [cm/s]','fontweight','bold')
title('Meridional Velocity off Abrolhos Bank - 18^\circ S','fontweight','bold')
gtext('Southward = 63%','fontweight','bold')
gtext('Mean Velocity = -1.8 cm/s','fontweight','bold','color','r')
datetick('x',4,'keeplimits');

subplot(212)
plot([tt(1) tt(end)],[0 0],'k'); hold on; grid on
plot(tt,Vcf,'k','linewidth',1);
plot(tt,mean(Vcf)*ones(size(tt)),'r','linewidth',2)
axis([tt(1) tt(end) -50 50])
xlabel('January 2007 - August 2009','fontweight','bold')
ylabel('V [cm/s]','fontweight','bold')
title('Meridional Velocity off Cabo de Sao Tome - 22^\circ S','fontweight','bold')
gtext('Southward = 89%','fontweight','bold')
gtext('Mean Velocity = -10 cm/s','fontweight','bold','color','r')
datetick('x',4,'keeplimits');
eval(['print -dpng ',figdir,'serie_aviso_omar09.png'])


















