%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      OBSERVED STREAM FUNCTION 
%      
%            Rafael Soutelino - March/2008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc;warning off

% configurations ----------------------------------------------------
pathname = ['/home/rafaelgs/doutorado/data/']
alis = 'n';
jan = 11;
cruz = 'leste2';
%  cruz = inputdlg('Pick the name of the dataset'); cruz = cell2mat(cruz);
CRUZ = 'Leste 2'
cc = 'n';
maperr = 'n';
PRINT = 'y';
lonlim = [-45 -30]; latlim = [-24 -9];
dx = 1; dy = 1; % subsampling vectors for plotting 
dec = 1; % sub-sampling factor for input adcp data into AOV
cs = 50; % depth of the dinamical boundary (continental shelf break)
lc = 1; 
E = 0.02;
s = 1; % vector size scale
lpsi = -2e5:2000:2e5;
%  load rainbow.mat
eval(['load ',pathname,cruz,'/grid_curv.mat']);
%-------------------------------------------------------------------
% setting map projection
m_proj('mercator','long',[lonlim],'lat',[latlim],'on');
[zb,xb,yb] = m_tbase([lonlim latlim]);
m_gshhs_l('save','costa.mat')

% loading data 
eval(['load ',pathname,cruz,'/adcp/',cruz,'_uv.mat;']);
eval(['load ',pathname,cruz,'/adcp/',cruz,'_xy.mat;']);
u = uv(:,1:2:end-1); v = uv(:,2:2:end); 
lon = xyt(1,:); lon = lon-360; lat = xyt(2,:); t = xyt(3,:);  

% setting the ploting level
p = input('Please chose a depth level (data from 32 m to 800 m, de 16/16 m): ')
n = near(zc,p,1);
u = u(n,:); v = v(n,:);
mod = sqrt(u.^2 + v.^2);

% removing NaNs and remaining spikes
thd = 1; 
f = find(u >= -thd & u <= thd & isnan(u)==0);
u=u(f);v=v(f);lon=lon(f);lat=lat(f);t=t(f);
f = find(v >= -thd & v <= thd & isnan(v)==0);
u=u(f);v=v(f);lon=lon(f);lat=lat(f);t=t(f);

%  if cruz == 'abrolhos2'
%     u = u([1:65 88:end-20]); 
%     v = v([1:65 88:end-20]);
%     lon = lon([1:65 88:end-20]);
%     lat = lat([1:65 88:end-20]);
%  end

% decimation
u = u(1:dec:end); v = v(1:dec:end);
lon = lon(1:dec:end); lat = lat(1:dec:end);

if alis == 's'
   u = weim(jan,'hann',u);
   v = weim(jan,'hann',v);
end

% preparing to the vectorial Objective Analysis
xgi = grid3(:,:,1); 
ygi = grid3(:,:,2);

[l1,l2] = size(xgi);

xg = reshape(xgi,1,l2*l1);
yg = reshape(ygi,1,l2*l1);

u = u*8.9993e-06; % converting m/s to degree/s 
v = v*8.9993e-06;

% calculating observed psi trough OA
[psiob] = vectoa(xg,yg,lon,lat,u,v,lc,E,0);

if maperr == 's'
  [ans,er] = scaloa(xg,yg,lon,lat,u,lc,E);
  er = 100*sqrt(er);
end

% using dirichlet boundary conditions
figure(1)
if p <= cs
   [c,h] = contour(xb,yb,zb,[-cs -cs]);
   pl = cs;
else 
   [c,h] = contour(xb,yb,zb,[-p -p]);
   pl = p;
end
xcont = get(h(1),'xdata'); xcont = xcont(find(not(isnan(xcont))));
ycont = get(h(1),'ydata'); ycont = ycont(find(not(isnan(ycont))));
%  xcont = weim(31,'hann',xcont);
%  ycont = weim(31,'hann',ycont);

xcont=xcont';
ycont=ycont';

close(1); 

psiob = psiob*(111120^2);
psiob = psiob-mean(psiob);

if cc == 's'
   xg2 = [xg xcont]; yg2 = [yg ycont];
   psiob = [psiob' zeros(size(xcont))];
   [psiob] = scaloa(xg,yg,xg2,yg2,psiob,lc,E);
else
   xg2 = xg'; yg2 = yg';
   psiob = psiob';
end

psiob = reshape(psiob,l1,l2);
[uo,vo] = psi2uv(xgi,ygi,psiob);

if maperr == 's'
   er = reshape(er,l1,l2);
end

u = u/8.9993e-06; % reconvertendo para plotar vetores brutos
v = v/8.9993e-06;

figure(1)
set(gcf,...
        'Color',[1 1 1],...
        'InvertHardcopy','on',...
        'PaperUnits','inches',...
        'Units','inches',...
        'PaperOrientation','portrait',...
        'PaperPosition',[0 0 8.5 11],...
        'PaperPositionMode','manual',...
        'PaperType','usletter',...
        'Position',[.2 .2 8 9],...
        'ShareColors','off',...
        'Clipping','on');
m_contourf(xgi,ygi,psiob,30);hold on; shading flat
cc = colorbar;
m_quiver(xgi(1:dx:end,1:dy:end),ygi(1:dx:end,1:dy:end),uo(1:dx:end,1:dy:end)*s,vo(1:dx:end,1:dy:end)*s,0,'k');
[c,h]=m_contourf(xb,yb,zb,[-pl -pl],'k');
set(h,'facecolor',[.7 .7 .7]);
[c,h]=m_contour(xb,yb,zb,[-pl -pl],'k');
set(h,'color',[.7 .7 .7]);
m_usercoast('costa.mat','patch',[0.542 0.422 0.000],'LineStyle','-');
m_grid('box','fancy','yaxislocation','left','xaxislocation','bottom');
m_quiver(lonlim(1)+0.5,latlim(2)-0.5,1*s,0*s,0,'k')
text = ['1 m s^{-1}'];
m_text(lonlim(1)+0.5,latlim(2)-0.7,'1 m s^{-1}','Fontsize',10,'FontWeight','Bold','Color','k');
tex = [num2str(p),' m'];
m_text(lonlim(1)+0.12,latlim(1)+0.21,tex,'fontsize',11,'FontWeight','bold','EdgeColor','k','BackgroundColor','w')
title(['\Psi Observada - ',CRUZ,' [m s^{-1}]'],'fontsize',10,'fontweight','bold');
cc = colorbar;
%  colormap(rainbow);
pos = get(cc,'Position');
set(cc,'Position',[pos(1)/1.15 pos(2)*4 pos(3)/3.3 pos(4)/3.5])
if PRINT == 'y'
   print(1,'-depsc',[pathname,'figures/psiob_',cruz,'_',num2str(p),'m']);
   eval(['!epstopdf ',pathname,'figures/psiob_',cruz,'_',num2str(p),'m.eps'])
end
drawnow


% figura com o mapa de erro da AO
if maperr == 'y'
figure(2)
set(gcf,...
        'Color',[1 1 1],...
        'InvertHardcopy','on',...
        'PaperUnits','inches',...
        'Units','inches',...
        'PaperOrientation','portrait',...
        'PaperPosition',[0 0 8.5 11],...
        'PaperPositionMode','manual',...
        'PaperType','usletter',...
        'Position',[.2 .2 8 9],...
        'ShareColors','off',...
        'Clipping','on');
m_contourf(xgi,ygi,er,50);hold on; shading flat
caxis([5 100])
[c,h]=m_contourf(xb,yb,zb,[-pl -pl],'k');
set(h,'facecolor',[.7 .7 .7]);
[c,h]=m_contour(xb,yb,zb,[-pl -pl],'k');
set(h,'color',[.7 .7 .7]);
m_usercoast('costa.mat','patch',[1 0.669 0.231],'LineStyle','-');
m_grid('box','fancy','yaxislocation','left','xaxislocation','bottom');
title('Erro de Interpolacao [%]','fontsize',10,'fontweight','bold');
tex = [num2str(p),' m'];
m_text(lonlim(1)+0.12,latlim(1)+0.21,tex,'fontsize',11,'FontWeight','bold','EdgeColor','k','BackgroundColor','w')
cc = colorbar;
pos = get(cc,'Position');
set(cc,'Position',[pos(1)/1.14 pos(2)*4.1 pos(3)/3 pos(4)/3])
if PRINT == 's'
   print(2,'-depsc',[pathname,'figures/errpsiob_',cruz,'_',num2str(p),'m']);
   eval(['!epstopdf ',pathname,'figuras/errpsiob_',cruz,'_',num2str(p),'m.eps'])
end
drawnow
end



