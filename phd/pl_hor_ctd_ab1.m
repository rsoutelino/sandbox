%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT HORIZONTAL HYDROGRAPHIC (T,S,rho) MAPS
% CALCULATES GEOSTROPHIC STREAMFUNCTION RELATIVE TO LEVEL OF NO MOTION
% Abrolhos1 - Sep/2004
% Rafael Soutelino - Jun, 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc;warning off

%%% BASIC SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CRUZ = 'abrolhos1'
PATHNAME = ['/Users/rsoutelino/rsoutelino/phd/data/abrolhos1/ctd'];
DATE = 'Sep, 2004';
eval(['load ',PATHNAME,'/costa.mat']); % coastline
eval(['load ',PATHNAME,'/etopo2_leste.mat']); % topography
zb2 = smoo2(zb,-9999,5,1);
eval(['load ',PATHNAME,'/posicoes.dat']); pb = posicoes; % positions file
eval(['load ',PATHNAME,'/seagrid.dat'])
eval(['load ',PATHNAME(1:28),'/mymatlab/lado/m_map/redblue']);
nest = pb(:,6); lat = -[pb(:,1) + pb(:,2)/60]; lon = -[pb(:,3) + pb(:,4)/60];
prof = pb(:,5); clear posicoes pb; % when available

pmin = 15; % minimum depth
nr = 1000; % level of no motion
lonlim = [-41 -32]; latlim = [-22 -9];
f0 = sw_f(mean(lat));    % coriolis parameter
fig1 = 'y';       % key for horizontal T,S plots - linear interpolation
fig2 = 'y';       % key for streamfunction plots - OA interpolation
fig3 = 'n';
mask = 100;       % mask for continental shelf
maxerr = 50;      % maximum allowed OA interpolation error (in %)
cc = 'n';                 % condicoes de contorno
n = input('Pick depth to plot:   ');
PRINT = 'n';

lT = 2:29; lS = 34:0.1:37.5; lD = 23.5:0.1:28; % color scales

xx = [7550:7610];  

%%% READING HYDROGRAPHIC DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lx=length(xx); 

c = 1;
for k=1:length(xx);
    XX = num2str(xx(k));
    eval(['load ',PATHNAME,'/abrolhos1_ctd_wb',XX,'.asc']);
    dados = eval(['abrolhos1_ctd_wb',XX]);
    eval(['clear abrolhos1_ctd',XX]);

    p = dados(:,1);          
    t = dados(:,2);          
    s = dados(:,5);          

    if fig3 == 'y'       
       pmax = nr;
    else   
       pmax = n;  
    end

    f = find(p >= pmin & p <= pmax);

  if max(p) >= pmax
     T(1:length(f),c) = t(f);
     S(1:length(f),c) = s(f);
     lats(c) = lat(find(nest==xx(k)));
     lons(c) = lon(find(nest==xx(k)));
     c = c+1;
     clear s t p XX dados 
  else
     disp(['This station does not reach required depth: ',num2str(pmax),' m'])
  end
end

clear c k f pos posicoes
clc

p = pmin:pmax; p = p';
D = sw_dens0(S,T); D = D - 1000*ones(size(D));

%%% COMPUTING STREAMFUNCTION POINT BY POINT, NIVEL BY NIVEL %%%%%%%%%%%%%%%%%%%%%%%

gpan = sw_gpan(S,T,p);
if fig3 == 'y'
   gpan = gpan - ones(size(p)) * gpan(nr-pmin,:);
end
psig = -[gpan./f0];

j    = n - pmin;
T    = T(j,:);
S    = S(j,:);
D    = D(j,:);
psig = psig(j,:);

% removing average to compute gradients later on
psig = psig - mean(psig);

%%% INTERPOLATION TO HORIZONTAL GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xg  = seagrid(:,9); yg = seagrid(:,10); lj = max(seagrid(:,1)); li = max(seagrid(:,2));         
xgi = reshape(xg,li,lj); ygi = reshape(yg,li,lj);

disp(' ')
disp('Interpolation of Topography to mask data .........')
disp(' ')
zbi = griddata(xb,yb,zb2,xgi,ygi);

if fig1 == 'y'
% linear interpolation
disp(' ')
disp('Linear Interpolation of T, S, Rho, Psi .........')
disp(' ')
Ti      = griddata(lons,lats,T,xgi,ygi);    disp('Temp')
Si      = griddata(lons,lats,S,xgi,ygi);    disp('Salt')
Di      = griddata(lons,lats,D,xgi,ygi);    disp('Dens')
psigi   = griddata(lons,lats,psig,xgi,ygi); disp('Psig')
[ui,vi] = psi2uv(xgi,ygi,psigi);

% masking filds for topography
Ti(find(zbi > -n)) = NaN; Si(find(zbi > -n)) = NaN; Di(find(zbi > -n)) = NaN;
psigi(find(zbi > -n)) = NaN; ui(find(zbi > -n)) = NaN; vi(find(zbi > -n)) = NaN;

% plotting raw fields ====================================================================
%  disp(' ')
%  disp('Plotting raw fields .........')
%  disp(' ')
%  figure; set(gcf,'color','w');
%  pcolor(xgi,ygi,Ti); shading flat; hold on;  
%  contour(xb,yb,zb,[-200 -1000],'k'); plot(ncst(:,1),ncst(:,2),'k');
%  plot(lons,lats,'k.',lons,lats,'ow')
%  caxis([min(min(Ti)) max(max(Ti))]);
%  axis('equal'); axis([lonlim latlim]); colorbar('east')
%  tit = ['Temp. Linear Int. - ',num2str(n),'m'];
%  text(lonlim(1)+0.5,latlim(2)-1,tit,'fontweight','bold')
%  text(lonlim(1)+0.5,latlim(2)-2,DATE,'fontweight','bold')
%  
%  figure; set(gcf,'color','w');
%  pcolor(xgi,ygi,Si); shading flat; hold on;  
%  contour(xb,yb,zb,[-200 -1000],'k'); plot(ncst(:,1),ncst(:,2),'k');
%  plot(lons,lats,'k.',lons,lats,'ow')
%  caxis([min(min(Si)) max(max(Si))]);
%  axis('equal'); axis([lonlim latlim]); colorbar('east')
%  tit = ['Sal. Linear Int. - ',num2str(n),'m'];
%  text(lonlim(1)+0.5,latlim(2)-1,tit,'fontweight','bold')
%  text(lonlim(1)+0.5,latlim(2)-2,DATE,'fontweight','bold')
%  
%  figure; set(gcf,'color','w');
%  pcolor(xgi,ygi,Di); shading flat; hold on;  
%  contour(xb,yb,zb,[-200 -1000],'k'); plot(ncst(:,1),ncst(:,2),'k');
%  plot(lons,lats,'k.',lons,lats,'ow')
%  caxis([min(min(Di)) max(max(Di))]);
%  axis('equal'); axis([lonlim latlim]); colorbar('east')
%  tit = ['\sigma_\theta. Linear Int. - ',num2str(n),'m'];
%  text(lonlim(1)+0.5,latlim(2)-1,tit,'fontweight','bold')
%  text(lonlim(1)+0.5,latlim(2)-2,DATE,'fontweight','bold')
end

if fig2 == 'y'
% objective analysis
disp(' ')
disp('OA Interpolation of T, S, Psi .......')
disp(' ')

% removing averages to perform OA
mT= mean(T); mS = mean(S); mD = mean(D);

xg=xg';yg=yg';
lc = 1.5;       
E = 0.02;           
[To,er]    = scaloa(xg,yg,lons,lats,T-mT,lc,E);    disp('Temp')
[So,er]    = scaloa(xg,yg,lons,lats,S-mS,lc,E);    disp('Salt')
[Do,er]    = scaloa(xg,yg,lons,lats,D-mD,lc,E);    disp('Dens')
[psigo,er] = scaloa(xg,yg,lons,lats,psig,lc,E); disp('Psig')
er = 100*sqrt(er);

% adding removed average
To = To + mT; So = So + mS; Do = Do + mD;

To = reshape(To,li,lj); So = reshape(So,li,lj);
Do = reshape(Do,li,lj); psigo = reshape(psigo,li,lj); er = reshape(er,li,lj);
psigo = psigo - mean(mean(psigo)); % removing streamfunction average for velocity computing
[uo,vo] = psi2uv(xgi,ygi,psigo);

% masking fields for topography
To(find(zbi > -n)) = NaN; So(find(zbi > -n)) = NaN; Do(find(zbi > -n)) = NaN;
psigo(find(zbi > -n)) = NaN; uo(find(zbi > -n)) = NaN; vo(find(zbi > -n)) = NaN;

% masking fields for OA interpolation error %
To(find(er > maxerr)) = NaN; So(find(er > maxerr)) = NaN; Do(find(er > maxerr)) = NaN;
psigo(find(er > maxerr)) = NaN; uo(find(er > maxerr)) = NaN; vo(find(er > maxerr)) = NaN;


% plotting OA fields ===========================================================================
disp(' ')
disp('Plotting OA Fields ...................')
disp(' ')

figure; set(gcf,'color','w');
pcolor(xgi,ygi,To); shading flat; hold on;  
contour(xb,yb,zb,[-200 -1000],'k'); plot(ncst(:,1),ncst(:,2),'k');
plot(lons,lats,'k.',lons,lats,'ow')
caxis([min(min(Ti)) max(max(Ti))]);
axis('equal'); axis([lonlim latlim]); colorbar('east')
tit = ['Temp. OA - ',num2str(n),'m'];
text(lonlim(1)+0.5,latlim(2)-1,tit,'fontweight','bold')
text(lonlim(1)+0.5,latlim(2)-2,DATE,'fontweight','bold')
saveas(gcf,[CRUZ,'_temp_',num2str(n),'m'],'fig');

figure; set(gcf,'color','w');
pcolor(xgi,ygi,So); shading flat; hold on;  
contour(xb,yb,zb,[-200 -1000],'k'); plot(ncst(:,1),ncst(:,2),'k');
plot(lons,lats,'k.',lons,lats,'ow')
caxis([min(min(Si)) max(max(Si))]);
axis('equal'); axis([lonlim latlim]); colorbar('east')
tit = ['Sal. OA. - ',num2str(n),'m'];
text(lonlim(1)+0.5,latlim(2)-1,tit,'fontweight','bold')
text(lonlim(1)+0.5,latlim(2)-2,DATE,'fontweight','bold')
saveas(gcf,[CRUZ,'_salt_',num2str(n),'m'],'fig');

figure; set(gcf,'color','w');
pcolor(xgi,ygi,Do); shading flat; hold on;  
contour(xb,yb,zb,[-200 -1000],'k'); plot(ncst(:,1),ncst(:,2),'k');
plot(lons,lats,'k.',lons,lats,'ow')
caxis([min(min(Di)) max(max(Di))]);
axis('equal'); axis([lonlim latlim]); colorbar('east')
tit = ['\sigma_\theta. OA. - ',num2str(n),'m'];
text(lonlim(1)+0.5,latlim(2)-1,tit,'fontweight','bold')
text(lonlim(1)+0.5,latlim(2)-2,DATE,'fontweight','bold')
saveas(gcf,[CRUZ,'_dens_',num2str(n),'m'],'fig');
end

