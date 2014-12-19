%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT SECTIONAL HYDROGRAPHIC (T,S,rho) MAPS
% CALCULATES GEOSTROPHIC VELOCITIES RELATIVE TO LEVEL OF NO MOTION
% Abrolhos 1 - Set/2004
% Rafael Soutelino - Jun, 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc

%%% BASIC SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PATHNAME = ['/Users/rsoutelino/rsoutelino/phd/data/abrolhos1/ctd'];
DATE = 'Set, 2004';
load costa.mat; % coastline
load etopo2_leste.mat; % topography
eval(['load ',PATHNAME,'/posicoes.dat']); pb = posicoes; % positions file
eval(['load ',PATHNAME(1:28),'/mymatlab/lado/m_map/redblue']);
nest = pb(:,6); lat = -[pb(:,1) + pb(:,2)/60]; lon = -[pb(:,3) + pb(:,4)/60]; 
prof = pb(:,5); clear posicoes pb; % when available

pmin = 15; % minimum depth
pplot = 1200; % max depth for plotting
nr = 1000; % level of no motion
lonlim = [-44 -31]; latlim = [-22 -7];

lT = 2:29; lS = 34:0.1:37.5; lD = 23.5:0.1:28; l = ['lTlSlD']; lV = -1:0.01:1; lV2 = [-0.5 0.5];
pl = [1]; % properties to be plotted (1 = T, 2 = S, 3 = sigma_theta)
EXTRAP = 'y';
PRINT = 'n';

figure(100);
plot(ncst(:,1),ncst(:,2),'k'); hold on; p1 = plot(lon,lat,'.r'); axis([lonlim latlim]); axis('equal');
contour(xb,yb,zb,[-200 -1000],'k')
text(lon(find(nest==7559))+0.2,lat(find(nest==7559)),'1'); 
text(lon(find(nest==7560))+0.2,lat(find(nest==7560)),'2');
text(lon(find(nest==7596))+0.2,lat(find(nest==7596)),'3');
text(lon(find(nest==7610))+0.2,lat(find(nest==7610)),'4');  

RR = menu('Pick Transect','# 1','# 2','# 3','# 4');
set(p1,'visible','off')

if     RR==1;  % perpendicular to the coast (W --> E)
xx = [7550:7559];
elseif RR==2;  % perpendicular to the coast (W --> E)
xx = [7573:-1:7560];
elseif RR==3;  % perpendicular to the coast (W --> E)
xx = [7574:7596];
elseif RR==4;  % perpendicular to the coast (W --> E)
xx = [7597:7610];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% READING HYDROGRAPHIC DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lx=length(xx); clear zb

c=1;
for k=1:lx;
  XX = num2str(xx(k));
  eval(['load ',PATHNAME,'/abrolhos1_ctd',num2str(RR),'_wb',XX,'.asc']);
  dados = eval(['abrolhos1_ctd',num2str(RR),'_wb',XX]);
  eval(['clear abrolhos1_ctd',num2str(RR),'_wb',XX]);

  p = dados(:,1); 
  t = dados(:,2);
  s = dados(:,5);
  f = find(p >= pmin & p <= pplot); % avoiding different minimum depths
  if max(p) >= 10
     T(1:length(f),c) = t(f);
     S(1:length(f),c) = s(f);
%       zb(c) = p(length(p)); % maximum depth is the last profiled depth
     zb(c) = prof(find(nest==xx(k))); % maximum depth is local depth (when available)
     lats(c) = lat(find(nest==xx(k)));
     lons(c) = lon(find(nest==xx(k)));
     c=c+1;
     clear s t p XX dados latg latm long lonm 
  else
     disp('This station is shallower than 10m')
  end
end

% plotting the transect on map
figure(100); plot(lons,lats,'r.');

f = find(T==0); T(f) = NaN; S(f) = NaN; % masking data
[lz,lx]=size(T); z = 15:lz+14; z = z'; zmax = max(z); % depth matrix
[dst,ang] = sw_dist(lats,lons,'km'); x = [0 cumsum(dst)]; xmax = max(x); % distance matrix 

if RR == 0
   ybi = linspace(min(lats),max(lats),100); zbi = interp1(lats,zb,ybi,'cubic'); % topography
   ybi = [min(ybi) ybi];
else
   xbi = linspace(min(lons),max(lons),100); zbi = interp1(lons,zb,xbi,'cubic'); % topography
   xbi = [min(xbi) xbi];
end
 zbi = [(max(zbi)) zbi];

clear ang f dst c k % making memory available

%%% REID-MANTYLLA EXTRAPOLATION TECHNIQUE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ti=T; Si=S; xi=x;

if EXTRAP == 'y';
disp('Lateral extrapolation of the fields towards to ocean floor ')
decT = 2; decS = 2; % gradient decay rate towards to the ocean floor (laterally)

for i = 1:length(z)   
   f = find(isnan(Ti(i,:))==1); % findind masked points on each line
   if isempty(f)==1; clc
   else
     ff = find(diff(f) > 1);
     if isempty(ff)==1  % searching for seamounts

        % extrapolation WITHOUT seamounts
        if f(end)+2 > lx % checking if there is any available value to calculate gradient
            for k=fliplr(f)
              Tult = Ti(i,k+1); Sult = Si(i,k+1);
              Ti(i,k) = Tult; Si(i,k) = Sult;
            end
        else
            dt = Ti(i,f(end)+2) - Ti(i,f(end)+1); % gradient   
            ds = Si(i,f(end)+2) - Si(i,f(end)+1);
            for k=fliplr(f)
              dt = dt/decT; ds = ds/decS;
              Tult = Ti(i,k+1); Sult = Si(i,k+1);                      
              Tex = Tult - dt; Sex = Sult - ds;        
              Ti(i,k) = Tex; Si(i,k) = Sex;              
            end
        end

     else

         % extrapolation WITH seamounts
         f1 = f(1:ff); % NaNs on the continental slope
         if isnan(Ti(i,ff+2)) == 1 % checking if there is any available value to calculate gradient
            for k=fliplr(f1)
              Tult = Ti(i,k+1); Sult = Si(i,k+1);
              Ti(i,k) = Tult; Si(i,k) = Sult;
            end
         else             
            dt = Ti(i,f1(end)+2) - Ti(i,f1(end)+1); % gradient  
            ds = Si(i,f1(end)+2) - Si(i,f1(end)+1);
            for k=fliplr(f1)
              dt = dt/decT; ds = ds/decS;
              Tult = Ti(i,k+1); Sult = Si(i,k+1);                      
              Tex = Tult - dt; Sex = Sult - ds;        
              Ti(i,k) = Tex; Si(i,k) = Sex;              
            end     
         end

         % left side of the seamount
         f2 = f(ff+1:end); % NaNs do monte
         
         for k=f2(1:ceil(length(f2)/2))
           Tult = Ti(i,k-1); Sult = Si(i,k-1);                             
           Ti(i,k) = Tult; Si(i,k) = Sult;              
         end   
         % right side of the seamount
         for k=fliplr(f2(ceil(length(f2)/2):end))
           Tult = Ti(i,k+1); Sult = Si(i,k+1);                              
           Ti(i,k) = Tult; Si(i,k) = Sult;              
         end  
               
      end
   end
end
end
clear i k f ff f1 f2 ds dt % making memory available

Di = sw_dens0(Si,Ti)-1000; % sigma_theta computation

%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P = ['TiSiDi'];
if RR == 0
   LON = num2str(abs(round(mean(lons))));
   t1 = ['TEMP: ',LON,'^\circ W'];
   t2 = ['SAL: ',LON,'^\circ W'];
   t3 = ['\sigma_\theta: ',LON,'^\circ W'];
else
   LAT = num2str(abs(round(mean(lats))));
   t1 = ['TEMP: ',LAT,'^\circ S'];
   t2 = ['SAL: ',LAT,'^\circ S'];
   t3 = ['\sigma_\theta: ',LAT,'^\circ S'];
end

tit = ['t1t2t3'];

c=0;
for k = pl
  figure(k)
  set(k,'Color','w'); 
  hold on;
  if RR == 0
     eval(['[c1,h1] = contourf(lats,z,',P(c+k:c+k+1),',',l(c+k:c+k+1),');']); 
     hold on; set(gca,'ydir','reverse');
     eval(['caxis([',l(c+k:c+k+1),'(1)  ',l(c+k:c+k+1),'(end)])'])
     axis([min(lats) max(lats) 0 pplot])
%    clabel(c1,h1,'labelspacing',500)
     fill(ybi,zbi,[.8 .8 .8])
     plot(lats,10*ones(size(lats)),'kv')
     eval(['text(min(lats)+0.01,pplot-120,',tit(c+k:c+k+1),',''fontweight'',''bold'')']);
     text(min(lats)+0.01,pplot-50,DATE,'fontweight','bold');
     xlabel('Latitude')
     ylabel('Depth')
  else
     eval(['[c1,h1] = contourf(lons,z,',P(c+k:c+k+1),',',l(c+k:c+k+1),');']); 
     hold on; set(gca,'ydir','reverse'); 
     eval(['caxis([',l(c+k:c+k+1),'(1)  ',l(c+k:c+k+1),'(end)])'])
     axis([min(lons) max(lons) 0 pplot])
%    clabel(c1,h1,'labelspacing',500)
     fill(xbi,zbi,[.8 .8 .8])
     plot(lons,10*ones(size(lons)),'kv')
     eval(['text(min(lons)+0.01,pplot-120,',tit(c+k:c+k+1),',''fontweight'',''bold'')']);
     text(min(lons)+0.01,pplot-50,DATE,'fontweight','bold');
     xlabel('Longitude')
     ylabel('Depth')
  end
 
  cc = colorbar('eastoutside');
  pbaspect([1 0.7 1]);
  hold off;

c=c+1; 
end
clear k c

%%% COMPUTING RELATIVE GEOSTROPHIC VELOCITIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if EXTRAP == 'y'
f0 = sw_f(mean(lats)); % coriolis parameter
gpan = sw_gpan(Si,Ti,z); % geopotential anomaly
dgpan = (diff(gpan'))'; % dGPAN
dx = diff(xi); dx = dx*1000; % dX 
dx = ones(size(z)) * dx; 

vg = -(1/f0) .* (dgpan./dx); % gvel equation
vg = vg - ones(size(z))*vg(nr,:);

if RR == 0
   for j=1:length(lats)-1
     latsv(j) = (lats(j) + lats(j+1))./2;
   end
titvg = ['G.VEL: ',LON,'^\circ W'];
else
   for j=1:length(lons)-1
     lonsv(j) = (lons(j)+lons(j+1))./2;
   end
titvg = ['G.VEL: ',LAT,'^\circ S'];
end

figure(4)
  set(4,'Color','w'); 
if RR == 0
  contourf(latsv,z,vg,lV,'k'); shading flat ; hold on
  caxis(lV2);
  colormap(flipud(redblue));
  set(gca,'ydir','reverse')
  axis([min(lats) max(lats) 0 pplot]) 
  fill(ybi,zbi,[.8 .8 .8])
  plot(lats,10*ones(size(lats)),'kv')
  plot(latsv,10*ones(size(latsv)),'k*')
  cc = colorbar('eastoutside');
  pbaspect([1 0.7 1])
  hold off;
  text(min(lats)+0.01,pplot-120,titvg,'fontweight','bold');
  text(min(lats)+0.01,pplot-50,DATE,'fontweight','bold');
  xlabel('Latitude')
  ylabel('Depth')  
else
  contourf(lonsv,z,vg,lV,'k'); shading flat ; hold on
  caxis(lV2);
  colormap(flipud(redblue));
  set(gca,'ydir','reverse')
  axis([min(lons) max(lons) 0 pplot]) 
  fill(xbi,zbi,[.8 .8 .8])
  plot(lons,10*ones(size(lons)),'kv')
  plot(lonsv,10*ones(size(lonsv)),'k*')
  cc = colorbar('eastoutside');
  pbaspect([1 0.7 1])
  hold off;
  text(min(lons)+0.01,pplot-120,titvg,'fontweight','bold');
  text(min(lons)+0.01,pplot-50,DATE,'fontweight','bold');
  xlabel('Longitude')
  ylabel('Depth')  
end

end
