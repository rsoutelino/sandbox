clear all; close all; clc;

ans = menu('Which of those would you like to plot?','TEMPERATURE','SALINITY','DENSITY');
if     ans == 1;
   PROP = 'temp';
elseif ans == 2;
   PROP = 'salt';
else
   PROP = 'dens';
end

ans = menu('Which depth would you like to plot?','50 m','100 m','200 m','500 m','1000 m','1500 m','2000 m');
if     ans == 1;
   DEPTH = '50m';
   if PROP == 'temp'; scale = [20 29]; elseif PROP == 'salt'; scale = [36.5 37.5]; else; scale = [23.5 25.5]; end
elseif ans == 2;
   DEPTH = '100m';
   if PROP == 'temp'; scale = [17 26]; elseif PROP == 'salt'; scale = [35.6 37.5]; else; scale = [24.5 26.3]; end
elseif ans == 3;
   DEPTH = '200m';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if PROP == 'temp'; scale = [13 22]; elseif PROP == 'salt'; scale = [35.3 36.7]; else; scale = [25.5 26.5]; end
elseif ans == 4;
   DEPTH = '500m';
   if PROP == 'temp'; scale = [6 10.5]; elseif PROP == 'salt'; scale = [34.5 34.9]; else; scale = [26.76 27.17]; end
elseif ans == 5;
   DEPTH = '1000m';
   if PROP == 'temp'; scale = [3.5 4.3]; elseif PROP == 'salt'; scale = [34.35 34.65]; else; scale = [27.26 27.52]; end
elseif ans == 6;
   DEPTH = '1500m';
   if PROP == 'temp'; scale = [3.6 4.4]; elseif PROP == 'salt'; scale = [34.7 35]; else; scale = [27.56 27.78]; end
else
   DEPTH = '2000m';
   if PROP == 'temp'; scale = [3 3.8]; elseif PROP == 'salt'; scale = [34.9 35.01]; else; scale = [27.775 27.84]; end
end

eval(['open figures/leste1_',PROP,'_',DEPTH,'.fig']);
caxis([scale])
eval(['open figures/leste2_',PROP,'_',DEPTH,'.fig']);
caxis([scale])
eval(['open figures/abrolhos2_',PROP,'_',DEPTH,'.fig']);
caxis([scale])
eval(['open figures/abrolhos1_',PROP,'_',DEPTH,'.fig']);
caxis([scale])
eval(['open figures/proab_',PROP,'_',DEPTH,'.fig']);
caxis([scale])
