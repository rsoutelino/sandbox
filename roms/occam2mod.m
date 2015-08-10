% lev2mod.m
%       reads in ASCII-data file for LEVITUS Summer Ddata and 
%       returns to MOD file.
%

%   1/12/04     Hyun-Sook Kim       
%---------------------------------------------------------------

%  z=[0.,    10.,   20.,   30.,   50.,   75.,   100.,  125.,  150.,  200., ...
%          250.,  300.,  400.,  500.,  600.,  700.,  800.,  900.,  1000., 1100., ...
%          1200., 1300., 1400., 1500.,1750., 2000., 2500., 3000., 3500., 4000., ...
%          4500., 5000., 5500];            % seasonal/annual mean Levitus hydro
%  %==== loading data ===========================================================
%  %tdat=load('Tsummer_0.25d_MB_lrg.dat');
%  %sdat=load('Ssummer_0.25d_MB_lrg.dat');
%  %=============================================================================
%  disp('      ****        Levitus ASCII data to MOD file')
%  tname=uigetfile('*.dat','   Pick a Levitus Temperature ASCII Data File');
%  tdat=load(tname);
%      k=findstr((tname),'T');
%      sname=tname;
%      sname(k)='S';
%  sdat=load(sname);

clear all; close all; clc

%%%% BASIC SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xx = 1:112; lx = length(xx);  % number of the CTD stations
datadir = ['/home/rafaelgs/mestrado/dados/leste2/ctd/filtrados/'];
eval(['load ',datadir,'posicoes_leste2.dat']); 
                              % lon, lat positions of the stations
pmax = 2000                   % maximum depth for hydrographic data
pmin = 10                     % minimum depth for hydrographic data             
depthmin = 100                % minimum total depth of interest       
dz = 10                        % depth interval between samples                       
z=[0.,    10.,   20.,   30.,   50.,   75.,   100.,  125.,  150.,  200., ...
        250.,  300.,  400.,  500.,  600.,  700.,  800.,  900.,  1000., 1100., ...
        1200., 1300., 1400., 1500.,1750., 2000., 2500., 3000., 3500., 4000., ...
        4500., 5000., 5500];            % seasonal/annual mean Levitus hydro
stime = 14833

%%% CTD data reading 

pb = posicoes_leste2;
profest=pb(:,6);
nest=pb(:,1);
latg=pb(:,2); latm=pb(:,3); lats=latg+latm/60;
long=pb(:,4); lonm=pb(:,5); lons=long+lonm/60;

clear  latg latm long lonm pb posicoes_leste2

c=1;
for k=1:lx;
  XX=num2str(xx(k));

  eval(['load ',datadir,'lesteII_ctd',XX,'.dat']);
  dados=eval(['lesteII_ctd',XX]);
  eval(['clear lesteII_ctd',XX]);

  p=dados(:,1); 
  t=dados(:,2);
  s=dados(:,3);

  F = [];
  for j = z
     f = find(p == j);
     F = [F f]; 
  end

  if max(p) >= depthmin
     tdat(1:length(F),c) = t(F);
     sdat(1:length(F),c) = s(F);
     lat(c) = -lats(find(nest==xx(k)));
     lon(c) = -lons(find(nest==xx(k)));
     c=c+1;
  else
     txt = ['Profile ',XX,' does not reach ',num2str(depthmin),' m'];
     disp(txt)
  end
end

clear  p XX dados latg latm long lonm posicoes_leste2 f k  t s lons lats 

f=find(tdat==0);

tdat(f)=NaN;
sdat(f)=NaN;
stop
%%% adding the surface level by repeating the 10 m field
tdat = [tdat(1,:); tdat]; sdat = [sdat(1,:); sdat];

tdat=tdat'; sdat=sdat'; lon=lon'; lat=lat';

%  %---
%  disp(' ')
%  txt=sprintf('   input T&S data filenames are : %s, %s',tname,sname);
%  disp(txt)
%  disp(' ')


k=find(~isnan(tdat(:,1)));
tdat=tdat(k,:);
sdat=sdat(k,:);
lat=lat(k); lon=lon(k);
clear k

[R,C]=size(tdat);

zdat=repmat(z,R,1);

%--- parameters for MOD-header
mxlon=max(lon); mnlon=min(lon);
mxlat=max(lat); mnlat=min(lat);

disp(' ')
disp('  ... Parameter Info for MOD-Header ...')

txt=sprintf(' [min max]-latitude = [%6.2f %6.2f]: ',mnlat,mxlat);
disp(txt)
txt=sprintf(' [min max]-longitude = [%6.2f %6.2f]: ',mnlon,mxlon);
disp(txt)

txt=sprintf('   A total # of Stations: %g',R);
disp(txt)

disp(' ')

hdr=whdr;           % writing MOD-header in an interactive mode
hdr=char(hdr);


%--- parameters for htype
htype=repmat(' ''CTD: z t s''',R,1);

%--- parameters for hinfo
%     hinfo   Cast header information matrix
%              hinfo(nsta,1):     nhvar   - number of variables
%              hinfo(nsta,2):     nhpts   - number of data points
%              hinfo(nsta,3):     castid  - cast identifier
%              hinfo(nsta,4):     hlng    - longitude
%              hinfo(nsta,5):     hlat    - latitude
%              hinfo(nsta,6):     hdpth   - maximum depth
%              hinfo(nsta,7):     htime   - time
%              hinfo(nsta,.):     hscl(1:nhvar) - data scales
%              hinfo(nsta,..):    hflag   - special flag

hinfo(:,1)=repmat(3,R,1);
hinfo(:,3)=[1:1:R]';
hinfo(:,4)=lon;
hinfo(:,5)=lat;
hinfo(:,7)=repmat(julian(2003,8,15,0)-2440000,R,1);
hinfo(:,8:10)=repmat([0.1 0.001 0.001],R,1);
hinfo(:,11)=zeros(R,1);

npts=[];
zmax=[];
for k=1:R
    a=find(~isnan(tdat(k,:)));
    if isempty(a)
        npts=[npts; length(z)];
        zmax=[zmax; z(end)];
    else
        npts=[npts; length(a)];
        zmax=[zmax; z(a(end))];
    end
end

hinfo(:,2)=npts(:);
hinfo(:,6)=zmax(:);

%===== writing output to a MOD file ============================================
%       modname='Levitus_Summer2.mod';
%===============================================================================
disp(' ')
modname=input('     *. Enter the MOD-output Filename w/ suffix:  ','s');
status=whydro(modname,hdr,hinfo,htype,zdat,tdat,sdat);
