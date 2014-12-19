%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   LE OS ARQUIVOS BINARIOS DA RDI E SALVA 
%      ARQUIVOS .MAT COM AS VARIAVEIS
%     Rafael Soutelino - Mestrado IOUSP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all; clc ; warning off

U=[]; V=[]; lon=[]; lat=[]; dt=[]; pg=[];
for j = [1:54]
    if j < 10
       eval(['adcp = rdradcp(''../../../dados/leste2/adcp/brutos/OE00',num2str(j),'_000000.LTA'',1)']);
    else
       eval(['adcp = rdradcp(''../../../dados/leste2/adcp/brutos/OE0',num2str(j),'_000000.LTA'',1)']);
    end
   
    U = [U adcp.east_vel(1:end-1)];
    V = [V adcp.north_vel(1:end-1)];
    lon = [lon adcp.nav_slongitude(1:end-1)];
    lat = [lat adcp.nav_slatitude(1:end-1)];
    dt = [dt adcp.mtime(1:end-1)];
    pg4 = adcp.perc_good(1:end-1);
    pg4 = squeeze(pg4(:,4,:));
    pg = [pg pg4];

end