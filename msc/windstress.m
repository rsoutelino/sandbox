%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    PROGRAMA PARA MANIPULACAO DOS DADOS DE
%     TENSAO DE CIZALHAMENTO DO VENTO PARA
%      O PERIODO DE REALIZACAO DO CRUZEIRO 
%              OCEANO LESTE II
%   Rafael Guarino Soutelino - mestrado IOUSP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc

lonlim = [-41 -33.7]; latlim = [-20.5 -10.2];
m_proj('mercator','long',[min(lonlim) max(lonlim)],'lat',[min(latlim) max(latlim)],'on');

load ../../common/etopo2_leste.mat
load /usr/local/matlab/toolbox/lado/m_map/ocean_colormap

% lendo os dados de TCV para o periodo do cruzeiro

U=[]; V=[]; c=1; 
for k = 32:92
   eval(['load tau0',num2str(k),'.mat;'])
   U(:,:,c) = u;
   V(:,:,c) = v;
   c=c+1;
end 

[l,c,g] = size(U);

% calculando campo medio de TAU para o periodo do cruzeiro

for k = 1:l
   um = U(k,:,:);
   vm = V(k,:,:);
   um = reshape(um,g,c);
   vm = reshape(vm,g,c);
      for i = 1:c
          taux(k,i) = nanmean(um(:,i));
          tauy(k,i) = nanmean(vm(:,i));
      end
end

% calculando transporte de ekman

rho = 1026.954;
f0 = sw_f(-15);

ue = tauy./(rho*f0);
ve = - taux./(rho*f0);