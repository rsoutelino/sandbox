%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   SOLUCAO DA EQUACAO DA ADVECCAO/DIFUSAO BI-DIMENSIONAL LINEAR 
%        POR ESQUEMA AVANCADO NO TEMPO E CENTRADO NO ESPACO, 
%        USANDO SPLITTING DE PROCESSOS DIFUSIVOS/ADVECTIVOS 
%   PROBLEMA: DISPERSAO DE UM CONTAMINANTE NA BAIA DE GUANABARA  
%     MODELAGEM NUMERICA DE PROCESSOS COSTEIROS E ESTUARINOS
%       LISTA DE EXERCICIOS 1 - RAFAEL GUARINO SOUTELINO
%
% Ultima Modificacao: 06/10/2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc

% carregando colorbar monotonica
load redblue

% carregando linha de costa
load costa.mat; lonlc = ncst(:,1); latlc = ncst(:,2);

% PARAMETROS DO MODELO
u = 0; v = 0.3; % corrente na baia
h = 200; % espaçamento de grade
T = 60; % passo de tempo
D = 20; % coeficiente de difusao
ci = 200; % concentracao inicial no centro da baia
mmax = 10000; % numero de passos de tempo

% declarando limites da baia de guanabara
lonlim = [-43.33 -43]; latlim = [-22.95 -22.6];
m_proj('mercator','lat',[latlim],'lon',[lonlim])

% estabelecendo grade do modelo
jend = sw_dist([latlim(1) latlim(1)],[lonlim(1) lonlim(2)],'km')*1000;
kend = sw_dist([latlim(1) latlim(2)],[lonlim(1) lonlim(1)],'km')*1000;
j = 0:h:jend; k = 0:h:kend; lk = length(k); lj = length(j);

% estabelecendo grade em lat lon para plotar
lon = linspace(lonlim(1),lonlim(end),length(j));
lat = linspace(latlim(1),latlim(end),length(k));
[J,K] = meshgrid(j,k);
[lon,lat] = meshgrid(lon,lat);

% definindo condicoes iniciais do modelo
catu = zeros(lk,lj); 
catu(ceil(lk/2),ceil(lj/2)) = ci;
crenadv = zeros(lk,lj); 
crendif = zeros(lk,lj);

% calculando constantes:
b = 2*T*D/h^2;
ax = u*T/h;
ay = v*T/h;

% LOOP NO TEMPO

for m = 2:mmax
   tempo = T*m;
   % calculando splitting para adveccao
   crenadv(2:lk-1,2:lj-1) = catu(2:lk-1,2:lj-1)...
                          - ax*(catu(2:lk-1,3:lj) - catu(2:lk-1,1:lj-2))...
                          - ay*(catu(3:lk,2:lj-1) - catu(1:lk-2,2:lj-1));
   f = find(crenadv < 0);
   crenadv(f)=0;

   % calculando splitting para difusao
   crendif(2:lk-1,2:lj-1) = catu(2:lk-1,2:lj-1)...
                          + b*(catu(2:lk-1,3:lj) -2*catu(2:lk-1,2:lj-1) + catu(2:lk-1,1:lj-2)...
                          +    catu(3:lk,2:lj-1) -2*catu(2:lk-1,2:lj-1) + catu(1:lk-2,2:lj-1));
   f = find(crendif < 0);
   crendif(f)=0;

   % solucao final
   cren = crenadv/2 + crendif/2;

   % plotando
   pcolor(lon,lat,cren); shading flat; caxis([0 ci/2]); colormap(redblue(32:end,:)); colorbar; hold on;
   axis('equal')
   plot(lonlc,latlc,'k');	
%     axis([-43.175 -43.15 -22.79 -22.76])
   title(['Dispersao de poluente - ',num2str(tempo/3600),' h'])
   xlabel('Longitude')
   ylabel('Latitude')
   hold off
   pause(0.001)
   
   % preparando para proximo passo de tempo
   catu = cren;
   catu(ceil(lk/2),ceil(lj/2)) = ci; % mantendo emissao maxima no centro da baia
end

% confeccao de figuras finais

plot(lat(:,ceil(lj/2)),cren(:,ceil(lj/2)),'linewidth',2);hold on
axis([latlim(1) latlim(2) -10 250])
xlabel('Latitude')
ylabel('Concentracao')
title(['Dispersao de poluente - ',num2str(round(tempo/3600)),' h'],'fontweight','bold')
grid on
pbaspect([1.5 0.5 1.5])

m_contourf(lon,lat,cren,50);shading flat; caxis([0 ci/2]); colormap(redblue(32:end,:));
cc = colorbar('horiz'); pos = get(cc,'position');
hold on;
m_usercoast('costa.mat','patch',[0.169 0.687 0.639])
m_grid('box','fancy')
set(cc,'position',[pos(1)*2.5 pos(2)*7.3 pos(3)/3 pos(4)/3])
m_text(-43.314,-22.6376,['Dispersao de poluente - ',num2str(round(tempo/3600)),' h'],'fontweight','bold')
m_text(-43.3082,-22.69,['C = [kg^2 m^{-1}]'],'fontweight','bold')












