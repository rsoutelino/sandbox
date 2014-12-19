%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLUCAO DA EQUACAO DA ADVECCAO - DIFUSAO BI-DIMENSIONAL LINEAR 
% PELO METODO DE LEAPFROG - CENTRADO NO TEMPO E NO ESPACO. 
% DISPERSAO DE UM CONTAMINANTE NA BAIA DE GUANABARA  
% MODELAGEM NUMERICA DE PROCESSOS COSTEIROS E ESTUARINOS
% LISTA DE EXERCICIOS 1 - RAFAEL GUARINO SOUTELINO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc

% carregando colorbar
load /usr/local/matlab2008b/toolbox/lado/m_map/redblue

% carregando linha de costa
load costa.mat; lonlc = ncst(:,1); latlc = ncst(:,2);

% CONSTANTES DO MODELO
u = 0; v = 0.3; % corrente na baia
h = 200; % resolucao espacial
T = 3600; % passo de tempo
D = 3.6e-5; % coeficiente de difusao
ci = 90; % concentracao inicial no centro da baia
% lembrando CE => 0 <= T <= h^2/4D

mmax = 10000; % numero de passos de tempo

% declarando limites da baia de guanabara
lonlim = [-43.33 -43]; latlim = [-22.95 -22.6];

% criando eixos de distancia em m para construir grade do modelo
jend = sw_dist([latlim(1) latlim(1)],[lonlim(1) lonlim(2)],'km')*1000;
kend = sw_dist([latlim(1) latlim(2)],[lonlim(1) lonlim(1)],'km')*1000;

j = 0:h:jend; k = 0:h:kend; lk = length(k); lj = length(j);
lon = linspace(lonlim(1),lonlim(end),length(j));
lat = linspace(latlim(1),latlim(end),length(k));
[J,K] = meshgrid(j,k);
[lon,lat] = meshgrid(lon,lat);

% definindo condicoes iniciais do modelo
catu = zeros(lk,lj); catu(ceil(lj/2),ceil(lk/2)) = ci;
cren = zeros(lk,lj);

% calculando constantes:
b = (T*D)/(h^2);
ax = (u*T)/(2*h);
ay = (v*T)/(2*h);

% LOOP NO TEMPO
% CONDICOES DE CONTORNO
% FORMULA DE RECORRENCIA
% PLOTAGEM (PRESSIONE ENTER PARA EVOLUIR NO TEMPO)
% EVOLUCAO NO TEMPO DAS VARIAVEIS
% FAZENDO SEM SPLITTING


for m = 2:mmax
   tempo = T*m;
   cren(2:lk-1,2:lj-1) = catu(2:lk-1,2:lj-1)...
                       - ax*(catu(2:lk-1,3:lj) - catu(2:lk-1,1:lj-2))... % adveccao zonal
                       - ay*(catu(3:lk,2:lj-1) - catu(1:lk-2,2:lj-1))... % adveccao meridional
                       + b*(catu(2:lk-1,3:lj) - 2*catu(2:lk-1,2:lj-1) + catu(2:lk-1,1:lj-2)... % difusao
                       +    catu(3:lk,2:lj-1) - 2*catu(2:lk-1,2:lj-1) + catu(1:lk-2,2:lj-1));  % difusao
   % plotando
   pcolor(lon,lat,cren); shading flat; caxis([-ci ci]); colormap(redblue); colorbar; hold on;
   plot(lonlc,latlc,'k');	
%     axis([-43.175 -43.15 -22.79 -22.76])
%     plot(lon(ceil(lj/2),:),cren(ceil(lj/2),:));
%     axis([lonlim(1) lonlim(2) 0 200])
   title(['Dispersao de poluente no centro da baia de guanabara - ',num2str(tempo),' s'])
   xlabel('Longitude')
   ylabel('Latitude')
   hold off
   pause
   
   % preparando para proximo passo de tempo
   catu=cren;
%     catu(ceil(lj/2),ceil(lk/2)) = ci; % mantendo emissao maxima no centro da baia
end
