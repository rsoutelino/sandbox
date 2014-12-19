function [U,V]=intdir2uv(int,dir,decl_mag,ang_rot)
%  [U,V]=intdir2uv(int,dir,decl_mag,ang_rot)
% FUNCAO INTDIR2UV.M Decompoem o vetor corrente definido pela intensidade
%                    e direcao (ref. Norte -> Este) considerando a 
%	             declinacao magnetica e a rotacao do sistema de coordenadas
% Uso: entre com intensidade, direcao,declinacao magnetica e orientacao do
% eixo Y, a partir do Norte (0 deg.). Por exemplo, Canal de Sao Sebastiao =
% 51 deg. A declinacao para oeste e' negativa, p.ex.: -18 deg. 
% 
% Se valor de ang_rot nao for suprido, e' assumido que nao ha' rotacao
% alem disso, se decl_mag tambem nao for suprido nao e' feito correcao
% magnetica


% Author: Roberto Fioravanti Carelli Fontes
% Depto. de Oceanografia Fisica IOUSP
% Laboratorio de Hidrodinamica Costeira (LHICO)

 if nargin < 4,
   ang_rot=0;
 if nargin == 2,
 decl_mag=0;
 end 
 end

% decompor o vetor
dir=dir+decl_mag;
dir=mod(dir,360);
dir=dir*pi/180;
% inlinacao da linha de costa
 
ang_rot=(ang_rot)*pi/180; % alinhar o sistema com a costa

u=int.*sin(dir); % aqui eu passo de NE para
v=int.*cos(dir); %                  XY

U = u*cos(ang_rot) - v*sin(ang_rot); % aqui eu faco a rotacao
V = u*sin(ang_rot) + v*cos(ang_rot); % segundo o alinhamento da costa



