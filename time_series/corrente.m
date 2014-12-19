%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Corrente %

% elaborado por Sandro Vianna Paixao (2006)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Descricao das colunas:

%AVN - velocidade media (1 min) - componente norte-sul (positivo para norte) - cm/s
%AVE - velocidade media (1 min) - componente leste-oeste (positivo para leste) - cm/s
%ASPD - velocidade media (1 min) - modulo - cm/s
%AVDIR - velocidade media (1 min) - direcao (graus, norte = 0.0)
%ATLT - inclinacao media (1 min) - graus, vertical = 0.0
%TIME - hora GMT
%DATE - data GMT
%HDNG - direcao media da bussola (graus) - nao utilizar !
%BATT - tensao media da bateria (V) 
%TX e TY - inclinacao media em eixos referenciados ao instrumento - nao utilizar !
%VN e VE - componentes N/S e E/W da velocidade instantanea (cm/s)
%STEMP - temperatura (graus Centigrados) 


% clear all remove todas as variaveis, funcoes, e arquivos M da memoria, deixando o "workspace" vazio.
% close all apaga todas as figuras.
clear all;close all;


% Leitura dos Dados do arquivo ss230306t.dat.
load ss230306t.dat 

reg =  ss230306t;

AVN=reg(:,1);
AVE=reg(:,2);
ASPD=reg(:,3);
AVDIR=reg(:,4);
ATLT=reg(:,5);
TIME=reg(:,6);
DATE=reg(:,7);
BATT=reg(:,9); 
VN_VE=reg(:,11);

Modulo=ASPD;
Direcao=AVDIR;

% Conversao dos dados da Corrente de graus para radianos.
THETA=Direcao*pi/180;
RHO=Modulo;

% Conversao dos dados de coordenada polar para cartesiana.
[u,v] = pol2cart(THETA,RHO);


% DATENUM([Y,M,D,H,MI,S])
% data inicial: 10-02-2006 11:01:00
datai=datenum(2006,02,10,11,01,0);
% data final: 23-03-2006 11:33:00
dataf=datenum(2006,03,23,11,33,0);
% Intervalo: 00:30:00
int=datenum(0,0,0,0,30,0);

DATA=[datai:int:dataf];

figure(1)
    plot(DATA)
    title('Registro de Corrente - Sao Sebastiao')
    xlabel('tempo em horas')
    ylabel('Amplitude') 
    datetick('x',7)

% visualiza de forma simples os graficos das series temporais
feather(u,v);


% save serie_corrente.dat % cria arquivo para ser lido pelo t_tide ...

