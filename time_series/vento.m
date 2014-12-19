%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vento %

% elaborado por Sandro Vianna Paixao (2006)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Descricao das colunas:

%col 1 - id. do sistema
%col 2 - ano
%col 3 - dia juliano 
%col 4 - hora e minuto 
%col 5 - segundo
%col 6 - temperatura interna (graus centigrados)
%col 7 - tensao da bateria (V)
%col 8 - temperatura do ar (graus centigrados)
%col 9 - umidade relativa (%)
%col 10 - chuva (mm) - ltimos 10 minutos) 
%col 11 - velocidade do vento (km/h)
%col 12 - direcao do vento (graus) - referenciada ao norte magnetico, sentido de onde o vento vem 
%col 13 - pressao barometrica (hPa)
%col 14 - CO2 (%)
%col 15 - UV (mW/cm2)
%col 16 - temperatura do sensor UV (graus centigrados)
%col 17 - piranometro (W/m2)
%col 18-36 - dados do ondografo (sensor inoperante)

%Datalogger Year	integer
%Datalogger Day	integer
%Datalogger Hour-Minute	integer	1st 2 digits = hour	2nd 2 digits = minute
%Datalogger Second	integer
%Datalogger Internal Temp	float	C
%Datalogger Input Voltage	float	Volts
%Air Temperature	float	C
%Relative Humidity	float	%
%Rainfall (over 10 min)	float	mm
%Wind Speed	float	Km/Hr
%Wind Direction	float	degrees
%Barometric Pressure	float	hPa
%Carbon Dioxide	float	% X100
%UV A & B Light	float	mW/cm2
%UV Sensor Temperature	float	C
%Quantum Light Sensor	float	W/m2




% clear all remove todas as variaveis, funcoes, e arquivos M da memoria, deixando o "workspace" vazio.
% close all apaga todas as figuras.

clear all;close all;


% load carrega os arquivos de dados e_0206.txt e e_0306.txt.
load e_0206.txt 
load e_0306.txt

reg1=e_0206;
reg2=e_0306;


% Concatenacao de reg1 e reg2 para criacao de uma unica matriz.

reg3=[reg1;reg2];


% Os dados de Intensidade do Vento estao na 11 coluna, em km/h.
% Os dados de Direcao do Vento estao na 12 coluna, em graus, e referenciada ao norte magnetico, sentido de onde o vento vem. 

Intensidade=reg3(:,11);
Direcao=reg3(:,12);


% Referenciar o vento ao Norte Verdadeiro.
Direcao1=Direcao+22;


% Conversao dos dados da Intensidade do Vento de graus para radianos.

THETA=Direcao1*pi/180;
RHO=Intensidade;


% Inversao das componentes de saida para obter u e v, corretamente.
% Conversao dos dados de coordenada polar para cartesiana.

[v,u] = pol2cart(THETA,RHO);

      

% Sublot: Visualizacao dos graficos em series temporais.
% Datetick: plota no eixo 'x', na configuracao de data n 7. 
% "grid on" insere grades nos graficos.


figure(1)
    feather(u,v);
    datetick('x',7);grid on
    axis equal;
    title('Intensidade e Direcao do Vento')
    xlabel('Dia Juliano')
    ylabel('Intensidade em [Km/h] e Direcao do Vento') 




% feather plota vetores de velocidade.

 
% Criacao de arquivo para ser lido pelo t_tide.
save serie_vento.dat 

