%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Converte o registro do maregrafo Anderaa em pressao %

% elaborado por Sandro Vianna Paixao (2006)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Utilizados os dados a partir de 10 de fevereiro de 2006

%- descri�o das colunas:

% coluna 1 - data
% coluna 2 - hora
% coluna 3 - contagem de Condutividade = N
% coluna 4 - contagem de Temperatura = N
% coluna 5 - contagem de Press� (high bits) = N3
% coluna 6 - contagem de Press� (low bits) = N4

% obs: as linhas onde a coluna 3 contem (0007) N� s� linhas de dados - elas indicam data e hora
    

% clear all remove todas as variaveis, funcoes, e arquivos M da memoria, deixando o "workspace" vazio.
% close all apaga todas as figuras.

clear all;close all;


% fopen abre um arquivo ou obtem informacao sobre arquivos abertos.
% fid = fopen(ceb0306.Asc) abre o arquivo ceb0306.Asc para leitura. 
% fid = fopen(ceb0306.Asc,'permissao') abre o arquivo ceb0306.Asc no modo especificado pela permissao.
% permissao = 'r' abre o arquivo para leitura (default).

fid=fopen('ceb0306.Asc','r');


% tline = fgetl(fid) le uma linha do texto do arquivo ceb0306.Asc e retorna o dado para tline.

tline=fgetl(fid);


% Contador igual a 1 enquanto tline for diferente de -1.

cont=1;
while tline ~= -1
    
    
% Se lenght(tline) for igual a 49 ou 48, entao, tline=fgetl(fid).
    
if length(tline)==49 | length(tline)==48

tline=fgetl(fid);
        
% Se lenght(tline) nao for igual a 49 ou 48, entao: 
 
% A funcao str2num converte "string" para numero.
% NN: os numeros que estao nas posicoes de 25 a 29 - valor da Temperatura.
% N3: os numeros que estao nas posicoes de 30 a 34 - valor da Pressao (high bits).
% N4: os numeros que estao nas posicoes de 35 a 39 - valor da Pressao (low bits).
% D:  os numeros que representam o dia.
% M:  os numeros que representam o mes.
% Y:  os numeros que representam o ano.
% MI: os numeros que representam os minutos.
    
    else
               
       NN(cont)=str2num(tline(25:29));
       N3(cont)=str2num(tline(30:34));
       N4(cont)=str2num(tline(35:39));
       D=str2num(tline(1:2));
       M=str2num(tline(4:5));
       Y=str2num(tline(7:10));
       MI=str2num(tline(14:15));
              

% Se length(MI) for igual a zero, entao, em MI os numeros que representam os minutos estarao nas posicoes 15 e 16 e
% os numeros que representam as horas estarao nas posicoes 12 e 13.
    
if length(MI)==0
MI=str2num(tline(15:16));
H=str2num(tline(12:13));

% Se length(MI) nao for igual a zero, entao, os numeros que representam as horas estarao nas posicoes 11 e 12.

else
H=str2num(tline(11:12));
end

% Armazenamento dos dados(cont) e retorna como um numero os valores do ano, mes, dia, hora e minuto, datenum(Y,M,D,H,MI,S).
% A funcao datenum converte "date strings" e "date vectors" em numeros seriais de data. 

       DATA(cont)=datenum([Y,M,D,H,MI,0]);
       cont=cont+1;
    end

% Realiza isso linha por linha. 
tline=fgetl(fid);
end


% Polinomio de conversao dos fatores em pressao (psia).

A= -2.337687e03;
B=  4.275187e-03;
C= -1.139844e-09;
D=  3.393918e-17;

N = (N3.*1024)+N4;

P = A+(B.*N)+(C.*(N.^2))+(D.*(N.^3));

pp = (P - 14.7).*0.6894759;



%Conversao de dbar em metros via seawater

% Introducao da latitude do maregrafo.
lat_maregrafo=-23.83;             


% Convercao para metros via rotinas sea water.
depth1=sw_pres(pp',lat_maregrafo); 


% Polinomio de conversao da temperatura (C)

AA=-3.028;
BB=3.602e-2;
CC=-7.117e-6;
DD=1.015e-8;

     T=AA + BB * NN + CC * (NN.^2) + DD * (NN.^3);


% Sublot: Visualizacao dos graficos em series temporais.
% Datetick: plota no eixo 'x', na configuracao de data n 7. 
% "grid on" insere grades nos graficos.


figure(1)
     subplot(311)
plot(DATA,pp)
datetick('x',1);grid on
title('Amplitude de Mare [dbar]')
xlabel('Data')
ylabel('Amplitude [dbar]') 

  
  subplot(312)
plot(DATA,depth1,'k')
datetick('x',1);grid on
title('Amplitude de Mare [cm]')
xlabel('Data')
ylabel('Amplitude [cm]') 


    subplot(313)
plot(DATA,T,'r')
datetick('x',1);grid on
title('Temperatura [^C]')
xlabel('Data')
ylabel('Temperatura [^C]')


% Criacao de arquivo para ser lido pelo t_tide.

save serie_mare.dat depth1 -ascii 