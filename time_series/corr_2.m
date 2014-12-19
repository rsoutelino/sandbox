%  - descrição das colunas:
%  
%   AVN - velocidade média (1 min) - componente norte-sul (positivo para norte) - cm/s
%   AVE - velocidade média (1 min) - componente leste-oeste (positivo para leste) - cm/s
%   ASPD - velocidade média (1 min) - módulo - cm/s
%   AVDIR - velocidade média (1 min) - direção (graus, norte = 0.0)
%   ATLT - inclinação média (1 min) - graus, vertical = 0.0
%   TIME - hora GMT
%   DATE - data GMT
%   HDNG - direção média da bússola (graus) - não utilizar !
%   BATT - tensão média da bateria (V) 
%   TX e TY - inclinação média em eixos referenciados ao instrumento - não utilizar !
%   VN e VE - componentes N/S e E/W da velocidade instantânea (cm/s)
%   STEMP - temperatura (graus Centígrados) 

clear all;close all;

%load ceb0306.Asc % Como sai da convesao de binario para ascii do aparelho

fid=fopen('ss230306t.Dat','r'); 
t = waitbar(0,'Relax...');

for i=1:12
linha=fgetl(fid);
end

cont=1;
while linha ~= -1

       AVN(cont)=str2num(linha(3:8));
       AVE(cont)=str2num(linha(13:18));

       DATA(cont,1)=datenum(linha(50:69));

waitbar(cont/5910)

cont=cont+1;
linha=fgetl(fid);
end

AVE=AVE(4:end-5)';
AVN=AVN(4:end-5)';
DATA=DATA(4:end-5);
close(t)

%Calculando as direçoes das velocidades horizontais

vel=sqrt(AVE.^2 + AVN.^2);

dir=atan2(AVE,AVN);
dir=dir.*180./pi;

decl=-22;
dir=decl+dir;

r = find(dir < 0);
dir(r) = 360 + dir(r);


r1 = find(dir > 360);
dir(r1) = dir(r1)-360;



% Plotagem
figure(1)
subplot(311)
plot(DATA,AVE./100)
datetick('x',1);grid on
title('Componente Zonal (E-W) - [m.s^{-1}]')
ylabel('[m.s^{-1}]') 

subplot(312)
plot(DATA,AVN./100)
datetick('x',1);grid on
title('Componente Meridional (N-S) - [m.s^{-1}]')
ylabel('[m.s^{-1}]') 

subplot(313)
estick(DATA,AVE./100,AVN./100,45,'m.s{-1}');
datetick('x',1);
title('Série Temporal Correntográfica - São Sebastião')


figure(2)
mypolar2(dir*pi/180,vel,'r.');









