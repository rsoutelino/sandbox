close all;clear all

%***************************************************************
% cria série sintética de onda horária com periódo de 12 horas
% série espressa em dias
%**************************************************************

t=0:1/24:7; t=t';

onda1= 4*sin(2*pi/.5*t);

onda2= 2*sin(2*pi/.5*t);

figure

plot(t,onda1,'r',t,onda2,'b')
     xlabel('tempo em dias')
     ylabel('amplitude em metros')

%*********************************************
% interpretando figuras de mérito estatísticas
%*********************************************

% calcula coeficiente de correlação

[r,p,rlo,rup]=corrcoef(onda1, onda2);

corr=r(1,2);
prob=p(1,2);

% calcula rms (ou remq)

rms= sqrt(sum( (onda1-onda2).^2)./sum(onda2.^2))

%******************************************
% fazendo ajustes por mínimos quadrados...
%******************************************
stop
% caso 1 -só onda

arg=sin(2*pi/.5*t);

amp=arg\onda2

% caso2 - onda+bias

onda3=1+onda2;

arg=[ ones(size(t)) sin(2*pi/.5*t)];

amp2=arg\onda3

% caso 3 - onda+bias+tendência linear

onda4=1+.5*t+onda2;

arg=[ ones(size(t)) t sin(2*pi/.5*t)];

amp3=arg\onda4

%**********************************************
% usem a rotina sinfitc para fazer estes ajustes
% escolhendo um período e divirtam-se...
%**********************************************

%**********************************************
% em termos de análise harmônica, que tal brincar
% com a sinfitc antes de rodar o t_tide?
%***********************************************

% para as componentes M2, S2, K2, N2, K1, O1

%    ampM2=sqrt((coef(2).^2)+(coef(3).^2));  % M2
%    phaseM2=atan2(coef(3),coef(2));
%    uM2=ampM2*sin((2*pi/12.421*t + phaseM2));

%    ampS2=sqrt((coef(6).^2)+(coef(7).^2));  % S2
%    phaseS2=atan2(coef(7),coef(6));
%    uS2=ampS2*sin((2*pi/12*t + phaseS2));

%    ampN2=sqrt((coef(10).^2)+(coef(11).^2));  % N2
%    phaseN2=atan2(coef(11),coef(10));
%    uN2=ampN2*sin((2*pi/12.658*t + phaseN2));

%    ampK1=sqrt((coef(4).^2)+(coef(5).^2));  % K1
%    phaseK1=atan2(coef(5),coef(4));
%    uK1=ampK1*sin((2*pi/23.395*t + phaseK1));

%    ampO1=sqrt((coef(8).^2)+(coef(9).^2));  % O1
%    phaseO1=atan2(coef(9),coef(8));
%    uO1=ampO1*sin((2*pi/25.819*t + phaseO1));