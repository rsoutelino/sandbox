

%   - descrição das colunas:
%  
%  col 1 - id. do sistema
%  col 2 - ano
%  col 3 - dia juliano 
%  col 4 - hora e minuto 
%  col 5 - segundo
%  col 6 - temperatura interna (graus centígrados)
%  col 7 - tensão da bateria (V)
%  col 8 - temperatura do ar (graus centígrados)
%  col 9 - umidade relativa (%)
%  col 10 - chuva (mm) - últimos 10 minutos) 
%  col 11 - velocidade do vento (km/h)
%  col 12 - direção do vento (graus) - referenciada ao norte magnético, sentido de onde o vento vem 
%  col 13 - pressão barométrica (hPa)
%  col 14 - CO2 (%)
%  col 15 - UV (mW/cm2)
%  col 16 - temperatura do sensor UV (graus centígrados)
%  col 17 - piranometro (W/m2)
%  col 18-36 - dados do ondógrafo (sensor inoperante)
%  
%  Datalogger Year	integer
%  Datalogger Day	integer
%  Datalogger Hour-Minute	integer	1st 2 digits = hour	2nd 2 digits = minute
%  Datalogger Second	integer
%  Datalogger Internal Temp	float	°C
%  Datalogger Input Voltage	float	Volts
%  Air Temperature	float	°C
%  Relative Humidity	float	%
%  Rainfall (over 10 min)	float	mm
%  Wind Speed	float	Km/Hr
%  Wind Direction	float	degrees
%  Barometric Pressure	float	hPa
%  Carbon Dioxide	float	% X100
%  UV A & B Light	float	mW/cm2
%  UV Sensor Temperature	float	°C
%  Quantum Light Sensor	float	W/m2


clear all;close all;

load e_0206.txt

T=e_0206(:,8);
UR=e_0206(:,9);
CV=e_0206(:,10);

vel=e_0206(:,11)./3.6;
dir=e_0206(:,12)+180;
pp=e_0206(:,13);

CO2=e_0206(:,14);
UV=e_0206(:,15);


D=e_0206(:,3)-31;
H=e_0206(:,4);
Y=e_0206(:,2);
M=2*ones(1,length(H));
hh=num2str(H);

h1=hh(:,1);
h2=hh(:,2);
m1=hh(:,3);
m2=hh(:,4);
t = waitbar(0,'AAAAhhhhhhhh...');
for i=1:length(hh)
waitbar(i/ 4032)
    a=str2num(h1(i));
    if length(a)==0
       H1(i)=0;
    else
       H1(i)=a;
    end
    b=str2num(h2(i));
    if length(b)==0
       H2(i)=0;
    else
       H2(i)=b;
    end
    c=str2num(m1(i));
    if length(c)==0
       M1(i)=0;
    else
       M1(i)=c;
    end
    d=str2num(m2(i));
    if length(d)==0
       M2(i)=0;
    else
       M2(i)=d;
    end
end


HORA=(H1*10)+H2;
MIN=(M1*10)+M2;


DATA=datenum(Y',M,D',HORA,MIN,0);

close(t)

%Calculando as direçoes das velocidades horizontais


decl=-22;
dir=decl+dir+180;

r = find(dir < 0);
dir(r) = 360 + dir(r);


r1 = find(dir > 360);
dir(r1) = dir(r1)-360;


D=dir.*pi./180;
[v,u] = pol2cart(D,vel);



% Plotagem
figure(1)
subplot(311)
plot(DATA,u./100)
datetick('x',19);grid on
title('Componente Zonal (E-W) - [m.s^{-1}]')
ylabel('[m.s^{-1}]') 

subplot(312)
plot(DATA,v./100)
datetick('x',19);grid on
title('Componente Meridional (N-S) - [m.s^{-1}]')
ylabel('[m.s^{-1}]') 

subplot(313)
estick(DATA(1:10:end),u(1:10:end),v(1:10:end),45,'m.s{-1}');
datetick('x',19);
title('Série Temporal de Vento - São Sebastião')


figure(2)
mypolar2(dir*pi/180,vel,'r.');








