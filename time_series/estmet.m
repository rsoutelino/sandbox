

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% FEVEREIRO%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load e_0206.txt

T=e_0206(:,8);
UR=e_0206(:,9);
CV=e_0206(:,10);

vel=e_0206(:,11)./3.6;
dir=e_0206(:,12);
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
t = waitbar(0,'Lendo Fevereiro...');
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


DATA1=datenum(Y',M,D',HORA,MIN,0);

close(t)

%Calculando as direçoes das velocidades horizontais



decl=-22;
dir=decl+dir+180;

r = find(dir < 0);
dir(r) = 360 + dir(r);


r1 = find(dir > 360);
dir(r1) = dir(r1)-360;

clear r r1

Dir1=dir.*pi./180;
[v1,u1] = pol2cart(Dir1,vel);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% MARÇO   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load e_0306.txt

%  T2=e_0306(:,8);
%  UR2=e_0306(:,9);
%  CV2=e_0306(:,10);

vel2=e_0306(:,11)./3.6;
dir2=e_0306(:,12);
%  pp2=e_0306(:,13);

%  CO22=e_0306(:,14);
%  UV2=e_0306(:,15);


D2=e_0306(:,3)-59;
H2=e_0306(:,4);
Y2=e_0306(:,2);
M2=3*ones(1,length(H2));
hh2=num2str(H2);

h12=hh2(:,1);
h22=hh2(:,2);
m12=hh2(:,3);
m22=hh2(:,4);

t = waitbar(0,'Lendo Marco...');

for i=1:length(hh2)
waitbar(i/ 3240)
    a2=str2num(h12(i));
    if length(a2)==0
       H12(i)=0;
    else
       H12(i)=a2;
    end
    b2=str2num(h22(i));
    if length(b2)==0
       H22(i)=0;
    else
       H22(i)=b2;
    end
    c2=str2num(m12(i));
    if length(c2)==0
       M12(i)=0;
    else
       M12(i)=c2;
    end
    d2=str2num(m22(i));
    if length(d2)==0
       M22(i)=0;
    else
       M22(i)=d2;
    end
end

close(t)


HORA2=(H12*10)+H22;
MIN2=(M12*10)+M22;


DATA2=datenum(Y2',M2,D2',HORA2,MIN2,0);

%Calculando as direçoes das velocidades horizontais


% corrigindo declinacao magnetica

dir2=decl+dir2+180;

r = find(dir2 < 0);
dir2(r) = 360 + dir2(r);


r1 = find(dir2 > 360);
dir2(r1) = dir2(r1)-360;


Dir2=dir2.*pi./180;
[v2,u2] = pol2cart(Dir2,vel2);

stop


dddd=DATA2(3030:3168);
udd=u2(3030:3168);
vdd=v2(3030:3168);

mu=mean(udd);
mv=mean(vdd);

estick(dddd,udd./100,vdd./100,1,'m.s{-1}');
datetick('x',13)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  JUNTANDO OS DADOS


DATA=[DATA1 DATA2];
u=[u1;u2];
v=[v1;v2];

%%% extraindo 3 ultimos dias

dia1=e_0206(:,3);
dia2=e_0306(:,3);
dia=[dia1;dia2];

f20=find(dia==79);
f21=find(dia==80);
f22=find(dia==81);

u20=u(f20); v20=v(f20);
u21=u(f21); v21=v(f21);
u22=u(f22); v22=v(f22);

% salvando o vento medio para o dia 22, do adcp

mu=mean(u22);
mv=mean(v22);

save 




%%%%%%%%%%%%%%% ACERTANDO AS DATAS %%%%%%%%%%%%%%%%%%%%%%

DI=find(DATA==datenum('17-Feb-2006 19:03'));
DT=find(DATA==datenum('23-Mar-2006 10:33'));
DATA=DATA(DI:3:DT);

u=u(DI:3:DT);
v=v(DI:3:DT);


%%%%%%%%%%%%%%%%%% PLOTAGEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
estick(DATA(1:10:end),u(1:10:end)./100,v(1:10:end)./100,45,'m.s{-1}');
datetick('x',19);
title('Série Temporal de Vento - São Sebastião')


figure(2)
mypolar2(dir*pi/180,vel./100,'r.');



estmet_v=v;
estmet_u=u;

save estmet estmet_v estmet_u
!cp estmet.mat ../.

