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
lat=-23.83;
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

close all

%ACERTANDO DADOS CORRENTE - FAZENDO MÉDIA PARA 3 VALORES

DATA2=DATA(1:3:end-1);
for i=1:ceil(length(AVN)/3)-1


V(i)=mean(AVN(3*i-2:i*3));
U(i)=mean(AVE(3*i-2:i*3));

end

vel2=sqrt(U.^2 + V.^2);

dir2=atan2(U,V);
dir2=dir2.*180./pi;

decl=-22;
dir2=decl+dir2;

r2 = find(dir2 < 0);
dir2(r2) = 360 + dir2(r2);


r3 = find(dir2 > 360);
dir2(r3) = dir2(r3)-360;



% Plotagem
figure
subplot(311)
plot(DATA2,U./100)
datetick('x',1);grid on
title('Componente Zonal (E-W) - [m.s^{-1}]')
ylabel('[m.s^{-1}]') 

subplot(312)
plot(DATA2,V./100)
datetick('x',1);grid on
title('Componente Meridional (N-S) - [m.s^{-1}]')
ylabel('[m.s^{-1}]') 

subplot(313)
estick(DATA2,U./100,V./100,45,'m.s{-1}');
datetick('x',1);
title('Série Temporal Correntográfica - São Sebastião')


figure
mypolar2(dir2*pi/180,vel2,'r.');



% Análise de Mare T_tide

path(path,'~/monografia/dados/mare/t_tide_v1.1')

diai = DATA2(1);
xinv=V./100;
xinu=U./100;

[tidestrucu,poutu]=t_tide(xinu','interval',.5,'start time',diai,'latitude',lat,'output',['rel_SS.txt'],'error','cboot','synthesis',1);
[tidestrucv,poutv]=t_tide(xinv','interval',.5,'start time',diai,'latitude',lat,'output',['rel_SS.txt'],'error','cboot','synthesis',1);

figure
plot(DATA2,poutu,'k')
datetick('x',7);grid on
title('Corrente de Mare Componente Zonal [m.s^{-1}]')
xlabel('Dia do Mes')
ylabel('Corrente de Mare [m.s^{-1}]') 

figure
plot(DATA2,poutv,'k')
datetick('x',7);grid on
title('Corrente de Mare Componente Meridional [m.s^{-1}]')
xlabel('Dia do Mes')
ylabel('Corrente de Mare [m.s^{-1}]') 

% Plotando Constituintes Significativas
figure
fsig=tidestrucu.tidecon(:,1)>tidestrucu.tidecon(:,2); % Picos significativos
semilogy([tidestrucu.freq(~fsig),tidestrucu.freq(~fsig)]',[.0005*ones(sum(~fsig),1),tidestrucu.tidecon(~fsig,1)]','color',([.5 .5 .5]));
line([tidestrucu.freq(fsig),tidestrucu.freq(fsig)]',[.0005*ones(sum(fsig),1),tidestrucu.tidecon(fsig,1)]','marker','.','color','k');
line(tidestrucu.freq,tidestrucu.tidecon(:,2),'linestyle',':','color',[.5 .5 .5]);
set(gca,'ylim',[.0005 1],'xlim',[0 .5]);
xlabel('frequencia (cph)');
text(tidestrucu.freq,tidestrucu.tidecon(:,1),tidestrucu.name,'rotation',45,'vertical','base');
ylabel('Amplitude (m.s{-1})');
text(.35,.2,'Constituintes significativas ','color','k');
text(.35,.1,'Constituintes nao-significantivas' ,'color',([.5 .5 .5]));
text(.35,.05,'95%','color',[.5 .5 .5]);



figure
fsig=tidestrucv.tidecon(:,1)>tidestrucv.tidecon(:,2); % Picos significativos
semilogy([tidestrucv.freq(~fsig),tidestrucv.freq(~fsig)]',[.0005*ones(sum(~fsig),1),tidestrucv.tidecon(~fsig,1)]','color',([.5 .5 .5]));
line([tidestrucv.freq(fsig),tidestrucv.freq(fsig)]',[.0005*ones(sum(fsig),1),tidestrucv.tidecon(fsig,1)]','marker','.','color','k');
line(tidestrucv.freq,tidestrucv.tidecon(:,2),'linestyle',':','color',[.5 .5 .5]);
set(gca,'ylim',[.0005 1],'xlim',[0 .5]);
xlabel('frequencia (cph)');
text(tidestrucv.freq,tidestrucv.tidecon(:,1),tidestrucv.name,'rotation',45,'vertical','base');
ylabel('Amplitude (m.s{-1})');
text(.35,.2,'Constituintes significativas ','color','k');
text(.35,.1,'Constituintes nao-significantivas' ,'color',([.5 .5 .5]));
text(.35,.05,'95%','color',[.5 .5 .5]);





