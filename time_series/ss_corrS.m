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

fid=fopen('ES1801.Dat','r');                          % abre o arquivo com os dados de corrente 
t = waitbar(0,'Relax...');
lat=-23.83;					      % insere a latitude do correntografo	
for i=1:12
linha=fgetl(fid);				      % lê da linha 1 até 12 do arquivo
end

cont=1;
while linha ~= -1

       AVN(cont)=str2num(linha(3:8));		      % armazena a componente meridional da velocidade
       AVE(cont)=str2num(linha(13:18));		      % armazena a componente zonal da velocidade

       DATA(cont,1)=datenum(linha(50:69));	      % monta o vetor de data

waitbar(cont/7711)

cont=cont+1;
linha=fgetl(fid);
end

AVE=AVE(4:end-5)';
AVN=AVN(4:end-5)';
DATA=DATA(4:end-5);
close(t)



%ACERTANDO DADOS CORRENTE - FAZENDO MÉDIA PARA 3 VALORES

DATA2=DATA(1:3:end-1);
for i=1:ceil(length(AVN)/3)-1


V(i)=mean(AVN(3*i-2:i*3));			      % calcula a média do componente meridional da velocidade
U(i)=mean(AVE(3*i-2:i*3));			      % calcula a média do componente zonal da velocidade

end

vel2=sqrt(U.^2 + V.^2);				      % monta o vetor velocidade resultante

dir2=atan2(U,V);				      % monta o vetor de direção
dir2=dir2.*180./pi;

decl=-22;					      % corrige a declinação
dir2=decl+dir2;

r2 = find(dir2 < 0);
dir2(r2) = 360 + dir2(r2);


r3 = find(dir2 > 360);
dir2(r3) = dir2(r3)-360;


D=dir2.*pi./180;
[v,u] = pol2cart(D,vel2);

DF=find(DATA2==datenum('30-Apr-2006 23:31'));	      % monta o vetor da data final
DATA2=DATA2(1:DF);
u=u(1:DF);					      
v=v(1:DF);					      
dir2=dir2(1:DF);
vel2=vel2(1:DF);

% Plotagem
figure;pause
subplot(311)
plot(DATA2,u./100)
datetick('x',19);grid on;axis tight
title('Componente Zonal (E-W) Superficie - [m.s^{-1}]')
ylabel('[m.s^{-1}]') 

subplot(312)
plot(DATA2,v./100)
datetick('x',19);grid on;axis tight
title('Componente Meridional (N-S) Superficie - [m.s^{-1}]')
ylabel('[m.s^{-1}]') 

subplot(325)
estick(DATA2,u./50,v./50,45,'m.s{-1}');
datetick('x',19);set(gca,'YColor',[1 1 1]);grid on;axis tight
title('Serie Temporal Correntografica Superficie')

subplot(326)
mypolar2(dir2*pi/180,vel2./100,'r.');

print -depsc corrS_brt
!epstopdf corrS_brt.eps


% Análise de Mare T_tide

path(path,'~/monografia/dados/mare/t_tide_v1.1')

diai = DATA2(1);
xinv=v./100;
xinu=u./100;

[tidestrucu,poutu]=t_tide(xinu','interval',.5,'start time',diai,'latitude',lat,'output',['relu_SS.txt'],'error','cboot','synthesis',1);
[tidestrucv,poutv]=t_tide(xinv','interval',.5,'start time',diai,'latitude',lat,'output',['relv_SS.txt'],'error','cboot','synthesis',1);


figure;pause
subplot(211)
plot(DATA2,(u./100),DATA2,poutu,'r',DATA2,zeros(size(DATA2)),'k')
datetick('x',19);grid on;axis tight
title('Corrente de Mare Componente Zonal [Superficie]')
ylabel('[m.s^{-1}]') 


subplot(212)
plot(DATA2,(v./100),DATA2,poutv,'r',DATA2,zeros(size(DATA2)),'k')
datetick('x',19);grid on;axis tight
title('Corrente de Mare Componente Meridional [Superficie]')
ylabel('[m.s^{-1}]')

print -depsc mareS
!epstopdf mareS.eps


% Plotando Constituintes Significativas
figure;pause
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
title('Constituintes de Mare - Componente Zonal [Superficie]')

print -depsc constSu
!epstopdf constSu.eps

figure;pause
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
title('Constituintes de Mare - Componente Meridional [Superficie]')

print -depsc constSv
!epstopdf constSv.eps



pout_corrS_u=poutu; 
pout_corrS_v=poutv; 
corrS_u=u; 
corrS_v=v; 
corrS_dir=dir2; 
corrS_vel=vel2;

save corrS pout_corrS_u pout_corrS_v corrS_u corrS_v corrS_dir corrS_vel

