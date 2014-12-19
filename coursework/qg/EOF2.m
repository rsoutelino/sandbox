%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%              Cálculo de EOFs            %%%
%%%               WESTRAX 2                 %%%
%%%     Rafael Soutelino - agosto / 2006    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc

%%% LOOPING PARA LEITURA DOS DADOS

lon=[]; lat=[]; nest=[]; k=0; 

for i=[22:37 39:50];

    eval(['load ../ctd/w2cpz0',num2str(i),'.pro;']);
    eval(['data = w2cpz0',num2str(i),';'])
    p=data(:,1); 
    if data(10,8)==-9999 % tirando as estacoes que nao tem perfilagem pegasus 
        disp(['Estacao ',num2str(i),' nao tem pegasus!'])    
    else 
        disp('.')
        f=find(p==3000);
        if isempty(f)==0;
          k=k+1;
          lon=[lon data(1,6)*-1];
          lat=[lat data(1,5)];
          Usec(5:3002,k) = data(5:3002,8);
          Vsec(5:3002,k) = data(5:3002,9);
          psec(5:3002,k) = data(5:3002,1);
          nest=[nest i];
        else
          disp(['Estacao ',num2str(i),' nao atinge profundidade desejada'])
        end
    end
    eval(['clear w2cpz0',num2str(i)])    
end

clear data i f k p
 
Usec=Usec'/100;
Vsec=Vsec'/100;


psi = [Usec;Vsec];


%  dist=sw_dist(lat,lon,'km');
%  dist=[0 cumsum(dist)];

% calculando as EOFs

disp(' ')
disp('Calculando EOFs ...')
[lam_eof,F_eof,A_eof]=eoft(psi);

m=length(F_eof);

leg4=['Modo 4 [',num2str(round(lam_eof(m-3)*100)),'%]'];
leg3=['Modo 3 [',num2str(round(lam_eof(m-2)*100)),'%]'];
leg2=['Modo 2 [',num2str(round(lam_eof(m-1)*100)),'%]'];
leg1=['Modo 1 [',num2str(round(lam_eof(m)*100)),'%]'];


prof=15:10:2995;

load ../lista_pr5/F.mat
dz=10;H=2980;


Amp_md1=A_eof(end,13:10:end-5);
Amp_md2=A_eof(end-1,13:10:end-5);
Amp_md3=A_eof(end-2,13:10:end-5);
Amp_md4=A_eof(end-3,13:10:end-5);


F1=F1(1:end-100);
F2=F2(1:end-100);
F3=F3(1:end-100);
F4=F4(1:end-100);


Vsec1=mean(Vsec(:,13:10:end-5));
Usec1=mean(Usec(:,13:10:end-5));

for i=1:4

%  eval(['a1(' num2str(i) ')=1/H*sum(Amp_md' num2str(i) '.*F' num2str(i) '''*dz);'])


eval(['a1(' num2str(i) ')=1/H*sum(Amp_md1.*F' num2str(i) '''*dz);'])
eval(['a2(' num2str(i) ')=1/H*sum(Amp_md2.*F' num2str(i) '''*dz);'])
eval(['rms(' num2str(i) ')= (round((sqrt(mean((Amp_md' num2str(i) '- a1(' num2str(i) ').*F' num2str(i) ''').^2)/sqrt(mean(Amp_md' num2str(i) '.^2))))*10000))/100;'])

end


%  RMS1_4= sqrt(mean((Amp_md1 - (a1(1).*F1' + a1(2).*F2' + a1(3).*F3' + a1(4).*F4')).^2)/sqrt(mean(Amp_md1.^2)));
RMS1_3= (round((sqrt(mean((Amp_md1 - (a1(1).*F1' + a1(2).*F2' + a1(3).*F3')).^2)/sqrt(mean(Amp_md1.^2))))*10000))/100;
RMS1_3=100-RMS1_3;
RMS1_2= (round((sqrt(mean((Amp_md1 - (a1(1).*F1' + a1(2).*F2')).^2)/sqrt(mean(Amp_md1.^2))))*10000))/100;
RMS1_2=100-RMS1_2;

%  RMS2_4= sqrt(mean((Amp_md2 - (a1(1).*F1' + a1(2).*F2' + a1(3).*F3' + a1(4).*F4')).^2)/sqrt(mean(Amp_md2.^2)));
RMS2_3 = (round((sqrt(mean((Amp_md2 - (a2(1).*F1' + a2(2).*F2' + a2(3).*F3')).^2)/sqrt(mean(Amp_md2.^2))))*10000))/100;
RMS2_3=100-RMS2_3;
RMS2_2 = (round((sqrt(mean((Amp_md2 - (a2(1).*F1' + a2(2).*F2')).^2)/sqrt(mean(Amp_md2.^2))))*10000))/100;
RMS2_2=100-RMS2_2;


for i=1:3
subplot(1,3,i)
eval(['plot(Amp_md' num2str(i) ',-prof)'])
hold on;
eval(['plot(a1(' num2str(i) ')*F' num2str(i) ',-prof,''r'',''linewidth'',2)'])
plot([0 0],[-3500 0],'k')
title(['Modo ' num2str(i) ': ' num2str(rms(i)) ' %' ],'fontweight','bold')
if i==1;ylabel('Profundidade [m]');elseif i==2;xlabel('Amplitude Modal');else; legend('Modo EOF','Modo Dinamico',4);end
end

set(gcf,'color','w')

print -depsc Mods_ind
!epstopdf Mods_ind.eps

figure
set(gcf,'color','w')
subplot(1,2,1)
plot(Amp_md1,-prof)
hold on;
plot(a1(1).*F1 + a1(2).*F2,-prof,'r','linewidth',2)
plot([0 0],[-3500 0],'k')
title(['2 Modos Dinamicos --> EOF 1: ' num2str(RMS1_2) ' %' ],'fontweight','bold')
ylabel('Profundidade [m]');
xlabel('Amplitude Modal');
legend('Modo EOF','Modo Dinamico',4)

subplot(1,2,2)
plot(Amp_md1,-prof)
hold on;
plot(a1(1).*F1 + a1(2).*F2 + a1(3).*F3,-prof,'r','linewidth',2)
plot([0 0],[-3500 0],'k')
title(['3 Modos Dinamicos --> EOF 1: ' num2str(RMS1_3) ' %' ],'fontweight','bold')
ylabel('Profundidade [m]');
xlabel('Amplitude Modal');
legend('Modo EOF','Modo Dinamico',4)

print -depsc Mods1_comb
!epstopdf Mods1_comb.eps

figure
set(gcf,'color','w')
subplot(1,2,1)
plot(Amp_md2,-prof)
hold on;
plot(a2(1).*F1 + a2(2).*F2,-prof,'r','linewidth',2)
plot([0 0],[-3500 0],'k')
title(['2 Modos Dinamicos --> EOF 2: ' num2str(RMS2_2) ' %' ],'fontweight','bold')
ylabel('Profundidade [m]');
xlabel('Amplitude Modal');
legend('Modo EOF','Modo Dinamico',4)

subplot(1,2,2)
plot(Amp_md2,-prof)
hold on;
plot(a2(1).*F1 + a2(2).*F2 + a2(3).*F3,-prof,'r','linewidth',2)
plot([0 0],[-3500 0],'k')
title(['3 Modos Dinamicos --> EOF 2: ' num2str(RMS2_3) ' %' ],'fontweight','bold')
ylabel('Profundidade [m]');
xlabel('Amplitude Modal');
legend('Modo EOF','Modo Dinamico',4)

print -depsc Mods2_comb
!epstopdf Mods2_comb.eps

%%%%%%%%%%%%%%%%%%
stop
% MONTE CARLO

[numlam,lamV,mlamc]=montecarlo(psi,lam_eof,1000);

figure(4)
set(gcf,'color','w')
print -depsc montcarl1
!epstopdf montcarl1.eps

figure(5)
set(gcf,'color','w')
print -depsc montcarl2
!epstopdf montcarl2.eps
