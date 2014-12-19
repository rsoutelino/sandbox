%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%        PROGRAMA PARA CALCULO DOS        %%%
%%%   PERFIS MEDIOS PARA POSTERIOR CALCULO  %%%
%%%         DOS MODOS DINAMICOS -           %%%
%%%               WESTRAX 2                 %%%
%%%     Rafael Soutelino - julho / 2006     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc

%%% LOOPING PARA LEITURA DOS DADOS

lon=[]; lat=[]; nest=[]; cont=1; Tsec=[]; Ssec=[];

for i=[35 36 37 39 40 41 42] 

    eval(['load ../ctd/w2cpz0',num2str(i),'.pro;']);
    eval(['data = w2cpz0',num2str(i),';'])
    
    p=data(:,1);
    f=find(p==4000);
    if isempty(f)==0;
       lon=[lon data(1,6)*-1];
       lat=[lat data(1,5)];
       eval(['p',num2str(i),' = data(5:4002,1);'])
       eval(['T',num2str(i),' = data(5:4002,3);'])
       eval(['S',num2str(i),' = data(5:4002,5);'])
       nest=[nest i];
    end

    eval(['clear w2cpz0',num2str(i)])    
end

clear data   

Tsec=[T35 T36 T39 T40 T41];
Ssec=[S35 S36 S39 S40 S41];

Tm=mean(Tsec'); Tm=Tm';
Sm=mean(Ssec'); Sm=Sm';
p=-p35;

%%% calculando frequencia de estratificacao

N=sw_bfrq(Sm,Tm,-p,5); N=[N;NaN];
N=weim(41,'hann',N);

%%% decimando N2 e p,t,s
 Dm = sw_dens0(Sm,Tm);
Dm2=Dm; Sm2=Sm; Tm2=Tm; p2=p';
Sm=[]; Tm=[]; Dm=[]; N2=[]; p=[];

for i=10:10:length(N)
    N2(i/10) = N(i);
    p(i/10) = p2(i);
    Sm(i/10) = Sm2(i);
    Tm(i/10) = Tm2(i);
    Dm(i/10) = Dm2(i);
end

%  save N2.mat N2 p2

figure

subplot(131)
plot(Tm,p,'r'); hold on
axis([0 30 -3000 0])
title('Perfil Medio de Temperatura','fontweight','bold')
xlabel('T (\circ C)')
ylabel('Profundidade (m)')
set(gca,'plotboxaspectratio',[.5 1 .5])
set(gcf,'color','w')

subplot(132)
plot(Sm,p,'g'); hold on
axis([34 37 -3000 0])
title('Perfil Medio de Salinidade','fontweight','bold')
set(gca,'plotboxaspectratio',[.5 1 .5])
set(gcf,'color','w')

subplot(133)
plot(N2,p); hold on
axis([-1e-4 5e-4 -3000 0])
title('Frequencia de Estratificacao','fontweight','bold')
xlabel('cph')
set(gca,'plotboxaspectratio',[.5 1 .5])
set(gcf,'color','w')

clc







