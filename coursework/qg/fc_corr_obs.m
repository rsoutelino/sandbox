%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%        PROGRAMA PARA CALCULO DE         %%%
%%%     FUNCAO DE CORRENTE OBSERVADA -    %%%
%%%               WESTRAX 2                 %%%
%%%     Rafael Soutelino - agosto / 2006     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc
tic

%%% LOOPING PARA LEITURA DOS DADOS

lon=[]; lat=[]; lonest=[]; latest=[];
lonpeg=[]; latpeg=[];
u50=[]; v50=[]; u500=[]; v500=[]; u2000=[]; v2000=[];

disp('Primeira Leitura')
disp(' ')

for i=[22:37 39:50]

    eval(['load ../ctd/w2cpz0',num2str(i),'.pro;']);
    eval(['data = w2cpz0',num2str(i),';'])
    lonest=[lonest ; data(1,6)*-1];
    latest=[latest ; data(1,5)];
    if data(10,8)==-9999 % tirando as estacoes que nao tem perfilagem pegasus
        x=-9999; % só pra constar no if
    else
        p=data(:,1);
        lonpeg=[lonpeg ; data(1,6)*-1];
        latpeg=[latpeg ; data(1,5)];
        if max(p) > 2000
           u2000=[u2000 ; data(2000,8)];
           v2000=[v2000 ; data(2000,9)];         
           u500=[u500 ; data(500,8)];
           v500=[v500 ; data(500,9)];          
           u50=[u50 ; data(50,8)];
           v50=[v50 ; data(50,9)];
           lon=[lon ; data(1,6)*-1];
           lat=[lat ; data(1,5)];
        else 
            disp('prof max < 2000');
        end
        
    end
    eval(['clear w2cpz0',num2str(i)])    
end


u50 = u50/100;
u500 = u500/100;
u2000 = u2000/100;
v50 = v50/100;
v500 = v500/100;
v2000 = v2000/100;


clear data p i x 

%%% montando imagens para a condicao de contorno

loni=[]; lati=[]; 
u50i=[]; v50i=[]; 
u500i=[]; v500i=[];
u2000i=[]; v2000i=[];

for i=1:length(lon)
   [loni2,lati2,u50i2,v50i2] = mirror(lon(i),lat(i),u50(i),v50(i));
   loni=[loni ; loni2]; lati=[lati ; lati2];
   u50i=[u50i ; u50i2]; v50i=[v50i ; v50i2];
   [loni2,lati2,u500i2,v500i2] = mirror(lon(i),lat(i),u500(i),v500(i));
   u500i=[u500i ; u500i2]; v500i=[v500i ; v500i2];
   [loni2,lati2,u2000i2,v2000i2] = mirror(lon(i),lat(i),u2000(i),v2000(i));
   u2000i=[u2000i ; u2000i2]; v2000i=[v2000i ; v2000i2];
end   

clear u50i2 v50i2 u500i2 v500i2 u2000i2 v2000i2 loni2 lati2

lon=[lon ; loni]; lat=[lat ; lati];
u50=[u50 ; u50i]; v50=[v50 ; v50i];
u500=[u500 ; u500i]; v500=[v500 ; v500i];
u2000=[u2000 ; u2000i]; v2000=[v2000 ; v2000i];

LON=lon;
LAT=lat;

% limites para plotagem

lonlim=[-53 -44]; latlim=[0 9]; % limites de silveira_etal2000

m_proj('mercator','long',[lonlim(1) lonlim(2)],'lat',[latlim(1) latlim(2)],'on');


% montando reta L200, para usar como espelho

xL200=lonlim(1):0.25:lonlim(2);
yL200=xL200*(-1) - 45.5;

[xL200,yL200]=m_ll2xy(xL200,yL200,'clip','patch');


% criando uma grade regular

xi=-53:0.25:-44; yi=0:0.25:9;
[xg,yg]=meshgrid(xi,yi);
[lj,li]=size(xg);


% ANALISE OBJETIVA -----------------------------------------

xg2=reshape(xg,1,lj*li);
yg2=reshape(yg,1,lj*li);

u50 = u50*8.9993e-06; % fator de conversao para grau/s (9e-6)
v50 = v50*8.9993e-06;
u500 = u500*8.9993e-06; % fator de conversao para grau/s (9e-6)
v500 = v500*8.9993e-06;
u2000 = u2000*8.9993e-06; % fator de conversao para grau/s (9e-6)
v2000 = v2000*8.9993e-06;


lc=4;
E=0.09^2;

[psiob50]=vectoa(xg2,yg2,lon,lat,u50,v50,lc,E,0);
[psiob500]=vectoa(xg2,yg2,lon,lat,u500,v500,lc,E,0);
[psiob2000]=vectoa(xg2,yg2,lon,lat,u2000,v2000,lc,E,0);

psiob50 = psiob50*(111120^2);
psiob500 = psiob500*(111120^2);
psiob2000 = psiob2000*(111120^2);

psiob50=reshape(psiob50,li,lj);
psiob500=reshape(psiob500,li,lj);
psiob2000=reshape(psiob2000,li,lj);
xg=reshape(xg2,li,lj);
yg=reshape(yg2,li,lj);


% calculando vetores de velocidade, a partir de psi

[V50,U50]=gradient(psiob50); U50=-U50;
[V500,U500]=gradient(psiob500); U500=-U500;
[V2000,U2000]=gradient(psiob2000); U2000=-U2000;



% PLOTANDO ==========================================================================================

lpsi=[-254:5:254];

load iso200.mat; 

%  % pegasus ---------------------------------------------------
%  
%  figure
%  m_contourf(xg,yg,psiob50,lpsi); shading flat; hold on; 
%  cc=colorbar;
%  m_quiver(xg,yg,U50,V50,2,'k'); hold on
%  %  fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]);
%  fill([xL200(1); xL200' ; xL200(end)],[yL200(end);yL200'; yL200(end)],[.7 .7 .7]);
%  m_usercoast('costa.mat','patch',[0 0 0])
%  m_plot(lonpeg,latpeg,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4);
%  m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',10);
%  title('Funcao de Corrente Observada \Psi: [m^2 s^{-1}] - 50 m','fontsize',12,'fontweigh','bold')
%  set(gcf,'color','w')
%  
%  
%  print -depsc fc_corr_obs_50m.eps
%  !epstopdf fc_corr_obs_50m.eps

%  lpsi=[-49:5:49];
%  
%  figure
%  m_contourf(xg,yg,psiob500,lpsi); shading flat; hold on; 
%  cc=colorbar;
%  m_quiver(xg,yg,U500,V500,2,'k'); hold on
%  fill([xL200(1); xL200' ; xL200(end)],[yL200(end);yL200'; yL200(end)],[.7 .7 .7]);
%  %  fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]);
%  m_usercoast('costa.mat','patch',[0 0 0])
%  m_plot(lonpeg,latpeg,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4);
%  m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',10);
%  title('Funcao de Corrente Observada \Psi: [m^2 s^{-1}] - 500 m','fontsize',12,'fontweigh','bold')
%  set(gcf,'color','w')
%  
%  print -depsc fc_corr_obs_500m.eps
%  !epstopdf fc_corr_obs_500m.eps
%  
%  lpsi=[-52:5:43];
%  
%  figure
%  m_contourf(xg,yg,psiob2000,lpsi); shading flat; hold on; 
%  cc=colorbar;
%  m_quiver(xg,yg,U2000,V2000,2,'k'); hold on
%  fill([xL200(1); xL200' ; xL200(end)],[yL200(end);yL200'; yL200(end)],[.7 .7 .7]);
%  %  fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]);
%  m_usercoast('costa.mat','patch',[0 0 0])
%  m_plot(lonpeg,latpeg,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4);
%  m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',10);
%  title('Funcao de Corrente Observada \Psi: [m^2 s^{-1}] - 2000 m','fontsize',12,'fontweigh','bold')
%  set(gcf,'color','w')
%  
%  print -depsc fc_corr_obs_2000m.eps
%  !epstopdf fc_corr_obs_2000m.eps
%  


%%% CALCULANDO A VORTICIDADE RELATIVA
stop
[zeta50]=vortoa(xg2,yg2,lon,lat,u50,v50,lc,E);
zeta50=reshape(zeta50,li,lj);
lzeta=[min(min(zeta50)):0.1:max(max(zeta50))];

u50 = u50/8.9993e-06; % fator de conversao para grau/s (9e-6)
v50 = v50/8.9993e-06;
u500 = u500/8.9993e-06; % fator de conversao para grau/s (9e-6)
v500 = v500/8.9993e-06;
u2000 = u2000/8.9993e-06; % fator de conversao para grau/s (9e-6)
v2000 = v2000/8.9993e-06;

figure

subplot(121)

m_contourf(xg,yg,zeta50,lzeta); shading flat; hold on; 
cc=colorbar('horiz');
fill([xL200(1); xL200' ; xL200(end)],[yL200(end);yL200'; yL200(end)],[.7 .7 .7]);
%  fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]);
m_usercoast('costa.mat','patch',[0 0 0])
%  m_plot(lonpeg,latpeg,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4);
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',10);
title('Vorticidade Relativa - 50 m','fontsize',12,'fontweigh','bold')
set(gcf,'color','w')



%%% MAPEANDO NUMERO DE ROSSBY

f0=sw_f(5);
Ro=abs(zeta50)/f0;
lRo=[min(min(Ro)):5000:max(max(Ro))];

subplot(122)

m_contourf(xg,yg,Ro,lRo); shading flat; hold on; 
cc=colorbar('horiz');
fill([xL200(1); xL200' ; xL200(end)],[yL200(end);yL200'; yL200(end)],[.7 .7 .7]);
%  fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]);
m_usercoast('costa.mat','patch',[0 0 0])
%  m_plot(lonpeg,latpeg,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4);
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',10);
title('Numero de Rossby Pontual','fontsize',12,'fontweigh','bold')
set(gcf,'color','w')

print -depsc zeta_Ro.eps
!epstopdf zeta_Ro.eps 


% valor medio para o Numero de Rossby

Romean=triu(Ro); 
f=find(Romean==0);
Romean(f)==NaN;
Romean=mean(mean(Romean));

%%% MAPEANDO AS AMPLITUDES MODAIS

load ../lista_pr2/F.mat;

% lendo novamente os dados pois precisamos de todos os perfis

u=[]; v=[]; lon=[]; lat=[];

disp(' ')
disp('Segunda Leitura')
disp(' ')


for i=[22:37 39:50]

    eval(['load ../ctd/w2cpz0',num2str(i),'.pro;']);
    eval(['data = w2cpz0',num2str(i),';']);
    if data(10,8)==-9999 % tirando as estacoes que nao tem perfilagem pegasus
        x=-9999; % só pra constar no if
    else
        p=data(:,1);
    
        if max(p) > 2000
           u=[u  data(15:2000,8)];
           v=[v  data(15:2000,9)];      
           lon=[lon ; data(1,6)*-1];
           lat=[lat ; data(1,5)];
        else 
            disp('prof max < 2000');
        end
        
    end
    eval(['clear w2cpz0',num2str(i)])    
end


% criando eixo de pressão

prof=15:2000;

% decimando todos de 10 em 10

U=[]; V=[]; p=[];

for i=10:10:(length(prof)-1)
    U(i/10,:) = u(i,:);
    V(i/10,:) = v(i,:);
    p(i/10,:) = prof(i);
end

f=near(po,max(p)); f=f-1;

F0=F0(1:f); F1=F1(1:f); F2=F2(1:f); F3=F3(1:f); F4=F4(1:f); F5=F5(1:f); 

% passando para o SI

U=U/100; V=V/100;

% calculando as aplitudes modais

for i=0:5
   for k=1:20
       eval(['Aumod',num2str(i),'(1,k) = 10/4000*sum(U(:,k).*(F',num2str(i),'));'])
       eval(['Avmod',num2str(i),'(1,k) = 10/4000*sum(V(:,k).*(F',num2str(i),'));'])
   end
end

for i=0:5
   for k=1:20
       eval(['Umod',num2str(i),'(:,k) = Aumod',num2str(i),'(1,k).*F',num2str(i),';'])
       eval(['Vmod',num2str(i),'(:,k) = Avmod',num2str(i),'(1,k).*F',num2str(i),';'])
   end
end

f50 = near(p,50);
f500 = near(p,500);
f2000 = near(p,2000);

for i=0:5
    eval(['U',num2str(i),'_50 = Umod',num2str(i),'(f50,:);']);
    eval(['U',num2str(i),'_500 = Umod',num2str(i),'(f500,:);']);
    eval(['U',num2str(i),'_2000 = Umod',num2str(i),'(f2000,:);']);
    eval(['V',num2str(i),'_50 = Vmod',num2str(i),'(f50,:);']);
    eval(['V',num2str(i),'_500 = Vmod',num2str(i),'(f500,:);']);
    eval(['V',num2str(i),'_2000 = Vmod',num2str(i),'(f2000,:);']);
end


% criando as imagens

for i=0:5
   eval(['U',num2str(i),'_50IM = [];']);
   eval(['U',num2str(i),'_500IM = [];']);
   eval(['U',num2str(i),'_2000IM = [];']);
   eval(['V',num2str(i),'_50IM = [];']);
   eval(['V',num2str(i),'_500IM = [];']);
   eval(['V',num2str(i),'_2000IM = [];']);
end

for i=0:5
   for k=1:20
     eval(['[loni,lati,U',num2str(i),'_50IM2,V',num2str(i),'_50IM2] = mirror(lon(k),lat(k),U',num2str(i),'_50(k),V',num2str(i),'_50(k));']); 
     eval(['[loni,lati,U',num2str(i),'_500IM2,V',num2str(i),'_500IM2] = mirror(lon(k),lat(k),U',num2str(i),'_500(k),V',num2str(i),'_500(k));']);
     eval(['[loni,lati,U',num2str(i),'_2000IM2,V',num2str(i),'_2000IM2] = mirror(lon(k),lat(k),U',num2str(i),'_2000(k),V',num2str(i),'_2000(k));']);

     eval(['U',num2str(i),'_50IM = [U',num2str(i),'_50IM U',num2str(i),'_50IM2];']);
     eval(['V',num2str(i),'_50IM = [V',num2str(i),'_50IM V',num2str(i),'_50IM2];']);

     eval(['U',num2str(i),'_500IM = [U',num2str(i),'_500IM U',num2str(i),'_500IM2];']);
     eval(['V',num2str(i),'_500IM = [V',num2str(i),'_500IM V',num2str(i),'_500IM2];']);

     eval(['U',num2str(i),'_2000IM = [U',num2str(i),'_2000IM U',num2str(i),'_2000IM2];']);
     eval(['V',num2str(i),'_2000IM = [V',num2str(i),'_2000IM V',num2str(i),'_2000IM2];']);
   end   
end

for i=0:5
    eval(['U',num2str(i),'_50 = [U',num2str(i),'_50  U',num2str(i),'_50IM];']);
    eval(['V',num2str(i),'_50 = [V',num2str(i),'_50  V',num2str(i),'_50IM];']);

    eval(['U',num2str(i),'_500 = [U',num2str(i),'_500  U',num2str(i),'_500IM];']);
    eval(['V',num2str(i),'_500 = [V',num2str(i),'_500  V',num2str(i),'_500IM];']);

    eval(['U',num2str(i),'_2000 = [U',num2str(i),'_2000  U',num2str(i),'_2000IM];']);
    eval(['V',num2str(i),'_2000 = [V',num2str(i),'_2000  V',num2str(i),'_2000IM];']);
end


% fazendo analise objetiva  (usar os vetores LAT LON)

for i=0:5
    eval(['[PSI_mod',num2str(i),'_50] = vectoa(xg2,yg2,LON,LAT , U',num2str(i),'_50 , V',num2str(i),'_50 ,lc,E,0);']);
    eval(['PSI_mod',num2str(i),'_50 = reshape(PSI_mod',num2str(i),'_50 ,li,lj);']);

    eval(['[PSI_mod',num2str(i),'_500] = vectoa(xg2,yg2,LON,LAT , U',num2str(i),'_500 , V',num2str(i),'_500 ,lc,E,0);']);
    eval(['PSI_mod',num2str(i),'_500 = reshape(PSI_mod',num2str(i),'_500 ,li,lj);']);

    eval(['[PSI_mod',num2str(i),'_2000] = vectoa(xg2,yg2,LON,LAT , U',num2str(i),'_2000 , V',num2str(i),'_2000 ,lc,E,0);']);
    eval(['PSI_mod',num2str(i),'_2000 = reshape(PSI_mod',num2str(i),'_2000 ,li,lj);']);

end

% calculando u e v das amplitudes modais
    
for i=0:5

    eval(['[V',num2str(i),'_50 , U',num2str(i),'_50]  = gradient(PSI_mod',num2str(i),'_50);']);
    eval(['U',num2str(i),'_50 = - U',num2str(i),'_50  ;']);
    
    eval(['[V',num2str(i),'_500 , U',num2str(i),'_500]  = gradient(PSI_mod',num2str(i),'_500);']);
    eval(['U',num2str(i),'_500 = - U',num2str(i),'_500  ;']);

    eval(['[V',num2str(i),'_2000 , U',num2str(i),'_2000]  = gradient(PSI_mod',num2str(i),'_2000);']);
    eval(['U',num2str(i),'_2000 = - U',num2str(i),'_2000  ;']);


end 


% plotando os modos dinamicos

figure % 50 m de profundidade

subplot (221) % pegasus ---------------------------------------------------
lpsi=-2.49:0.1:2.5;

m_contourf(xg,yg,psiob50/100,lpsi); shading flat; hold on; 
%  cc=colorbar;
m_quiver(xg,yg,U50,V50,2,'k'); hold on
%  fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]);
fill([xL200(1); xL200' ; xL200(end)],[yL200(end);yL200'; yL200(end)],[.7 .7 .7]);
m_usercoast('costa.mat','patch',[0 0 0])
%  m_plot(lonpeg,latpeg,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4);
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',10);
title('PEGASUS \Psi: [m^2 s^{-1}] - 50 m','fontsize',10,'fontweigh','bold')
set(gcf,'color','w')
set(findobj('tag','m_grid_color'),'facecolor','none')

subplot (222) % primeiro modo baroclinico + modo barotropico --------------

m_contourf(xg,yg,PSI_mod0_50+PSI_mod1_50,lpsi); shading flat; hold on; 
%  cc=colorbar;
m_quiver(xg,yg,U0_50+U1_50,V0_50+V1_50,2,'k'); hold on
%  fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]);
fill([xL200(1); xL200' ; xL200(end)],[yL200(end);yL200'; yL200(end)],[.7 .7 .7]);
m_usercoast('costa.mat','patch',[0 0 0])
%  m_plot(lonpeg,latpeg,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4);
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',10);
title('\Psi: 2 modos [m^2 s^{-1}] - 50 m','fontsize',10,'fontweigh','bold')
set(gcf,'color','w')
set(findobj('tag','m_grid_color'),'facecolor','none')

subplot (223) % 3 modos --------------

m_contourf(xg,yg,PSI_mod0_50+PSI_mod1_50+PSI_mod2_50,lpsi); shading flat; hold on; 
%  cc=colorbar;
m_quiver(xg,yg,U0_50+U1_50+U2_50,V0_50+V1_50+V2_50,2,'k'); hold on
%  fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]);
fill([xL200(1); xL200' ; xL200(end)],[yL200(end);yL200'; yL200(end)],[.7 .7 .7]);
m_usercoast('costa.mat','patch',[0 0 0])
%  m_plot(lonpeg,latpeg,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4);
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',10);
title('\Psi: 3 modos [m^2 s^{-1}] - 50 m','fontsize',10,'fontweigh','bold')
set(gcf,'color','w')
set(findobj('tag','m_grid_color'),'facecolor','none')

subplot (224) % 6 modos --------------

m_contourf(xg,yg,PSI_mod0_50+PSI_mod1_50+PSI_mod2_50+PSI_mod3_50+PSI_mod4_50+PSI_mod5_50,lpsi); shading flat; hold on; 
%  cc=colorbar;
m_quiver(xg,yg,U0_50+U1_50+U2_50+U3_50+U4_50+U5_50,V0_50+V1_50+V2_50+V3_50+V4_50+V5_50,2,'k'); hold on
%  fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]);
fill([xL200(1); xL200' ; xL200(end)],[yL200(end);yL200'; yL200(end)],[.7 .7 .7]);
m_usercoast('costa.mat','patch',[0 0 0])
%  m_plot(lonpeg,latpeg,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4);
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',10);
title('\Psi: 6 modos [m^2 s^{-1}] - 50 m','fontsize',10,'fontweigh','bold')
set(gcf,'color','w')
set(findobj('tag','m_grid_color'),'facecolor','none')

print -depsc psi_50m.eps
!epstopdf psi_50m.eps 


figure % 500 m de profundidade

subplot (221) % pegasus ---------------------------------------------------
lpsi=-.6:0.05:.6;

m_contourf(xg,yg,psiob500/100,lpsi); shading flat; hold on; 
%  cc=colorbar;
m_quiver(xg,yg,U500,V500,2,'k'); hold on
%  fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]);
fill([xL200(1); xL200' ; xL200(end)],[yL200(end);yL200'; yL200(end)],[.7 .7 .7]);
m_usercoast('costa.mat','patch',[0 0 0])
%  m_plot(lonpeg,latpeg,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4);
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',10);
title('PEGASUS \Psi: [m^2 s^{-1}] - 500 m','fontsize',10,'fontweigh','bold')
set(gcf,'color','w')
set(findobj('tag','m_grid_color'),'facecolor','none')

subplot (222) % primeiro modo baroclinico + modo barotropico --------------

m_contourf(xg,yg,PSI_mod0_500+PSI_mod1_500,lpsi); shading flat; hold on; 
%  cc=colorbar;
m_quiver(xg,yg,U0_500+U1_500,V0_500+V1_500,2,'k'); hold on
%  fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]);
fill([xL200(1); xL200' ; xL200(end)],[yL200(end);yL200'; yL200(end)],[.7 .7 .7]);
m_usercoast('costa.mat','patch',[0 0 0])
%  m_plot(lonpeg,latpeg,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4);
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',10);
title('\Psi: 2 modos [m^2 s^{-1}] - 500 m','fontsize',10,'fontweigh','bold')
set(gcf,'color','w')
set(findobj('tag','m_grid_color'),'facecolor','none')

subplot (223) % 3 modos  --------------

m_contourf(xg,yg,PSI_mod0_500+PSI_mod1_500+PSI_mod2_500,lpsi); shading flat; hold on; 
%  cc=colorbar;
m_quiver(xg,yg,U0_500+U1_500+U2_500,V0_500+V1_500+V2_500,2,'k'); hold on
%  fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]);
fill([xL200(1); xL200' ; xL200(end)],[yL200(end);yL200'; yL200(end)],[.7 .7 .7]);
m_usercoast('costa.mat','patch',[0 0 0])
%  m_plot(lonpeg,latpeg,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4);
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',10);
title('\Psi: 3 modos [m^2 s^{-1}] - 500 m','fontsize',10,'fontweigh','bold')
set(gcf,'color','w')
set(findobj('tag','m_grid_color'),'facecolor','none')

subplot (224) % 6 modos --------------

m_contourf(xg,yg,PSI_mod0_500+PSI_mod1_500+PSI_mod2_500+PSI_mod3_500+PSI_mod4_500+PSI_mod5_500,lpsi); shading flat; hold on; 
%  cc=colorbar;
m_quiver(xg,yg,U0_500+U1_500+U2_500+U3_500+U4_500+U5_500,V0_500+V1_500+V2_500+V3_500+V4_500+V5_500,2,'k'); hold on
%  fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]);
fill([xL200(1); xL200' ; xL200(end)],[yL200(end);yL200'; yL200(end)],[.7 .7 .7]);
m_usercoast('costa.mat','patch',[0 0 0])
%  m_plot(lonpeg,latpeg,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4);
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',10);
title('\Psi: 6 modos [m^2 s^{-1}] - 500 m','fontsize',10,'fontweigh','bold')
set(gcf,'color','w')
set(findobj('tag','m_grid_color'),'facecolor','none')

print -depsc psi_500m.eps
!epstopdf psi_500m.eps 


figure % 2000 m de profundidade

subplot (221) % pegasus ---------------------------------------------------
lpsi=-.51:.05:.51;

m_contourf(xg,yg,psiob2000/100,lpsi); shading flat; hold on; 
%  cc=colorbar;
m_quiver(xg,yg,U2000,V2000,2,'k'); hold on
%  fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]);
fill([xL200(1); xL200' ; xL200(end)],[yL200(end);yL200'; yL200(end)],[.7 .7 .7]);
m_usercoast('costa.mat','patch',[0 0 0])
%  m_plot(lonpeg,latpeg,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4);
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',10);
title('PEGASUS \Psi: [m^2 s^{-1}] - 2000 m','fontsize',10,'fontweigh','bold')
set(gcf,'color','w')
set(findobj('tag','m_grid_color'),'facecolor','none')

subplot (222) % primeiro modo baroclinico + modo barotropico --------------

m_contourf(xg,yg,PSI_mod0_2000+PSI_mod1_2000,lpsi); shading flat; hold on; 
%  cc=colorbar;
m_quiver(xg,yg,U0_2000+U1_2000,V0_2000+V1_2000,2,'k'); hold on
%  fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]);
fill([xL200(1); xL200' ; xL200(end)],[yL200(end);yL200'; yL200(end)],[.7 .7 .7]);
m_usercoast('costa.mat','patch',[0 0 0])
%  %  m_plot(lonpeg,latpeg,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4);
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',10);
title('\Psi: 2 modos [m^2 s^{-1}] - 2000 m','fontsize',10,'fontweigh','bold')
set(gcf,'color','w')
set(findobj('tag','m_grid_color'),'facecolor','none')

subplot (223) % 3 modos --------------

m_contourf(xg,yg,PSI_mod0_2000+PSI_mod1_2000+PSI_mod2_2000,lpsi); shading flat; hold on; 
%  cc=colorbar;
m_quiver(xg,yg,U0_2000+U1_2000+U2_2000,V0_2000+V1_2000+V2_2000,2,'k'); hold on
%  fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]);
fill([xL200(1); xL200' ; xL200(end)],[yL200(end);yL200'; yL200(end)],[.7 .7 .7]);
m_usercoast('costa.mat','patch',[0 0 0])
%  m_plot(lonpeg,latpeg,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4);
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',10);
title('\Psi: 3 modos [m^2 s^{-1}] - 2000 m','fontsize',10,'fontweigh','bold')
set(gcf,'color','w')
set(findobj('tag','m_grid_color'),'facecolor','none')

subplot (224) % 6 modos --------------

m_contourf(xg,yg,PSI_mod0_2000+PSI_mod1_2000+PSI_mod2_2000+PSI_mod3_2000+PSI_mod4_2000+PSI_mod5_2000,lpsi); shading flat; hold on; 
%  cc=colorbar;
m_quiver(xg,yg,U0_2000+U1_2000+U2_2000+U3_2000+U4_2000+U5_2000,V0_2000+V1_2000+V2_2000+V3_2000+V4_2000+V5_2000,2,'k'); hold on
%  fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]);
fill([xL200(1); xL200' ; xL200(end)],[yL200(end);yL200'; yL200(end)],[.7 .7 .7]);
m_usercoast('costa.mat','patch',[0 0 0])
%  m_plot(lonpeg,latpeg,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4);
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',10);
title('\Psi: 6 modos [m^2 s^{-1}] - 2000 m','fontsize',10,'fontweigh','bold')
set(gcf,'color','w')
set(findobj('tag','m_grid_color'),'facecolor','none')

print -depsc psi_2000m.eps
!epstopdf psi_2000m.eps 

clc





toc