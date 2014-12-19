%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%            Modelo de Camadas            %%%
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
        f=find(p==1000);
        if isempty(f)==0;
          k=k+1;
          lon=[lon data(1,6)*-1];
          lat=[lat data(1,5)];
          Usec(5:1002,k) = data(5:1002,8);
          Vsec(5:1002,k) = data(5:1002,9);
          psec(5:1002,k) = data(5:1002,1);
          nest=[nest i];
        else
          disp(['Estacao ',num2str(i),' nao atinge profundidade desejada'])
        end
    end
    eval(['clear w2cpz0',num2str(i)])    
end

clear data i f  p

%%% dando os valores dos raios de deformacao

Rd_0=15583.3942*1e3;
Rd_1=196.132203*1e3;
Rd_2=124.124119*1e3;
Rd_3=76.5253272*1e3;
Rd_4=56.8374664*1e3;
Rd_5=44.8154785*1e3;

load F.mat

f0 = sw_f(5);
H = 4000; 
g = 9.8;
p=psec(:,1);

met = 2;% input('Escolha o metodo de calibracao (1=prioriza efeitos de Ekman; 2=prioriza efeitos nao-lineares): ')

%%% aplicando primeiro metodo de calibracao (efeitos de ekaman)
if met==1

   e = ( Rd_1^2 * f0^2 * H * (F1(1)^2+1)^2 ) / (H^2 * F1(1)^2)
   H1 = H / (F1(1)^2 +1) 
%%% aplicando segundo metodo de calibracao (efeitos nao-lineares)
else
   ksi = 1/H.*sum(F1.^3).*10;
   H1 = ( H/4 * (sqrt(ksi^2+4) - ksi)^2) / ( 1/4 * ((sqrt(ksi^2+4) - ksi))^2 + 1 )
   e = ( Rd_1^2 * f0^2 * H * (F1(1)^2+1)^2 ) / (H1*(H-H1))
end
%%% calculando os modos

H2 = H - H1;
delta = H1/H2;

f1 = sqrt(1/delta);
f2 = -sqrt(delta);
  
%  figure
%  set(gcf,'color','w')
%  plot([1 1],[-4000 0],'r');hold on
%  plot([f1 f1 f2 f2 f2],[0 -H1 -H1 -H2 -H]); 
%  plot([0 0],[-4000 0],'k--')
%  axis([-4 5 -4000 0])
%  legend('Barotropico','Baroclinico',0)
%  legend boxoff
%  ylabel('Profundidade (m)')
%  xlabel('Modos')
%  tit=['Metodo ',num2str(met),' de calibracao'];
%  title(tit,'fontsize',12,'fontweight','bold')
%  set(gca,'plotboxaspectratio',[0.5 1 0.5])
%  
%  eval(['print -depsc modos_met',num2str(met),'.eps'])
%  eval(['!epstopdf modos_met',num2str(met),'.eps'])
%  eval(['!rm -rf modos_met',num2str(met),'.eps'])

if met==1, break,end

%%% Mapeando funcao de corrente para as camadas

%%% calculando amplitudes modais
U1 = []; U2 = []; V1 = []; V2 = [];
% achando indice para H1
m=near(p,H1,1);

for i=1:k
    ampU1 = 1/H * sum(Usec(1:m,i))*f1*10;
    ampU2 = 1/H * sum(Usec(m:end,i))*f2*10;
    U1 = [U1 ampU1];
    U2 = [U2 ampU2];
    
    ampV1 = 1/H * sum(Vsec(1:m,i))*f1*10;
    ampV2 = 1/H * sum(Vsec(m:end,i))*f2*10;
    V1 = [V1 ampV1];
    V2 = [V2 ampV2];
end

u1 = U1.*f1;
u2 = U2.*f2; 
v1 = V1.*f1;
v2 = V2.*f2; 

% limites para plotagem

lonlim=[-53 -44]; latlim=[0 9]; % limites de silveira_etal2000

m_proj('mercator','long',[lonlim(1) lonlim(2)],'lat',[latlim(1) latlim(2)],'on');


% montando reta L200, para usar como espelho

xL200=lonlim(1):0.25:lonlim(2);
yL200=xL200*(-1) - 45.5;

[xL200,yL200]=m_ll2xy(xL200,yL200,'clip','patch');

% calculando as imagens
for i=1:length(lon)
  [lon2(i),lat2(i),u1_2(i),v1_2(i)]=mirror(lon(i),lat(i),u1(i),v1(i));
  [lon2(i),lat2(i),u2_2(i),v2_2(i)]=mirror(lon(i),lat(i),u2(i),v2(i));
end
% incorporando as imagens
lonm=[lon lon2];
latm=[lat lat2];
u1m=[u1  u1_2];
u2m=[u2  u2_2];
v1m=[v1  v1_2];
v2m=[v2  v2_2];
% criando uma grade regular

xi=-53:0.25:-44; yi=0:0.25:9;
[xg,yg]=meshgrid(xi,yi);
[lj,li]=size(xg);

% ANALISE OBJETIVA -----------------------------------------

xg2=reshape(xg,1,lj*li);
yg2=reshape(yg,1,lj*li);

lc=4;
E=0.09^2;

[psi1]=vectoa(xg2,yg2,lonm,latm,u1m,v1m,lc,E,0);
[psi2]=vectoa(xg2,yg2,lonm,latm,u2m,v2m,lc,E,0);

psi1=reshape(psi1,li,lj);
psi2=reshape(psi2,li,lj);
xg=reshape(xg2,li,lj);
yg=reshape(yg2,li,lj);

% calculando vetores de velocidade, a partir de psi

[v1psi,u1psi]=gradient(psi1); u1psi=-u1psi;
[v2psi,u2psi]=gradient(psi2); u2psi=-u2psi;

% PLOTANDO ==========================================================================================

lpsi1=[min(min(psi1)):100:max(max(psi1))];

load iso200.mat; 

% psi na camada superior ---------------------------------------------------
figure(1)
set(gcf,...
        'Color',[1 1 1],...
        'InvertHardcopy','on',...
        'PaperUnits','inches',...
        'Units','inches',...
        'PaperOrientation','portrait',...
        'PaperPosition',[0 0 8.5 11],...
        'PaperPositionMode','manual',...
        'PaperType','usletter',...
        'Position',[.2 .2 8.5 11],...
        'ShareColors','off',...
        'Clipping','on');

subplot(321)

m_contourf(xg,yg,psi1,lpsi1); shading flat; hold on; 
cc=colorbar;
m_quiver(xg,yg,u1psi,v1psi,2,'k'); hold on
%  fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]);
fill([xL200(1); xL200' ; xL200(end)],[yL200(end);yL200'; yL200(end)],[.7 .7 .7]);
m_usercoast('costa.mat','patch',[0 0 0])
%  m_plot(lon,lat,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4);
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',6);
title('\Psi - Camada 1 (m^2 s^{-1})','fontweigh','bold')
set(gcf,'color','w')


% psi na camada inferior ---------------------------------------------------
lpsi=[min(min(psi2)):0.5:max(max(psi2))];

figure(2)
set(gcf,...
        'Color',[1 1 1],...
        'InvertHardcopy','on',...
        'PaperUnits','inches',...
        'Units','inches',...
        'PaperOrientation','portrait',...
        'PaperPosition',[0 0 8.5 11],...
        'PaperPositionMode','manual',...
        'PaperType','usletter',...
        'Position',[.2 .2 8.5 11],...
        'ShareColors','off',...
        'Clipping','on');

subplot(321)

m_contourf(xg,yg,psi2,lpsi); shading flat; hold on; 
cc=colorbar;
m_quiver(xg,yg,u2psi,v2psi,2,'k'); hold on
%  fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]);
fill([xL200(1); xL200' ; xL200(end)],[yL200(end);yL200'; yL200(end)],[.7 .7 .7]);
m_usercoast('costa.mat','patch',[0 0 0])
%  m_plot(lon,lat,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4);
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',6);
title('\Psi - Camada 2 (m^2 s^{-1})','fontweigh','bold')
set(gcf,'color','w')


%%% calculando Vorticidade Potencial

% vorticidade relativa
vr1 = vortoa(xg,yg,lonm,latm,u1m,v1m,lc,E);
vr2 = vortoa(xg,yg,lonm,latm,u2m,v2m,lc,E);
vr1=reshape(vr1,li,lj);
vr2=reshape(vr2,li,lj);

% vorticidade de estiramento
ve1 = (f0.^2./(e.*g.*H1)) .* (psi2 - psi1);
ve2 = (f0.^2./(e.*g.*H2)) .* (psi1 - psi2);

% vorticidade planetaria
y = [0 ; cumsum(sw_dist(yg(:,1),xg(:,1),'km'))*1000];
y = y - y((length(y)-1)./2);
[ans,y] = meshgrid(y);
vpL = sw_b(5).*y; 

% vorticidade potencial
VP1 = vr1 + vpL + ve1;
VP2 = vr2 + vpL + ve2;

% vorticidade relativa na camada superior ---------------------------------------------------
lvr=[min(min(vr1)):100:max(max(vr1))];

figure(1)
subplot(322)

m_contourf(xg,yg,vr1,lvr); hold on; 
cc=colorbar;
%  fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]);
fill([xL200(1); xL200' ; xL200(end)],[yL200(end);yL200'; yL200(end)],[.7 .7 .7]);
m_usercoast('costa.mat','patch',[0 0 0])
%  m_plot(lon,lat,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4);
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',6);
title('V.Rel - Camada 1 (s^{-1})','fontweigh','bold')
set(gcf,'color','w')

% vorticidade relativa na camada inferior ---------------------------------------------------
lvr=[min(min(vr2)):0.5:max(max(vr2))];

figure(2)
subplot(322)

m_contourf(xg,yg,vr2,lvr); hold on; 
cc=colorbar;
%  fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]);
fill([xL200(1); xL200' ; xL200(end)],[yL200(end);yL200'; yL200(end)],[.7 .7 .7]);
m_usercoast('costa.mat','patch',[0 0 0])
%  m_plot(lon,lat,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4);
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',6);
title('V.Rel - Camada 2 (s^{-1})','fontweigh','bold')
set(gcf,'color','w')

% vorticidade de estiramento na camada superior ---------------------------------------------------
lve=[min(min(ve1)):7e-13:max(max(ve1))];

figure(1)
subplot(323)

m_contourf(xg,yg,ve1,lve); hold on; 
cc=colorbar;
%  fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]);
fill([xL200(1); xL200' ; xL200(end)],[yL200(end);yL200'; yL200(end)],[.7 .7 .7]);
m_usercoast('costa.mat','patch',[0 0 0])
%  m_plot(lon,lat,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4);
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',6);
title('V.Est - Camada 1 (s^{-1})','fontweigh','bold')
set(gcf,'color','w')

% vorticidade de estiramento na camada inferior ---------------------------------------------------
lve=[min(min(ve2)):2e-13:max(max(ve2))];

figure(2)
subplot(323)

m_contourf(xg,yg,ve2,lve); shading flat; hold on; 
m_contour(xg,yg,ve2,lve,'k');
cc=colorbar;
%  fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]);
fill([xL200(1); xL200' ; xL200(end)],[yL200(end);yL200'; yL200(end)],[.7 .7 .7]);
m_usercoast('costa.mat','patch',[0 0 0])
%  m_plot(lon,lat,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4);
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',6);
title('V.Est - Camada 2 (s^{-1})','fontweigh','bold')
set(gcf,'color','w')

% vorticidade planetaria ---------------------------------------------------
lvpL=[min(min(vpL)):1e-6:max(max(vpL))];

figure(1)
subplot(324)

m_contourf(xg,yg,vpL,lvpL); shading flat; hold on; 
cc=colorbar;
%  fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]);
fill([xL200(1); xL200' ; xL200(end)],[yL200(end);yL200'; yL200(end)],[.7 .7 .7]);
m_usercoast('costa.mat','patch',[0 0 0])
%  m_plot(lon,lat,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4);
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',6);
title('V. Planetaria (s^{-1})','fontweigh','bold')
set(gcf,'color','w')

% vorticidade planetaria ---------------------------------------------------
lvpL=[min(min(vpL)):1e-6:max(max(vpL))];

figure(2)
subplot(324)

m_contourf(xg,yg,vpL,lvpL); shading flat; hold on; 
cc=colorbar;
%  fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]);
fill([xL200(1); xL200' ; xL200(end)],[yL200(end);yL200'; yL200(end)],[.7 .7 .7]);
m_usercoast('costa.mat','patch',[0 0 0])
%  m_plot(lon,lat,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4);
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',6);
title('V. Planetaria (s^{-1})','fontweigh','bold')
set(gcf,'color','w')


% vorticidade potencial na camada superior ---------------------------------------------------
lvp=[min(min(VP1)):100:max(max(VP1))];

figure(1)
subplot(325)

m_contourf(xg,yg,VP1,lvp); shading flat; hold on; 
m_contour(xg,yg,VP1,lvp,'k');
cc=colorbar;
%  fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]);
fill([xL200(1); xL200' ; xL200(end)],[yL200(end);yL200'; yL200(end)],[.7 .7 .7]);
m_usercoast('costa.mat','patch',[0 0 0])
%  m_plot(lon,lat,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4);
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',6);
title('V. Pot - Camada 1 (s^{-1})','fontweight','bold')
set(gcf,'color','w')

% vorticidade potencial na camada inferior ---------------------------------------------------
lvp=[min(min(VP2)):.5:max(max(VP2))];

figure(2)
subplot(325)

m_contourf(xg,yg,VP2,lvp); shading flat; hold on; 
m_contour(xg,yg,VP2,lvp,'k'); 
cc=colorbar;
%  fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]);
fill([xL200(1); xL200' ; xL200(end)],[yL200(end);yL200'; yL200(end)],[.7 .7 .7]);
m_usercoast('costa.mat','patch',[0 0 0])
%  m_plot(lon,lat,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4);
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',6);
title('V. Pot - Camada 2 (s^{-1})','fontweigh','bold')
set(gcf,'color','w')


% vorticidade potencial + psi - camada superior ---------------------------------------------------
lvp=[min(min(VP1)):100:max(max(VP1))];

figure(1)
subplot(326)

%  m_contourf(xg,yg,VP1,lvp); shading flat; hold on; 
m_contour(xg,yg,VP1,lvp,'b');hold on
cc=colorbar;
m_quiver(xg,yg,u1psi,v1psi,2,'k'); hold on
%  fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]);
fill([xL200(1); xL200' ; xL200(end)],[yL200(end);yL200'; yL200(end)],[.7 .7 .7]);
m_usercoast('costa.mat','patch',[0 0 0])
%  m_plot(lon,lat,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4);
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',6);
title('VP x \Psi - Camada 1','fontweigh','bold')
set(gcf,'color','w')


% vorticidade potencial + psi - camada inferior ---------------------------------------------------
lvp=[min(min(VP2)):0.5:max(max(VP2))];

figure(2)
subplot(326)

%  m_contourf(xg,yg,VP2,lvp); shading flat; hold on; 
m_contour(xg,yg,VP2,lvp,'b');hold on
cc=colorbar;
m_quiver(xg,yg,u2psi,v2psi,2,'k'); hold on
%  fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]);
fill([xL200(1); xL200' ; xL200(end)],[yL200(end);yL200'; yL200(end)],[.7 .7 .7]);
m_usercoast('costa.mat','patch',[0 0 0])
%  m_plot(lon,lat,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4);
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',6);
title('VP x \Psi - Camada 1','fontweigh','bold')
set(gcf,'color','w')


print(1,'-depsc','camada1.eps')
!epstopdf camada1.eps
!rm -rf camada1.eps

print(2,'-depsc','camada2.eps')
!epstopdf camada2.eps
!rm -rf camada2.eps













