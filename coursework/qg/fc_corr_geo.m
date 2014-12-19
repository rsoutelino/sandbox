%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%        PROGRAMA PARA CALCULO DE         %%%
%%%     FUNCAO DE CORRENTE GEOSTROFICA -    %%%
%%%               WESTRAX 2                 %%%
%%%     Rafael Soutelino - junho / 2006     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc

%%% LOOPING PARA LEITURA DOS DADOS

lon=[]; lat=[]; nest=[];

for i=[22:37 39:50]

    eval(['load ../ctd/w2cpz0',num2str(i),'.pro;']);
    eval(['data = w2cpz0',num2str(i),';'])
    if data(2,8)==9999 % tirando as estacoes que nao tem perfilagem pegasus
        x=0; % só pra constar no if
    else
        p=data(:,1);
        f=find(p==1500);
        if isempty(f)==0;
           lon=[lon data(1,6)*-1];
           lat=[lat data(1,5)];
           eval(['p',num2str(i),' = data(15:end,1);'])
           eval(['T',num2str(i),' = data(15:end,3);'])
           eval(['S',num2str(i),' = data(15:end,5);'])
           nest=[nest i];
        end
    end
    eval(['clear w2cpz0',num2str(i)])    
end

clear data

%%% calculando anomalia do geopotencial

for i=nest
    eval(['gpan',num2str(i),' = sw_gpan(S',num2str(i),',T',num2str(i),',p',num2str(i),');'])
end

f0=sw_f(5);

gpan100=[]; gpan1500=[];

for i=nest
    eval(['gpan100=[gpan100 gpan',num2str(i),'(100)];'])
    eval(['gpan1500=[gpan1500 gpan',num2str(i),'(1500)];'])
end

for i=nest
    eval(['clear T',num2str(i),])
    eval(['clear S',num2str(i),])
    eval(['clear p',num2str(i),])
    eval(['clear gpan',num2str(i),])
end

psig=(gpan100-gpan1500)/f0; psig=-psig;
 
% extraindo a anomalia de psig

psig=psig-mean(psig);

save geostream.mat lon lat psig

% aplicando condicoes de contorno

load coast.dat; xcoast=coast(:,2); ycoast=coast(:,1);
load m200.dat; 
x200=m200(:,2); 
y200=m200(:,1);	
x200=x200'; 
y200=y200';

zero=zeros(size(x200));

% livre escorregamento
lon2=[lon  x200];
lat2=[lat  y200];
psig2=[psig  zero];

save geostream_bc.mat lon2 lat2 psig2


% limites para plotagem

lonlim=[-53.5 -43]; latlim=[-0.5 10]; % limites de silveira_etal2000

m_proj('mercator','long',[lonlim(1) lonlim(2)],'lat',[latlim(1) latlim(2)],'on');


% montando reta L200, para usar como espelho

xL200=lonlim(1):0.25:lonlim(2);
yL200=xL200*(-1) - 45.5;

[xL200,yL200]=m_ll2xy(xL200,yL200,'clip','patch');


% interpolando pra uma grade regular

xi=-53:0.25:-44; yi=0:0.25:9;
[xg,yg]=meshgrid(xi,yi);
[lj,li]=size(xg);

lpsi=[min(psig):10000:max(psig)];

% interpolacao v4
psigI=griddata(lon2,lat2,psig2,xg,yg,'v4');


% ANALISE OBJETIVA -----------------------------------------

xg2=reshape(xg,1,lj*li);
yg2=reshape(yg,1,lj*li);

lc=4;
E=0.09^2;

[psigO,er]=scaloa(xg2,yg2,lon2,lat2,psig2,lc,E);

er=100*(sqrt(er));
psigO=reshape(psigO,li,lj);
er=reshape(er,li,lj);
xg=reshape(xg2,li,lj);
yg=reshape(yg2,li,lj);

% calculando vetores de velocidade

[v,u]=gradient(psigI); u=-u;
[vo,uo]=gradient(psigO); uo=-uo;


% preenchendo poligonos

xcoast=[min(xcoast); xcoast]; ycoast=[min(ycoast); ycoast]; 
x200=[min(x200) x200]; y200=[min(y200) y200];


% PLOTANDO ==========================================================================================


load iso200.mat; 

% v4 ---------------------------------------------------

figure(1)
m_contourf(xg,yg,psigI,lpsi); shading flat; hold on; 
cc=colorbar;
m_quiver(xg,yg,u,v,2,'w'); hold on
plot(xL200,yL200,'w','linewidth',4);
%  fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]);
m_usercoast('costa.mat','patch',[0 0 0])
m_plot(lon,lat,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',6);
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',12);
title('Funcao de Corrente Geostrofica \Psi: [m^2 s^{-1}] - V4 - 100 m','fontsize',12,'fontweigh','bold')
set(gcf,'color','w')

print -depsc fc_corr_geo.eps
!epstopdf fc_corr_geo.eps


% AO ---------------------------------------------------

figure(2)
m_contourf(xg,yg,psigO,lpsi); shading flat; hold on; 
cc=colorbar;
m_quiver(xg,yg,uo,vo,2,'w'); hold on
plot(xL200,yL200,'w','linewidth',3);
plot(x1,y1,'k','linewidth',2)
%  fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]);
m_usercoast('costa.mat','patch',[0 0 0])
m_plot(lon,lat,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',6);
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',12);
title('\Psi: [m^2 s^{-1}] - AO - 100 m - cc = zeros','fontsize',12,'fontweigh','bold')
set(gcf,'color','w')

print -depsc fc_corr_geo_AO.eps
!epstopdf fc_corr_geo_AO.eps

lE=min(min(er)):4:max(max(er));
lE=round(lE);

figure(3)
[c,h]=m_contour(xg,yg,er,lE,'k'); hold on
clabel(c,h);
%  m_contourf(xg,yg,er,lE); shading flat; hold on; 
%  caxis([0 20])
%  cc=colorbar;
fill([x1';x1(end);x1(1)],[y1';y1(1);y1(1)],[.7 .7 .7]); hold on
m_usercoast('costa.mat','patch',[0 0 0])
m_plot(lon,lat,'ko','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',6);
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',12);
title('Erro de Interpolacao (%) - AO - cc = zeros','fontsize',12,'fontweigh','bold')
set(gcf,'color','w')

print -depsc erro_AO.eps
!epstopdf erro_AO.eps

%%% APLICANDO CONDICAO DE CONTORNO ATRAVES DO METODO DAS IMAGENS ================================================================
clc
% matrizes originais = lon, lat, psig

% calculando as imagens
[lon2,lat2,psig2]=mirror_sc(lon,lat,psig);

% incorporando as imagens
lon3=[lon lon2];
lat3=[lat lat2];
psig3=[psig psig2];

% analise objetiva

[psigO2,er2] = scaloa(xg2,yg2,lon3,lat3,psig3,lc,E);

% redimensionando para plotar
er2=100*(sqrt(er2));
psigO2=reshape(psigO2,li,lj);
er2=reshape(er2,li,lj);

% calculando vetores de velocidade

[vo2,uo2]=gradient(psigO2); uo2=-uo2;


figure(4)
m_contourf(xg,yg,psigO2,lpsi); shading flat; hold on; 
cc=colorbar;
m_quiver(xg,yg,uo2,vo2,2,'w');hold on
plot(xL200,yL200,'w','linewidth',3);
plot(x1,y1,'k','linewidth',2)
%  fill([xL200(1) xL200],[yL200(end) yL200],[.7 .7 .7]);
m_usercoast('costa.mat','patch',[0 0 0])
m_plot(lon,lat,'ko','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',6);
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',12);
title('\Psi: [m^2 s^{-1}] - AO - 100 m - cc = imagens','fontsize',12,'fontweigh','bold')
set(gcf,'color','w')

print -depsc fc_corr_geo_AO_cc2.eps
!epstopdf fc_corr_geo_AO_cc2.eps
clc

figure(5)
[c,h]=m_contour(xg,yg,er2,lE,'k'); hold on
clabel(c,h);
%  m_contourf(xg,yg,er,lE); shading flat; hold on; 
%  caxis([0 20])
%  cc=colorbar;
fill([xL200(1) xL200],[yL200(end) yL200],[.7 .7 .7]); hold on
m_usercoast('costa.mat','patch',[0 0 0])
m_plot(lon,lat,'ko','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',6);
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',12);
title('Erro de Interpolacao (%) - AO - cc = imagens ','fontsize',12,'fontweigh','bold')
set(gcf,'color','w')

print -depsc erro_AO_cc2.eps
!epstopdf erro_AO_cc2.eps

clc