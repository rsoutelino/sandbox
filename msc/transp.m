%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   CALCULO DE TRANSPORTE EM SECOES VERTICAIS
%          DE VELOCIDADE GEOSTROFICA 
%               Oceano Leste 2
%         out/2006 - Mestrado - IOUSP
%           Rafael Guarino Soutelino
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
clear all;close all;clc;warning off

load ../../common/costa_leste.mat
lon = ncst(:,1); lat = ncst(:,2);
load ../../common/seagrid_leste2.dat
xg = seagrid_leste2(:,9);
yg = seagrid_leste2(:,10);
lj = max(seagrid_leste2(:,1));
li = max(seagrid_leste2(:,2));
xg = reshape(xg,li,lj);
yg = reshape(yg,li,lj);

figure(1)
plot(lon,lat,'k');hold on
axis([-42 -32 -21 -9])
for k = 1:3:length(xg)
    plot(xg(k,:),yg(k,:))
    text(xg(k,end),yg(k,end),num2str(k),'color','r') 
end

sec = input('Escolha a secao para o calculo do transporte: ')
close(1)

eval(['load ../mat/lesteII_vgao_sec',num2str(sec),'.mat'])

%%% CALCULO DE TRANSPORTE

%%% CALCULA OS VALORES DE VELOCIDADE NO CENTRO DA CELULA DE GRADE %%%

vgaom = 0.5*(vgao(:,2:end)+vgao(:,1:end-1));
vgaom = 0.5*(vgaom(2:end,:)+vgaom(1:end-1,:));

%%% MATRIZ DE dx's E dz's %%%

[xx,zz] = meshgrid(dist,prof);

dx = xx(:,2:end)-xx(:,1:end-1); dx = dx*1000; % passando para m
dz = zz(2:end,:)-zz(1:end-1,:); dz=-dz;

xm = 0.5*(xx(:,2:end)+xx(:,1:end-1));
xm = 0.5*(xm(2:end,:)+xm(1:end-1,:));
zm = 0.5*(zz(:,2:end)+zz(:,1:end-1));
zm = 0.5*(zm(2:end,:)+zm(1:end-1,:));

%%% TRANSPORTE INDIVIDUAL DAS CORRENTES %%%

dxbt = sw_dist(ybt,xbt,'km');
%  dxbt = sw_dist(zeros(size(xbt)),xbt,'km');
dxbt = [0 cumsum(dxbt)];

figure
set(gcf,'Color',[1 1 1]);
lv = [-0.5:0.01:0.5];
contourf(dist,-prof,-vgao,lv); shading flat; hold on; caxis([lv(1) lv(end)]);
set(gca,'plotboxaspectratio',[1 .7 1])
fill([min(dxbt) dxbt],[min(zbt) zbt],[.7 .7 .7])
cc = colorbar;
pos = get(cc,'Position');
set(cc,'Position',[pos(1) pos(2)*1.9 pos(3)/2.5 pos(4)/1.33])
pos=get(cc,'Position');
cc1=-str2num(get(cc,'YTickLabel'));cc1(find(cc1==0))=0;
set(cc,'YTickLabel',cc1)


%%% CELULAS DE VELOCIDADE CENTRAL SUPERIOR A 3 cm/s %%%
disp(' ')
disp('Selecione os limites laterais e verticais do fluxo para sul que deseja medir: ')
disp(' ')
sul = ginput(4);
fsul = find(vgaom<0 & xm>=sul(1,1) & xm<=sul(2,1) & -zm<=sul(3,2) & -zm>=sul(4,2));
plot(xm(fsul),-zm(fsul),'k.','markersize',3)

disp(' ')
disp('Selecione os limites laterais e verticais do fluxo para norte que deseja medir: ')
disp(' ')
norte = ginput(4);
fnorte = find(vgaom>0 & xm>=norte(1,1) & xm<=norte(2,1) & -zm<=norte(3,2) & -zm>=norte(4,2));
plot(xm(fnorte),-zm(fnorte),'w.','markersize',3)

pause 
close(1)

%%% TRANSPORTE DE VOLUME %%%

Tsul = sum(vgaom(fsul).*dx(fsul).*dz(fsul))*1e-6
Tnorte = sum(vgaom(fnorte).*dx(fnorte).*dz(fnorte))*1e-6
 
%%% VISUALIZACAO

latsec = -round(mean(ysec));
dxbt = sw_dist(ybt,xbt,'km');
%  dxbt = sw_dist(zeros(size(xbt)),xbt,'km');
dxbt = [0 cumsum(dxbt)];


figure
set(gcf,'Color',[1 1 1]);
lv = [-0.5:0.01:0.5];
contourf(dist,-prof,-vgao,lv); shading flat; hold on; caxis([lv(1) lv(end)]);
axis([0 dist(end) -max(prof) 0])
set(gca,'plotboxaspectratio',[1 .7 1])
plot(xm(fsul),-zm(fsul),'k.','markersize',3)
plot(xm(fnorte),-zm(fnorte),'w.','markersize',3)
fill([min(dxbt) dxbt],[min(zbt) zbt],[.7 .7 .7])
tit = [' ',num2str(latsec),'^\circS - Tv Sul = ',num2str(abs(round(Tsul))),' Sv,  Tv Norte = ',num2str(abs(round(Tnorte))),' Sv'];
t = title(tit,'fontweight','bold');
pos = get(t,'Position');
set(t,'Position',[pos(1) pos(2)*3 pos(3)])
xlabel('[^\circ W] / [km]')
ylabel('Profundidade [m]')
cc = colorbar;
pos = get(cc,'Position');
set(cc,'Position',[pos(1) pos(2)*1.9 pos(3)/2.5 pos(4)/1.33])
%  pos=get(cc,'Position');
cc1=str2num(get(cc,'YTickLabel')); %cc1(find(cc1==0))=0;
set(cc,'YTickLabel',-cc1)
ax1 = gca;
% eixo 2
ax2 = axes('Position',get(ax1,'Position'));
p=plot(xsec,xsec); hold on
set(gca,'color','none','XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k')
set(gca,'YTickLabel',[]);
set(p,'linestyle','none')
set(gca,'plotboxaspectratio',[1 .7 1])

   print(1,'-depsc',['../figuras/transp_OEII_',num2str(latsec),'S']);
   eval(['!epstopdf ../figuras/transp_OEII_',num2str(latsec),'S.eps'])

