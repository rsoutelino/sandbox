%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROGRAMA PARA CONSTRUIR UM MODELO DE STOMMEL
%   PARA OS GIROS OCEANICOS , USANDO MEDIAS 
%     MENSAIS DE VENTO DO QUICKSCAT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc; tic

load qk00_05_v2.mat;

Ujun=U(:,:,6);
Vjun=V(:,:,6);

Ujul=U(:,:,7);
Vjul=V(:,:,7);

Uago=U(:,:,8);
Vago=V(:,:,8);

u = Ujun + Ujul + Uago; u = u/3;
v = Vjun + Vjul + Vago; v = v/3;

udec = u(1:5:end,1:5:end);  vdec = v(1:5:end,1:5:end);
v0 = zeros(size(vdec)); 

[lon,lat]=meshgrid(lon,lat);
londec = lon(1:5:end,1:5:end);  latdec = lat(1:5:end,1:5:end);

mod = sqrt(udec.^2+vdec.^2); lmod = min(min(mod)):0.3:max(max(mod));

lonlim=[-100 50];
latlim=[-65 65];

m_proj('mercator','long',[lonlim(1) lonlim(2)],'lat',[latlim(1) latlim(2)],'on');

% plotando ventos para toda a bacia do atlantico
%  
%  figure
%  m_contourf(lon,lat,sqrt(u.^2+v.^2),lmod);shading flat; hold on
%  m_quiver(londec,latdec,udec,vdec,2,'k');
%  m_usercoast('costa_atl.mat','patch',[0.5 0.5 0.5]);
%  m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','right','xaxislocation','top','fontsize',10);
%  m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',10);
%  set(gcf,'color','w')
%  %  title('Vento Medio para o verao do H.S nas Bacias do Oceano Atlantico')
%  
%  print -depsc vento_atl.eps
%  !epstopdf vento_atl.eps
%  print -dpng vento_atl.png
%  
%  clear Ujun Ujul Uago Vjun Vjul Vago U V
%  
%  % plotando os limites geograficos para o modelo
%  % definindo limites:
%  
%  lonbox2=[-33 -18 -18 -33];
%  latbox2=[-60 -60 60 60];
%  [lonbox,latbox] = m_ll2xy(lonbox2,latbox2);
%  
%  figure
%  fill(lonbox,latbox,[0.7 0.7 0.7])
%  m_usercoast('costa_atl.mat','patch',[0 0 0]);hold on
%  m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',10);
%  title('Limites Geograficos para o Modelo','fontweight','bold')
%  set(gcf,'color','w')
%  m_text(-46,-62,'33 W')
%  m_text(-20,-62,'18 W')
%  print -depsc box.eps
%  !epstopdf box.eps
%  print -dpng box.png

% truncando os ventos para os limites do modelo

flon=find(lon>=-33 & lon<=-18);
lon=lon(flon);lat=lat(flon);u=u(flon);v=v(flon);
flat=find(lat>=-60 & lat<=60);
lon=lon(flat);lat=lat(flat);u=u(flat);v=v(flat);

clear flat flon
% convertendo vento para tensao de cizalhamento
mod = sqrt(u.^2+v.^2);
rho_ar = 1.0235; % kg/m^3 
Cd = (1.1 + 0.0536.*mod).*0.001; % adimensional
Tx = rho_ar.*Cd.*mod.*u; % Pa = kg m^-1 s^-2
Ty = rho_ar.*Cd.*mod.*v;

% interpolando para area do modelo

xi=[-33:0.25:-18];
yi=[-60:0.25:60];
[xg,yg]=meshgrid(xi,yi);
%  xg=flipud(xg); yg=flipud(yg);
Tx = griddata(lon,lat,Tx,xg,yg);
Ty = griddata(lon,lat,Ty,xg,yg);
[li,lj]=size(Tx);

clear xi yi
% calculando a media zonal da TCV

Txm=nanmean(Tx'); Txm=Txm';
Tym=nanmean(Ty'); Tym=Tym';

for k=1:lj;
   Txm2(:,k)=Txm;
   Tym2(:,k)=Tym;
end

Txm=Txm2; T0=zeros(size(Txm));
Tym=Tym2; 
modm=sqrt(Txm.^2+Tym.^2);
resm=abs(Txm)+abs(Tym);

clear Txm2 Tym2

% plotando a media zonal para TCV e comparando com media meridional

load ygiro.mat

figure

set(gcf,'Color',[1 1 1]);
subplot(121)
x0 = [min(min(Txm))-0.05 :mean(mean(diff(Txm))): max(max(Txm))+0.05];
for i=1:length(ygiro) 
   plot(x0,ygiro(i)*ones(size(x0)),'g');hold on
end
plot(T0(:,1),yg(:,1),'k'); hold on
plot(Txm(:,1),yg(:,1),'r','linewidth',2);
plot(yg(:,1),T0(:,1),'--k')
axis([min(min(Txm))-0.05 max(max(Txm))+0.05 -61 61])
grid
title('Media Zonal da TCV (Pa)','fontweight','bold')
xlabel('TCV')
ylabel('Latitude')
set(gca,'plotboxaspectratio',[.5 .8 .5])

subplot(122)

plot(abs(Txm(:,1)),yg(:,1),'r','linewidth',2); hold on
plot(abs(Tym(:,1)),yg(:,1),'b');
%  plot(resm(:,1),yg(:,1),'k');
plot(T0(:,1),yg(:,1),'k'); 
plot(yg(:,1),T0(:,1),'--k')
legend('Zonal','Meridional')
title('Magnitudes da TCV','fontweight','bold')
xlabel('TCV')
ylabel('Latitude')
axis([min(min(Txm))-0.05 max(max(Txm))+0.05 -61 61])
grid
for i=1:length(ygiro) 
   plot(x0,ygiro(i)*ones(size(x0)),'g');hold on
end
set(gca,'plotboxaspectratio',[.5 .8 .5])

%  print -depsc perfis_medios.eps
%  print -dpng perfis_medios.png
%  !epstopdf perfis_medios.eps


% calculando os transportes de EKMAN, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% negligenciando as regioes entre 5 S e 5 N

rho_agua = 1026.954; % kg/m^3
f = sw_f(yg);
fnan = find(yg >=-1 & yg <=1);
% colocando NaN no valor de f0 para as latitudes entre 5 S e 5 N
f(fnan)=NaN;

Ue = (1/rho_agua).*(Tym./f);
Ve = (-1/rho_agua).*(Txm./f);

figure
set(gcf,'color','w')
subplot(121)
x0 = [min(min(Ue))-10 :nanmean(nanmean(abs(diff(Ue)))): max(max(Ue))+10];
for i=1:length(ygiro) 
   plot(x0,ygiro(i)*ones(size(x0)),'g');hold on
end
plot(Ue(:,1),yg(:,1),'r','linewidth',2); hold on
plot(T0(:,1),yg(:,1),'k'); 
plot(yg(:,1),T0(:,1),'--k')
title('Transporte de Ekman Zonal','fontweight','bold')
xlabel('m^2 s^{-1}')
ylabel('Latitude')
axis([-2 2 -61 61])
grid
set(gca,'plotboxaspectratio',[.5 .8 .5])

subplot(122)
x0 = [min(min(Ve))-10 :nanmean(nanmean(abs(diff(Ve)))): max(max(Ve))+10];
for i=1:length(ygiro) 
   plot(x0,ygiro(i)*ones(size(x0)),'g');hold on
end
plot(Ve(:,1),yg(:,1),'b','linewidth',2); hold on
plot(T0(:,1),yg(:,1),'k'); 
plot(yg(:,1),T0(:,1),'--k')
title('Transporte de Ekman Meridional','fontweight','bold')
xlabel('m^2 s^{-1}')
ylabel('Latitude')
axis([-2 2 -61 61])
grid
set(gca,'plotboxaspectratio',[.5 .8 .5])

print -depsc perfis_ekman.eps
print -dpng perfis_ekman.png
!epstopdf perfis_ekman.eps


figure
set(gcf,'color','w')
subplot(121)
x0 = [-40:0.25:-10];
for i=1:length(ygiro) 
   plot(x0,ygiro(i)*ones(size(x0)),'g');hold on
end
quiver(xg(1:15:end,1:10:end),yg(1:15:end,1:10:end),T0(1:15:end,1:10:end),Ve(1:15:end,1:10:end),3,'k')
hold on
plot(T0(:,1),yg(:,1),'k'); 
plot(yg(:,1),T0(:,1),'--k')
axis([-35 -17  -61 60])
grid
title('Transporte de Ekman Meridional para a Area do Modelo','fontweight','bold')
xlabel('Longitude')
ylabel('Latitude')

subplot(122)
x0 = [-40:0.25:-10];
for i=1:length(ygiro) 
   plot(x0,ygiro(i)*ones(size(x0)),'g');hold on
end
quiver(xg(1:10:end,37),yg(1:10:end,37),Txm(1:10:end,37),T0(1:10:end,37),.5)
hold on
plot(T0(:,1),yg(:,1),'k'); 
plot(yg(:,1),T0(:,1),'--k')
axis([-36 -17  -61 60])
grid
title('Vento Zonal Medio Para a Area do Modelo','fontweight','bold')
xlabel('Longitude')

print -depsc distr_hor_ekman.eps
print -dpng distr_hor_ekman.png
!epstopdf distr_hor_ekman.eps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculando a componente vertical do rotacional do vento, considerando 
% o vento total

[ans,dTx_dy] = gradient(Tx,27780,27780);
[dTy_dx,ans] = gradient(Ty,27780,27780);

RotT = dTy_dx - dTx_dy;
%  RotT = curl(xg,yg,Tx,Ty);
lrot=[min(min(RotT)):0.00000003:max(max(RotT))];


figure

set(gcf,'color','w')

subplot(133)
contourf(xg,yg,RotT,lrot);shading flat;hold on;%caxis([-5e-3 4e-3]);
colorbar
[c,h]=contour(xg,yg,RotT,[0 0],'k');
%  quiver(xg(1:10:end),yg(1:10:end),Tx(1:10:end),Ty(1:10:end),'k')
plot(yg(:,1),T0(:,1),'--k','linewidth',2)
axis([-33 -18  -61 60])
title('Rotacional da TCV','fontweight','bold')
xlabel('Pa m^{-1}')
%  set(gca,'plotboxaspectratio',[.5 .8 .5])

subplot(131)
quiver(xg(1:10:end,37),yg(1:10:end,37),Txm(1:10:end,37),T0(1:10:end,37),.5)
hold on
plot(T0(:,1),yg(:,1),'k'); 
plot(yg(:,1),T0(:,1),'--k')
axis([-36 -17  -61 60])
grid
title('Vento Zonal Medio','fontweight','bold')
ylabel('Latitude')
%  set(gca,'plotboxaspectratio',[.5 .8 .5])

subplot(132)
plot(nanmean(RotT'),yg(:,1),'b','linewidth',2); hold on
plot(T0(:,1),yg(:,1),'k'); 
plot(yg(:,1),T0(:,1),'--k')
grid
axis([-2.2e-7 2.2e-7 -61 61])
title('Perfil Medio - Rotacional TCV','fontweight','bold')
xlabel('Pa m^{-1}')
%  set(gca,'plotboxaspectratio',[.5 .8 .5])

print -depsc rot_vento_total.eps
print -dpng rot_vento_total.png
!epstopdf rot_vento_total.eps

% calculando a componente vertical do rotacional do vento, considerando 
% apenas o vento zonal
[ans,dTx_dy] = gradient(Txm,27780,27780);
RotTx = -dTx_dy;
lrotx=[min(min(RotTx)):0.000000003:max(max(RotTx))];

figure

set(gcf,'color','w')

subplot(133)
contourf(xg,yg,RotTx,lrot);shading flat;hold on;colorbar
[c,h]=contour(xg,yg,RotTx,[0 0],'k');hold on
%  quiver(xg(1:10:end),yg(1:10:end),Tx(1:10:end),Ty(1:10:end),'k')
plot(yg(:,1),T0(:,1),'--k','linewidth',2)
axis([-33 -18  -61 60])
%  title('Rotacional da componente zonal da TCV','fontweight','bold')
%  set(gca,'plotboxaspectratio',[.5 .8 .5])

subplot(131)
x0 = [-40:.25:-10];
for i=1:length(ygiro) 
   plot(x0,ygiro(i)*ones(size(x0)),'g');hold on
end
quiver(xg(1:10:end,37),yg(1:10:end,37),Txm(1:10:end,37),T0(1:10:end,37),.5)
hold on
plot(T0(:,1),yg(:,1),'k'); 
plot(yg(:,1),T0(:,1),'--k')
axis([-36 -17  -61 60])
grid
%  title('Vento Zonal Medio','fontweight','bold')
ylabel('Latitude')
%  set(gca,'plotboxaspectratio',[.5 .8 .5])

subplot(132)
x0 = [min(min(RotTx)):nanmean(nanmean(abs(diff(RotTx)))): max(max(RotTx))];
for i=1:length(ygiro) 
   plot(x0,ygiro(i)*ones(size(x0)),'g');hold on
end
plot(nanmean(RotTx'),yg(:,1),'b','linewidth',2); hold on
plot(T0(:,1),yg(:,1),'k'); 
plot(yg(:,1),T0(:,1),'--k')
grid
axis([min(min(RotTx)) max(max(RotTx)) -61 61])
title('Rotacional da componente zonal media da TCV','fontweight','bold')
xlabel('Longitude')
%  set(gca,'plotboxaspectratio',[.5 .8 .5])

print -depsc rot_vento_zonal.eps
print -dpng rot_vento_zonal.png
!epstopdf rot_vento_zonal.eps

% calculando o bombeamento de Ekman

We=(rho_agua.*flipud(f)).^(-1).*RotTx.*(-1);

figure

set(gcf,'color','w')
x0 = [min(min(We)):nanmean(nanmean(abs(diff(We)))): max(max(We))];
for i=1:length(ygiro) 
   plot(x0,ygiro(i)*ones(size(x0)),'g');hold on
end
plot(nanmean(We'),yg(:,1),'b','linewidth',2); hold on
plot(T0(:,1),yg(:,1),'k'); 
plot(yg(:,1),T0(:,1),'--k')
grid
axis([-2.5e-6 2.5e-6 -61 61])
title('Perfil Medio do Bombeamento de Ekman','fontweight','bold')
xlabel('m s^{-1}')
set(gca,'plotboxaspectratio',[.5 .8 .5])

print -depsc bomb.eps
print -dpng bomb.png
!epstopdf bomb.eps


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% APLICANDO A SOLUCAO DE STOMMEL PARA CRIAR UM CAMPO DE 
%%%      FUNCAO DE CORRENTE NO BOX MODEL CRIADO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% alisando o vento para melhor a derivada segunda

Txm = smoo2(Txm,-9999,15,.1);
Txm = smoo2(Txm,-9999,15,.1);

% CALCULANDO PARAMETROS CONSTANTES PARA TODOS OS GIROS

L1 = sw_dist([0 0],[-33 -18],'km')*1000;
rho = rho_agua;
H = 3000;
d = 10;

%% separando as regiões dos giros, usando isolinhas de rotacional = 0
% calculando latitude no centro de cada giro para calcular o 
% f0 para cada um

Ygiro=[];
for i=1:5
    eval(['ygiro' num2str(i) '=yg(near(yg,ygiro(' num2str(i) '),1));'])

if i >=2
   eval(['Cgiro(' num2str(i-1) ')=(ygiro' num2str(i) ' + ygiro' num2str(i-1) ')/2;']);
end

end

f0 = sw_f(Cgiro);
beta = sw_b(Cgiro);
a=6378000;% raio da terra

% truncando as matrizes de lat e lon para cada giro

%%% GIRO SUB-TROPICAL DO ATLANTICO NORTE

g1N = near(yg(:,1),ygiro(1),1);
g1S = near(yg(:,1),ygiro(2),1);

xg_AN = xg(g1S:g1N,:);
yg_AN = yg(g1S:g1N,:);
Txm_AN = Txm(g1S:g1N,:);

%  %************************************
%  y_AN = 0:500000:3500000;
%  x_AN = 0:500000:5000000;
%  [x_AN,y_AN] = meshgrid(x_AN,y_AN);
%  y_AN=flipud(y_AN);
%  Txm_AN=-.1*cos(pi.*y_AN/3500000);
%  L1=5000000;
%  %************************************
[li,lj] = size(xg_AN);

for i=1:li
   x_AN(i,:) = [0 sw_dist(zeros(1,lj),xg_AN(i,:),'km')]; 
   x_AN(i,:) = [cumsum(x_AN(i,:))];
end

x_AN = x_AN*1000;

for i=1:lj
   y_AN(:,i) = [0 ; sw_dist(yg_AN(:,i),zeros(1,li),'km')]; 
   y_AN(:,i) = [cumsum(y_AN(:,i))];
end

y_AN = y_AN*1000;

[ans,dTx_dy_AN] = gradient(Txm_AN,27780,27780); 
[ans,d2Tx_dy2_AN] = gradient(dTx_dy_AN,27780,27780);

parmtU_AN = -L1/(beta(1)*H*rho);
parmtV_AN = 1/(beta(1)*H*rho);
parmtP_AN = (f0(1)*L1) / (beta(1)*H);
parmtV2_AN = (2*beta(1)*H*L1) / (f0(1)*d);
Ls = exp( - ( ((2*beta(1)*H).*x_AN) ./ (f0(1)*d) ));

U_AN = parmtU_AN.*d2Tx_dy2_AN.* ( 1 - x_AN./L1 - Ls ); 
V_AN = parmtV_AN.*dTx_dy_AN.* (-1 + parmtV2_AN.*Ls);
P_AN = parmtP_AN.*dTx_dy_AN.* ( 1 - x_AN./L1 - Ls );

psi_AN = P_AN./(rho*f0(1));

%  figure;
%  contourf(x_AN,y_AN,psi_AN)

%%% GIRO TROPICAL

g2N = near(yg(:,1),ygiro(2),1);
g2S = near(yg(:,1),ygiro(3),1);

xg_TR = xg(g2S+1:g2N,:);
yg_TR = yg(g2S+1:g2N,:);
Txm_TR = Txm(g2S+1:g2N,:);

[li,lj] = size(xg_TR);

for i=1:li
   x_TR(i,:) = [0 sw_dist(zeros(1,lj),xg_TR(i,:),'km')]; 
   x_TR(i,:) = [cumsum(x_TR(i,:))];
end

x_TR = x_TR*1000;

for i=1:lj
   y_TR(:,i) = [0 ; sw_dist(yg_TR(:,i),zeros(1,li),'km')]; 
   y_TR(:,i) = [cumsum(y_TR(:,i))];
end

y_TR = y_TR*1000;


[ans,dTx_dy_TR] = gradient(Txm_TR,27780,27780); 
[ans,d2Tx_dy2_TR] = gradient(dTx_dy_TR,27780,27780);

parmtU_TR = -L1/(beta(2)*H*rho);
parmtV_TR = 1/(beta(2)*H*rho);
parmtP_TR = (f0(2)*L1) / (beta(2)*H);
parmtV2_TR = (2*beta(2)*H*L1) / (f0(2)*d);
Ls = exp( - ( ((2*beta(2)*H).*x_TR) ./ (f0(2)*d) ));

U_TR = parmtU_TR.*d2Tx_dy2_TR.* ( 1 - x_TR./L1 - Ls ); 
V_TR = parmtV_TR.*dTx_dy_TR.* (-1 + parmtV2_TR.*Ls);
P_TR = parmtP_TR.*dTx_dy_TR.* ( 1 - x_TR./L1 - Ls );

psi_TR = P_TR./(rho*f0(2));
%  
%  figure;
%  contourf(x_TR,y_TR,psi_TR)

%%% GIRO EQUATORIAL

g3N = near(yg(:,1),ygiro(3),1);
g3S = near(yg(:,1),ygiro(4),1);

xg_EQ = xg(g3S:g3N+1,:);
yg_EQ = yg(g3S:g3N+1,:);
Txm_EQ = Txm(g3S:g3N+1,:);

[li,lj] = size(xg_EQ);

for i=1:li
   x_EQ(i,:) = [0 sw_dist(zeros(1,lj),xg_EQ(i,:),'km')]; 
   x_EQ(i,:) = [cumsum(x_EQ(i,:))];
end

x_EQ = x_EQ*1000;

for i=1:lj
   y_EQ(:,i) = [0 ; sw_dist(yg_EQ(:,i),zeros(1,li),'km')]; 
   y_EQ(:,i) = [cumsum(y_EQ(:,i))];
end

y_EQ = y_EQ*1000;

[ans,dTx_dy_EQ] = gradient(Txm_EQ,27780,27780); 
[ans,d2Tx_dy2_EQ] = gradient(dTx_dy_EQ,27780,27780);

parmtU_EQ = -L1/(beta(3)*H*rho);
parmtV_EQ = 1/(beta(3)*H*rho);
parmtP_EQ = (f0(3)*L1) / (beta(3)*H);
parmtV2_EQ = (2*beta(3)*H*L1) / (f0(3)*d);
Ls = exp( - ( ((2*beta(3)*H).*x_EQ) ./ (abs(f0(3))*d) ));

U_EQ = parmtU_EQ.*d2Tx_dy2_EQ.* ( 1 - x_EQ./L1 - Ls ); 
V_EQ = parmtV_EQ.*dTx_dy_EQ.* (-1 + parmtV2_EQ.*Ls);
P_EQ = parmtP_EQ.*dTx_dy_EQ.* ( 1 - x_EQ./L1 - Ls );

psi_EQ = P_EQ./(rho*f0(3));

%  figure;
%  contourf(x_EQ,y_EQ,psi_EQ)

%%% GIRO SUBTROPICAL DO ATLANTICO SUL

g4N = near(yg(:,1),ygiro(4),1);
g4S = near(yg(:,1),ygiro(5),1);

xg_AS = xg(g4S:g4N,:);
yg_AS = yg(g4S:g4N,:);
Txm_AS = Txm(g4S:g4N,:);

[li,lj] = size(xg_AS);

for i=1:li
   x_AS(i,:) = [0 sw_dist(zeros(1,lj),xg_AS(i,:),'km')]; 
   x_AS(i,:) = [cumsum(x_AS(i,:))];
end

x_AS = x_AS*1000;

for i=1:lj
   y_AS(:,i) = [0 ; sw_dist(yg_AS(:,i),zeros(1,li),'km')]; 
   y_AS(:,i) = [cumsum(y_AS(:,i))];
end

y_AS = y_AS*1000;

[ans,dTx_dy_AS] = gradient(Txm_AS,27780,27780); 
[ans,d2Tx_dy2_AS] = gradient(dTx_dy_AS,27780,27780);

parmtU_AS = -L1/(beta(4)*H*rho);
parmtV_AS = 1/(beta(4)*H*rho);
parmtP_AS = (f0(4)*L1) / (beta(4)*H);
parmtV2_AS = (2*beta(4)*H*L1) / (f0(4)*d);
Ls = exp( - ( ((2*beta(4)*H).*x_AS) ./ (abs(f0(4))*d) ));

U_AS = parmtU_AS.*d2Tx_dy2_AS.* ( 1 - x_AS./L1 - Ls ); 
V_AS = parmtV_AS.*dTx_dy_AS.* (-1 + parmtV2_AS.*Ls);
P_AS = parmtP_AS.*dTx_dy_AS.* ( 1 - x_AS./L1 - Ls );

psi_AS = P_AS./(rho*f0(4));

%  figure;
%  contourf(x_AS,y_AS,psi_AS)

% juntando os giros

lpsian = min(min(psi_AN)):200:max(max(psi_AN));
lpsitr = min(min(psi_TR)):200:max(max(psi_TR));
lpsieq = min(min(psi_EQ)):200:max(max(psi_EQ));
lpsias = min(min(psi_AS)):200:max(max(psi_AS));

lonlim=[-33.1 -18];
latlim=[yg(g4S) yg(g1N)];

m_proj('mercator','long',[lonlim(1) lonlim(2)],'lat',[latlim(1) latlim(2)],'on');

figure
 m_contourf(xg_AN,yg_AN,psi_AN,lpsian);hold on
 m_contourf(xg_TR,yg_TR,psi_TR,lpsitr);
 m_contourf(xg_EQ,yg_EQ,psi_EQ,lpsieq);
 m_contourf(xg_AS,yg_AS,psi_AS,lpsias);shading flat
colorbar
%  caxis([min(min(psi_TR))  max(max(psi_AN))])
m_plot(xg_AN(1,:),yg_AN(1,:),'w','linewidth',3)
m_plot(xg_AN(end,:),yg_AN(end,:),'k','linewidth',3)
m_plot(xg_TR(end,:),yg_TR(end,:),'w','linewidth',3)
m_plot(xg_EQ(end,:),yg_EQ(end,:),'w','linewidth',3)
m_plot(xg_AS(end,:),yg_AS(end,:),'w','linewidth',3)
m_plot(xg_AS(1,:),yg_AS(1,:),'k','linewidth',3)
%  cc=colorbar;
m_grid('xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',10);
title('Campo de \Psi - Modelo de Stommel (m^2 s^{-1})','fontweight','bold')
set(gcf,'color','w')
daspect([.2 1 .2])

print -depsc stommel.eps
print -djpeg100 stommel.jpg
!epstopdf stommel.eps

%%% plotando os giros semparados

lonlim=[-33.4 -18];
latlim=[yg(g1S) yg(g1N)];

m_proj('mercator','long',[lonlim(1) lonlim(2)],'lat',[latlim(1) latlim(2)],'on');

lpsian = min(min(psi_AN)):200:max(max(psi_AN));

figure
set(gcf,'color','w')
m_contourf(xg_AN,yg_AN,psi_AN,lpsian);colorbar
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',10);
daspect([0.35 1 0.35])
title('Giro Sub-tropical do Atlantico Norte - Modelo de Stommel -\Psi (m^2 s^{-1})','fontweight','bold')

print -depsc giro_AN.eps
!epstopdf giro_AN.eps
print -djpeg100 giro_AN.jpg

%*******************************************************************************

lonlim=[-33.2 -18];
latlim=[yg(g2S) yg(g2N)];

m_proj('mercator','long',[lonlim(1) lonlim(2)],'lat',[latlim(1) latlim(2)],'on');

lpsitr = min(min(psi_TR)):200:max(max(psi_TR));

figure
set(gcf,'color','w')
m_contourf(xg_TR,yg_TR,psi_TR,lpsitr);colorbar('horiz')
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',10);
daspect([0.5 1 0.5])
title('Giro Tropical - Modelo de Stommel - \Psi (m^2 s^{-1})','fontweight','bold')

print -depsc giro_TR.eps
print -djpeg100 giro_TR.jpg
!epstopdf giro_TR.eps

% ******************************************************************************

lonlim=[-33.2 -18];
latlim=[yg(g3S) yg(g3N)];

m_proj('mercator','long',[lonlim(1) lonlim(2)],'lat',[latlim(1) latlim(2)],'on');

lpsieq = min(min(psi_EQ)):200:max(max(psi_EQ));

figure
set(gcf,'color','w')
m_contourf(xg_EQ,yg_EQ,psi_EQ,lpsieq);colorbar('horiz')
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',10);
daspect([0.4 1 0.4])
title('Giro Equatorial - Modelo de Stommel - \Psi (m^2 s^{-1})','fontweight','bold')

print -depsc giro_EQ.eps
!epstopdf giro_EQ.eps
print -djpeg100 giro_EQ.jpg

%********************************************************************************
 
lonlim=[-33.5 -18];
latlim=[yg(g4S) yg(g4N)];

m_proj('mercator','long',[lonlim(1) lonlim(2)],'lat',[latlim(1) latlim(2)],'on');

lpsias = min(min(psi_AS)):200:max(max(psi_AS));

figure
set(gcf,'color','w')
m_contourf(xg_AS,yg_AS,psi_AS,lpsias); colorbar
m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',10);
daspect([0.3 1 0.3])
title('Giro Sub-tropical do Atlantico Sul - Modelo de Stommel - \Psi (m^2 s^{-1})','fontweight','bold')

print -depsc giro_AS.eps
print -djpeg100 giro_AS.jpg
!epstopdf giro_AS.eps



%%% CALCULANDO OS TRANSPORTES DE CG E CB

Tcg = max(max(abs(V_AN))) * H * 100000;
Tcg = Tcg / 1e6; 

Tcb = max(max(abs(V_AS))) * H * 100000;
Tcb = Tcb / 1e6; 

disp(['Transporte da Corrente do Golfo: ' num2str(round(Tcg)) ' Sv'])
disp(' ')
disp(['Transporte da Corrente do Brasil: ' num2str(round(Tcb)) ' Sv'])

!rm -rf *.png
!rm -rf *.jpg

%%% figura extra

RotTx = smoo2(RotTx,-9999,15,.1);

figure
set(gcf,'color','w')

subplot(131)
x0 = [min(min(RotTx)):nanmean(nanmean(abs(diff(RotTx)))): max(max(RotTx))];
for i=1:length(ygiro) 
   plot(x0,ygiro(i)*ones(size(x0)),'g');hold on
end
plot(nanmean(RotTx'),yg(:,1),'b','linewidth',2); hold on
plot(T0(:,1),yg(:,1),'k'); 
plot(yg(:,1),T0(:,1),'--k')
grid
axis([-9e-8 9e-8 -51 51])
title('Rotacional da componente zonal media da TCV','fontweight','bold')
xlabel('Longitude')
%  set(gca,'plotboxaspectratio',[.25 1 .25])

subplot(1,1.5,1.5)


lonlim=[-33.1 -18];
latlim=[yg(g4S) yg(g1N)];

m_proj('mercator','long',[lonlim(1) lonlim(2)],'lat',[latlim(1) latlim(2)],'on');

 m_contourf(xg_AN,yg_AN,psi_AN,lpsian);hold on
 m_contourf(xg_TR,yg_TR,psi_TR,lpsitr);
 m_contourf(xg_EQ,yg_EQ,psi_EQ,lpsieq);
 m_contourf(xg_AS,yg_AS,psi_AS,lpsias);shading flat
%  colorbar
%  caxis([min(min(psi_TR))  max(max(psi_AN))])
m_plot(xg_AN(1,:),yg_AN(1,:),'w','linewidth',3)
m_plot(xg_AN(end,:),yg_AN(end,:),'k','linewidth',3)
m_plot(xg_TR(end,:),yg_TR(end,:),'w','linewidth',3)
m_plot(xg_EQ(end,:),yg_EQ(end,:),'w','linewidth',3)
m_plot(xg_AS(end,:),yg_AS(end,:),'w','linewidth',3)
m_plot(xg_AS(1,:),yg_AS(1,:),'k','linewidth',3)
%  cc=colorbar;
m_grid('xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',10);
title('Campo de \Psi - Modelo de Stommel (m^2 s^{-1})','fontweight','bold')
set(gcf,'color','w')
%  daspect([.7 2.5 .7])


print -depsc giros_rot.eps
!epstopdf giros_rot.eps


 toc








