%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   PROGRAMA PARA CALCULAR VELOCIDADES GEOSTROFICAS
%   RELATIVAS COM BASE EM MEDIDAS DIRETAS, OU SEJA
%      REFERENCIAMENTO DO CALCULO GEOSTROFICO
%
%        Rafael Soutelino - Mestrado IOUSP
%                   jul/2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% este programa faz o procedimento atraves do uso 
% da interpolacao linear das vgeo relativas a 
% superficie e das velocidades geostroficas
% oriundas da filtragem dos dados de ADCP
clear all;close all;clc;warning off

% declarando os limites da area, projecao, etc...
m_proj('mercator','long',[-41 -33.7],'lat',[-20.5 -10.2],'on');
load ../common/etopo2_leste.mat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      P A R T E   I      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A PRIMEIRA PARTE DO PROGRAMA SE DESTINA A 
% GRADEAR OS DADOS DE ADCP PARA A GRADE CURVILINEA
% DA OEII, INTERPOLANDO LINEARMENTE

% carregando .mat com u,v,lon,lat dos dados de ADCP
load ../adcp/mat/adcp_filt.mat;

% carregando a grade
load ../common/seagrid_leste2.dat;
grid = seagrid_leste2;
xg = grid(:,9); 
yg = grid(:,10);

lj = max(grid(:,1));
li = max(grid(:,2));

xgi = reshape(xg,li,lj);
ygi = reshape(yg,li,lj);

clear seagrid_leste2 grid

% interpolando linearmente as velocidades do ADCP
disp('Interpolacao Linear.........')
Uadcp = griddata(lonadcp,latadcp,Uadcp,xgi,ygi);
Vadcp = griddata(lonadcp,latadcp,Vadcp,xgi,ygi);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     P A R T E   II      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% A SEGUNDA PARTE DO PROGRAMA SE DESTINA A 
% CALCULAR E GRADEAR OS CAMPOS DE VGEO CALCULADAS 
% COM BASE NOS DADOS TERMOHALINOS E REFERENCIADAS NA 
% SUPERFICIE

% carregando os .mat com as anomalias do geopotencial 
% previamente calculadas:

load posicoes_leste2.dat; pb=posicoes_leste2;
latg=pb(:,2); latm=pb(:,3); latctd=latg+latm/60; latctd=-latctd; latctd=latctd'; 
long=pb(:,4); lonm=pb(:,5); lonctd=long+lonm/60; lonctd=-lonctd; lonctd=lonctd'; 

clear latg latm long lonm posicoes_leste2 pb

% como as velocidades estao relativas a 1000 dbar, 
% e queremos que ela fique relativa a superficie,
% todos os niveis terao de ser subtraidos da superficie
% antes de serem referenciados pelo ADCP

% ---------------------------------------------------------
% buscando entao as velocidades em superficie

% carregando matrizes com os dados hidrograficos
gpan=[];T=[];S=[];lon=[];lat=[];
for k=1:12
   eval(['load ../hidrografia/mat/lesteII_agp_rad',num2str(k),'.mat']);
   gpan = [gpan agp(1,:)];
   lon = [lon lons];
   lat = [lat lats];  
end

clear agp Ti Si lons lats

%%% calculando PSI
f0=sw_f(-15);
psig=gpan/f0; psig=-psig;
psig=psig-mean(psig);

% interpolando linearmente PSI
disp('Interpolacao Linear.........')
psig = griddata(lon,lat,psig,xgi,ygi);

% calculando as velocidades geostroficas
[Usupf,Vsupf] = psi2uv(xgi,ygi,psig);

clear psig gpan f0 lon lat lons lats k 
% ----------------------------------------------------------


% ------------------------------------------------------
% calculando a velocidade o nivel desejado:

pp = input('Escolha a profundidade: '); 

% carregando matrizes com os dados hidrograficos
gpan=[];T=[];S=[];lon=[];lat=[];
for k=1:12
   eval(['load ../hidrografia/mat/lesteII_agp_rad',num2str(k),'.mat']);
   gpan = [gpan agp(pp,:)];
   lon = [lon lons];
   lat = [lat lats];  
end

clear agp Ti Si lons lats

%%% calculando PSI
f0=sw_f(-15);
psig=gpan/f0; psig=-psig;
psig=psig-mean(psig);

% interpolando linearmente PSI
disp('Interpolacao Linear.........')
psig = griddata(lon,lat,psig,xgi,ygi);

% calculando as velocidades geostroficas
[Uctd,Vctd] = psi2uv(xgi,ygi,psig);

clear gpan f0 lon lat lons lats k

%  figure(2)
%  subplot(141)
%  quiver(xgi,ygi,Uctd,Vctd,0)
%  title('Velocidade no nivel escolhido')

% referenciando na superficie:
Uctd = Uctd - Usupf;
Vctd = Vctd - Vsupf;
%  
%  subplot(142)
%  quiver(xgi,ygi,Usupf,Vsupf,0)
%  title('Velocidade na superficie')

%  subplot(143)
%  quiver(xgi,ygi,Uctd,Vctd,0)
%  title('Velocidade no nivel desejado referenciada na superficie')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     P A R T E   III     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% REFERENCIANDO AS VELOCIDADES A PARTIR DAS VELOCIDADES DE ADCP E 
% APLICANDO ANALISE OBJETIVA VETORIAL PARA OBTER A FUNCAO DE CORRENTE


% referenciando a partir dos dados de ADCP
U = Uadcp + Uctd;
V = Vadcp + Vctd;

%  subplot(144)
%  quiver(xgi,ygi,U,V,0)
%  title('Velocidade no nivel desejado referenciada pelo ADCP')

% aplicando analise objetiva vetorial para calcular o campo final
%%% INTERPOLACAO POR ANALISE OBJETIVA

% preparando campos
U = reshape(U,1,li*lj);
V = reshape(V,1,li*lj);

f = find(not(isnan(U)));
U = U(f); V = V(f); xg2 = xg(f); yg2 = yg(f);

lc = 1; % em graus, entao, precisa passar as velocidades para grau/s
E = 0.02;

U = U*9e-6; % fator de conversao para grau/s (9e-6)
V = V*9e-6;

% analise objetiva vetorial para obter psi ABSOLUTO
disp('ANALISE OBJETIVA VETORIAL')
[psi] = vectoa(xg',yg',xg2',yg2',U,V,lc,E,0);

% passando psi para m^2/s
psi = psi*111120; psi = psi*1e5;

%%% APLICANDO CONDICOES DE CONTORNO PARA PSI ********************************
%%% NAO ESCORREGAMENTO

% eliminando valor medio de psi
psi = psi-mean(psi);

% buscando valor de lon e lat das isobatas de ?? m
figure(1)
if pp <= 100
   [c,h] = contour(xb,yb,zb,[-100 -100]);
else 
   [c,h] = contour(xb,yb,zb,[-pp -pp]);
end
xcont = get(h(1),'xdata'); xcont = xcont(find(not(isnan(xcont))));
ycont = get(h(1),'ydata'); ycont = ycont(find(not(isnan(ycont))));
%  xcont = weim(31,'hann',xcont);
%  ycont = weim(31,'hann',ycont);
close(1); 

xf = [min(xcont); xcont];
yf = [max(ycont); ycont];
[xf,yf] = m_ll2xy(xf,yf,'patch');

%% aplicando condicao de contorno em psi
xg2 = [xg; xcont]; yg2 = [yg; ycont];
psi = [psi; zeros(size(xcont))];

% analise objetiva escalar em psiob para aplicar cond de contorno
disp('ANALISE OBJETIVA ESCALAR')
[psi] = scaloa(xg',yg',xg2',yg2',psi',lc*1.5,E);
psi = reshape(psi,li,lj);
[U,V] = psi2uv(xgi,ygi,psi);

% preparando para plotar
lpsi = -70000:800:70000;
inc = ( max(max(psi)) - min(min(psi)) ) / 20;
lpsi2 = -70000:inc:70000;
int = num2str(100*round(inc/100));
fc = 1;

figure(6)
set(gcf,...
        'Color',[1 1 1],...
        'InvertHardcopy','on',...
        'PaperUnits','inches',...
        'Units','inches',...
        'PaperOrientation','portrait',...
        'PaperPosition',[0 0 8.5 11],...
        'PaperPositionMode','manual',...
        'PaperType','usletter',...
        'Position',[.2 .2 8 9],...
        'ShareColors','off',...
        'Clipping','on');

m_contourf(xgi,ygi,psi,lpsi);shading flat;hold on
m_contour(xgi,ygi,psi,lpsi2,'w')
m_quiver(xgi,ygi,U*fc,V*fc,0,'k')
fill(xf,yf,[0.8 0.8 0.8]);
m_usercoast('costa_leste.mat','patch',[0 0 0])
m_quiver(-40,-11,0.5*fc,0,0,'w')
text = ['50 cm s^{-1}'];
m_text(-40,-11.3,text,'color','w','fontsize',10,'fontweight','bold')
text = [num2str(pp),' m'];
m_text(-34.8,-20.1,text,'color','r','fontweight','bold')
title('\Psi Geostrofico Referenciado por ADCP - OEII [m^2 s^{-1}]','fontsize',10,'fontweight','bold');
text = ['Contornos: ',int,' m^2s^{-1}'];
m_text(-40.5,-11.8,text,'color','w','fontsize',8,'fontweight','bold')
m_grid('box','fancy','yaxislocation','left','xaxislocation','bottom','fontsize',10);
cc = colorbar;
pos = get(cc,'Position');
set(cc,'Position',[pos(1) pos(2)*1.3 pos(3)/2.5 pos(4)/1.1])
   print(6,'-depsc',['figuras/psi_ref_adcp_',num2str(pp),'m']);
   eval(['!epstopdf figuras/psi_ref_adcp_',num2str(pp),'m.eps'])
drawnow


















