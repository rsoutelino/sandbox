%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   WORLD OCEAN ATLAS
% Plota seções verticais de velocidade via metodo dinamico
% e salva matrizes de gpan para usar posteriormente nos 
%                    calculos de psi
%      Rafael Guarino Soutelino - Mestrado IOUSP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc;

%%% LIMITES OESTE-LESTE e SUL-NORTE DA AREA DE INTERESSE %%%
lonlim = [-51 -33];
latlim = [-33 -20];
m_proj('mercator','long',[min(lonlim) max(lonlim)],'lat',[min(latlim) max(latlim)],'on');
%%% PERIODO DESEJADO %%%

disp(' ');
disp(['  ### Escolha do Periodo ###']);
disp(' ');
disp(['  - Para o campo anual: [anual]']);
disp(['  - Para o campo mensal: [jan] ou [fev] ou ...']);
disp(['  - Para media dos campos mensais: [jan,fev,mar,...]']);
disp(' ');
time = input(['  Periodo: '],'s');

time([1 end]) = [];

ft = find(time==',');

if isempty(ft) == 1
  Time = cellstr(time);
else
  for i = 1:length(ft)
    Time{i} = time(ft(i)-3:ft(i)-1);
  end
  Time{i+1} = time(ft(end)+1:ft(end)+3);
end


%%% LEITURA DA CLIMATOLOGIA DECLARADA %%%

for i = 1:length(Time)

  eval(['load ../../../dados/WOA2001/',Time{i},'/Twoa_',Time{i},'_atlantic.dat;']);
  eval(['load ../../../dados/WOA2001/',Time{i},'/Swoa_',Time{i},'_atlantic.dat;']);

  eval(['Twoa = Twoa_',Time{i},'_atlantic;']);
  eval(['Swoa = Swoa_',Time{i},'_atlantic;']);
  
  eval(['clear Twoa_',Time{i},'_atlantic;']);
  eval(['clear Swoa_',Time{i},'_atlantic;']);
  
  lat = Twoa(:,1)';
  lon = Twoa(:,2)';
  
  flat = find(lat >= min(latlim) & lat <= max(latlim));
  flon = find(lon >= min(lonlim) & lon <= max(lonlim));


  ii = 1;
  for jj = 1:length(flon)
    if isempty(find(flat == flon(jj))) == 0
      q(ii) = flon(jj);
      ii = ii+1;
    end
  end

  eval(['T',Time{i},' = Twoa(q,3:end);']);
  eval(['S',Time{i},' = Swoa(q,3:end);']);

end

%%% MEDIA DOS CAMPOS DESEJADOS %%%

strT = [];
strS = [];
for i = 1:length(Time)
  strT = [strT,'T',Time{i},'+'];
  strS = [strS,'S',Time{i},'+'];
end
strT(end) = [];
strS(end) = [];

eval(['Twoam = ((',strT,')/',num2str(i),')'';']);
eval(['Swoam = ((',strS,')/',num2str(i),')'';']);

latw = lat(q);
lonw = lon(q);

%%% MANIPULACAO DOS CAMPOS %%%

if size(time,2) == 5
  Pw = [0;10;20;30;50;75;100;125;150;200;250;300;400;500;600;700;...
        800;900;1000;1100;1200;1300;1400;1500;1750;2000;2500;3000;...
        3500;4000;4500;5000;5500];
else
  Pw = [0;10;20;30;50;75;100;125;150;200;250;300;400;500;600;700;...
        800;900;1000;1100;1200;1300;1400;1500];
end

pplot = 3000; % profmax
pf = near(Pw,pplot,1);
pi = input('Escolha a profundidade: ') % profundidade de interesse para plotar psi
pfi = near(Pw,pi,1);
z = Pw(1:pf);
nr = 560;
fnr = near(z,nr,1);


clear Twoa Swoa Sanual Tanual f ft jj ii q strT strS


%%% MONTANDO MATRIZES DE LAT LON
% descobrindo tamanho da grade

f = find(diff(lonw)~=0);
li = f(2)-f(1);
lj = length(f)+1;
lon = reshape(lonw,li,lj);
lat = reshape(latw,li,lj);

gpan = []; LON = []; LAT = []; Ti = []; Si = [];

%começa looping para as radiais WOA
for m=1:li	

	for j=1:lj 
           la = lat(m,1);
           fla = find(latw==la);
           T(:,j) = Twoam(1:pf,fla(j));
	   S(:,j) = Swoam(1:pf,fla(j));
 	   
	end
        lons = lon(1,:);
        lats = lat(m,:);
        Ti = []; Si = [];
        
 	% retirando os perfis que nao atingem o NR
        for b = 1:length(T(1,:))
           t = T(:,b);
           s = S(:,b);
           if pfi < fnr               
	       fn = find(isnan(t(1:fnr))==1);
           else
               fn = find(isnan(t(1:pfi))==1);
           end
           if isempty(fn) == 1
              Ti = [Ti t];
              Si = [Si s];
  	      LON = [LON lons(b)]; 
              LAT = [LAT lats(b)]; 
	   else
%  	      disp('Perfil nao atinge NR')
              clear t s
           end	
	end
        %%% CALCULANDO ANOMALIA DO GEOPOTENCIAL
	% anomalia do geopotencial relativa a superficie
	agp = sw_gpan(Si,Ti,z);
	
	% anomalia do geopotencial relativa ao NR
	agp = agp - ones(size(z)) * agp(fnr,:);
	
	gpan = [gpan agp(pfi,:)];
        clear Ti Si lons lats
        
end

lon = LON;
lat = LAT;
f0=sw_f(-25);
psig=gpan/f0; psig=-psig;
psig=psig-mean(psig);
if pi > nr
   psig = -psig;
end

%%% usando grade curvilinea do seagrid

load seagrid_ceres_woa.grd;

xg = seagrid_ceres_woa(:,9);
yg = seagrid_ceres_woa(:,10);

lj = max(seagrid_ceres_woa(:,1));
li = max(seagrid_ceres_woa(:,2));

xgi = reshape(xg,li,lj);
ygi = reshape(yg,li,lj);

%%% criando grade para interpolar via AO

%  load ../levitus_larga.grd;
%  grid = levitus_larga;
%  
%  xg = lonlim(1):0.3:lonlim(2);
%  yg = latlim(1):0.3:latlim(2);
%  
%  [xgi,ygi] = meshgrid(xg,yg);
%  
%  [li,lj] = size(xgi);
%  
%  xg = reshape(xgi,1,li*lj);
%  yg = reshape(ygi,1,li*lj);

%%% interpolacao linear com grade curvilinear %%%%%%%%%%%%%

disp('Interpolacao Linear.........')
psigi = griddata(lon,lat,psig,xgi,ygi);
[ui,vi]=psi2uv(xgi,ygi,psigi);

%%% ANALISE OBJETIVA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xg=xg';yg=yg';
corrlen = 4.3; 
err = 0.005; 

%%% APLICANDO CONDICOES DE CONTORNO PARA PSI ********************************

% carregando a isobata de 1000 do tbase
[zb,xb,yb] = m_tbase([lonlim(1) lonlim(2) latlim(1) latlim(2)]);
[c,h] = contour(xb,yb,zb,[-1000 -1000]);

xcont = get(h(1),'xdata'); xcont = xcont(find(not(isnan(xcont))));
ycont = get(h(1),'ydata'); ycont = ycont(find(not(isnan(ycont))));
xcont = weim(31,'hann',xcont);
ycont = weim(31,'hann',ycont);
%  xcont=xcont';
%  ycont=ycont';

close(1); 

xf=[min(xcont)-10 xcont min(xcont)-10];
yf=[min(ycont) ycont max(ycont)];
[xf,yf] = m_ll2xy(xf,yf,'patch');

 % enxertando valores da isobata nos vetores
lon = [lon xcont]; lat = [lat ycont];
psig = [psig zeros(size(xcont))];

disp('Interpolacao Objetiva de Psi')
[psigo,er] = scaloa(xg,yg,lon(1:10:end),lat(1:10:end),psig(1:10:end),corrlen,err);
er=100*sqrt(er); 

psigo=reshape(psigo,li,lj);
er=reshape(er,li,lj);

% calculando componentes u e v
[uo,vo]=psi2uv(xgi,ygi,psigo);

     
%%% PLOTANDO OS RESULTADOS **********************************************************
isob=[-1000 -1000];
% Distribuicoes termohalinas horizontais:------------------------------------------------------

P = ['ToSoDo'];

lT=26.8:0.05:28.4; %% para superficie %%
lS=36.6:0.03:37.4;  %% para superficie %%
lD=23.3:0.03:24.3; %% para superficie %%
l = ['lTlSlD'];

%  t1 = ['Temperatura ( \circ C) - OEII - ',nn,' m'];
%  t2 = ['Salinidade - OEII - ',nn,' m'];
%  t3 = ['Densidade Potencial (kg m^{-3}) - OEII - ',nn,' m'];
%  tit = ['t1t2t3'];
%  
%  c=0;
%  for k = 1:3
%    figure(k)
%    set(k,'Color','w');
%    hold on;
%    eval(['[c1,h1] = m_contourf(xgi,ygi,',P(c+k:c+k+1),',',l(c+k:c+k+1),');']); 
%    shading flat; hold on
%    %  clabel(c1,h1); % para aparecer os numeros no mapa %
%  %    [c2,h2] = m_contour(xb,yb,zb,isob,'k');
%  %    clabel(c2,h2,'labelspacing',500)
%    fill(xf,yf,[0.8 0.8 0.8]);
%    m_usercoast('costa_leste.mat','patch',[0 0 0],'LineStyle','-');
%    %  m_grid('box','fancy','yaxislocation','right','xaxislocation','top','fontsize',10);
%    m_grid('box','fancy','yaxislocation','left','xaxislocation','bottom','fontsize',8);
%    cc = colorbar;
%    pos = get(cc,'Position');
%    set(cc,'Position',[pos(1)/1.2 pos(2) pos(3)/2.5 pos(4)])
%    hold off;
%    eval(['title(',tit(c+k:c+k+1),',''fontsize'',12,''fontweight'',''bold'')']);
%  drawnow
%       print(k,'-depsc',['figuras/',P(c+k),'_OEII_',nn,'m']);
%       eval(['!epstopdf figuras/',P(c+k),'_OEII_',nn,'m.eps'])
%       eval(['!rm -rf figuras/',P(c+k),'_OEII_',nn,'m.eps'])
%  c=c+1; 
%  end
%  clear k c

% ***********************************************************************************************

% carregando a grade ceres I para plotar aqui

    load posicoes_ceres1_nav.dat
    latg = posicoes_ceres1_nav(:,1);
    latm = posicoes_ceres1_nav(:,2);
    long = posicoes_ceres1_nav(:,3);
    lonm = posicoes_ceres1_nav(:,4);

    lonc = long + lonm/60; lonc = -lonc;
    latc = latg + latm/60; latc = -latc;

fc = 10; % fator de escala dos vetores

% construindo escalas de contorno para psi

lpsi = -50000:800:50000;
inc = ( max(max(psigo)) - min(min(psigo)) ) / 50;
lpsi2 = -50000:inc:50000;
int = num2str(100*round(inc/100));

figure(4)
set(gcf,'color','w')
m_contourf(xgi,ygi,psigo,lpsi);hold on;shading flat;
caxis([-50000 50000])
[c,h] = m_contour(xgi,ygi,psigo,lpsi2,'k');
set(h,'color',[.7 .7 .7])
m_quiver(xgi,ygi,uo*fc,vo*fc,0,'k')
fill(xf,yf,[0.8 0.8 0.8]);
m_plot(lonc,latc,'ko','markersize',3.5,'markerfacecolor','r','markeredgecolor','y')
m_usercoast('../../common/costa_larga.mat','patch',[0 0 0])
vmax = round(max(max(sqrt((uo*100).^2+(vo*100).^2))));
m_quiver(-48,-21,0.2*fc,0,0,'w')
text = ['20 cm s^{-1}'];
m_text(-48,-21.5,text,'color','w','fontsize',10,'fontweight','bold')
text = [num2str(pi),' m'];
m_text(-48,-22,text,'color','y','fontweight','bold')
text = ['Contornos: ',int,' m^2s^{-1}'];
m_text(-49,-23,text,'color','w','fontsize',8,'fontweight','bold')
m_grid('box','fancy','yaxislocation','left','xaxislocation','bottom','fontsize',8);
tit = ['\Psi - WOA - ',time,' [m^2 s^{-1}]'];
title(tit,'fontsize',10,'fontweight','bold');
cc = colorbar;
pos = get(cc,'Position');
set(cc,'Position',[pos(1) pos(2)*1.6 pos(3)/2.5 pos(4)/1.2])
   print(4,'-depsc',['psi_ceres_woa_',time,'_',num2str(pi),'m']);
   eval(['!epstopdf psi_ceres_woa_',time,'_',num2str(pi),'m.eps'])
drawnow







