%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   WORLD OCEAN ATLAS
% Plota seções verticais de velocidade via metodo dinamico
% e salva matrizes de gpan para usar posteriormente nos 
%                    calculos de psi
%      Rafael Guarino Soutelino - Mestrado IOUSP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc;

%%% LIMITES OESTE-LESTE e SUL-NORTE DA AREA DE INTERESSE %%%
lonlim = [-41 -33.7]; 
latlim = [-20.5 -10.2];
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

  eval(['load ../../../dados/WOA2001/',Time{i},'/Twoa_',Time{i},'_brazilcoast.dat;']);
  eval(['load ../../../dados/WOA2001/',Time{i},'/Swoa_',Time{i},'_brazilcoast.dat;']);

  eval(['Twoa = Twoa_',Time{i},'_brazilcoast;']);
  eval(['Swoa = Swoa_',Time{i},'_brazilcoast;']);
  
  eval(['clear Twoa_',Time{i},'_brazilcoast;']);
  eval(['clear Swoa_',Time{i},'_brazilcoast;']);
  
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
nr = 1000;
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
	      disp('Perfil nao atinge NR')
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
f0=sw_f(-15);
psig=gpan/f0; psig=-psig;
psig=psig-mean(psig);

% salvando campos brutos para calcular a funcao de correlacao
save psi_woa.mat psig lon lat

%%% carregando a grade construida no seagrid

load ../../common/seagrid_leste2.dat;
grid = seagrid_leste2;
xg = grid(:,9); 
yg = grid(:,10);

lj = max(grid(:,1));
li = max(grid(:,2));

xgi = reshape(xg,li,lj);
ygi = reshape(yg,li,lj);

%%% interpolacao linear com grade curvilinear %%%%%%%%%%%%%

disp('Interpolacao Linear.........')
psigi = griddata(lon,lat,psig,xgi,ygi);
[ui,vi]=psi2uv(xgi,ygi,psigi);

%%% ANALISE OBJETIVA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xg=xg';yg=yg';
corrlen = 3; 
err = 0.005; 

%%% APLICANDO CONDICOES DE CONTORNO PARA PSI ********************************

load ../../common/etopo2_leste.mat;

% buscando valor de lon e lat das isobatas de ?? m
if pi <= 1000
   [c,h] = contour(xb,yb,zb,[-1000 -1000]);
else 
   [c,h] = contour(xb,yb,zb,[-pi -pi]);
end
xcont = get(h(1),'xdata'); xcont = xcont(find(not(isnan(xcont))));
ycont = get(h(1),'ydata'); ycont = ycont(find(not(isnan(ycont))));
%  xcont = weim(31,'hann',xcont);
%  ycont = weim(31,'hann',ycont);
xcont=xcont';
ycont=ycont';

close(1); 

xf=[min(xcont) xcont];
yf=[max(ycont) ycont];
[xf,yf] = m_ll2xy(xf,yf,'patch');

 % enxertando valores da isobata nos vetores
lon = [lon xcont]; lat = [lat ycont];
psig = [psig zeros(size(xcont))];

disp('Interpolacao Objetiva de Psi')
[psigo,er] = scaloa(xg,yg,lon,lat,psig,corrlen,err);
er=100*sqrt(er); 

psigo=reshape(psigo,li,lj);
er=reshape(er,li,lj);

% calculando componentes u e v
[uo,vo]=psi2uv(xgi,ygi,psigo);

     
%%% PLOTANDO OS RESULTADOS **********************************************************
isob=[-100 -100];
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

fc = 2; % fator de escala dos vetores

% construindo escalas de contorno para psi

lpsi = -40000:800:30000;
inc = ( max(max(psigo)) - min(min(psigo)) ) / 17;
lpsi2 = -30000:inc:30000;
int = num2str(100*round(inc/100));



figure(4)
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

m_contourf(xgi,ygi,psigo,lpsi);hold on;shading flat;
caxis([-30000 30000])
m_contour(xgi,ygi,psigo,lpsi2,'w')
m_quiver(xgi,ygi,uo*fc,vo*fc,0,'k')
fill(xf,yf,[0.8 0.8 0.8]);
%  plot(xf,yf,'color',[.8 .8 .8])
%  m_plot(lonctd,latctd,'w.','markersize',6)
m_usercoast('costa_leste.mat','patch',[0 0 0])
vmax = round(max(max(sqrt((uo*100).^2+(vo*100).^2))));
m_quiver(-40,-11,0.2*fc,0,0,'w')
text = ['20 cm s^{-1}'];
m_text(-40,-11.3,text,'color','w','fontsize',10,'fontweight','bold')
text = [num2str(pi),' m'];
m_text(-34.8,-20.1,text,'color','r','fontweight','bold')
text = ['Contornos: ',int,' m^2s^{-1}'];
m_text(-40.5,-11.8,text,'color','w','fontsize',8,'fontweight','bold')
m_grid('box','fancy','yaxislocation','left','xaxislocation','bottom','fontsize',8);
tit = ['\Psi - WOA - ',time,' [m^2 s^{-1}]'];
title(tit,'fontsize',10,'fontweight','bold');
cc = colorbar;
pos = get(cc,'Position');
set(cc,'Position',[pos(1) pos(2)*1.3 pos(3)/2.5 pos(4)/1.1])
   print(4,'-depsc',['../figuras/psi_woa_',time,'_',num2str(pi),'m']);
   eval(['!epstopdf ../figuras/psi_woa_',time,'_',num2str(pi),'m.eps'])
drawnow







