%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plota seções verticais de velocidade via metodo dinamico
% e salva matrizes de gpan para usar posteriormente nos 
% calculos de psi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc;

%%% LIMITES OESTE-LESTE e SUL-NORTE DA AREA DE INTERESSE %%%

lonlim = [-41 -33];
latlim = [-20 -10];

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

%%% PROFUNDIDADE DE PLOTAGEM [m] %%%

pi = 100; 

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


	Tw = Twoam(near(Pw,pi,1),:);
	Sw = Swoam(near(Pw,pi,1),:);

f=find(diff(lonw)>0);

lonwr = reshape(lonw,f(1),length(lonw)/f(1));
latwr = reshape(latw,f(1),length(lonw)/f(1));
Tw = reshape(Tw,f(1),length(lonw)/f(1));
Sw = reshape(Sw,f(1),length(lonw)/f(1));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLOTAGEM CAMPOS TEMPERATURA E SALINIDADE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m_proj('mercator','long',[min(lonlim) max(lonlim)],'lat',[min(latlim) max(latlim)],'on');

%%% usando batimetria do etopo

load ../../common/etopo2_E_751x901.xyz;
etopo = etopo2_E_751x901;
xbat = etopo(:,1);
ybat = etopo(:,2);
zbat = etopo(:,3);

xbat = reshape(xbat,751,901);
ybat = reshape(ybat,751,901);
zbat = reshape(zbat,751,901);

fx = find(xbat(:,1) >= lonlim(1) & xbat(:,1) <= lonlim(2));
fy = find(ybat(1,:) >= latlim(1) & ybat(1,:) <= latlim(2));

xbat = xbat(fx,fy);
ybat = ybat(fx,fy);
zbat = zbat(fx,fy);

isob=[-200 -200];


% PLOTANDO AS DISTRIBUICOES TERMOHALINAS HORIZONTAIS ********************************************

Dw=sw_pden(Sw,Tw,pi,0)-1000;
P = ['TwSwDw'];

lT=21:0.05:26; %% para superficie %%
lS=36.6:0.03:37.4;  %% para superficie %%
lD=24.3:0.03:25.9; %% para superficie %%
l = ['lTlSlD'];

t1 = ['Temperatura ( \circ C) - Levitus - ',num2str(Pw(near(Pw,pi,1))),' m'];
t2 = ['Salinidade - Levitus - ',num2str(Pw(near(Pw,pi,1))),' m'];
t3 = ['Densidade Potencial (kg m^{-3}) - Levitus - ',num2str(Pw(near(Pw,pi,1))),' m'];
tit = ['t1t2t3'];

c=0;
for k = 1:3
  figure(k)
  set(k,'Color','w');
  hold on;
  eval(['[c1,h1] = m_contourf(lonwr,latwr,',P(c+k:c+k+1),',',l(c+k:c+k+1),');']); 
  shading flat; hold on
  %  clabel(c1,h1); % para aparecer os numeros no mapa %
  [c2,h2] = m_contour(xbat,ybat,zbat,isob,'k');
  clabel(c2,h2,'labelspacing',500)
  m_usercoast('costa_leste.mat','patch',[.7 .7 .7],'LineStyle','-');
  %  m_grid('box','fancy','yaxislocation','right','xaxislocation','top','fontsize',10);
  m_grid('box','fancy','yaxislocation','left','xaxislocation','bottom','fontsize',10);
  cc = colorbar;
  pos = get(cc,'Position');
  set(cc,'Position',[pos(1)/1.1 pos(2) pos(3)/2 pos(4)])
  hold off;
  eval(['title(',tit(c+k:c+k+1),',''fontsize'',12,''fontweight'',''bold'')']);

     print(k,'-depsc',['../figuras/',P(c+k),'_levitus_',num2str(pi),'m']);
     eval(['!epstopdf ../figuras/',P(c+k),'_levitus_',num2str(pi),'m.eps'])
%       eval(['!rm -rf ../figuras/',P(c+k),'_levitus_',num2str(pi),'m.eps'])
c=c+1; 
end
clear k c

% ***********************************************************************************************





