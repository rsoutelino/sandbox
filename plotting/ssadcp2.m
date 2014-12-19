% ADCP
%  - para visualização e conversão dos dados de ADCP, instalar o programa VIEW
%    . para transformar em dados ASCII, usar a opção File / Export
%  
%  - tipos de arquivos:
%  
%  *.adp - arquivos binário contendo os dados de corrente - devem ser abertos pelo programa VIEW
%  *.gga - dados de tempo e posição fornecidos pelo GPS
%  *.gps - total dos dados de gps (não deve ser utilizado - informações redundantes)

%ve --> vel. zonal 
%vn--> vel. merdional

clear all;close all;

% declinação
decl=-22;
% celulas de profundidade vai de 4 até 80
cel_prof=4:4:80;cel_prof=cel_prof';


% Carregando Dados

for i=1:6    % 6 radiais

% comando eval valida o q está escrito
eval(['load export/CURS' num2str(i) '.ve']);
eval(['Rad' num2str(i) '_ve= CURS' num2str(i) '(:,2:end);'])
eval(['a=find(Rad' num2str(i) '_ve==3276.7);'])
eval(['Rad' num2str(i) '_ve(a)= NaN;'])

eval(['clear a']);
eval(['clear CURS' num2str(i)]);


eval(['load export/CURS' num2str(i) '.vn']);
eval(['Rad' num2str(i) '_vn= CURS' num2str(i) '(:,2:end);'])
eval(['a=find(Rad' num2str(i) '_vn==3276.7);'])
eval(['Rad' num2str(i) '_vn(a)= NaN;'])

eval(['clear a']);
eval(['clear CURS' num2str(i)]);


%Calculando as direçoes das velocidades horizontais

eval(['vel' num2str(i) '=sqrt(Rad' num2str(i) '_ve.^2 + Rad' num2str(i) '_vn.^2);'])
eval(['dir' num2str(i) '=atan2(Rad' num2str(i) '_ve,Rad' num2str(i) '_vn);'])
eval(['dir' num2str(i) '=dir' num2str(i) '.*180/pi;'])

% corrige a declinaçao
eval(['dir' num2str(i) '=dir' num2str(i) '+decl;'])


eval(['r' num2str(i) '=find(dir' num2str(i) ' < 0);'])
eval(['dir' num2str(i) '(r' num2str(i) ')=dir' num2str(i) '(r' num2str(i) ') + 360;'])

eval(['rb' num2str(i) '=find(dir' num2str(i) ' > 360);'])
eval(['dir' num2str(i) '(rb' num2str(i) ')=dir' num2str(i) '(rb' num2str(i) ') - 360;'])

eval(['clear rb' num2str(i) ' r' num2str(i)])

eval(['dir' num2str(i) '= dir' num2str(i) '.*pi./180;'])
eval(['[RadC' num2str(i) '_vn,RadC' num2str(i) '_ve] = pol2cart(dir' num2str(i) ',vel' num2str(i) ');'])

eval(['load export/CURS' num2str(i) '.gps']);
eval(['Pos_rad' num2str(i) '= [CURS' num2str(i) '(:,3:4) CURS' num2str(i) '(:,6:7)];'])

eval(['Lat_rad' num2str(i) '= (Pos_rad' num2str(i) '(:,1) + Pos_rad' num2str(i) '(:,3))./2;'])
eval(['Lon_rad' num2str(i) '= (Pos_rad' num2str(i) '(:,2) + Pos_rad' num2str(i) '(:,4))./2;'])

eval(['clear CURS' num2str(i)]);

%  eval(['[dist_rad' num2str(i) ',phase_rad' num2str(i) ']=sw_dist(Lat_rad' num2str(i) ',Lon_rad' num2str(i) ', ''km'');'])
%  eval(['dist_rad' num2str(i) '= cumsum(dist_rad' num2str(i) ');'])
%  
%  eval(['[x' num2str(i) ',y' num2str(i) ']=meshgrid([0;dist_rad' num2str(i) '],cel_prof);']) 
%  eval(['intv' num2str(i) '=griddata([0;dist_rad' num2str(i) '],cel_prof,Rad' num2str(i) '_vn'',x' num2str(i) ',y' num2str(i) ',''linear'');'])
%  eval(['intu' num2str(i) '=griddata([0;dist_rad' num2str(i) '],cel_prof,Rad' num2str(i) '_ve'',x' num2str(i) ',y' num2str(i) ',''linear'');'])




end




%  Plotagem do Stickplots
for i=1:6
  figure
for j=1:4
eval(['Teste_e = weim(7,''hann'',RadC' num2str(i) '_ve(:,' num2str(j) '));'])
eval(['Teste_n = weim(7,''hann'',RadC' num2str(i) '_vn(:,' num2str(j) '));'])

% plota pra cada radial os vetores de velocidade
subplot(4,1,j)
if i==1
estick(1:length(Teste_e),Teste_e./100,Teste_n./100,120,'m.s{-1}');
elseif i==2
estick(1:length(Teste_e),Teste_e./100,Teste_n./100,110,'m.s{-1}');
elseif i==3
estick(1:length(Teste_e),Teste_e./100,Teste_n./100,40,'m.s{-1}');
elseif i==4
estick(1:length(Teste_e),Teste_e./100,Teste_n./100,40,'m.s{-1}');
else
estick(1:length(Teste_e),Teste_e./100,Teste_n./100,60,'m.s{-1}');
end

title(['Radial ' num2str(i) ' - Camada ' num2str(j) ': ' num2str(cel_prof(j)) ' m'])

end
end


%  PLOTANDO NO MAPA


format bank

load estacoes.dat;

% separa as estacoes
xest=estacoes(:,1);
yest=estacoes(:,2);

load batimetria_GMS.dat;

ssb = batimetria_GMS;
lon=ssb(:,1);
lat=ssb(:,2);
z=ssb(:,3);  % profundidade

% limites
lonlim=[min(lon) max(lon)];
latlim=[min(lat) max(lat)];

% faz a projecao a partir dos limites de lat e long 
% (nessa latitude é melhor a projeção de mercator--> cilindrica)
m_proj('mercator','long',[lonlim(1) lonlim(2)],'lat',[latlim(1) latlim(2)],'on');

%%% CRIANDO GRADE PARA INTERPOLAR %%%%%%%%
% varia de 0.0003°
x=min(min(lon)):0.0003:max(max(lon));
y=min(min(lat)):0.0003:max(max(lat));

% meshgrid--> monta a grade (ver help)
[xc,yc]=meshgrid(x,y);
[x,y,zc]=griddata(lon,lat,z,xc,yc); % interpolando as posições
[xc,yc]=m_ll2xy(xc,yc); % muda o formato

figure
% contorno preenchido
  contourf(xc,yc,zc,min(min(zc)):max(max(zc)));shading flat;  % shading flat-->forma de apresentaçao da grade
  caxis([-100 200])     % arruma a escala de cores
  hold on
%    m_usercoast('coast_ssebastiao.mat','patch',[.36 .43 .36],'LineStyle','-');
  m_gshhs_f('patch',[.6 .6 .6],'LineStyle','-');
  m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','right','xaxislocation','top','fontsize',12);
  m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',12);



  hold on
%cria um conjunto vazio
latf=[];E=[];V=[];lonf=[];
j=1;

for i=1:6
% monta so um vetor pra u, v, lat e long pra todas as radiais
eval(['lat2=Lat_rad' num2str(i) ';'])
eval(['latf=[latf  lat2''];'])
clear lat2
eval(['lon2=Lon_rad' num2str(i) ';'])
eval(['lonf=[lonf  lon2''];'])
clear lon2

% volta pra u e v
eval(['E=[E RadC' num2str(i) '_ve(:,' num2str(j) ')''];'])
eval(['V=[V RadC' num2str(i) '_vn(:,' num2str(j) ')''];'])

end

  hold on



% plotagem dos vetores
m_vec(200,lonf,latf,E,V,'k','shaftwidth',.5,'headangle',15)



%  for i=1:6
%  eval(['m_vec(200,Lon_rad' num2str(i) '(1:3:end),Lat_rad' num2str(i) '(1:3:end),RadC' num2str(i) '_ve(1:3:end,' num2str(j) '),RadC' num2str(i) '_vn(1:3:end,' num2str(j) '),''k'',''shaftwidth'',.5,''headangle'',15)'])
%  end


