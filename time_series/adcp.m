%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     PROGRAMA PARA PROCESSAMENTO E VISUALIZACAO 
%         DOS DADOS DE ADCP - SAO SEBASTIAO
%            maio/2006 - Observacional
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  - para visualização e conversão dos dados de ADCP, instalar o programa VIEW
%    . para transformar em dados ASCII, usar a opção File / Export
%  
%  - tipos de arquivos:
%  
%  *.adp - arquivos binário contendo os dados de corrente - devem ser abertos pelo programa VIEW
%  *.gga - dados de tempo e posição fornecidos pelo GPS
%  *.gps - total dos dados de gps (não deve ser utilizado - informações redundantes)
%   ve --> vel. zonal 
%   vn--> vel. merdional
%          ------ PRODUTOS ------
%
% - plota vetores de velocidade para cada radial e cada camada
% - plota mapas com os vetores em cada radial para cada camada individuamente
% - plota mapas filtrados, interpolados para uma grade curvilinear atraves de AO, para cada camada
%%

clear all;close all;clc

% declinação magnetica
decl=-22;
% celulas de profundidade vai de 4 até 80
cel_prof=4:4:80;cel_prof=cel_prof';


%%%%  LEITURA DOS DADOS - o "eval" é usado para usar variaveis-coringa para os nomes.

for i=1:6    % 6 radiais

	% comando eval executa o q está escrito
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
	
	% corrige a declinaçao magnetica
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




%%% Plotagem dos Stickplots com cada camada em cada radial

for i=1:6
    figure
    for j=1:4
	eval(['Teste_e = weim(7,''hann'',RadC' num2str(i) '_ve(:,' num2str(j) '));'])
	eval(['Teste_n = weim(7,''hann'',RadC' num2str(i) '_vn(:,' num2str(j) '));'])
	
	% plota pra cada radial os vetores de velocidade
	subplot(4,1,j)
	if i==1
	   estick(1:length(Teste_e),Teste_e./100,Teste_n./100,120,'m.s^{-1}');
	elseif i==2
	   estick(1:length(Teste_e),Teste_e./100,Teste_n./100,110,'m.s^{-1}');
	elseif i==3
	   estick(1:length(Teste_e),Teste_e./100,Teste_n./100,40,'m.s^{-1}');
	elseif i==4
	   estick(1:length(Teste_e),Teste_e./100,Teste_n./100,40,'m.s^{-1}');
	else
	   estick(1:length(Teste_e),Teste_e./100,Teste_n./100,60,'m.s^{-1}');
	end
	
	title(['Radial ' num2str(i) ' - Camada ' num2str(j) ': ' num2str(cel_prof(j)) ' m'])

    end
end


%  PLOTANDO MAPAS


format bank

load ../batimetria/estacoes.dat; 

% separa as estacoes
xest=estacoes(:,1);
yest=estacoes(:,2);

load ../batimetria/batimetria_GMS.dat; % carrega batimetria para compor o fundo

ssb = batimetria_GMS;
lon=ssb(:,1);
lat=ssb(:,2);
z=ssb(:,3);  % profundidade

% limites para a projecao
lonlim=[min(lon) max(lon)];
latlim=[min(lat) max(lat)];

% monta a projecao a partir dos limites de lat e long 
% (nessa latitude é melhor a projeção de mercator--> cilindrica)

m_proj('mercator','long',[lonlim(1) lonlim(2)],'lat',[latlim(1) latlim(2)],'on');

%%% CRIANDO GRADE PARA INTERPOLAR BATIMETRIA %%%%%%%%
% varia de 0.0003°
x=min(min(lon)):0.0003:max(max(lon));
y=min(min(lat)):0.0003:max(max(lat));

% meshgrid--> monta a grade (ver help)
[xc,yc]=meshgrid(x,y);
[x,y,zc]=griddata(lon,lat,z,xc,yc); % interpolando as posições
[xc,yc]=m_ll2xy(xc,yc); % muda o formato das corrdenadas para plotar


%cria um conjunto vazio
latf=[];U=[];V=[];lonf=[];

j=1; % numero do bin (1=superficie, 4=fundo)

%%% looping para montar os vetores u e v usado nos mapas
for i=1:6
	% monta so um vetor pra u, v, lat e long pra todas as radiais
	eval(['lat2=Lat_rad' num2str(i) ';'])
	eval(['latf=[latf  lat2''];'])
	clear lat2
	eval(['lon2=Lon_rad' num2str(i) ';'])
	eval(['lonf=[lonf  lon2''];'])
	clear lon2
	
	% volta pra u e v
	eval(['U=[U RadC' num2str(i) '_ve(:,' num2str(j) ')''];'])
	eval(['V=[V RadC' num2str(i) '_vn(:,' num2str(j) ')''];'])

end

%%% escala para os vetores

% calcula-se o modulo maximo entre todos os vetores e adiciona-se ao final dos mesmos
% para plotar o quiver com a mesma escala

mod_bruto=sqrt(U.^2+V.^2);
max_bruto=max(mod_bruto);
U2=[U max_bruto]; V2=[V 0]; lonf2=[lonf -45.5201282365647]; latf2=[latf -23.8861931150716];

format short g

esc1=([num2str(round(max_bruto)/100),' m s^{-1}']); % strig para a escala

% carregando vento medio no dia do ADCP (mu e mv) - nao foi feito, usamos o quickscat
%  
%  load vento_medio_adcp.mat;  
%  
%  xmu=0.00091961;
%  ymu=-0.42906;
%  
%  [xmu,ymu]=m_xy2ll(xmu,ymu);

% carregando serie de vento de 2 dias antes e do dia do ADCP

%  load serie_vento_adcp.mat;


figure

  contourf(xc,yc,zc,min(min(zc)):max(max(zc)));shading flat;  % shading flat-->forma de apresentaçao da batimetria
  caxis([-100 200])     % arruma a escala de cores
  hold on
  m_quiver(lonf2,latf2,U2,V2,'y') 
  m_usercoast('../ctd/costa_ssebastiao.mat','patch',[.6 .6 .6],'LineStyle','-');
%    m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','right','xaxislocation','top','fontsize',12);
  m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',12);
  title('Mapeamento - ADCP','fontsize',12,'fontweight','bold')
  m_text(-45.5201282365647,-23.8785101740827,esc1)
%    m_quiver(xmu,ymu,mu,mv,0.007,'w')
%    m_text(-45.406,-23.871,'Direcao do vento medio','color','y') 
%    m_text(-45.402,-23.878,'durante a perfilagem','color','y')
  set(gcf,'color','w')

print -depsc adcp_bruto.eps
!epstopdf adcp_bruto.eps

%%% montando retas com o trajeto do ADCP, usando ginput, para ficarem RETAS e nao CURVAS ---------------------------

x1_r1= -0.00137177810330082;   y1_r1= -0.428745387810814;
x2_r1= -0.000559098870610355;  y2_r1= -0.429533625562822;

x1_r2= -0.000583540351292775;  y1_r2= -0.428494862633819;
x2_r2= -0.000161924809521032;  y2_r2= -0.429154782612245;

x1_r3= 0.000204697400715267;   y1_r3= -0.428427648561943;
x2_r3= 0.000375787765492206;   y2_r3= -0.428678173738938;

x1_r4= 0.000840175898458184;   y1_r4= -0.427755507843176;
x2_r4= 0.00112736329647662;    y2_r4= -0.427786059694029;

x1_r5= 0.00100515589306452;    y1_r5= -0.427217795268163;
x2_r5= 0.001353446992789;      y2_r5= -0.427211684897992;

x1_r6= 0.00118846699818267;    y1_r6= -0.426423447145984;
x2_r6= 0.00154286846807776;    y2_r6= -0.426496771588032;

%%% INTERPOLANDO OS DADOS DE ADCP PARA A GRADE CURVILINEAR USANDO INTERPOLACAO LINEAR
load iso8.mat  % lat lon das isobatas de 8 m 

% convertendo unidades para plotar
[xx1,yy1]=m_ll2xy(x1,y1);
[xx2,yy2]=m_ll2xy(x2,y2);

%%% carregando a grade construida no seagrid -  mesmo procedimento de pl_hor.m

load grade_ssebastiao_adcp.txt;
grid=grade_ssebastiao_adcp;
xg=grid(:,9); xg=xg';
yg=grid(:,10); yg=yg';
li=max(grid(:,1));
lj=max(grid(:,2));

%%% tirando NaN

f=find(not(isnan(U)));

U=U(f); V=V(f); latf=latf(f); lonf=lonf(f);

clc % limpa a tela de comando do matlab
disp('Interpolacao Linear.........')
Ui=griddata(lonf,latf,U,xg,yg);
Vi=griddata(lonf,latf,V,xg,yg);

%%% escala para os vetores - mesmo procedimento dos brutos

mod_linear=sqrt(Ui.^2+Vi.^2);
max_linear=max(mod_linear);
Ui2=[Ui max_linear]; Vi2=[Vi 0]; xg2=[xg -45.5201282365647]; yg2=[yg -23.8861931150716];

format short g

esc2=([num2str(round(max_linear)/100),' m s^{-1}']);


figure 

  contourf(xc,yc,zc,min(min(zc)):max(max(zc)));shading flat;  % shading flat-->forma de apresentaçao da grade
  caxis([-100 200])     % arruma a escala de cores
  hold on
  plot([x1_r1 x2_r1],[y1_r1 y2_r1],'y',[x1_r2 x2_r2],[y1_r2 y2_r2],'y',[x1_r3 x2_r3],[y1_r3 y2_r3],'y',[x1_r4 x2_r4],[y1_r4 y2_r4],'y',[x1_r5 x2_r5],[y1_r5 y2_r5],'y',[x1_r6 x2_r6],[y1_r6 y2_r6],'y')
  m_quiver(xg2,yg2,Ui2,Vi2,'k')
  fill(xx1,yy1,[.8 .8 .8]);
  fill(xx2,yy2,[.8 .8 .8]);
  m_usercoast('../ctd/costa_ssebastiao.mat','patch',[.4 .4 .4],'LineStyle','-');
%    m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','right','xaxislocation','top','fontsize',12);
  m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',12);
  title('Mapeamento - ADCP - Interpolacao Linear','fontsize',12,'fontweight','bold')
  plot([-0.00137177810330082 -0.00109070107545299],[-0.426289019002231 -0.426289019002231],'y','linewidth',2)
  m_text(-45.4847682959388,-23.713,'Trajetoria','color','y')
  m_text(-45.4847682959388,-23.72,'do ADCP','color','y')
  m_text(-45.5201282365647,-23.8785101740827,esc2) 
%    m_quiver(xmu,ymu,mu,mv,0.007,'w')
%    m_text(-45.406,-23.871,'Direcao do vento medio','color','y') 
%    m_text(-45.402,-23.878,'durante a perfilagem','color','y')
  set(gcf,'color','w')

print -depsc adcp_linear.eps
!epstopdf adcp_linear.eps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INTERPOLANDO OS DADOS DE ADCP PARA A GRADE CURVILINEAR USANDO ANALISE OBJETIVA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% aplicando condicao de contorno de NAO-escorregamento -------------------------

% subamostrando (decimando) o contorno, tornando-o compatível a matriz de amostragem do ADCP

xcon1=[]; ycon1=[];
for i=1:30:length(x1) % o intervalo de decimacao foi escolhido por tentativa e erro
    xcon1=[xcon1 x1(i)];
    ycon1=[ycon1 y1(i)];
end

xcon2=[]; ycon2=[];
for j=1:30:length(x2)
    xcon2=[xcon2 x2(j)];
    ycon2=[ycon2 y2(j)];
end

% acrescentando o contorno nas matrizes de U e V 
% (adiciona-se as corrdenadas da vel=0 ao final dos vetores, antes de interpolar!!!)

lonf=[lonf xcon1 xcon2];
latf=[latf ycon1 ycon2];
U=[U zeros(size(xcon1)) zeros(size(xcon2))]; 
V=[V zeros(size(ycon1)) zeros(size(ycon2))];

%%% entrando com o comprimento de correlacao e variancia do erro aleatorio
%%% lembrando que para obter o mapa de erro REAL, deve-se rodar a AO sem 
%%% os valores do contorno --> procedimento já realizado.

%%% nesse caso, foram encontrados os valores "a priori" - tentativa e erro,
%%% sempre comparando com o campo da int. linear
corrlen = 0.2; % este no caso, é o comprimento da grade usada
err = 0.01^2; % valor que proporcionou o melhor resultado

% lembrando que o mapa de erro independe da propriedade
% depende apenas da relacao "pontos amostrais x pontos de grade"

clc
disp('Interpolacao Objetiva.........')
[Uo,er] = scaloa(xg,yg,lonf,latf,U,corrlen,err);
[Vo,er] = scaloa(xg,yg,lonf,latf,V,corrlen,err);

%%% escala para os vetores - da mesma forma que anteriormente

mod_AO=sqrt(Uo.^2+Vo.^2);
max_AO=max(mod_AO);
Uo2=[Uo max_AO]; Vo2=[Vo 0]; 

format short g

esc3=([num2str(round(max_AO)/100),' m s^{-1}']);


figure 

  contourf(xc,yc,zc,min(min(zc)):max(max(zc)));shading flat;  % shading flat-->forma de apresentaçao da grade
  caxis([-100 200])     % arruma a escala de cores
  hold on
  plot([x1_r1 x2_r1],[y1_r1 y2_r1],'y',[x1_r2 x2_r2],[y1_r2 y2_r2],'y',[x1_r3 x2_r3],[y1_r3 y2_r3],'y',[x1_r4 x2_r4],[y1_r4 y2_r4],'y',[x1_r5 x2_r5],[y1_r5 y2_r5],'y',[x1_r6 x2_r6],[y1_r6 y2_r6],'y')
  m_quiver(xg2,yg2,Uo2,Vo2,'k');
  fill(xx1,yy1,[.8 .8 .8]);
  fill(xx2,yy2,[.8 .8 .8]);
  m_usercoast('../ctd/costa_ssebastiao.mat','patch',[.4 .4 .4],'LineStyle','-');
%  m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','right','xaxislocation','top','fontsize',12);
  m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',12);
  title('Mapeamento - ADCP - Analise Objetiva','fontsize',12,'fontweight','bold')
  plot([-0.00137177810330082 -0.00109070107545299],[-0.426289019002231 -0.426289019002231],'y','linewidth',2)
  m_text(-45.4847682959388,-23.713,'Trajetoria','color','y')
  m_text(-45.4847682959388,-23.72,'do ADCP','color','y')
  m_text(-45.5201282365647,-23.8785101740827,esc3)
%    m_quiver(xmu,ymu,mu,mv,0.007,'w')
%    m_text(-45.406,-23.871,'Direcao do vento medio','color','y') 
%    m_text(-45.402,-23.878,'durante a perfilagem','color','y')
  set(gcf,'color','w')

print -depsc adcp_AO.eps
!epstopdf adcp_AO.eps


%%% plotando series de vento nos dias proximos ao adcp - nao foi feito, usamos o quickscat

%  hora20=0:0.167:24; hora21=hora20;
%  y20=0*ones(size(hora20)); y21=y20;
%  
%  hora22=0:0.336:24;
%  y22=0*ones(size(hora22));

%  figure
%  
%  subplot(311); quiver(hora20,y20,u20,v20); 
%  title('Direcao e intensidade do vento em 20/03')
%  subplot(312); quiver(hora21,y21,u21,v21);
%  ylabel('Intensidade do vento (m s^{-1})')
%  title('Direcao e intensidade do vento em 21/03')
%  subplot(313); quiver(hora22,y22,u22,v22);
%  title('Direcao e intensidade do vento em 22/03 (dia da perfilagem)')
%  xlabel('Tempo (horas)')
%  set(gcf,'color','w')
%  
%  print -depsc vento_3dias.eps
%  !epstopdf vento_3dias.eps
 




