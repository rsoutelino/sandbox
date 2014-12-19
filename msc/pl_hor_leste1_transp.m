%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    PROGRAMA PARA PLOTAR AS DISTRIBUICOES
%   VERTICAIS E CALCULAR TRANSPORTE A PARTIR
%   DA JUNCAO DOS CAMPOS HORIZONTAIS AO
%               Oceano Leste 2
%         out/2006 - Mestrado - IOUSP
%           Rafael Guarino Soutelino
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
clear all;close all;clc;warning off
tic

for sec = 1:49

% definindo os limites da area
lonlim = [-41 -33.7]; latlim = [-20.5 -10.2];
m_proj('mercator','long',[min(lonlim) max(lonlim)],'lat',[min(latlim) max(latlim)],'on');
%  m_gshhs_f('save','costa_leste.mat')

% lendo as posicoes das estacoes

load ../posicoes_leste1.dat; pb=posicoes_leste1;
latctd = pb(:,2); 
lonctd = pb(:,3);

clear posicoes_leste1

%%% carregando os dados hidrograficos %%%%%

%  n=input('Escolha a profundidade: '); 
n=1; nn=num2str(n);

% carregando matrizes com os dados hidrograficos
gpan=[];T=[];S=[];lon=[];lat=[];
for k=1:13
   eval(['load ../mat/lesteI_agp_rad',num2str(k),'.mat']);
   eval(['load ../mat/lesteI_T_rad',num2str(k),'.mat']);
   eval(['load ../mat/lesteI_S_rad',num2str(k),'.mat']);

   gpan = [gpan agp(n,:)];
   lon = [lon -lons];
   lat = [lat -lats];  
   T = [T Ti(n,:)];
   S = [S Si(n,:)];
end

clear agp Ti Si lons lats

%%% calculando PSI
f0=sw_f(-15);
psig=gpan/f0; psig=-psig;
psig=psig-mean(psig);
psigctd = psig;

%%% carregando a grade construida no seagrid

load ../../common/seagrid_leste2.dat;
grid = seagrid_leste2;
xg = grid(:,9); 
yg = grid(:,10);

lj = max(grid(:,1));
li = max(grid(:,2));

xgi = reshape(xg,li,lj);
ygi = reshape(yg,li,lj);

%calculando a ortonogonalidade da grade
[ort,xo,yo]=orthog(xgi,ygi);

%%% interpolacao linear com grade curvilinear %%%%%%%%%%%%%

disp('Interpolacao Linear.........')
Ti = griddata(lon,lat,T,xgi,ygi);
Si = griddata(lon,lat,S,xgi,ygi);
psigi = griddata(lon,lat,psig,xgi,ygi);
[ui,vi]=psi2uv(xgi,ygi,psigi);

%%% ANALISE OBJETIVA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp('Interpolacao Objetiva das variaveis termohalinas')
disp(' ')
xg=xg';yg=yg';
corrlen = 1; 
err = 0.02; 
[To,er] = scaloa(xg,yg,lon,lat,T,corrlen,err); 
[So,er] = scaloa(xg,yg,lon,lat,S,corrlen,err);
Do=sw_dens0(So,To); Do=Do-1000*ones(size(Do));

%%% APLICANDO CONDICOES DE CONTORNO PARA PSI ********************************
%%% DIRICHILET

%  load ../../levitus/etopo2_E_751x901.xyz;
%  x=etopo2_E_751x901(:,1);
%  y=etopo2_E_751x901(:,2);
%  z=etopo2_E_751x901(:,3);
%  clear etopo2_E_751x901
%  x = reshape(x,751,901);
%  y = reshape(y,751,901);
%  z = reshape(z,751,901);
%  fx = find(x(:,1) >= -43  & x(:,1) <= -30);
%  fy = find(y(1,:) >= -23 & y(1,:) <= -10);
%  xb = x(fx,fy);
%  yb = y(fx,fy);
%  zb = z(fx,fy);
%  save etopo2_leste2.mat xb yb zb

load ../../common/etopo2_leste.mat;

% buscando valor de lon e lat das isobatas de ?? m
figure(1)
if n <= 100
   [c,h] = contour(xb,yb,zb,[-100 -100]);
else 
   [c,h] = contour(xb,yb,zb,[-n -n]);
end
xcont = get(h(1),'xdata'); xcont = xcont(find(not(isnan(xcont))));
ycont = get(h(1),'ydata'); ycont = ycont(find(not(isnan(ycont))));
xcont = weim(31,'hann',xcont);
ycont = weim(31,'hann',ycont);
%  xcont=xcont';
%  ycont=ycont';

close(1); 

xf=[min(xcont) xcont];
yf=[max(ycont) ycont];
[xf,yf] = m_ll2xy(xf,yf,'patch');

% enxertando valores da interpolacao linear
xadd=[]; yadd=[]; psiadd=[];
%  for k = 40:45
%     xadd = [xadd xgi(k,:)];
%     yadd = [yadd ygi(k,:)];
%     psiadd = [psiadd psigi(k,:)];
%  end

% enxertando valores da isobata nos vetores
%  lon = [lon xcont xadd]; lat = [lat ycont yadd];
%  psig = [psig zeros(size(xcont)) psiadd];

% retirando eventuais NaNs

psig = psig(find(not(isnan(psig))));
lon = lon(find(not(isnan(psig))));
lat = lat(find(not(isnan(psig))));

disp('Interpolacao Objetiva de Psi')
[psigo,er] = scaloa(xg,yg,lon,lat,psig,corrlen,err);
er=100*sqrt(er); 

To=reshape(To,li,lj);
So=reshape(So,li,lj);
Do=reshape(Do,li,lj);
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

t1 = ['Temperatura ( \circ C) - OEI - ',nn,' m'];
t2 = ['Salinidade - OEI - ',nn,' m'];
t3 = ['Densidade Potencial (kg m^{-3}) - OEI - ',nn,' m'];
tit = ['t1t2t3'];
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
%       print(k,'-depsc',['../figuras/',P(c+k),'_OEI_',nn,'m']);
%       eval(['!epstopdf ../figuras/',P(c+k),'_OEI_',nn,'m.eps'])
%  c=c+1; 
%  end
%  clear k c

% ***********************************************************************************************
% fator de escala dos vetores
fc = 1; 

% construindo escalas de contorno para psi

lpsi = -30000:800:30000;
inc = ( max(max(psigo)) - min(min(psigo)) ) / 17;
lpsi2 = -30000:inc:30000;
int = num2str(100*round(inc/100));

%  figure(4)
%  set(gcf,...
%          'Color',[1 1 1],...
%          'InvertHardcopy','on',...
%          'PaperUnits','inches',...
%          'Units','inches',...
%          'PaperOrientation','portrait',...
%          'PaperPosition',[0 0 8.5 11],...
%          'PaperPositionMode','manual',...
%          'PaperType','usletter',...
%          'Position',[.2 .2 8 9],...
%          'ShareColors','off',...
%          'Clipping','on');
%  
%  m_contourf(xgi,ygi,psigo,lpsi);hold on;shading flat;
%  caxis([-30000 30000])
%  m_contour(xgi,ygi,psigo,lpsi2,'w')
%  q = m_quiver(xgi,ygi,uo*fc,vo*fc,0,'k');
%  %  set(q,'color',[.3 .3 .3])
%  fill(xf,yf,[0.8 0.8 0.8]);
%  %  plot(xf,yf,'color',[.8 .8 .8])
%  m_plot(lonctd,latctd,'w.','markersize',6)
%  m_usercoast('costa_leste.mat','patch',[0 0 0])
%  vmax = round(max(max(sqrt((uo*100).^2+(vo*100).^2))));
%  m_quiver(-40,-11,0.2*fc,0,0,'w')
%  text = ['20 cm s^{-1}'];
%  m_text(-40,-11.3,text,'color','w','fontsize',10,'fontweight','bold')
%  text = [nn,' m'];
%  m_text(-34.8,-20.1,text,'color','r','fontweight','bold')
%  text = ['Contornos: ',int,' m^2s^{-1}'];
%  m_text(-40.5,-11.8,text,'color','w','fontsize',8,'fontweight','bold')
%  m_grid('box','fancy','yaxislocation','left','xaxislocation','bottom','fontsize',8);
%  title('\Psi Geostrofico OEII [m^2 s^{-1}]','fontsize',10,'fontweight','bold');
%  for  l = 1:3:length(xgi) 
%      text = ['Secao: ',num2str(l)];
%      m_text(xgi(l,end),ygi(l,end),text,'color','k','fontweight','bold')
%  end
%  cc = colorbar;
%  pos = get(cc,'Position');
%  set(cc,'Position',[pos(1) pos(2) pos(3)/2.5 pos(4)])



%%% COMECA AQUI A SELECAO DA SECAO PARA CALCULO DE TRANSPORTE
%  disp(' ')
%  sec = input('Selecione a secao para efetuar o calculo do transporte: ');
%  disp(' ')

%%% construindo o eixo de distancia da costa:
xsec = xgi(sec,:); ysec = ygi(sec,:);
dist = [sw_dist(ysec,xsec,'km')];
dist = [0 cumsum(dist)];

% encontrando coordenada do ponto da linha de costa alinhado com a secao,
% para atribuir valor zero a ele:
load ../../common/costa_leste.mat
xlc = ncst(:,1);
ylc = ncst(:,2);

disp(' ')
disp('Achando o ponto mais proximo a linha de costa ...........')
disp(' ')

dd = [];
for t = 1:length(xsec)
    for u = 1:length(xlc)
        dd(u,t) = sw_dist([ysec(t) ylc(u)],[xsec(t) xlc(u)],'km'); 
    end
end    

[i,j] = find(dd == min(min(dd)));
dist = dist - dist(j)*ones(size(dist));

if sec > 1 
close(2) 
end

%%% criando topografia para a secao a partir do etopo2

latbt = mean(ysec);
blat = yb(1,:);
b = near(blat,latbt,1);
xbt = xb(:,b);
ybt = yb(:,b);
zbt = zb(:,b);

zbt = interp1(xbt,zbt,[xsec(1):0.01:xsec(end)],'cubic');
ybt = interp1(xbt,ybt,[xsec(1):0.01:xsec(end)],'linear');
xbt = [xsec(1):0.01:xsec(end)];

%%% construindo a secao a partir das distribuicoes horizontais
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('Construindo secao ... vai demorar....')
disp(' ')

vgao = [];
cont = 1;

load ../posicoes_leste1.dat; pb=posicoes_leste1;
latctd = pb(:,2); 
lonctd = pb(:,3);

clear posicoes_leste1

%%% carregando os dados hidrograficos %%%%%
inc = 10; % resolucao vertical
for n = 1:inc:1000

	nn=num2str(n);
	disp(' ')
	disp(['Calculando Velocidade para ',num2str(n),' m ........'])
	disp(' ')

	% carregando matrizes com os dados hidrograficos
	gpan=[];T=[];S=[];lon=[];lat=[];
	for k=1:12
		eval(['load ../mat/lesteI_agp_rad',num2str(k),'.mat']);
		eval(['load ../mat/lesteI_T_rad',num2str(k),'.mat']);
		eval(['load ../mat/lesteI_S_rad',num2str(k),'.mat']);
		
		gpan = [gpan agp(n,:)];
		lon = [lon -lons];
		lat = [lat -lats];  
		T = [T Ti(n,:)];
		S = [S Si(n,:)];
	end
	
	clear agp Ti Si lons lats
	
	%%% calculando PSI
	f0=sw_f(-15);
	psig=gpan/f0; psig=-psig;
	psig=psig-mean(psig);
	psigctd = psig;
	
	%%% carregando a grade construida no seagrid
	
	load ../../common/seagrid_leste2.dat;
	grid = seagrid_leste2;
	xg = grid(:,9); 
	yg = grid(:,10);
	
	lj = max(grid(:,1));
	li = max(grid(:,2));
	
	xgi = reshape(xg,li,lj);
	ygi = reshape(yg,li,lj);
	
	%calculando a ortonogonalidade da grade
	[ort,xo,yo]=orthog(xgi,ygi);
	
	%%% interpolacao linear com grade curvilinear %%%%%%%%%%%%%
	Ti = griddata(lon,lat,T,xgi,ygi);
	Si = griddata(lon,lat,S,xgi,ygi);
	psigi = griddata(lon,lat,psig,xgi,ygi);
	[ui,vi]=psi2uv(xgi,ygi,psigi);
	
	%%% ANALISE OBJETIVA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	xg=xg';yg=yg';
	corrlen = 1; 
	err = 0.02; 
	[To,er] = scaloa(xg,yg,lon,lat,T,corrlen,err); 
	[So,er] = scaloa(xg,yg,lon,lat,S,corrlen,err);
	Do=sw_dens0(So,To); Do=Do-1000*ones(size(Do));
	
	%%% APLICANDO CONDICOES DE CONTORNO PARA PSI ********************************
	%%% DIRICHILET
	
	load ../../common/etopo2_leste.mat;
	
	% buscando valor de lon e lat das isobatas de ?? m
        figure(5)
	if n <= 100
	    [c,h] = contour(xb,yb,zb,[-100 -100]);
	else 
	    [c,h] = contour(xb,yb,zb,[-n -n]);
	end
	xcont = get(h(1),'xdata'); xcont = xcont(find(not(isnan(xcont))));
	ycont = get(h(1),'ydata'); ycont = ycont(find(not(isnan(ycont))));
	xcont = weim(31,'hann',xcont);
	ycont = weim(31,'hann',ycont);
	%  xcont=xcont';
	%  ycont=ycont';
	
	close(5); 
	
	xf=[min(xcont) xcont];
	yf=[max(ycont) ycont];
	[xf,yf] = m_ll2xy(xf,yf,'patch');
	
	% enxertando valores da isobata nos vetores
%  	lon = [lon xcont xadd]; lat = [lat ycont yadd];
%  	psig = [psig zeros(size(xcont)) psiadd];
	
	% retirando eventuais NaNs
	
	psig = psig(find(not(isnan(psig))));
	lon = lon(find(not(isnan(psig))));
	lat = lat(find(not(isnan(psig))));

	[psigo,er] = scaloa(xg,yg,lon,lat,psig,corrlen,err);
	er=100*sqrt(er); 
	
	To=reshape(To,li,lj);
	So=reshape(So,li,lj);
	Do=reshape(Do,li,lj);
	psigo=reshape(psigo,li,lj);
	er=reshape(er,li,lj);
	
	% calculando componentes u e v
	[uo,vo]=psi2uv(xgi,ygi,psigo);

        vgao(cont,:) = vo(sec,:); 
        cont = cont + 1;

end

% construindo eixo vertical:
prof = 1:inc:n;

eval(['save ../mat/lesteI_vgao_sec',num2str(sec),'.mat  vgao dist prof ysec xsec xbt ybt zbt'])

figure(2)
set(gcf,'Color',[1 1 1]);
lv = [-0.5:0.01:0.5];
contourf(xsec,-prof,-vgao,lv); shading flat; hold on; caxis([lv(1) lv(end)]);
if max(max(abs(vgao))) <= 0.2
   [c,h]=contour(xsec,-prof,vgao,[-0.8:0.1:0.8],'k');
else 
   [c,h]=contour(xsec,-prof,vgao,[-0.8:0.2:0.8],'k');
end
clabel(c,h,'labelspacing',500);
set(gca,'plotboxaspectratio',[1 .7 1])
fill([min(xbt) xbt],[min(zbt) zbt],[.7 .7 .7])
latsec = -(round(ysec(1)*10))./10;
tit = ['Velocidade Geostrofica - OEI - ',num2str(latsec),'^\circS [m s^{-1}]'];
t = title(tit,'fontweight','bold');
pos = get(t,'Position');
set(t,'Position',[pos(1) pos(2)*3 pos(3)])
xlabel('[^\circ W] / [km]')
ylabel('Profundidade [m]')
cc = colorbar;
pos = get(cc,'Position');
set(cc,'Position',[pos(1) pos(2)*1.9 pos(3)/2.5 pos(4)/1.33])
pos=get(cc,'Position');
cc1=-str2num(get(cc,'YTickLabel'));cc1(find(cc1==0))=0;
set(cc,'YTickLabel',cc1)
ax1 = gca;
% eixo 2
ax2 = axes('Position',get(ax1,'Position'));
p=plot(dist,dist); hold on
set(gca,'color','none','XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k')
set(gca,'YTickLabel',[]);
set(p,'linestyle','none')
set(gca,'plotboxaspectratio',[1 .7 1])

   eval(['print -depsc ../figuras/vgAO_leste1_sec',num2str(sec),'.eps'])
   eval(['!epstopdf ../figuras/vgAO_leste1_sec',num2str(sec),'.eps'])


% criacao de video: eval(['!mencoder mf://figs/horiz_',varnam{nvar},'*_sig',num2str(sig),'.png -mf w=1200:h=901:fps=5:type=png -ovc raw -oac copy -o figs/horiz_',varnam{nvar},'_sig',num2str(sig),'.avi']);

clear Do       b        f0       k        lon      p        vgao     xgi      yg P        blat     fc       l        lonctd   pb       t        vi       xi       ygi S        c        gpan     lD       lonlim   pos      t1       vmax     xlc      ylc Si       cc       grid     lS       lpsi     prof     t2       vo       xo       yo So       cc1      h        lT       lpsi2    psiadd   t3       xadd     xsec     ysec T        cont     i        lat      lv       psig     text     xb       yadd     z Ti       corrlen  inc      latbt    n        psigctd  tit      xbi      yb       zb To       dd       int      latctd   ncst     psigi    u        xbt      ybt      zbi ax1      dist     isob     latlim   nn       psigo    ui       xcont    ycont    zbt ax2      er       j        li       ort      q        uo       xf       yf

end

toc

