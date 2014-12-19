%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    PROGRAMA PARA PLOTAR O CAMPO MEDIO DE 
%      PSI NA CAMADA SUPERIOR DO OCEANO
%               Oceano Leste 2
%         JUL/2007 - Mestrado - IOUSP
%           Rafael Guarino Soutelino
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
clear all;close all;clc;warning off
tic
% definindo os limites da area
lonlim = [-41 -33.7]; latlim = [-20.5 -10.2];
m_proj('mercator','long',[min(lonlim) max(lonlim)],'lat',[min(latlim) max(latlim)],'on');
%  m_gshhs_f('save','costa_leste.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Construindo campo medio de psi para a camada superficial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PSI = zeros(50,30);
cont = 1;

load posicoes_leste2.dat; pb=posicoes_leste2;
latg=pb(:,2); latm=pb(:,3); latctd=latg+latm/60; latctd=-latctd; latctd=latctd'; 
long=pb(:,4); lonm=pb(:,5); lonctd=long+lonm/60; lonctd=-lonctd; lonctd=lonctd'; 

clear latg latm long lonm posicoes_leste2

%%% carregando os dados hidrograficos %%%%%
inc = 1; % resolucao vertical
pp = 20;

for n = 10:30

	nn=num2str(n);
	disp(' ')
	disp(['Calculando Velocidade para ',num2str(n),' m ........'])
	disp(' ')

	% carregando matrizes com os dados hidrograficos
	gpan=[];T=[];S=[];lon=[];lat=[];
	for k=1:12
		eval(['load ../mat/lesteII_agp_rad',num2str(k),'.mat']);
		eval(['load ../mat/lesteII_T_rad',num2str(k),'.mat']);
		eval(['load ../mat/lesteII_S_rad',num2str(k),'.mat']);
		
		gpan = [gpan agp(n,:)];
		lon = [lon lons];
		lat = [lat lats];  
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
%  	xcont = weim(31,'hann',xcont);
%  	ycont = weim(31,'hann',ycont);
	xcont=xcont';
	ycont=ycont';
	
	close(5); 
	
	xf=[min(xcont) xcont];
	yf=[max(ycont) ycont];
	[xf,yf] = m_ll2xy(xf,yf,'patch');
	
	% enxertando valores da isobata nos vetores
%  	lon = [lon xcont]; lat = [lat ycont];
%  	psig = [psig zeros(size(xcont))];
%  	
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
	
        PSI = PSI + psigo;
        
        cont = cont + 1;

end

PSI = PSI./pp;

[U,V]=psi2uv(xgi,ygi,PSI);


% fator de escala dos vetores
fc = 1; 

% construindo escalas de contorno para psi

lpsi = -40000:800:30000;
inc = ( max(max(PSI)) - min(min(PSI)) ) / 17;
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

m_contourf(xgi,ygi,PSI,lpsi);hold on;shading flat;
caxis([-30000 30000])
m_contour(xgi,ygi,PSI,lpsi2,'w')
q = m_quiver(xgi,ygi,U*fc,V*fc,0,'k');
%  set(q,'color',[.3 .3 .3])
fill(xf,yf,[0.8 0.8 0.8]);
%  plot(xf,yf,'color',[.8 .8 .8])
m_plot(lonctd,latctd,'w.','markersize',6)
m_usercoast('costa_leste.mat','patch',[0 0 0])
vmax = round(max(max(sqrt((U*100).^2+(V*100).^2))));
m_quiver(-40,-11,0.2*fc,0,0,'w')
text = ['20 cm s^{-1}'];
m_text(-40,-11.3,text,'color','w','fontsize',10,'fontweight','bold')
text = ['1-20 m'];
m_text(-34.8,-20.1,text,'color','r','fontweight','bold')
text = ['Contornos: ',int,' m^2s^{-1}'];
m_text(-40.5,-11.8,text,'color','w','fontsize',8,'fontweight','bold')
m_grid('box','fancy','yaxislocation','left','xaxislocation','bottom','fontsize',8);
title('\Psi Geostrofico OEII [m^2 s^{-1}]','fontsize',10,'fontweight','bold');
cc = colorbar;
pos = get(cc,'Position');
set(cc,'Position',[pos(1) pos(2)*1.3 pos(3)/2.5 pos(4)/1.1])
   print(4,'-depsc',['../figuras/psi_OEII_med_',num2str(pp),'m']);
   eval(['!epstopdf ../figuras/psi_OEII_med_',num2str(pp),'m.eps'])
drawnow

toc

