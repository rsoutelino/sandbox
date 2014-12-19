%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   PROGRAMA PARA CALCULAR VELOCIDADES GEOSTROFICAS
%   RELATIVAS COM BASE EM MEDIDAS DIRETAS, OU SEJA
%      REFERENCIAMENTO DO CALCULO GEOSTROFICO
%
%        Rafael Soutelino - Mestrado IOUSP
%                   jul/2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% este programa faz o procedimento atraves do uso 
% da das vgeo relativas a superficie e das velocidades geostroficas
% oriundas da filtragem dos dados de ADCP, ambas oriundas da AO
clear all;close all;clc;warning off

% configurações
cc = 'n';

% declarando os limites da area, projecao, etc...
m_proj('mercator','long',[-41 -33.7],'lat',[-20.5 -10.2],'on');
load ../common/etopo2_leste.mat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      P A R T E   I      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A PRIMEIRA PARTE DO PROGRAMA SE DESTINA A 
% CARREGAR OS DADOS DE ADCP EM 150 m NA GRADE CURVILINEA
% DA OEII, INTERPOLADOS PELA AO VETORIAL

% carregando .mat com u,v,lon,lat dos dados de ADCP
load ../adcp/mat/psi_obs_150m.mat;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     P A R T E   II      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% A SEGUNDA PARTE DO PROGRAMA SE DESTINA A 
% CALCULAR E GRADEAR OS CAMPOS DE VGEO CALCULADAS 
% COM BASE NOS DADOS TERMOHALINOS E REFERENCIADAS NA 
% SUPERFICIE

% carregando os .mat com as velocidades geostroficas 
% previamente calculadas:

% ---------------------------------------------------------
% buscando as velocidades em 150 m, ou seja, abaixo 
% da camada de ekman, onde as velocidade do adcp estao livres
% da ageostrofia
%  load ../hidrografia/mat/uv_psi_AO_OEII_150m.mat;
%  ugo150 = ugo; vgo150 = vgo;

% ------------------------------------------------------
% buscando as velocidades no nivel desejado:

pp = input('Escolha a profundidade:  ');

        disp(['Calculando velocidades referenciadas por ADCP, em ',num2str(pp),' m...'])

	eval(['load ../hidrografia/mat/uv_psig_AO_OEII_',num2str(pp),'m.mat']);
	
	% referenciando em 150 m:
	Uctd = ugo;
	Vctd = vgo;
       
        if pp == 150;
           Uctd = 0*Uctd;
           Vctd = 0*Vctd;
        end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%     P A R T E   III     %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% REFERENCIANDO AS VELOCIDADES A PARTIR DAS VELOCIDADES DE ADCP E 
	% APLICANDO ANALISE OBJETIVA VETORIAL PARA OBTER A FUNCAO DE CORRENTE
	
	
	% referenciando a partir dos dados de ADCP
        
	U = Uadcp + Uctd;
        V = Vadcp + Vctd;
 
	
	% aplicando analise objetiva vetorial para calcular o campo final
	%%% INTERPOLACAO POR ANALISE OBJETIVA
	
	% preparando campos
	U = reshape(U,1,li*lj);
	V = reshape(V,1,li*lj);
	
	lc = 1; % em graus, entao, precisa passar as velocidades para grau/s
	E = 0.02;
	
	U = U*9e-6; % fator de conversao para grau/s (9e-6)
	V = V*9e-6;
	
	% analise objetiva vetorial para obter psi ABSOLUTO
	disp('ANALISE OBJETIVA VETORIAL:  u,v --> psi')
	[psi] = vectoa(xg,yg,xg,yg,U,V,lc,E,0);
	
	% passando psi para m^2/s
	psi = psi.*(111120^2); 
	
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
	
	xg2 = xg; yg2 = yg;
	
	%% aplicando condicao de contorno em psi
        if cc == 's'
	   xg2 = [xg; xcont]; yg2 = [yg; ycont];
	   psi = [psi; zeros(size(xcont))];
	
	% analise objetiva escalar em psiob para aplicar cond de contorno
	disp('ANALISE OBJETIVA ESCALAR ... para acertar cond de contorno')
	[psi] = scaloa(xg',yg',xg2',yg2',psi',lc,E);
	end

	% preparando para plotar
        psi = reshape(psi,li,lj);
	[U,V] = psi2uv(xgi,ygi,psi);
	lpsi = -60000:500:40000;
	inc = ( max(max(psi)) - min(min(psi)) ) / 40;
	lpsi2 = -70000:inc:70000;
	int = num2str(100*round(inc/100));
	fc = 1;


%  %%% FIGURA DE APOIO TEMPORARIA
%  
%  figure(10)
%  	set(gcf,...
%  		'Color',[1 1 1],...
%  		'InvertHardcopy','on',...
%  		'PaperUnits','inches',...
%  		'Units','inches',...
%  		'PaperOrientation','portrait',...
%  		'PaperPosition',[0 0 8.5 11],...
%  		'PaperPositionMode','manual',...
%  		'PaperType','usletter',...
%  		'Position',[.2 .2 8 9],...
%  		'ShareColors','off',...
%  		'Clipping','on');
%  	
%  	[c,h] = m_contourf(xgi,ygi,sqrt(U.^2 + V.^2),0:0.05:0.8);shading flat;hold on
%          caxis([0 0.6])
%  %  	m_contour(xgi,ygi,psi,lpsi2,'w')
%  	m_quiver(xgi,ygi,U*fc,V*fc,0,'w')
%          clabel(c,h)
%          if pp < 100
%             pl = 100;
%          else
%             pl = pp;
%          end
%  	[c,h]=m_contourf(xb,yb,zb,[-pl -pl],'k');
%          set(h,'facecolor',[.7 .7 .7]);
%          [c,h]=m_contour(xb,yb,zb,[-pl -pl],'k');
%          set(h,'color',[.7 .7 .7]);
%  	m_usercoast('../common/costa_leste.mat','patch',[0 0 0])
%  	m_quiver(-40,-11,0.5*fc,0,0,'w')
%  	text = ['50 cm s^{-1}'];
%  	m_text(-40,-11.3,text,'color','w','fontsize',10,'fontweight','bold')
%  	text = [num2str(pp),' m'];
%  	m_text(-34.8,-20.1,text,'color','r','fontweight','bold')
%  	title('\Psi Geostrofico Referenciado por ADCP - OEII [m^2 s^{-1}]','fontsize',10,'fontweight','bold');
%  	text = ['Contornos: ',int,' m^2s^{-1}'];
%  	m_text(-40.5,-11.8,text,'color','w','fontsize',8,'fontweight','bold')
%  	m_grid('box','fancy','yaxislocation','left','xaxislocation','bottom','fontsize',10);
%  	cc = colorbar;
%  	pos = get(cc,'Position');
%  	set(cc,'Position',[pos(1) pos(2)*1.3 pos(3)/2.5 pos(4)/1.1])
%          drawnow
%  

	figure(6)
        hold on
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
        caxis([-60000 40000])
%  	m_contour(xgi,ygi,psi,lpsi2,'w')
	m_quiver(xgi,ygi,U*fc,V*fc,0,'k')
        if pp < 100
           pl = 100;
        else
           pl = pp;
        end
	[c,h]=m_contourf(xb,yb,zb,[-pl -pl],'k');
        set(h,'facecolor',[.7 .7 .7]);
        [c,h]=m_contour(xb,yb,zb,[-pl -pl],'k');
        set(h,'color',[.7 .7 .7]);
	m_usercoast('../common/costa_leste.mat','patch',[0 0 0])
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
        drawnow
        pause(2)
        eval(['save mat/uv_psi_ref_OEII_',num2str(pp),'m.mat psi U V'])
        
        eval(['print -depsc figuras/psi_ref_adcp_OEII_',num2str(pp),'m.eps'])
        eval(['!epstopdf figuras/psi_ref_adcp_OEII_',num2str(pp),'m.eps'])

      
%  eval(['save mat/psi_ref_adcp_',num2str(pp),'m.mat psi'])















