%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Calcula Campos de Funcao de Corrente Geostrofica 
%                 Referenciada por ADCP
%            Rafael Soutelino - janeiro/20010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc;warning off

%%% CONFIGURACOES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p0 = 1000;         % nivel de referencia
pmin = 15;        % profundidade minima de amostragem
f0 = sw_f(-16);    % parametro de coriolis

fig1 = 'n';       % plota distr. horizontais termohalinas
fig2 = 's';       % plota psi geostrofico relativo (escolher p0 adequado!!!)
fig3 = 'n';       % plota psi observado puro
fig4 = 'n';       % plota figuras finais com psi referenciado

mask = 100;       % mascara de plataforma continental

Cc = 's';                 % condicoes de contorno
alis = 'n';  jan = 11;    % janela nos dados de adcp
dec = 'n'; dectx = 10;    % subamostragem dos dados de adcp
cruz = ''             % nome do diretorio CODAS
CRUZ = 'OEII'              % para constar em titulos de figuras
n = 50; % input('Escolha a profundidade para plotar psi:   ');

% DEFININDO OS LIMITES DA AREA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lonlim = [-41 -33.7]; latlim = [-22.5 -12];
m_proj('mercator','long',[lonlim(1) lonlim(2)],'lat',[latlim(1) latlim(2)],'on');

%%% CARREGANDO OS DADOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eval(['load ../posicoes_leste2.dat;']) 
eval(['pos = posicoes_leste2;']);   
load ../../common/etopo2_leste.mat;
load ../../common/seagrid_leste2.dat;
%  eval(['load ',cruz,'/contour/',cruz,'_uv.mat;']);
%  eval(['load ',cruz,'/contour/',cruz,'_xy.mat;']);

% preparando variaveis 
nest = pos(:,1);                       
lonctd2 = [pos(:,4) + pos(:,5)/60]*-1;   
latctd2 = [pos(:,2) + pos(:,3)/60]*-1;   

xx = [1:112];  

% CARREGANDO OS DADOS HIDROGRAFICOS JA TRATADOS

c = 1;
for k=1:length(xx);
    XX = num2str(xx(k));
    eval(['load  ~/mestrado/dados/leste2/ctd/filtrados/lesteII_ctd',XX,'.dat']);
    dados = eval(['lesteII_ctd',XX]);
    eval(['clear lesteII_ctd',XX]);

    p = dados(:,1);          
    t = dados(:,2);          
    s = dados(:,3);          

    if n < p0
       pmax = p0;
    else
       pmax = n;
    end

    f = find(p >= pmin & p <= pmax);

  if max(p) >= pmax
     T(1:length(f),c) = t(f);
     S(1:length(f),c) = s(f);
     latctd(c) = latctd2(find(nest==xx(k)));
     lonctd(c) = lonctd2(find(nest==xx(k)));
     c = c+1;
     clear s t p XX dados 
  else
     disp(['Estacao nao atinge ',num2str(pmax),' m'])
  end
end

clear c k f pos posicoesCTD_NEII xx
clc

p = pmin:pmax; p = p';        % REMONTANDO VETOR DE p

%%% CALCULANDO PSI PONTO A PONTO, NIVEL A NIVEL %%%%%%%%%%%%%%%%%%%%%%%

agp = sw_gpan(S,T,p);
agp = agp - ones(size(p)) * agp(p0-pmin,:);
psig = agp./f0; psig = -psig;

% separando o nivel para conduzir os calculos
j = n - pmin;
T = T(j,:);
S = S(j,:);
psig = psig(j,:);

% removendo a media de psi para calcular gradientes
psig = psig - mean(psig);

%%% INTERPOLACAO DE PSI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grid = seagrid_leste2;
xg = grid(:,9);              
yg = grid(:,10);             
lj = max(grid(:,1));         
li = max(grid(:,2));         
xgi = reshape(xg,li,lj);     
ygi = reshape(yg,li,lj);     
[ort,xo,yo]=orthog(xgi,ygi);
clear seagrid grid

%%% INTERPOLACAO LINEAR %%%%%%%%%%%%%%
disp('Interpolacao Linear.........')

Ti = griddata(lonctd,latctd,T,xgi,ygi);
Si = griddata(lonctd,latctd,S,xgi,ygi);
psigi = griddata(lonctd,latctd,psig,xgi,ygi);
[ui,vi] = psi2uv(xgi,ygi,psigi);

%%% INTERPOLACAO POR ANALISE OBJETIVA %%%%%%%%%
disp(' ')
disp('Analise Objetiva .......')
disp(' ')

xg=xg';yg=yg';
lc = 1;       
E = 0.02;           
[To,er] = scaloa(xg,yg,lonctd,latctd,T,lc,E);
[So,er] = scaloa(xg,yg,lonctd,latctd,S,lc,E);
Do=sw_dens0(So,To); Do=Do-1000*ones(size(Do));
[psigo,er] = scaloa(xg,yg,lonctd,latctd,psig,lc,E);
er=100*sqrt(er);

To=reshape(To,li,lj);
So=reshape(So,li,lj);
Do=reshape(Do,li,lj);
psigo=reshape(psigo,li,lj);
er=reshape(er,li,lj);
psigo = psigo - mean(mean(psigo)); % removendo a media

% condicao de contorno de dirichilet, caso tenha sido prescrito no inicio do programa
if Cc == 's'
figure(50)
if n <= 100
   [c,h] = contour(xb,yb,zb,[-100 -100]);
else 
   [c,h] = contour(xb,yb,zb,[-n -n]);
end
xcont = get(h(1),'xdata'); xcont = xcont(find(not(isnan(xcont))));
ycont = get(h(1),'ydata'); ycont = ycont(find(not(isnan(ycont))));
close(50)
xcont=xcont';
ycont=ycont';
xg2 = [xg xcont]; yg2 = [yg ycont];
psigo = reshape(psigo,1,li*lj);
psigo = [psigo zeros(size(xcont))];
psigo = psigo(find(not(isnan(psigo))));
xg2 = xg2(find(not(isnan(psigo))));
yg2 = yg2(find(not(isnan(psigo))));

disp('Aplicando condicao de contorno de dirichlet...........')
psigo = scaloa(xg,yg,xg2,yg2,psigo,lc,E);
psigo = reshape(psigo,li,lj);

end

[uctd,vctd] = psi2uv(xgi,ygi,psigo);

%% colocando mascara para os erros grandes
%  f = find(er > 40);
%  psigo(f) = nan; uctd(f) = nan; vctd(f) = nan;
s=1;

%% eliminando parte superior da grade
f2 = [42:50];
uctd(f2,:) = nan; vctd(f2,:) = nan;
psigo(f2,:) = nan;

figure(1)
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
hold on
%  m_contourf(xgi,ygi,psiob,lpsi);hold on; shading flat
%  caxis([-30000 30000])
%  m_contour(xgi,ygi,psiob,lpsi2,'k'); 
m_quiver(xgi(1:1:end),ygi(1:1:end),uctd(1:1:end)*s,vctd(1:1:end)*s,0,'k');
%  m_plot(xgi(f),ygi(f),'.w');
m_plot(xgi(f2,:),ygi(f2,:),'.w');
[c,h]=m_contourf(xb,yb,zb,[-mask -mask],'k');
set(h,'facecolor',[.7 .7 .7]);
[c,h]=m_contour(xb,yb,zb,[-mask -mask],'k');
set(h,'color',[.7 .7 .7]);
m_usercoast('../../common/costa_leste.mat','patch',[0 0 0])
m_grid('box','fancy','yaxislocation','left','xaxislocation','bottom');
m_quiver(-40,-11,0.5*s,0,0,'w')
text = ['50 cm s^{-1}'];
m_text(-40,-11.3,text,'color','w','fontsize',10,'fontweight','bold')
text = [num2str(n),' m'];
m_text(-34.8,-20.1,text,'color','k','fontweight','bold')
title('Geostrophic \Psi - OEII [m^2 s^{-1}]','fontsize',10,'fontweight','bold');
%  cc = colorbar;
%  pos = get(cc,'Position');
%  set(cc,'Position',[pos(1) pos(2)*2.2 pos(3)/2.5 pos(4)/1.5])
%     print(1,'-depsc',['../figuras/psiob_OEII_',num2str(p),'m']);
%     eval(['!epstopdf ../figuras/psiob_OEII_',num2str(p),'m.eps'])
drawnow

%% colored figure
figure(2)
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
hold on
m_contourf(xgi,ygi,psigo,50);hold on; shading flat
caxis([-30000 30000])
%  m_contour(xgi,ygi,psiob,lpsi2,'k'); 
m_quiver(xgi(1:1:end),ygi(1:1:end),uctd(1:1:end)*s,vctd(1:1:end)*s,0,'k');
%  m_plot(xgi(f),ygi(f),'.w');
m_plot(xgi(f2,:),ygi(f2,:),'.w');
[c,h]=m_contourf(xb,yb,zb,[-mask -mask],'k');
set(h,'facecolor',[.7 .7 .7]);
[c,h]=m_contour(xb,yb,zb,[-mask -mask],'k');
set(h,'color',[.7 .7 .7]);
m_usercoast('../../common/costa_leste.mat','patch',[0 0 0])
m_grid('box','fancy','yaxislocation','left','xaxislocation','bottom');
m_quiver(-40,-11,0.5*s,0,0,'w')
text = ['50 cm s^{-1}'];
m_text(-40,-11.3,text,'color','w','fontsize',10,'fontweight','bold')
text = [num2str(n),' m'];
m_text(-34.8,-20.1,text,'color','k','fontweight','bold')
title('Geostrophic \Psi - OEII [m^2 s^{-1}]','fontsize',10,'fontweight','bold');
cc = colorbar;
pos = get(cc,'Position');
set(cc,'Position',[pos(1) pos(2)*2.2 pos(3)/2.5 pos(4)/1.5])
%     print(1,'-depsc',['../figuras/psiob_OEII_',num2str(p),'m']);
%     eval(['!epstopdf ../figuras/psiob_OEII_',num2str(p),'m.eps'])
drawnow
stop

%%% DISTRIBUICAO HORIZONTAL DOS DADOS TERMOHALINOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if fig1 == 's'
zb = smoo2(zb,-9999,5,1);
figure
m_contourf(xgi,ygi,Ti,20); hold on; caxis([min(min(Ti)) max(max(Ti))]) ;colorbar
[c,h] = m_contourf(xb,yb,zb,[-mask -mask]);
set(h,'facecolor',[.5 .5 .5])
m_usercoast('costa.mat','patch',[0 0 0])
title(['T Int. Linear - ',num2str(n),' m - ',CRUZ])
m_grid

%  eval(['print -depsc figuras/T_linear_',num2str(n),'m_',cruz,'.eps'])
%  eval(['!epstopd T_linear_',num2str(n),'m_',cruz,'.eps'])

figure
m_contourf(xgi,ygi,To,20); hold on; caxis([min(min(To)) max(max(To))]) ;colorbar
[c,h] = m_contourf(xb,yb,zb,[-mask -mask]);
set(h,'facecolor',[.5 .5 .5])
m_usercoast('costa.mat','patch',[0 0 0])
title(['T Int. AO - ',num2str(n),' m - ',CRUZ])
m_grid

%  eval(['print -depsc figuras/T_AO_',num2str(n),'m_',cruz,'.eps'])
%  eval(['!epstopdf figuras/T_AO_',num2str(n),'m_',cruz,'.eps'])

figure
m_contourf(xgi,ygi,Si,20); hold on; caxis([min(min(Si)) max(max(Si))]) ;colorbar
[c,h] = m_contourf(xb,yb,zb,[-mask -mask]);
set(h,'facecolor',[.5 .5 .5])
m_usercoast('costa.mat','patch',[0 0 0])
title(['S Int. Linear - ',num2str(n),' m - ',CRUZ])
m_grid

%  eval(['print -depsc figuras/S_Linear_prof',num2str(n),'m_',Cruz,'.eps'])
%  eval(['!epstopdf figuras/S_Linear_prof',num2str(n),'m_',Cruz,'.eps'])

figure
m_contourf(xgi,ygi,So,20); hold on; caxis([min(min(So)) max(max(So))]) ;colorbar
[c,h] = m_contourf(xb,yb,zb,[-mask -mask]);
set(h,'facecolor',[.5 .5 .5])
m_usercoast('costa.mat','patch',[0 0 0])
title(['S Int. AO - ',num2str(n),' m - ',CRUZ])
m_grid

%  eval(['print -depsc figuras/S_AO_prof',num2str(n),'m_',Cruz,'.eps'])
%  eval(['!epstopdf figuras/S_AO_prof',num2str(n),'m_',Cruz,'.eps'])

end

%%% CAMPOS DE PSI SEM REFERENCIAMENTO POR ADCP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if fig2 == 's'

figure
dec = 1;
s = 1;
m_contourf(xgi,ygi,psigi,100); shading flat; hold on;
m_quiver(xgi(1:dec:end,1:dec:end),ygi(1:dec:end,1:dec:end),ui(1:dec:end,1:dec:end)*s,vi(1:dec:end,1:dec:end)*s,0,'k');
caxis([min(min(psigi)) max(max(psigi))]) ;cc = colorbar;
[c,h] = m_contourf(xb,yb,zb,[-mask -mask]);
set(h,'facecolor',[.8 .8 .8])
[c,h] = m_contour(xb,yb,zb,[-1000 -1000],'r'); % ISOBATA DE 100 m
m_usercoast('costa.mat','patch',[0.5 0.5 0.5])
t = [num2str(n),' m'];
m_text(max(lonlim) -1.8,min(latlim)+0.5,t,'fontsize',12,'fontweight','bold')
m_quiver(max(lonlim) -1.8,max(latlim)-0.5,1*s,0,0,'k');
t = ['1 m . s^{-1}'];
m_text(max(lonlim) -1.8,max(latlim)-0.9,t,'fontsize',10,'fontweight','bold')
title(['\Psi geo. [m^{2} . s^{-1}] - Int. Linear - NR ',num2str(p0),' m - prof ',num2str(n),' m - ',CRUZ])
pos = get(cc,'Position');
set(cc,'Position',[pos(1)/3.1 pos(2)*1.5 pos(3)/3 pos(4)/2])
m_grid('box','fancy')
%  m_plot(lonctd,latctd,'w.','markersize',6)
m_plot(lonctd,latctd,'Color','w','LineStyle','none','LineWidth',...
       1,'marker','o','markersize',4,'MarkerFaceColor',...
   'w','MarkerEdgeColor','k')
%  eval(['print -depsc figuras/PsiGeo_Linear_NR',num2str(p0),'m_prof',num2str(n),'m_',Cruz,'.eps'])
%  eval(['!epstopdf figuras/PsiGeo_Linear_NR',num2str(p0),'m_prof',num2str(n),'m_',Cruz,'.eps'])

%%%
% mascarando o erro

%  eval(['print -depsc figuras/PsiGeo_AO_NR',num2str(p0),'m_prof',num2str(n),'m_',Cruz,'.eps'])
%  eval(['!epstopdf figuras/PsiGeo_AO_NR',num2str(p0),'m_prof',num2str(n),'m_',Cruz,'.eps'])

end
stop
%%% REFERENCIAMENTO POR ADCP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% criando variaveis ADCP
u = uv(:,1:2:end-1);
v = uv(:,2:2:end);
lonadcp = xyt(1,:)-360; latadcp = xyt(2,:);
t = xyt(3,:);

% separando os dados no p0 estabelecido
n2 = near(zc,p0,1);
u = u(n2,:); v = v(n2,:);

% removendo eventuais dados espurios de ADCP
mod = sqrt(u.^2+v.^2);
f = find(mod > 2);
u(f)=nan; v(f)=nan; mod(f)=nan;
f = find(isnan(u)==0);
u=u(f);v=v(f);lonadcp=lonadcp(f);latadcp=latadcp(f);

% subamostrando os dados de ADCP
if dec == 's'
   u = u(1:dectx:end); v = v(1:dectx:end); 
   lonadcp = lonadcp(1:dectx:end); latadcp = latadcp(1:dectx:end); 
end

% alisamento dos dados de ADCP
if alis == 's'
   u = weim(jan,'hann',u);
   v = weim(jan,'hann',v);
end

%%% ANALISE OBJETIVA VETORIAL PARA CALCULAR PSI OBSERVADO

% convertendo unidade de velocidade m/s --> grau/s
u = u * 9e-6;
v = v * 9e-6;

disp(' ')
disp('Analise Objetiva Vetorial.......')
disp(' ')
psiobs = vectoa(xg,yg,lonadcp,latadcp,u,v,lc,E,0);
psiobs = reshape(psiobs,li,lj);

% re-convertendo para o SI
psiobs = psiobs * 111120.^2;
psiobs = psiobs - mean(mean(psiobs)); % removendo a media
[uadcp,vadcp] = psi2uv(xgi,ygi,psiobs); % u e v... 

%%% PLOTANDO FIGURA COM PSIOBS NO NIVEL p0 ESCOLHIDO

if fig3 == 's'

figure
dec = 2;
s = 1;
m_contourf(xgi,ygi,psiobs,100); shading flat; hold on;
m_quiver(xgi(1:dec:end,1:dec:end),ygi(1:dec:end,1:dec:end),uadcp(1:dec:end,1:dec:end)*s,vadcp(1:dec:end,1:dec:end)*s,0,'k');
caxis([min(min(psiobs)) max(max(psiobs))]) ;cc = colorbar;
[c,h] = m_contourf(xb,yb,zb,[-mask -mask]);
set(h,'facecolor',[.8 .8 .8])
[c,h] = m_contour(xb,yb,zb,[-1000 -1000],'r'); % ISOBATA DE 100 m
m_usercoast('costa.mat','patch',[.7 .7 .7])
t = [num2str(zc(n2)),' m'];
m_text(max(lonlim) -1.8,min(latlim)+0.5,t,'fontsize',12,'fontweight','bold')
m_quiver(max(lonlim) -1.8,max(latlim)-0.5,1*s,0,0,'k');
t = ['1 m . s^{-1}'];
m_text(max(lonlim) -1.8,max(latlim)-0.9,t,'fontsize',10,'fontweight','bold')
title(['\Psi Obs [m^{2} . s^{-1}] - Int. AOV - prof ',num2str(zc(n2)),' m - ',CRUZ])
pos = get(cc,'Position');
set(cc,'Position',[pos(1)/3.1 pos(2)*1.5 pos(3)/3 pos(4)/2])
m_grid('box','fancy')
m_plot(lonctd,latctd,'Color','w','LineStyle','none','LineWidth',...
       1,'marker','o','markersize',4,'MarkerFaceColor',...
   'w','MarkerEdgeColor','k');
%  m_text(-38.5124,-12.9704,'\bf Salvador','fontsize',8, 'HorizontalAlignment', 'right')
%  eval(['print -depsc figuras/PsiObs_prof',num2str(n),'m_',Cruz,'.eps'])
%  eval(['!epstopdf figuras/PsiObs_prof',num2str(n),'m_',Cruz,'.eps'])

end

%%% REFERENCIAMENTO PROPRIAMENTE DITO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psitot = psigo + psiobs;

%%% calculando finalmente as velocidades absolutas
[utot,vtot] = psi2uv(xgi,ygi,psitot);

%%% DISTRIBUICAO HORIZONTAL DA FUNCAO DE CORRENTE GEOSTROFICA ABSOULTA %%%%%%%%%%%%%%%%%%%%%%

if  fig4 == 's'

figure
dec = 2;
s = 1;
m_contourf(xgi,ygi,psitot,[-1e5:1000:1e5]); shading flat; hold on;
m_quiver(xgi(1:dec:end,1:dec:end),ygi(1:dec:end,1:dec:end),utot(1:dec:end,1:dec:end)*s,vtot(1:dec:end,1:dec:end)*s,0,'k');
caxis([min(min(psitot)) max(max(psitot))]) ;cc = colorbar;
[c,h] = m_contourf(xb,yb,zb,[-mask -mask],'k');
set(h,'facecolor',[.7 .7 .7]);
[c,h] = m_contour(xb,yb,zb,[-mask -mask],'k');
set(h,'color',[.7 .7 .7]);
m_usercoast('costa.mat','patch',[.5 .5 .5])
t = [num2str(n),' m'];
m_text(max(lonlim) -1.8,min(latlim)+0.5,t,'fontsize',12,'fontweight','bold')
m_quiver(max(lonlim) -1.8,max(latlim)-0.5,1*s,0,0,'k');
t = ['1 m . s^{-1}'];
m_text(max(lonlim) -1.8,max(latlim)-0.9,t,'fontsize',10,'fontweight','bold')
%  title(['\Psi Geostrofico Absoluto [m^{2} . s^{-1}] - ',num2str(n),' m - ',CRUZ])
pos = get(cc,'Position');
set(cc,'Position',[pos(1)/3.1 pos(2)*1.5 pos(3)/3 pos(4)/2])
m_grid
%  m_plot(lonctd,latctd,'Color','w','LineStyle','none','LineWidth',...
%         1,'marker','o','markersize',4,'MarkerFaceColor',...
%     'w','MarkerEdgeColor','k');
eval(['print -depsc figuras/psiref_',num2str(n),'m_',cruz,'.eps'])
eval(['!epstopdf figuras/psiref_',num2str(n),'m_',cruz,'.eps'])

end