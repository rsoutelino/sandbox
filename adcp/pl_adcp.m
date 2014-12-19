%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Lê, Processa e plota os dados de ADCP da     %
%          Comissão Oceano Leste II               %
%       Rafael Soutelino - Mestrado IOUSP         % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc; warning off

lonlim = [-41 -33.7]; latlim = [-20.5 -10.2];
m_proj('mercator','long',[min(lonlim) max(lonlim)],'lat',[min(latlim) max(latlim)],'on');
load ../../common/etopo2_leste.mat;

% leitura dos dados
U=[]; V=[]; lon=[]; lat=[]; mm=[]; dd=[]; hh=[]; mn=[]; pg=[];
for k = [1:54]
    eval(['load ../../../dados/leste2/adcp/exportados/adcp',num2str(k),'.mat'])
    U = [U SerEmmpersec'];
    V = [V SerNmmpersec'];
    lon = [lon  AnFLonDeg'];
    lat = [lat  AnFLatDeg'];
    mm = [mm SerMon'];
    dd = [dd SerDay'];	
    hh = [hh SerHour'];
    mn = [mn SerMin'];
    eval(['load ../../../dados/leste2/adcp/exportados/adcp_OEII_',num2str(k),'.mat'])
    pg = [pg SerPG4'];
end
clear Ser* RDI* An* k lonlim latlim

% separando os trajetos apenas nas radiais, de forma interativa
lon = [lon(413:2165) lon(2330:3750)];
lat = [lat(413:2165) lat(2330:3750)];
mm = [mm(413:2165) mm(2330:3750)];
dd = [dd(413:2165) dd(2330:3750)];
hh = [hh(413:2165) hh(2330:3750)];
mn = [mn(413:2165) mn(2330:3750)];
U = [U(:,413:2165) U(:,2330:3750)];
V = [V(:,413:2165) V(:,2330:3750)];
pg = [pg(:,413:2165) pg(:,2330:3750)];

% colocando NaN na flag de spike
f = find(U == -32768); 
U(f) = nan; V(f) = nan; 

% convertendo para o SI
U = U/1000; V = V/1000;

% escolhendo profundidade (bin) (cada bin corresponde a 8 m, com 8 m de
% blanking distance, ou seja, o primeiro bin corresponde a 8 m de prof.)

pp = input('Escolha a prof de interesse (no maximo 200m): ')
bin = ceil(pp/8)

u = U(bin,:); v = V(bin,:); pg = pg(bin,:);

% removendo nan
f = find(isnan(u)==0);
u = u(f); v = v(f); pg = pg(f); lon = lon(f); lat = lat(f); 
mm = mm(f); dd = dd(f); hh = hh(f); mn = mn(f); 

% eliminando medidas concomitantes a realização de estações de CTD
dist = sw_dist(lat,lon,'km');
f = find(dist >= 1);
u = u(f); v = v(f); pg = pg(f); lon = lon(f); lat = lat(f); 
mm = mm(f); dd = dd(f); hh = hh(f); mn = mn(f); 

% eliminando valores de pgood menores que 80 %
f = find(pg >= 80);
u = u(f); v = v(f); pg = pg(f); lon = lon(f); lat = lat(f); 
mm = mm(f); dd = dd(f); hh = hh(f); mn = mn(f);

% removendo spike: criterio do desvio padrao
f = find(abs(u) < (mean(u) + 8*std(u)));
u = u(f); v = v(f); pg = pg(f); lon = lon(f); lat = lat(f);
mm = mm(f); dd = dd(f); hh = hh(f); mn = mn(f);

f = find(abs(v) < (mean(v) + 8*std(v)));
Uadcp = u(f)'; Vadcp = v(f)'; pg = pg(f); lonadcp = lon(f)'; latadcp = lat(f)';
mm = mm(f)'; dd = dd(f)'; hh = hh(f)'; mn = mn(f)';

DT = datenum(2005,mm,dd,hh,mn,0); 

% vizualizando os vetores após tratamento básico:

fc = .5;
figure(2)
set(gcf,'color','w')
m_contourf(xb,yb,zb,[-6000:800:0]);shading flat;hold on
load /usr/local/matlab/toolbox/lado/m_map/ocean_colormap
colormap(ocean_colormap)
[c,h] = m_contour(xb,yb,zb,[-200 -1000],'k'); 
set(h,'color',[.3 .3 .3])
c=clabel(c,h,'labelspacing',500); set(c,'color',[.5 .5 .5])
m_usercoast('costa_leste.mat','patch',[0 0 0])
m_quiver(lon,lat,u*fc,v*fc,0,'y')
m_quiver(-40,-11,0.5*fc,0,0,'w')
text = ['50 cm s^{-1}'];
m_text(-40,-11.3,text,'color','w','fontsize',10,'fontweight','bold')
text = [num2str(pp),' m'];
m_text(-34.8,-20.1,text,'color','k','fontweight','bold')
title('ADCP - bruto - OEII','fontsize',10,'fontweight','bold');
m_grid('box','fancy','yaxislocation','left','xaxislocation','bottom','fontsize',8);
%     print(1,'-depsc',['../figuras/adcp_bruto_',num2str(pp),'m']);
%     eval(['!epstopdf ../figuras/adcp_bruto_',num2str(pp),'m.eps'])
drawnow



%% achando valores de velocidade entre cada par de estacao
load ../../hidrografia/posicoes_leste2.dat
pb = posicoes_leste2;
lonctd = (-(pb(:,5)/60 + pb(:,4)));
latctd = (-(pb(:,3)/60 + pb(:,2)));
  
dctd = [0; cumsum(sw_dist(latctd,lonctd,'km'))];
dadcp = [0; cumsum(sw_dist(latadcp,lonadcp,'km'))]; 

cont = 1;
for k = 0:30:dadcp(end)
    f = find(dadcp >= k & dadcp <= k+30); 
    lonadcp2(cont) = nanmean(lonadcp(f));
    latadcp2(cont) = nanmean(latadcp(f));
    Uadcp2(cont) = nanmean(Uadcp(f));
    Vadcp2(cont) = nanmean(Vadcp(f));
    pg2(cont) = nanmean(pg(f));
    mm2(cont) = round(nanmean(mm(f)));
    dd2(cont) = round(nanmean(dd(f)));
    hh2(cont) = round(nanmean(hh(f)));
    mn2(cont) = round(nanmean(mn(f)));
    cont = cont+1;
end

lonadcp = lonadcp2'; latadcp = latadcp2'; Uadcp = Uadcp2'; Vadcp = Vadcp2';
mm = mm2'; dd = dd2'; hh = hh2'; mn = mn2'; pg = pg2';

f = find(not(isnan(mm)));
lonadcp = lonadcp(f); latadcp = latadcp(f); Uadcp = Uadcp(f); Vadcp = Vadcp(f);
mm = mm(f); dd = dd(f); hh = hh(f); mn = mn(f); pg = pg(f);

clear c h fc u v lon lat U V dist f *2 cont k f

%%% gradeando o parametro de qualidade "percent good"

load ../../common/seagrid_leste2.dat;
grid = seagrid_leste2;
xg = grid(:,9); 
yg = grid(:,10);

l1 = max(grid(:,1));
l2 = max(grid(:,2));

xgi = reshape(xg,l2,l1);
ygi = reshape(yg,l2,l1);

pgi = griddata(lonadcp,latadcp,pg,xgi,ygi,'v4');


%  figure(3)
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
%  fc = 1;
%  m_contourf(xgi,ygi,pgi,[0:3:100]);shading flat;hold on
%  cc = colorbar;
%  caxis([0 120])
%  pos = get(cc,'Position');
%  m_usercoast('costa_leste.mat','patch',[0 0 0]);hold on
%  m_quiver(lonadcp,latadcp,Uadcp*fc,Vadcp*fc,0,'k')
%  m_quiver(-40,-11,0.5*fc,0,0,'w')
%  text = ['50 cm s^{-1}'];
%  m_text(-40,-11.3,text,'color','w','fontsize',10,'fontweight','bold')
%  text = [num2str(pp),' m'];
%  m_text(-34.8,-20.1,text,'color','k','fontweight','bold')
%  title('Velocidade Total - ADCP - OEII','fontsize',10,'fontweight','bold');
%  m_grid('box','fancy','yaxislocation','left','xaxislocation','bottom','fontsize',8);
%  set(cc,'Position',[pos(1) pos(2)*1.3 pos(3)/2.5 pos(4)/1.1])
%     print(3,'-depsc',['../figuras/adcp_trat_',num2str(pp),'m']);
%     eval(['!epstopdf ../figuras/adcp_trat_',num2str(pp),'m.eps']);


%%% REMOCAO DO TRANSPORTE DE EKMAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load ../mat/jj.dat % arquivo com calendario juliano para jan, fev, mar
Uek=[]; Vek=[]; HEk = []; 

for m = [2 3]; % intervalo temporal dos dados de ADCP (fev e mar)
    if m == 2 % looping para o mes de fev
       for d = 1:28 
           f1 = find(mm == m & dd == d); % achando os dados de adcp para o dia em questao
           if isempty(f1) == 0;
              mlat = mean(latadcp(f1)); % calculando lat media para procurar depois em tau
              f2 = find(jj(:,1) == m & jj(:,2) == d); % achando o dia juliano 
              jul = jj(f2,3); % dia juliano para carregar arquivo com tau
              cont = 0; w = 1;
              while w == 1
                 eval(['load ../mat/tau0',num2str(jul+cont),'.mat'])% carregando os dados de vento
                 [lj,li] = size(u);
                 Utau = reshape(u,li*lj,1); Vtau = reshape(v,li*lj,1);
                 lontau = reshape(lon,li*lj,1); lattau = reshape(lat,li*lj,1);
                 clear lon lat u v
                 f3 = find(lattau >= mlat-0.5 & lattau <= mlat+0.5); % achando tau para os pontos
                 Utau = nanmean(Utau(f3)); Vtau = nanmean(Vtau(f3)); % esta sera tau usada 
                 w = isnan(Utau);% checando se há alguma tau
                 cont = cont+1;
              end
              f0 = sw_f(mlat);
              modtau = sqrt(Utau^2+Vtau^2);
              velfric = sqrt(abs(modtau/1024));
              Hek = abs((0.4/f0)*velfric); HEk = [HEk Hek];
              if pp < Hek
                 uek = Vtau/(1024*f0); uek = uek/10;
                 vek = -Utau/(1024*f0); vek = vek/10;
                 Uadcp(f1) = Uadcp(f1) - uek;
                 Vadcp(f1) = Vadcp(f1) - vek;
              end
              Uek = [Uek; uek*ones(size(f1))];
              Vek = [Vek; vek*ones(size(f1))];
           end 
          clear f1 mlat f2 jul cont w Utau Vtau lontau lattau f3 f0 modtau velfric Hek
       end

    elseif m == 3 % looping para o mes de marco
       for d = 1:31 
           f1 = find(mm == m & dd == d); % achando os dados de adcp para o dia em questao
           if isempty(f1) == 0;
              mlat = mean(latadcp(f1)); % calculando lat media para procurar depois em tau
              f2 = find(jj(:,1) == m & jj(:,2) == d); % achando o dia juliano 
              jul = jj(f2,3); % dia juliano para carregar arquivo com tau
              cont = 0; w = 1;
              while w == 1
                 eval(['load ../mat/tau0',num2str(jul+cont),'.mat'])% carregando os dados de vento
                 [lj,li] = size(u);
                 Utau = reshape(u,li*lj,1); Vtau = reshape(v,li*lj,1);
                 lontau = reshape(lon,li*lj,1); lattau = reshape(lat,li*lj,1);
                 clear lon lat u v
                 f3 = find(lattau >= mlat-0.5 & lattau <= mlat+0.5); % achando tau para os pontos
                 Utau = nanmean(Utau(f3)); Vtau = nanmean(Vtau(f3)); % esta sera tau usada 
                 w = isnan(Utau);% checando se há alguma tau
                 cont = cont+1;
              end
              f0 = sw_f(mlat);
              modtau = sqrt(Utau^2+Vtau^2);
              velfric = sqrt(abs(modtau/1024));
              Hek = abs((0.4/f0)*velfric); HEk = [HEk Hek];
              if pp < Hek
                 uek = Vtau/(1024*f0); uek = uek/10;
                 vek = -Utau/(1024*f0); vek = vek/10;
                 Uadcp(f1) = Uadcp(f1); - uek;
                 Vadcp(f1) = Vadcp(f1); - vek;
              end
              Uek = [Uek; uek*ones(size(f1))];
              Vek = [Vek; vek*ones(size(f1))];
           end 
           clear f1 mlat f2 jul cont w Utau Vtau lontau lattau f3 f0 modtau velfric Hek
       end


    end

end
stop
% salvando um .mat com os dados filtrados para usar posteriormente 
% no programa de referenciamento da velocidade geostrofica
save ../mat/adcp_filt.mat Uadcp Vadcp lonadcp latadcp


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
%  m_usercoast('costa_leste.mat','patch',[0 0 0]);hold on
%  m_quiver(lonadcp,latadcp,Uek*fc,Vek*fc,0,'r')
%  m_quiver(-40,-11,0.5*fc,0,0,'w')
%  text = ['50 cm s^{-1}'];
%  m_text(-40,-11.3,text,'color','w','fontsize',10,'fontweight','bold')
%  text = [num2str(pp),' m'];
%  m_text(-34.8,-20.1,text,'color','k','fontweight','bold')
%  title('Velocidade Media na Camada de Ekman - OEII','fontsize',10,'fontweight','bold');
%  m_grid('box','fancy','yaxislocation','left','xaxislocation','bottom','fontsize',8);
%     print(4,'-depsc',['../figuras/adcp_ek_',num2str(pp),'m']);
%     eval(['!epstopdf ../figuras/adcp_ek_',num2str(pp),'m.eps']);
%  drawnow
%  
%  
%  figure(5)
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
%  m_contourf(xgi,ygi,pgi,[0:3:100]);shading flat;hold on
%  cc = colorbar;
%  caxis([0 120])
%  pos = get(cc,'Position');
%  m_usercoast('costa_leste.mat','patch',[0 0 0]);hold on
%  m_quiver(lonadcp,latadcp,Uadcp*fc,Vadcp*fc,0,'b')
%  m_quiver(-40,-11,0.5*fc,0,0,'w')
%  text = ['50 cm s^{-1}'];
%  m_text(-40,-11.3,text,'color','w','fontsize',10,'fontweight','bold')
%  text = [num2str(pp),' m'];
%  m_text(-34.8,-20.1,text,'color','k','fontweight','bold')
%  title('Velocidade sem Deriva de Ekman - ADCP - OEII','fontsize',10,'fontweight','bold');
%  m_grid('box','fancy','yaxislocation','left','xaxislocation','bottom','fontsize',8);
%  set(cc,'Position',[pos(1) pos(2)*1.3 pos(3)/2.5 pos(4)/1.1])
%     print(5,'-depsc',['../figuras/adcp_',num2str(pp),'m']);
%     eval(['!epstopdf ../figuras/adcp_',num2str(pp),'m.eps']);
%  drawnow


%%% INTERPOLACAO POR ANALISE OBJETIVA

lc = 1; % em graus, entao, precisa passar as velocidades para grau/s
E = 0.02;

Uadcp = Uadcp*9e-6; % fator de conversao para grau/s (9e-6)
Vadcp = Vadcp*9e-6;

%%% APLICANDO CONDICOES DE CONTORNO PARA PSI ********************************
%%% NAO ESCORREGAMENTO

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
xcont=xcont';
ycont=ycont';

close(1); 

xf=[min(xcont) xcont];
yf=[max(ycont) ycont];
[xf,yf] = m_ll2xy(xf,yf,'patch');

% analise objetiva vetorial para obter psi observado
[psiob] = vectoa(xg,yg,lonadcp,latadcp,Uadcp,Vadcp,lc,E,0);

% passando psi para m^2/s
psiob = psiob*111120;

% eliminando valor medio de psi
psiob = psiob-mean(psiob);

%% aplicando condicao de contorno em psiob
%  xg2 = [xg' xcont]; yg2 = [yg' ycont];
%  psiob = [psiob zeros(size(xcont))];
xg2 = xg; yg2 = yg;

% analise objetiva escalar em psiob para aplicar cond de contorno
[psiob] = scaloa(xg',yg',xg2,yg2,psiob,lc*1.5,E);
psiob = reshape(psiob,l2,l1);
[uo,vo] = psi2uv(xgi,ygi,psiob);

lpsi = -30000:800:30000;
inc = ( max(max(psiob)) - min(min(psiob)) ) / 40;
lpsi2 = -30000:inc:30000;
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

m_contourf(xgi,ygi,psiob,lpsi);shading flat;hold on
%  caxis([-30000 30000])
m_contour(xgi,ygi,psiob,lpsi2,'w')
m_quiver(xgi,ygi,uo*fc,vo*fc,0,'k')
fill(xf,yf,[0.8 0.8 0.8]);
m_usercoast('costa_leste.mat','patch',[0 0 0])
m_quiver(-40,-11,0.5*fc,0,0,'w')
text = ['50 cm s^{-1}'];
m_text(-40,-11.3,text,'color','w','fontsize',10,'fontweight','bold')
text = ['1-20 m'];
m_text(-34.8,-20.1,text,'color','r','fontweight','bold')
title('\Psi Observado - OEII [m^2 s^{-1}]','fontsize',10,'fontweight','bold');
text = ['Contornos: ',int,' m^2s^{-1}'];
m_text(-40.5,-11.8,text,'color','w','fontsize',8,'fontweight','bold')
m_grid('box','fancy','yaxislocation','left','xaxislocation','bottom','fontsize',10);
cc = colorbar;
pos = get(cc,'Position');
set(cc,'Position',[pos(1) pos(2)*1.3 pos(3)/2.5 pos(4)/1.1])
   print(6,'-depsc',['../figuras/adcp_AO_',num2str(pp),'m']);
   eval(['!epstopdf ../figuras/adcp_AO_',num2str(pp),'m.eps'])
drawnow

Uadcp = uo; Vadcp = vo;
save ../mat/adcp_filt_AO.mat Uadcp Vadcp





















