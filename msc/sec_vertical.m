%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	PROGRAMA PARA PLOTAR AS SECOES
%       HIDROGRAFICAS - Oceano Leste 2
%        mai/2008 - Mestrado - IOUSP
%          Rafael Guarino Soutelino
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Calcula Secao Vertical Velocidade Geostrofica via Metodo Dinamico Referenciado

clc;clear all;close all;warning off

% configuracoes---------------------------------
alis = 's';
jan = 11;
cruz = 'oeii'
CRUZ = 'OEII'
transp = 's';
cc = 's';
maperr = 'n';
%-----------------------------------------------

% carregando os dados de ADCP
% carregando os dados e tomando variaveis de lon,lat,u,v,t
eval(['load ../adcp/exper1/contour/',cruz,'_uv.mat']);
eval(['load ../adcp/exper1/contour/',cruz,'_xy.mat']);
u = uv(:,1:2:end-1); v = uv(:,2:2:end); 
lonadcp = xyt(1,:); lonadcp = lonadcp-360; latadcp = xyt(2,:); t = xyt(3,:); 
mod = sqrt(u.^2+v.^2);

% removendo spikes remanescentes
f = find(mod >= 1.5);
u(f)=nan;v(f)=nan;

% salvando as posicoes que limitam as radiais
rad1 = [213:278]; rad2 = [371:-1:303]; rad3 = [379:442]; rad4 = [526:-1:457];
rad5 = [528:607]; rad6 = [725:-1:628]; rad7 = [827:894];
rad8 = [977:-1:915]; rad9 = [978:1039]; rad10 = [1112:-1:1054]; 
rad11 = [1113:1176]; rad12 = [1295:-1:1196];

%%% CARREGANDO OS DADOS TERMOHALINOS

r = input('Escolha a(s) radial(ais) a serem processadas:  ');
for RR = r

%%% looping para escolher a Radial a ser plotada

%  RR = input('Entre com a radial desejada (Ex: 2):  '); 

if RR==1;
    xx=[1 2 3 4 5 6 7];
elseif RR==2;
    xx=16:-1:8;
elseif RR==3;
    xx=17:23;
elseif RR==4;
    xx=33:-1:24;
elseif RR==5;
    xx=34:42;
elseif RR==6;
    xx=51:-1:43;
elseif RR==7;
    xx=52:60;
elseif RR==8;
    xx=68:-1:61;
elseif RR==9;
    xx=69:78;
elseif RR==10;
    xx=90:-1:79;
elseif RR==11;
    xx=91:102;
else 
    xx=112:-1:103;
end


RR=num2str(RR);
lx=length(xx);

RADIAL=(['Radial ',RR]);
radial=(['rad',RR]);


%%% Leitura dos dados

load ../hidrografia/posicoes_leste2.dat; % arquivo com as posicoes das estacoes
pb = posicoes_leste2;
profest=pb(:,6);
nest=pb(:,1);

% convertendo para grau e decimo de grau
latg=pb(:,2); latm=pb(:,3); lat=latg+latm/60;
long=pb(:,4); lonm=pb(:,5); lon=long+lonm/60;

% prof max de plotagem
%  pplot = 1015;
pplot = 2000;

clear  latg latm long lonm
c=1;
for k=1:lx;
  XX=num2str(xx(k));

  eval(['load ../../dados/leste2/ctd/filtrados/lesteII_ctd',XX,'.dat']);
  dados=eval(['lesteII_ctd',XX]);
  eval(['clear lesteII_ctd',XX]);

  p=dados(:,1); % motando vetores com os perfis das propriedades
  t=dados(:,2);
  s=dados(:,3);
  f=find(p>=15 & p<=pplot); % colocando todas as estacoes iniciando da mesma profundidade (15m)
  if max(p) >= 20
     T(1:length(f),c)=t(f);
     S(1:length(f),c)=s(f);
     zb(c)=p(length(p)); % a profundidade maxima é cosiderada a ultima perfilada
%    pmax(k)=profest(find(nest==xx(k))); % prof. max. é a profundidade local

     lats(c)=-lat(find(nest==xx(k)));
     lons(c)=-lon(find(nest==xx(k)));
     c=c+1;
     clear s t p XX dados latg latm long lonm 
  else
     disp('Estacao nao atinge 100m')
  end
end


clear  p XX dados latg latm long lonm posicoes_leste2 f k lat lon

f=find(T==0); % substituindo 0 por NaN

T(f)=NaN;
S(f)=NaN;

[lz,lx]=size(T);
z=15:lz+14; z=z';  % motando um vetor com os valores de pressao
zmax=max(z);

[dst,ang]=sw_dist(lats,lons,'km'); % calculando distancia entre as estacoes
x=[0 cumsum(dst)];
xmax=max(x);

clear ang f dst

%%% TOPOGRAFIA %%%%%%%%%%%%%%%

xbi=0:0.1:(x(end));
zbi=interp1(x,zb,xbi,'cubic'); % suavizando a topografia

% alimentando os vetores pra fechar o poligono
xbi=[min(xbi) xbi];
zbi=[(max(zbi)) zbi];

%%% INTERPOLANDO CAMPOS TERMOHALINOS PARA ENRIQUECER A EXTRAPOLACAO
%  disp('Interpolando campos')
%  xi = 0:10:x(end);
%  Ti = griddata(x,z,T,xi,z);
%  Si = griddata(x,z,S,xi,z);

% no caso de nao interpolar (melhor opcao), muda-se o nome das variaveis
% antes da extrapolacao para nao perder os campos brutos
Ti=T;
Si=S;
xi=x;


%%% EXTRAPOLANDO CAMPOS TERMOHALINOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*****************************************************************************

disp('Extrapolando campos......')
decT = 2; % taxa de decaimento do gradiente
decS = 2;

for i = 1:length(z)   
   f = find(isnan(Ti(i,:))==1); % descobrindo os nan da linha
   if isempty(f)==1; clc
   else
     ff = find(diff(f) > 1);
     if isempty(ff)==1  % teste para saber se ha montes submarinos

        % extrapolando o talude num perfil sem monte
        if f(end)+2 > lx % checando se ha alguem para calcular o gradiente
            for k=fliplr(f)
              Tult = Ti(i,k+1); Sult = Si(i,k+1);
              Ti(i,k) = Tult; Si(i,k) = Sult;
            end
        else
            dt = Ti(i,f(end)+2) - Ti(i,f(end)+1); % grad entre as 2 ult. estacoes   
            ds = Si(i,f(end)+2) - Si(i,f(end)+1);
            for k=fliplr(f)
              dt = dt/decT; ds = ds/decS;
              Tult = Ti(i,k+1); Sult = Si(i,k+1);                      
              Tex = Tult - dt; Sex = Sult - ds;        
              Ti(i,k) = Tex; Si(i,k) = Sex;              
            end
        end

     else

         % extrapolando o talude num perfil com monte
         f1 = f(1:ff); % NaNs do talude
         if isnan(Ti(i,ff+2)) == 1 % checando se ha alguem para calcular o gradiente
            for k=fliplr(f1)
              Tult = Ti(i,k+1); Sult = Si(i,k+1);
              Ti(i,k) = Tult; Si(i,k) = Sult;
            end
         else             
            dt = Ti(i,f1(end)+2) - Ti(i,f1(end)+1); % grad entre as 2 ult. estacoes   
            ds = Si(i,f1(end)+2) - Si(i,f1(end)+1);
            for k=fliplr(f1)
              dt = dt/decT; ds = ds/decS;
              Tult = Ti(i,k+1); Sult = Si(i,k+1);                      
              Tex = Tult - dt; Sex = Sult - ds;        
              Ti(i,k) = Tex; Si(i,k) = Sex;              
            end     
         end

         % extrapolando lado esquerdo do monte
         f2 = f(ff+1:end); % NaNs do monte
         
         for k=f2(1:ceil(length(f2)/2))
           Tult = Ti(i,k-1); Sult = Si(i,k-1);                             
           Ti(i,k) = Tult; Si(i,k) = Sult;              
         end   
         % extrapolando lado direito do monte
         for k=fliplr(f2(ceil(length(f2)/2):end))
           Tult = Ti(i,k+1); Sult = Si(i,k+1);                              
           Ti(i,k) = Tult; Si(i,k) = Sult;              
         end  
               
      end
   end
end

clear i k f ff f1 f2 ds dt 
% *****************************************************************************

Di=sw_dens0(Si,Ti)-1000;

% ******************************************************************************

%%% CALCULANDO VELOCIDADE GEOSTROFICA VIA METODO DINAMICO CLASSICO

nr = 150;
f0 = sw_f(-16); % latitude central da grade OEII
% anomalia do geopotencial relativa a superficie
agp = sw_gpan(Si,Ti,z);

% anomalia do geopotencial relativa ao NR
agp = agp - ones(size(z)) * agp(nr,:);


% calculando velocidade puramente geostrofica

dagp = (diff(agp'))';
dx = diff(xi); dx = dx*1000; % passando para o SI
dx = ones(size(z)) * dx; 

vg = (1/f0) .* (dagp./dx); 

% criando eixo de distancia para plotar vg

for j=1:length(xi)-1
  xiv(j) = (xi(j)+xi(j+1))./2;
end

%%% REFERENCIANDO O METODO DINAMICO A PARTIR DO ADCP

% carregando os dados de ADCP da OEII para a radial em questao
eval(['u = u(:,rad',num2str(RR),');'])
eval(['v = v(:,rad',num2str(RR),');'])
eval(['x = lonadcp(rad',num2str(RR),');'])
eval(['y = latadcp(rad',num2str(RR),');'])

% escolhendo o nível de referencia
n = near(zc,nr,1);

if alis == 's'
   u = u(n,:); v = v(n,:);
   % removendo NaNs e spikes remanescentes
   f = find(isnan(u)==0);
   u = u(f); v = v(f); x = x(f); y = y(f);
   u=weim(jan,'hann',u);
   v=weim(jan,'hann',v);
else
   u = u(n,:); v = v(n,:);
end

% estabelecendo a velocidade normal a secao
% rotacionando os eixos para obter componente normal
[ans,ang] = sw_dist([y(1) y(end)],[x(1) x(end)],'km');
[mod,dir] = uv2intdir(u,v,0,ang);
[usec,vsec] = intdir2uv(mod,dir,0,0);
% daqui pra frente usaremos apenas vsec

% criando eixo de distancias para o ADCP
% este eixo tem que ser compativel ao de CTD. 
% o eixo do CTD sera considerado o "zero"

% calculando a distancia entre os pontos iniciais de ctd e adcp
d = sw_dist([lats(1) y(1)],[lons(1) x(1)]);

% calculando eixo de distancias do ADCP
xadcp = sw_dist(y,x,'km');
xadcp = [0 cumsum(xadcp)] + d;

%%% INTERPOLACAO PARA IGUALAR OS EIXOS

% criando grade nova
xi = min([xiv xadcp]):5:max([xiv xadcp]);
zi = min(z):5:max(z);
[xg,zg] = meshgrid(xi,zi);
[L1,L2] = size(xg);

% interpolando linearmente
vi = griddata(xiv,z,vg,xg,zg,'linear',{'QJ'});

% interp1 nos dados de ADCP no nivel de referencia
vnr = interp1(xadcp,vsec,xi);

% referenciando aos dados de ADCP
vi = vi - (ones(L1,1)*vnr);

% redimensionando grade nova para entrar na objmap.m
xg2 = reshape(xg,1,L1*L2);
zg2 = reshape(zg,1,L1*L2);

% redimensionando grade CTD antiga e removendo NaNs
x2 = reshape(xg,1,L1*L2);
z2 = reshape(zg,1,L1*L2);
vg2 = reshape(vi,1,L1*L2);
f = find(isnan(vg2)==0);
x2 = x2(f); z2 = z2(f); vg2 = vg2(f);

% escolhendo parametros de interpolacao
lx = 100;
lz = 50; 
E = 0.02;

% implementando condicoes de contorno
% usando topografia oriunda das estacoes CTD

if cc == 's'
   x2 = [x2 xbi];
   z2 = [z2 zbi];
   vg2 = [vg2 zeros(size(xbi))];
end

% interpolando
[vo,er] = objmap(x2',z2',vg2',xg2',zg2',[lx lz],E); 
vo = reshape(vo,L1,L2);
er = reshape(er,L1,L2);
er = 100*sqrt(er);
vo = -vo;% para acertar a barra de cores

% plotando a secao

lv1 = -0.8:0.005:0.5;
lv2 = -1:0.1:1;
jet2 = jet;
jet2 = flipud(jet);
close(1);

if transp == 'n'

figure(1);
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

subplot(211)
contourf(xg,-zg,vo,lv1);shading flat; 
hold on
[c,h] = contour(xg,-zg,vo,lv2,'k');
clabel(c,h);
fill(xbi,-zbi,[0 0 0])
caxis([lv1(1) lv1(end)]); cc = colorbar;
  colormap(jet2);
  pos = get(cc,'Position');
  set(cc,'Position',[pos(1) pos(2) pos(3)/2 pos(4)])
axis([min(min(xg)) max(max(xg)) -pplot 0])
xlabel('Distancia ao longo da Radial [km]')
ylabel('Profundidade [m]')
tit = ['Radial ',num2str(RR),' - V. Geostroficas Ref. por ADCP'];
title(tit,'fontweight','bold')
pbaspect([1 0.6 1])

zb2 = zb;
load ../common/etopo2_leste.mat;
m_proj('mercator','long',[-41 -33.7],'lat',[-20.5 -10.2],'on');

subplot(212)
[c,h] = m_contour(xb,yb,zb,[-200 -1000],'k'); hold on
clabel(c,h,'labelspacing',500);
p = m_plot(lonadcp,latadcp,'.k','markersize',6);
set(p,'color',[.5 .5 .5])
eval(['m_plot(lonadcp(rad',num2str(RR),'),latadcp(rad',num2str(RR),'),''.r'',''markersize'',6)'])
m_usercoast('../common/costa_leste.mat','patch',[0.542 0.422 0.000])
m_grid

eval(['print -depsc figuras/sec_adcp_objmap_rad',num2str(RR),'.eps'])
eval(['!epstopdf figuras/sec_adcp_objmap_rad',num2str(RR),'.eps'])

end

if maperr == 's'
figure(2);
set(gcf,'Color',[1 1 1])
contourf(xg,-zg,er,[1:0.1:30]);shading flat; 
hold on
%  [c,h] = contour(xg,-zg,er,[0:1:20],'k');
%  clabel(c,h);
fill(xbi,-zbi,[0 0 0])
%  fill(bx/1000,bz,[0 0 0])
caxis([1 30]); cc = colorbar;
  pos = get(cc,'Position');
  set(cc,'Position',[pos(1) pos(2)*2.3 pos(3)/2 pos(4)/1.5])
axis([min(min(xg)) max(max(xg)) -300 0])
xlabel('Distancia ao longo da radial [km]')
ylabel('Depth [m]')
tit = ['Radial ',num2str(RR),' - Erro de Interpolacao'];
title(tit,'fontweight','bold')
pbaspect([1 0.6 1])

eval(['print -depsc ../figuras/err_sec_adcp_rad',num2str(RR),'.eps'])
eval(['!epstopdf ../figuras/err_sec_adcp_rad',num2str(RR),'.eps'])
end

end

if transp == 's'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                VOLUME TRANSPORT CALCULATION                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CALCULA OS VALORES DE VELOCIDADE NO CENTRO DA CELULA DE GRADE %%%
xg2 = xg;

vom=0.5*(vo(:,2:end)+vo(:,1:end-1));
%  vom=0.5*(vo(2:end,:)+vo(1:end-1,:));

%%% MATRIZ DE dx's E dz's %%%

dx=xg2(:,2:end)-xg2(:,1:end-1); 
dz=zg(2:end,:)-zg(1:end-1,:); dz=-dz;

xm=0.5*(xg2(:,2:end)+xg2(:,1:end-1));
xm=0.5*(xm(2:end,:)+xm(1:end-1,:));
zm=0.5*(zg(:,2:end)+zg(:,1:end-1)); zm = -zm;
zm=0.5*(zm(2:end,:)+zm(1:end-1,:)); 

%%% ELIMINANDO 10 ULTIMOS NIVEIS SIGMA %%%
xm = xm(1:end-0,:);
zm = zm(1:end-0,:);
vom = vom(1:end-1,:);

figure(5)
set(gcf,'color','w')
contourf(xg,-zg,vo,lv1);shading flat; 
hold on
[c,h] = contour(xg,-zg,vo,[0 0],'k');
clabel(c,h);
%  fill([min(xi) xi max(xi)],-[max(zbi)+100 zbi max(zbi)+100],[1 1 1])
fill(xbi,-zbi,[0 0 0])
caxis([lv1(1) lv1(end)]); % cc = colorbar;
  colormap(jet2);
%    pos = get(cc,'Position');
%    set(cc,'Position',[pos(1) pos(2) pos(3)/2 pos(4)])
axis([min(min(xg)) max(max(xg)) -pplot 0])
xlabel('Along-track distance [km]')
ylabel('Depth [m]')
tit = ['Radial ',num2str(RR),' - V. Geostroficas Ref. por ADCP'];
title(tit,'fontweight','bold')
pbaspect([1 0.6 1])


zm = -zm;
% selecting the currents
disp(' ')
disp('Selecione os limites laterais e verticais da CB para sul que deseja medir: ')
disp(' ')
sul1 = ginput(4);
fsul1 = find(vom<0 & xm>=sul1(1,1) & xm<=sul1(2,1) & -zm<=sul1(3,2) & -zm>=sul1(4,2));
p = plot(xm(fsul1),-zm(fsul1),'k.','markersize',3)
set(p,'color',[.7 .7 .7])

disp(' ')
disp('Selecione os limites laterais e verticais do 2o fluxo para sul que deseja medir: ')
disp(' ')
sul2 = ginput(4);
fsul2 = find(vom<0 & xm>=sul2(1,1) & xm<=sul2(2,1) & -zm<=sul2(3,2) & -zm>=sul2(4,2));
p = plot(xm(fsul2),-zm(fsul2),'k.','markersize',3)
set(p,'color',[.7 .7 .7])

disp(' ')
disp('Selecione os limites laterais e verticais do 3o fluxo para sul que deseja medir: ')
disp(' ')
sul3 = ginput(4);
fsul3 = find(vom<0 & xm>=sul3(1,1) & xm<=sul3(2,1) & -zm<=sul3(3,2) & -zm>=sul3(4,2));
p = plot(xm(fsul3),-zm(fsul3),'k.','markersize',3)
set(p,'color',[.7 .7 .7])

disp(' ')
disp('Selecione os limites laterais e verticais da CCI para norte que deseja medir: ')
disp(' ')
nor1 = ginput(4);
fnor1 = find(vom>0 & xm>=nor1(1,1) & xm<=nor1(2,1) & -zm<=nor1(3,2) & -zm>=nor1(4,2));
plot(xm(fnor1),-zm(fnor1),'w.','markersize',3)

disp(' ')
disp('Selecione os limites laterais e verticais do 2o fluxo para norte que deseja medir: ')
disp(' ')
nor2 = ginput(4);
fnor2 = find(vom>0 & xm>=nor2(1,1) & xm<=nor2(2,1) & -zm<=nor2(3,2) & -zm>=nor2(4,2));
plot(xm(fnor2),-zm(fnor2),'w.','markersize',3)

disp(' ')
disp('Selecione os limites laterais e verticais do 3o fluxo para norte que deseja medir: ')
disp(' ')
nor3 = ginput(4);
fnor3 = find(vom>0 & xm>=nor3(1,1) & xm<=nor3(2,1) & -zm<=nor3(3,2) & -zm>=nor3(4,2));
plot(xm(fnor3),-zm(fnor3),'w.','markersize',3)

%%% TRANSPORTE DE VOLUME %%%

Tsul1 = sum(vom(fsul1).*dx(fsul1).*dz(fsul1))*1e-6*1e3;
Tsul2 = sum(vom(fsul2).*dx(fsul2).*dz(fsul2))*1e-6*1e3;
Tsul3 = sum(vom(fsul3).*dx(fsul3).*dz(fsul3))*1e-6*1e3;

Tnorte1 = sum(vom(fnor1).*dx(fnor1).*dz(fnor1))*1e-6*1e3;
Tnorte2 = sum(vom(fnor2).*dx(fnor2).*dz(fnor2))*1e-6*1e3;
Tnorte3 = sum(vom(fnor3).*dx(fnor3).*dz(fnor3))*1e-6*1e3;

%%% Criando meios de plotar os valores na propria figura
clear text
format bank

xsul1 = mean(mean(xm(fsul1)));
zsul1 = mean(mean(zm(fsul1)));
Tsul1 = (round(Tsul1*10))/10;
tsul1 = [num2str(abs(Tsul1)),' Sv'];
text(xsul1,-zsul1,tsul1,'fontweight','bold')
vmaxsul1 = round((min(min(vom(fsul1))))*100);
tsul1 = [num2str(vmaxsul1),' cm s^{-1}'];
text(xsul1,-(zsul1+100),tsul1,'fontweight','bold')

xsul2 = mean(mean(xm(fsul2)));
zsul2 = mean(mean(zm(fsul2)));
Tsul2 = (round(Tsul2*10))/10;
tsul2 = [num2str(abs(Tsul2)),' Sv'];
text(xsul2,-zsul2,tsul2,'fontweight','bold')
vmaxsul2 = round((min(min(vom(fsul2))))*100);
tsul2 = [num2str(vmaxsul2),' cm s^{-1}'];
text(xsul2,-(zsul2+100),tsul2,'fontweight','bold')

xsul3 = mean(mean(xm(fsul3)));
zsul3 = mean(mean(zm(fsul3)));
Tsul3 = (round(Tsul3*10))/10;
tsul3 = [num2str(abs(Tsul3)),' Sv'];
text(xsul3,-zsul3,tsul3,'fontweight','bold')
vmaxsul3 = round((min(min(vom(fsul3))))*100);
tsul3 = [num2str(vmaxsul3),' cm s^{-1}'];
text(xsul3,-(zsul3+100),tsul3,'fontweight','bold')

xnorte1 = mean(mean(xm(fnor1)));
znorte1 = mean(mean(zm(fnor1)));
Tnorte1 = (round(Tnorte1*10))/10;
tnorte1 = [num2str(abs(Tnorte1)),' Sv'];
text(xnorte1,-znorte1,tnorte1,'fontweight','bold')
vmaxnorte1 = round((max(max(vom(fnor1))))*100);
tnorte1 = [num2str(vmaxnorte1),' cm s^{-1}'];
text(xnorte1,-(znorte1+100),tnorte1,'fontweight','bold')

xnorte2 = mean(mean(xm(fnor2)));
znorte2 = mean(mean(zm(fnor2)));
Tnorte2= (round(Tnorte2*10))/10;
tnorte2 = [num2str(abs(Tnorte2)),' Sv'];
text(xnorte2,-znorte2,tnorte2,'fontweight','bold')
vmaxnorte2 = round((max(max(vom(fnor2))))*100);
tnorte2 = [num2str(vmaxnorte2),' cm s^{-1}'];
text(xnorte2,-(znorte2+100),tnorte2,'fontweight','bold')

xnorte3 = mean(mean(xm(fnor3)));
znorte3 = mean(mean(zm(fnor3)));
Tnorte3 = (round(Tnorte3*10))/10;
tnorte3 = [num2str(abs(Tnorte3)),' Sv'];
text(xnorte3,-znorte3,tnorte3,'fontweight','bold')
vmaxnorte3 = round((max(max(vom(fnor3))))*100);
tnorte3 = [num2str(vmaxnorte3),' cm s^{-1}'];
text(xnorte3,-(znorte3+100),tnorte3,'fontweight','bold')

fill(xbi,-zbi,[0 0 0])

print(5,'-dpng',['figuras/transp_adcp_rad',num2str(RR)])
print(5,'-depsc',['figuras/transp_adcp_rad',num2str(RR)])

end
