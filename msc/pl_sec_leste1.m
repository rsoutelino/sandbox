%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	PROGRAMA PARA PLOTAR AS SECOES
%       HIDROGRAFICAS - Oceano Leste 1
%        out/2006 - Mestrado - IOUSP
%          Rafael Guarino Soutelino
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Plota distribuições verticais de T e S e sigma theta, 
%   Extrapola campos de T e S
%   Calcula Velocidade Geostrofica via Metodo Dinamico



%  clc;clear all;close all

%  for RR = 1:13
for RR = 13

%%% looping para escolher a Radial a ser plotada

%  RR = input('Entre com a radial desejada (Ex: 2):  '); 

if RR==1;
    xx=1:11;
elseif RR==2;
    xx=12:18;
elseif RR==3;
    xx=26:-1:19;
elseif RR==4;
    xx=27:32;
elseif RR==5;
    xx=40:-1:33;
elseif RR==6;
    xx=41:48;
elseif RR==7;
    xx=[57:-1:53 51:-1:49];
elseif RR==8;
    xx=[58:63 65 67];
elseif RR==9;
    xx=[74:-1:70 68];
elseif RR==10;
    xx=[75:80 82];
elseif RR==11;
    xx=92:-1:83;
elseif RR==12;
    xx=93:101;
elseif RR==13; 
    xx=113:-1:102;
end

RR=num2str(RR);
lx=length(xx);

RADIAL=(['Radial ',RR]);
radial=(['rad',RR]);


%%% Leitura dos dados

load posicoes_leste1.dat; % arquivo com as posicoes das estacoes
pb = posicoes_leste1;
nest=pb(:,1);
lat=pb(:,2);
lon=pb(:,3);

% prof max de plotagem
pplot = 2000;

c=1;
for k=1:lx;
  XX=num2str(xx(k));

  eval(['load  ../../../dados/leste1/filtrados/lesteI_ctd',XX,'.dat']);
  dados=eval(['lesteI_ctd',XX]);
  eval(['clear lesteI_ctd',XX]);

  p=dados(:,1); % motando vetores com os perfis das propriedades
  t=dados(:,2);
  s=dados(:,3);
  f=find(p>=15 & p<=pplot); % colocando todas as estacoes iniciando da mesma profundidade (15m)
if max(p) >= 100
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

disp('Extrapolando campos')
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

% PLOTANDO AS DISTRIBUICOES TERMOHALINAS VERTICAIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

P = ['TiSiDi'];

lT=2:30; 
lS=34:0.1:37.5;  
lD=23.3:0.1:28; 
l = ['lTlSlD'];

t1 = ['Temperatura ( \circ C) - OEI - ',RADIAL];
t2 = ['Salinidade - OEI - ',RADIAL];
t3 = ['Densidade Potencial (kg m^{-3}) - OEI - ',RADIAL];
tit = ['t1t2t3'];

c=0;
for k = 1:3
  figure(k)
  set(k,'Color','w'); 
  hold on;
  eval(['[c1,h1] = contourf(xi,z,',P(c+k:c+k+1),',',l(c+k:c+k+1),');']); 
  shading flat; hold on
  set(gca,'ydir','reverse')
%    clabel(c1,h1,'labelspacing',500)
  axis([0 max(xi) 0 pplot])
  fill(xbi,zbi,[.8 .8 .8])
  plot(x,10*ones(size(x)),'kv')
  cc = colorbar;
  pos = get(cc,'Position');
  pbaspect([1 0.7 1])
  set(cc,'Position',[pos(1) pos(2)/.55 pos(3)/2.5 pos(4)/1.3])
  hold off;
  eval(['title(',tit(c+k:c+k+1),',''fontsize'',12,''fontweight'',''bold'')']); 
  xlabel('Distancia Aproximada da Costa (km)')
  ylabel('Profundidade (m)')
  drawnow
     print(k,'-depsc',['../figuras/',P(c+k),'_OEI_',radial,'.eps']);
     eval(['!epstopdf ../figuras/',P(c+k),'_OEI_',radial,'.eps'])
%       eval(['!rm -rf figuras/',P(c+k),'_OEI_',radial,'.eps'])
c=c+1; 
end
clear k c

% ******************************************************************************

%%% CALCULANDO VELOCIDADE GEOSTROFICA VIA METODO DINAMICO

nr = 1000;
f0 = sw_f(-16); % latitude central da grade OEII
% anomalia do geopotencial relativa a superficie
agp = sw_gpan(Si,Ti,z);

% anomalia do geopotencial relativa ao NR
agp = agp - ones(size(z)) * agp(nr,:);


% calculando velocidade geostrofica

dagp = (diff(agp'))';
dx = diff(xi); dx = dx*1000; % passando para o SI
dx = ones(size(z)) * dx; 

vg = (1/f0) .* (dagp./dx); 
vg(nr+1:end,:) = -vg(nr+1:end,:);

% criando eixo de distancia para plotar vg

for j=1:length(xi)-1
  xiv(j) = (xi(j)+xi(j+1))./2;
end

% plotando

titvg = ['Velocidade Geostrofica OEI - ',RADIAL,' (m s^{-1})'];
lV = -2:0.05:2;

figure(4)
  set(4,'Color','w'); 
  contourf(xiv,z,vg,lV); shading flat ; hold on
  caxis([-.8 .8])
  set(gca,'ydir','reverse')
  axis([0 max(xi) 0 pplot])
  fill(xbi,zbi,[.8 .8 .8])
  plot(x,10*ones(size(x)),'kv')
  plot(xiv,10*ones(size(xiv)),'k*')
  cc = colorbar;
  pos = get(cc,'Position');
  pbaspect([1 0.7 1])
  set(cc,'Position',[pos(1) pos(2)/.55 pos(3)/2.5 pos(4)/1.3])
  pos=get(cc,'Position');
  cc1=-str2num(get(cc,'YTickLabel'));cc1(find(cc1==0))=0;
  set(cc,'YTickLabel',cc1)
  hold off;
  title(titvg,'fontsize',12,'fontweight','bold'); 
  xlabel('Distancia aproximada da Costa (km)')
  ylabel('Profundidade (m)')
  drawnow  
     print(4,'-depsc',['../figuras/vg_OEI_',radial,'.eps']);
     eval(['!epstopdf ../figuras/vg_OEI_',radial,'.eps'])
%       eval(['!rm -rf figuras/vg_OEI_',radial,'.eps'])


% salvando matrizes com os dados extrapolados
 
eval(['save ../mat/lesteI_agp_',radial,'.mat agp lons lats'])
eval(['save ../mat/lesteI_T_',radial,'.mat Ti'])
eval(['save ../mat/lesteI_S_',radial,'.mat Si'])

clear 
close all
end
