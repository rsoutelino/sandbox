%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MODELO NUMERICO HIDRODINAMICO 3D LINEAR PARA O ATLANTICO SUDOESTE
%           COM FORCANTE BASEADA EM VENTO SINOTICO NCEP 6h-6h
%           LISTA DE EXERCICIOS 3 - RAFAEL GUARINO SOUTELINO
%
%   Ultima Modificacao: 8/12/2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Configuracao da grade (pontos tipo eta, u, v)
%NORTE (k)
%7 EUEUEUE
%6 V V V V
%5 EUEUEUE
%4 V V V V
%3 EUEUEUE
%2 V V V V
%1 EUEUEUE
%  1234567 LESTE (l)

clear all;close all;clc

% carregando colorbar monotonica
load redblue

% carregando linha de costa
load costa.mat; lonlc = ncst(:,1); latlc = ncst(:,2);

% carregando dados de vento do NCEP espacados de 6 em 6 h (0h 01/07/2009 --> 0h 01/08/2009)
ncw = netcdf('X189.62.144.10.318.18.3.28.nc');
lonw = ncw{'lon'}(:); lonw = lonw - 360;
latw = ncw{'lat'}(:);
uw = ncw{'uwnd'}(:);
ncw = netcdf('X189.62.144.10.318.18.4.14.nc');
vw = ncw{'vwnd'}(:);

% corrigindo os dados
scfac  = ncw{'vwnd'}.scale_factor(:);
offset = ncw{'vwnd'}.add_offset(:);
uw = uw*scfac + offset;
vw = vw*scfac + offset;

% PARAMETROS DO MODELO
dx = 30000; dy = dx; % espaçamento de grade
dz = 100; % espacamento vertical
dt = 60; % passo de tempo
ttot = (30*24*3600)/dt; % 30 dias

% declarando limites do modelo
lonlim = [-60 -30]; latlim = [-40 -10];
m_proj('mercator','lat',[latlim],'lon',[lonlim])

% estabelecendo grade do modelo
jend = sw_dist([latlim(1) latlim(1)],[lonlim(1) lonlim(2)],'km')*1000;
kend = sw_dist([latlim(1) latlim(2)],[lonlim(1) lonlim(1)],'km')*1000;
l = 1:dz:6000;
j = 0:dx:jend; k = 0:dy:kend; lk = length(k); lj = length(j); ll = length(l);

% estabelecendo grade em lat lon para plotar
lon = linspace(lonlim(1),lonlim(end),length(j));
lat = linspace(latlim(1),latlim(end),length(k));
[J,K] = meshgrid(j,k);
[lon,lat] = meshgrid(lon,lat);
clear j k

% interpolando batimetria do ETOPO 5 para a grade do modelo
disp('Interpolando Batimetria......')
[zb,xb,yb] = m_tbase([lonlim latlim]);
disp('Suavizando Batimetria......')
zb = smoo2(zb,-9999,21,0.5);
H = interp2(xb,yb,zb,lon,lat); H = -H;

% criando mascara para pontos oceanicos e continentais
mask = zeros(size(lon));
mask(H <= 0) = 0;
mask(H > 0) = 1;

% criando chave para definir pontos imediatamente adjacentes ao contorno
% para evitarmos o calculo dos termos nao-lineares nestas localidades
kint = zeros(size(mask));
for j = 2:lj-1
   for k = 2:lk-1
      if mask(k,j) == 1;
         if mask(k+1,j) == 0 | mask(k-1,j) == 0 | mask(k,j+1) == 0 |  mask(k,j-1) == 0
            kint(k,j) = 0;
         else
            kint(k,j) = 1;
         end
      end   
   end
end

kint(1:2,:) = 0;
kint(lk-1:lk,:) = 0;
kint(:,1:2) = 0;
kint(:,lj-1:lj) = 0;

% interpolando vento para a grade do modelo
disp('Interpolando Vento......')
lw = length(uw); Uw = []; Vw = [];
[lonw,latw] = meshgrid(lonw,latw);

for i = 1:lw
   Uw(i,:,:) = interp2(lonw,latw,squeeze(uw(i,:,:)),lon,lat);
   Vw(i,:,:) = interp2(lonw,latw,squeeze(vw(i,:,:)),lon,lat);
end

% plotando batimetria e vento para ilustrar  
figure
s = 0.1;
for i = 1:lw
   pcolor(lon,lat,H); shading flat; hold on; colorbar
   quiver(lonw,latw,squeeze(uw(i,:,:))*s,squeeze(vw(i,:,:))*s,0,'w')
   f = find(mask == 0);
   plot(lon(f),lat(f),'.k')
   axis('equal'); axis([lonlim latlim])
   title('Vento de 6/6h sobreposto a batimetria')
   xlabel('Longitude [W]')
   ylabel('Latitude [S]')
   hold off
   pause(0.01)
end

%%%%%%%%%  AUX  %%%%%%%%%%%%%
figure
f = find(mask == 1);
plot(lon(f),lat(f),'*g'); hold on
f = find(kint == 0);
plot(lon(f),lat(f),'*m');
f = find(mask == 0);
plot(lon(f),lat(f),'*k');

% calculando tensao de cisalhamento
rhoar = 1.025;
fric = 2.6e-3;
modw = sqrt(Uw.^2 + Vw.^2);
taux = fric.*rhoar.*Uw.*modw;
tauy = fric.*rhoar.*Vw.*modw;

% definicao de parametros
g = 9.8;
rho = 1024;
dt2 = dt*2; dx2 = dx*2; dy2 = dy*2; dx4 = dx*4; dy4 = dx*4; rho2 = rho*2;
rhoH = rho.*H;
rfric = 0.02; r = 1 + rfric.*dt2;
f0 = sw_f(lat);
difz = 0.001;    % coeficiente de difusao na vertical
difzdz=difz/dz;

% definindo condicoes iniciais do modelo --> repouso
eta1 = zeros(lk,lj);
u1 = zeros(lk,lj,ll);
v1 = zeros(lk,lj,ll);
eta2 = zeros(lk,lj);
u2 = zeros(lk,lj,ll);
v2 = zeros(lk,lj,ll);


% LOOP NO TEMPO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ETA = []; U = []; V = []; EK = [];
disp('  ')
disp('  ')
disp('========= Rodando o Modelo ============')
disp('  ')

for m = 1:ttot
   tempo = dt*m;
   disp(['Dia: ',num2str(tempo/86400)])
   % CALCULANDO CAMPO DE VENTO PARA CADA PASSO DE TEMPO
   % como o NCEP dispoe de resolucao de 6 h, e o passo de tempo eh bem menor
   % devemos interpolar linearmente o vento entre a dois intervalos consecutivos para estimar
   % o campo para este exato passo de tempo
   lw = size(taux); lw = lw(1);
   tw = 0:6*3600:[125*6*3600]-1;
   n = near(tw,tempo,2);
   
   Tx = []; Ty = [];
   for j = 1:lj
      for k = 1:lk
         Tx(k,j) = interp1(tw(n),squeeze(taux(n,k,j)),tempo);
         Ty(k,j) = interp1(tw(n),squeeze(tauy(n,k,j)),tempo); 
      end  
   end
  
   % EQUACAO DA CONTINUIDADE
   for j = 3:2:lj-2
      for k = 3:2:lk-2
         if mask(k,j) > 0;
           forcx = 0;
           forcy = 0;
           for l = 2:ll-1     
              forcx = (u1(k,j+1,l) - u1(k,j-1,l)).*dz/dx2 + forcx;
              forcy = (u1(k+1,j,l) - u1(k-1,j,l)).*dz/dy2 + forcy;
           end
           eta2(k,j) = eta1(k,j) - dt.*(forcx + forcy);      
         end
      end
   end

   % COMPONENTE ZONAL DA EQUACAO DO MOVIMENTO (x)   
    for j = 2:2:lj-1
        for k = 3:2:lk-2
            if (mask(k,j).*mask(k,j+1).*mask(k,j-1)) > 0
            for l = 2:ll-1
            dudzinf = (u1(k,j,l+1) - u1(k,j,l))*difzdz;
            if (l==2) dudzsup = -taux(k,j)/rho; end
            if (l~=2) dudzsup = (u1(k,j,l) - u1(k,j,l-1))*difzdz; end 
            vmed = (v1(k+1,j+1,l) + v1(k+1,j-1,l) + v1(k-1,j+1,l) + v1(k-1,j-1,l))/4;
            forc = f0(k,j)*vmed - g.*(eta2(k,j+1) - eta2(k,j-1))./dx2...
               +(dudzinf - dudzsup)/dz;
            u2(k,j,l) = u1(k,j,l) + forc*dt;
            end
            end
         end
    end

   % COMPONENTE MERIDIONAL DA EQUACAO DO MOVIMENTO (y)   
    for j = 3:2:lj-2
        for k = 2:2:lk-1
            if (mask(k,j).*mask(k,j+1).*mask(k,j-1)) > 0
            for l = 2:ll-1
            dvdzinf = (v1(k,j,l+1) - v1(k,j,l))*difzdz;
            if (l==2) dvdzsup = -tauy(k,j)/rho; end
            if (l~=2) dvdzsup = (v1(k,j,l) - v1(k,j,l-1))*difzdz; end 
            umed = (u1(k+1,j+1,l) + u1(k+1,j-1,l) + u1(k-1,j+1,l) + u1(k-1,j-1,l))/4;
            forc = -f0(k,j)*umed - g.*(eta2(k+1,j) - eta2(k-1,j))./dy2...
               +(dvdzinf - dvdzsup)/dz;
            v2(k,j,l) = v1(k,j,l) + forc*dt;
            end
            end
         end
    end

   %  Renovando variaveis
   eta1 = eta2; u1 = u2; v1 = v2;   
    
   % Preparando matrizes para plotar
   % Definindo u, v e eta nos pontos vazios (consequencia da grade alternada)
   for j = 2:2:lj-1
      for k = 2:2:lk-1    
         if mask(k,j) > 0
            for l = 2:ll-1
               u2(k,j,l) = (u2(k-1,j,l) + u2(k+1,j,l))/2; 
               v2(k,j,l) = (v2(k,j-1,l) + v2(k,j+1,l))/2;
            end
            eta2(k,j) = (eta2(k-1,j-1) + eta2(k-1,j+1) + eta2(k+1,j-1) + eta2(k+1,j+1))/4;
         end 
      end
   end
   
   % Definindo u e v nos pontos tipo eta
   for j = 3:2:lj-2
      for k = 3:2:lk-2    
         if mask(k,j) > 0
            for l = 2:ll-1
               u2(k,j,l) = (u2(k,j+1,l) + u2(k,j-1,l))/2;
               v2(k,j,l) = (v2(k+1,j,l) + v2(k-1,j,l))/2;
            end        
         end 
      end
   end

   % Definindo v e eta nos pontos tipo u
   for j = 2:2:lj-1
      for k = 3:2:lk-2    
         if mask(k,j) > 0
            for l = 2:ll-1
               v2(k,j,l) = (v2(k,j+1,l) + v2(k,j-1,l) + v2(k+1,j,l) + v2(k-1,j,l))/4; 
            end
            eta2(k,j) = (eta2(k+1,j) + eta2(k-1,j))/2;
         end 
      end
   end
  
   % Definindo u e eta nos pontos tipo v
   for j = 3:2:lj-2
      for k = 2:2:lk-1    
         if mask(k,j) > 0
            for l = 2:ll-1
               u2(k,j,l) = (u2(k,j+1,l) + u2(k,j-1,l) + u2(k+1,j,l) + u2(k-1,j,l))/4;  
            end
            eta2(k,j) = (eta2(k+1,j) + eta2(k-1,j))/2;
         end 
      end
   end

   % Implementando condicoes de contorno de nao-gradiente
   for j = 1:lj
      eta2(1,j) = eta2(2,j).*mask(1,j);  
      eta2(lk,j) = eta2(lk-1,j).*mask(lk,j);
      for l = 2:ll-1
         u2(1,j,l) = u2(2,j,l).*mask(1,j);  
         u2(lk,j,l) = u2(lk-1,j,l).*mask(lk,j);
         v2(1,j,l) = v2(2,j,l).*mask(1,j);  
         v2(lk,j,l) = v2(lk-1,j,l).*mask(lk,j);
      end
   end

   for k = 1:lk
      eta2(k,1) = eta2(k,2).*mask(k,1);  
      eta2(k,lj) = eta2(k,lj-1).*mask(k,lj);
      for l = 2:ll-1
         u2(k,1,l) = u2(k,2,l).*mask(k,1);  
         u2(k,lj,l) = u2(k,lj-1,l).*mask(k,lj);
         v2(k,1,l) = u2(k,2,l).*mask(k,1);  
         v2(k,lj,l) = v2(k,lj-1,l).*mask(k,lj);
      end
   end

   if rem(tempo,86400) == 0
   % Plotando resultados
   h = figure('visible','off');
   set(gcf,...
        'Color',[1 1 1],...
        'InvertHardcopy','on',...
        'PaperUnits','inches',...
        'Units','inches',...
        'PaperOrientation','portrait',...
        'PaperPosition',[0 0 18 9],...
        'PaperPositionMode','manual',...
        'PaperType','usletter',...
        'Position',[.2 .2 18 9],...
        'ShareColors','off',...
        'Clipping','on');
   
   subplot(122); s1 = 1000;
   pcolor(lon,lat,eta2);shading flat; hold on
   caxis([-0.01 0.01])
   colormap(redblue); colorbar('horiz')
   quiver(lon,lat,squeeze(u2(:,:,1))*s1,squeeze(v2(:,:,1))*s1,0,'k')
   axis('equal')
   axis([-60 -30 -40 -10])
   f = find(mask == 0);
   plot(lon(f),lat(f),'*k')
   title(['Elevacao [m] => ',num2str(tempo/86400),' dias'])
   xlabel('Longitude [W]')
   ylabel('Latitude [S]')
%     pause(0.002)
   hold off
  
   subplot(121); s2 = 2;
   pcolor(lon,lat,sqrt(Tx.^2 + Ty.^2)); shading flat; hold on; 
   caxis([-1.5 1.5]); colorbar('horiz'); 
   quiver(lon(1:3:end,1:3:end),lat(1:3:end,1:3:end),Tx(1:3:end,1:3:end)*s2,Ty(1:3:end,1:3:end)*s2,0,'k'); hold on
   axis('equal')
   axis([-60 -30 -40 -10])
   f = find(mask == 0);
   plot(lon(f),lat(f),'*k')
   title('Tensao de Cisalhamento do Vento [Pa]')
   xlabel('Longitude [W]')
   ylabel('Latitude [S]')
   hold off
   eval(['print -dpng figuras_lista3/eta_dia_',num2str(tempo/86400),'.png'])  
 
   % gravando outputs
   ETA = [ETA; eta2];
   U = [U; u2];
   V = [V; v2];
   
   end

   % calculando energia cinetica media
   ek = 0.5*rho*(mean(mean(u2.^2 + v2.^2)));
   EK = [EK ek];

end





