clear all;close all

load ../analises/mat/lon_lat_psi_OE2.mat
load ../proc/common/etopo2_leste.mat

% o programa é montado para fazer a cada radial

% carregando a isobata que sera usada como espelho

[c,h] = contour(xb,yb,zb,[-100 -100]);
xiso = get(h(1),'xdata'); xiso = xiso(find(not(isnan(xiso))));
yiso = get(h(1),'ydata'); yiso = yiso(find(not(isnan(yiso))));
xiso = weim(21,'hann',xiso);
yiso = weim(21,'hann',yiso);
close(1)

dd = sw_dist(lat,lon,'km');
fd = find(dd > 200);

c = 0; lonest=[]; latest=[]; lonim=[]; latim=[]; psiest=[]; psiIM=[];
for k = 1:length(fd)+1
   
    % separando a radial
    if k == length(fd)+1;
        xrad = lon(fd(end)+1:end);
        yrad = lat(fd(end)+1:end);
        psirad = psig(fd(end)+1:end);
    else
        xrad = lon((c+1):fd(k));
        yrad = lat((c+1):fd(k));
        psirad = psig((c+1):fd(k));
        c = fd(k);
    end
    
    % verificando se algum ponto da radial esta alem da isobata
    f = near(yiso,yrad(1),1);
    f2 = find(xrad >= xiso(f));
    xrad = xrad(f2); yrad = yrad(f2); psirad = psirad(f2);

    % colecionando pontos da isobata proximos a radial, 
    % para estabelecer uma reta
    xr = xiso(f-2:f+2); yr = yiso(f-2:f+2);
    p = polyfit(xr,yr,1); m1 = p(1); h1 = p(2); % eq da reta: y = m1.x + h1
    xr1 = min(lon):0.1:max(lon);
    yr1 = m1.*xr1 + h1;
    % espelhando ponto a ponto
    xim=[]; yim=[]; psiim=[];
    for q = 1:length(xrad)
        m2 = -1/m1; % montando reta perpendicular: y = m2.x + h2
        h2 = xrad(q)/m1 + yrad(q);
        xi = (h2-h1)/(m1-m2); yi = m1*((h2-h1)/(m1-m2)) + h1; % intersecao entre as retas
        u = xrad(q) - xi; v = yrad(q) - yi;
        x2 = xrad(q) - 2*u; y2 = yrad(q) - 2*v;

        % tirando os pontos imagem que ficaram dentro do dominio real
        fn = near(yiso,y2,1);
        if x2 < xiso(fn) 
            xim = [xim x2]; yim = [yim y2]; psiim = [psiim -psirad(q)];
        end


figure(1)
plot(xiso,yiso,'k');hold on
plot(lon,lat,'*b');hold on
plot(xrad,yrad,'r*');hold on
axis([-42 -33 -21 -5])
axis('equal')
plot(x2,y2,'*k')

       

    end 

   lonest = [lonest xrad]; latest = [latest yrad];
   lonim = [lonim xim]; latim = [latim yim];
   psiest = [psiest psirad]; psiIM = [psiIM psiim];

end

figure(2)
plot(xiso,yiso,'k') ; hold on
plot(lonest,latest,'k*')
plot(lonim,latim,'r*')
