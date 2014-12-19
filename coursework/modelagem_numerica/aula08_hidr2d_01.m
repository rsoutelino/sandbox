%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modelo Numerico hidrodinamico 2D linear (decaimento explicito) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%Configuracao da grade (pontos tipo eta, u, v)
%NORTE (k)
%7 EUEUEUE
%6 V V V V
%5 EUEUEUE
%4 V V V V
%3 EUEUEUE
%2 V V V V
%1 EUEUEUE
%  1234567 LESTE (j)

%Sistema Internacional de unidades (SI)

%tamanho da grade e batimetria (chaves para pontos maritimos)
jmax=11;
kmax=11;
bat(1:kmax,1:jmax)=10.;
bat(1,:)=0.;
bat(kmax,:)=0.;
bat(:,1)=0.;
bat(:,jmax)=0.;
kmar=bat*0;
kmar(bat>0)=1;
X=[1:1:jmax];
Y=[1:1:kmax];
    
%plotando batimetria
figure(1)
contourf(X,Y,bat);
colorbar;
title(['Batimetria da regiao modelada (m)'])
axis([1 jmax 1 kmax])
xlabel('EW')
ylabel('NS')
%print -djpeg fig_batim

% Constantes iniciais
g=9.8;         % aceleracao da gravidade
mmax=1000;     % niveis de tempo
dx=1000;       % espacamento da grade em x
dy=1000;       % espacamento da grade em y
dt=30;         % passo de tempo
rfric=0.02;    % coeficiente de friccao no fundo
freqplot=20;   % frequencia de plotagem
dens=1024;     % densidade media da agua do mar
latid=25*pi/180;            % latitude (em graus, para rad)
fco=2*7.292E-5*sin(latid);  % parametro de Coriolis

dx2=dx*2;
dy2=dy*2;

dens2=dens*2;
denbat=dens.*bat;

% Condicoes iniciais de repouso (1 - valores atuais, 2 - renovados)
eta1=zeros(kmax,jmax);
u1=zeros(kmax,jmax);
v1=zeros(kmax,jmax);
eta2=zeros(kmax,jmax);
u2=zeros(kmax,jmax);
v2=zeros(kmax,jmax);

%Condicoes de vento e calculo das tensoes de cisalhamento na superficie
dens_ar=1.025;
fric=2.6*1E-3;

uwind(1:kmax,1:jmax)=7.0711;
vwind(1:kmax,1:jmax)=7.0711;
wwind=sqrt(uwind.^2+vwind.^2);
taux=fric*dens_ar.*uwind.*wwind;
tauy=fric*dens_ar.*vwind.*wwind;

vento=sqrt(uwind.^2+vwind.^2);
ventoma=max(vento);
ventomax=max(ventoma);

%plotando os ventos
figure(2)
quiver(uwind,vwind)
title(['Vento - intensidade maxima ',num2str(ventomax),' m/s'])
axis([1 jmax 1 kmax])
xlabel('EW')
ylabel('NS')
%print -djpeg fig_vento

%Loop no tempo
kplot=0;
kfig=0;
for m=1:mmax
   tempo=m*dt;
   kplot=kplot+1;
    
    %Eq da Continuidade
    for j=3:2:jmax-2
        for k=3:2:kmax-2
           if kmar(k,j)>0
           forcx=(bat(k,j+1).*u1(k,j+1)-bat(k,j-1).*u1(k,j-1))/dx2; 
           forcy=(bat(k+1,j).*v1(k+1,j)-bat(k-1,j).*v1(k-1,j))/dy2;
           eta2(k,j)=eta1(k,j)-dt*(forcx+forcy);
           end
        end
     end
  
  %Eq. do movimento em x
    for j=2:2:jmax-1
        for k=3:2:kmax-2
            if (kmar(k,j)*kmar(k,j+1)*kmar(k,j-1))>0
            vmed=(v1(k+1,j+1)+v1(k+1,j-1)+v1(k-1,j+1)+v1(k-1,j-1))/4;
            forc=fco*vmed-g.*(eta2(k,j+1)-eta2(k,j-1))./dx2...
               +taux(k,j)./denbat(k,j)-rfric*u1(k,j);
            u2(k,j)=u1(k,j)+forc*dt;
            end
         end
    end
    
    %Eq. do movimento em y
    for j=3:2:jmax-2
        for k=2:2:kmax-1
            if (kmar(k,j)*kmar(k+1,j))*kmar(k-1,j)>0
            umed=(u2(k+1,j+1)+u2(k+1,j-1)+u2(k-1,j+1)+u2(k-1,j-1))/4;
            forc=-fco*umed-g.*(eta2(k+1,j)-eta2(k-1,j))./dy2...
               +tauy(k,j)./denbat(k,j)-rfric*v1(k,j);
            v2(k,j)=v1(k,j)+forc*dt;
            end
        end
    end
    
% Renovando as variaveis no tempo
eta1=eta2;
u1=u2;
v1=v2;
     
% Plotagem de resultados     
   if(kplot==freqplot)
   kplot=0;
   kfig=kfig+1;
   hold
   %definir u,v,eta nos pontos "vazios"
    for j=2:2:jmax-1
        for k=2:2:kmax-1
           if kmar(k,j)>0
              u2(k,j)=(u2(k-1,j)+u2(k+1,j))/2;
              v2(k,j)=(v2(k,j-1)+v2(k,j+1))/2;
              eta2(k,j)=(eta2(k-1,j-1)+eta2(k-1,j+1)+eta2(k+1,j-1)+eta2(k+1,j+1))/4;
           end
        end
     end
   %definir u e v nos pontos tipo eta 
    for j=3:2:jmax-2
        for k=3:2:kmax-2
           if kmar(k,j)>0
              u2(k,j)=(u2(k,j+1)+u2(k,j-1))/2;
              v2(k,j)=(v2(k+1,j)+v2(k-1,j))/2;
           end
        end
     end
    %definir v e eta nos pontos tipo u
    for j=2:2:jmax-1
        for k=3:2:kmax-2
            if kmar(k,j)>0
            v2(k,j)=(v2(k,j+1)+v2(k,j-1)+v2(k-1,j)+v2(k+1,j))/4;
            eta2(k,j)=(eta2(k,j+1)+eta2(k,j-1))/2;
            end
         end
    end
    %definir u e eta nos pontos tipo v      
    for j=3:2:jmax-2
        for k=2:2:kmax-1
            if kmar(k,j)>0
            u2(k,j)=(u2(k,j+1)+u2(k,j-1)+u2(k-1,j)+u2(k+1,j))/4;
            eta2(k,j)=(eta2(k-1,j)+eta2(k+1,j))/2;
            end
        end
    end

% Condicoes de contorno nao gradiente
for j=1:jmax
    eta2(1,j)=eta2(2,j)*kmar(1,j);
    eta2(kmax,j)=eta2(kmax-1,j)*kmar(kmax,j);
    u2(1,j)=u2(2,j)*kmar(1,j);
    v2(1,j)=v2(2,j)*kmar(1,j);
    u2(kmax,j)=u2(kmax-1,j)*kmar(kmax,j);
    v2(kmax,j)=v2(kmax-1,j)*kmar(kmax,j);
end
for k=1:kmax
    eta2(k,1)=eta2(k,2)*kmar(k,1);
    eta2(k,jmax)=eta2(k,jmax-1)*kmar(k,jmax);
    u2(k,1)=u2(k,2)*kmar(k,1);
    v2(k,1)=v2(k,2)*kmar(k,1);
    u2(k,jmax)=u2(k,jmax-1)*kmar(k,jmax);
    v2(k,jmax)=v2(k,jmax-1)*kmar(k,jmax);
end
   
% Plotando elevacoes e correntes
velo=sqrt(u2.^2+v2.^2);
veloma=max(velo);
velomax=max(veloma);
etama=max(eta2(:,:));
etami=min(eta2(:,:));
etamax=max(etama);
etamin=min(etami);
%figure (kfig)
figure(3)
subplot(2,1,1)
contourf(X,Y,eta2);
colorbar;
title(['Elevacoes (m) - tempo ',num2str(tempo/60),...
      ' minutos. Limites ',num2str(etamin),' a ',num2str(etamax),' m'])
axis([1 jmax 1 kmax])
subplot(2,1,2)
quiver(u2,v2);
title(['Correntes (m/s) - intensidade maxima ',num2str(velomax),' m/s'])
axis([1 jmax 1 kmax])
xlabel('longitude')
ylabel('latitude')
% print -djpeg fig_resul
   end

pause(1)

end
   