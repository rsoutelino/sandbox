%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modelo Numerico hidrodinamico 3D linear %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%Sistema Internacional de Unidades - SI

%tamanho da grade e batimetria (chaves para pontos maritimos e internos)
jmax=11;
kmax=11;
lmax=10;
bat(1:kmax,1:jmax)=10.;
bat(1,:)=0.;
bat(kmax,:)=0.;
bat(:,1)=0.;
bat(:,jmax)=0.;
kmar=bat*0;
kmar(bat>0)=1;
X=[1:1:jmax];
Y=[1:1:kmax];
[X3,Y3,Z3] = meshgrid(1:jmax,1:kmax,-1:-1:-lmax);
    
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
dz=1;          % espacamento de grade em z
dt=60;         % passo de tempo
difz=0.001;    % coeficiente de difusao na vertical
freqplot=3;    % frequencia de plotagem
dens=1024;     % densidade media da agua do mar
latid=25*pi/180;            % latitude
fco=2*7.292E-5*sin(latid);  % parametro de Coriolis

dx2=dx*2;
dy2=dy*2;
difzdz=difz/dz;

% Condicoes iniciais de repouso (1 - valores atuais, 2 - renovados)
eta1=zeros(kmax,jmax);
u1=zeros(kmax,jmax,lmax);
v1=zeros(kmax,jmax,lmax);
eta2=zeros(kmax,jmax);
u2=zeros(kmax,jmax,lmax);
v2=zeros(kmax,jmax,lmax);
w2=zeros(kmax,jmax,lmax);

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
               forcx=0;
               forcy=0;
               for l=2:lmax-1
               forcx=(u1(k,j+1,l)-u1(k,j-1,l))*dz/dx2+forcx; 
               forcy=(v1(k+1,j,l)-v1(k-1,j,l))*dz/dy2+forcy;
               end
           eta2(k,j)=eta1(k,j)-dt*(forcx+forcy);
           end
        end
     end
  
  %Eq. do movimento em x
    for j=2:2:jmax-1
        for k=3:2:kmax-2
            if (kmar(k,j)*kmar(k,j+1)*kmar(k,j-1))>0
            for l=2:lmax-1
            dudzinf=(u1(k,j,l+1)-u1(k,j,l))*difzdz;
            if (l==2) dudzsup=-taux(k,j)/dens; end
            if (l~=2) dudzsup=(u1(k,j,l)-u1(k,j,l-1))*difzdz; end 
            vmed=(v1(k+1,j+1,l)+v1(k+1,j-1,l)+v1(k-1,j+1,l)+v1(k-1,j-1,l))/4;
            forc=fco*vmed-g.*(eta2(k,j+1)-eta2(k,j-1))./dx2...
               +(dudzinf-dudzsup)/dz;
            u2(k,j,l)=u1(k,j,l)+forc*dt;
            end
            end
         end
    end
    
    %Eq. do movimento em y
    for j=3:2:jmax-2
        for k=2:2:kmax-1
            if (kmar(k,j)*kmar(k+1,j))*kmar(k-1,j)>0
            for l=2:lmax-1
            dvdzinf=(v1(k,j,l+1)-v1(k,j,l))*difzdz;
            if (l==2) dvdzsup=-tauy(k,j)/dens; end
            if (l~=2) dvdzsup=(v1(k,j,l)-v1(k,j,l-1))*difzdz; end  
            umed=(u2(k+1,j+1,l)+u2(k+1,j-1,l)+u2(k-1,j+1,l)+u2(k-1,j-1,l))/4;
            forc=-fco*umed-g.*(eta2(k+1,j)-eta2(k-1,j))./dy2...
               +(dvdzinf-dvdzsup)/dz;
            v2(k,j,l)=v1(k,j,l)+forc*dt;
            end
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
              for l=2:lmax-1
              u2(k,j,l)=(u2(k-1,j,l)+u2(k+1,j,l))/2;
              v2(k,j,l)=(v2(k,j-1,l)+v2(k,j+1,l))/2;
              end
              eta2(k,j)=(eta2(k-1,j-1)+eta2(k-1,j+1)+eta2(k+1,j-1)+eta2(k+1,j+1))/4;
          end
        end
     end
   %definir u e v nos pontos tipo eta 
    for j=3:2:jmax-2
        for k=3:2:kmax-2
           if kmar(k,j)>0
              for l=2:lmax-1
              u2(k,j,l)=(u2(k,j+1,l)+u2(k,j-1,l))/2;
              v2(k,j,l)=(v2(k+1,j,l)+v2(k-1,j,l))/2;
              end
           end
        end
     end
    %definir v e eta nos pontos tipo u
    for j=2:2:jmax-1
        for k=3:2:kmax-2
            if kmar(k,j)>0
            for l=2:lmax-1
            v2(k,j,l)=(v2(k,j+1,l)+v2(k,j-1,l)+v2(k-1,j,l)+v2(k+1,j,l))/4;
            end
            eta2(k,j)=(eta2(k,j+1)+eta2(k,j-1))/2;
            end
         end
    end
    %definir u e eta nos pontos tipo v      
    for j=3:2:jmax-2
        for k=2:2:kmax-1
            if kmar(k,j)>0
            for l=2:lmax-1
            u2(k,j,l)=(u2(k,j+1,l)+u2(k,j-1,l)+u2(k-1,j,l)+u2(k+1,j,l))/4;
            end
            eta2(k,j)=(eta2(k-1,j)+eta2(k+1,j))/2;
            end
        end
    end

% Condicoes de contorno nao gradiente
for j=1:jmax
    eta2(1,j)=eta2(2,j)*kmar(1,j);
    eta2(kmax,j)=eta2(kmax-1,j)*kmar(kmax,j);
    for l=2:lmax-1
    u2(1,j,l)=u2(2,j,l)*kmar(1,j);
    v2(1,j,l)=v2(2,j,l)*kmar(1,j);
    u2(kmax,j,l)=u2(kmax-1,j,l)*kmar(kmax,j);
    v2(kmax,j,l)=v2(kmax-1,j,l)*kmar(kmax,j);
    end
end
for k=1:kmax
    eta2(k,1)=eta2(k,2)*kmar(k,1);
    eta2(k,jmax)=eta2(k,jmax-1)*kmar(k,jmax);
    for l=2:lmax-1
    u2(k,1,l)=u2(k,2,l)*kmar(k,1);
    v2(k,1,l)=v2(k,2,l)*kmar(k,1);
    u2(k,jmax,l)=u2(k,jmax-1,l)*kmar(k,jmax);
    v2(k,jmax,l)=v2(k,jmax-1,l)*kmar(k,jmax);
    end
end

% Plotando elevacoes e correntes
etamax=max(max(eta2));
etamin=min(min(eta2));
velo=sqrt(u2.^2+v2.^2);
velomax=max(max(max(velo(:,:,:))));

%figure (kfig)
figure(3)
subplot(1,2,1)
contourf(X,Y,eta2);
colorbar;
title(['Elev ',num2str(etamin),' a ',num2str(etamax),' m, ',num2str(tempo/60),' min'])
axis equal
axis([1 jmax 1 kmax])
xlabel('longitude')
ylabel('latitude')
subplot(1,2,2)
quiver3(X3,Y3,Z3,u2,v2,w2)
title(['Veloc - max ',num2str(velomax),' m/s'])
axis([1 jmax 1 kmax -lmax -1])
xlabel('longitude')
ylabel('latitude')
zlabel('profundidade')
% print -djpeg fig_resul
   end

pause(1)

end
   