%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                          %
%    Lista Pratica # 7 - Problemas 2 e 3                                   %
%    Prof. Ilson Carlos Almeida da Silveira                                %
%                                                                          %
%                                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all


close all

% carrega arquivo com matriz das velocidades de camada calculadas

load w4_2ly

% cria vetores individuais para x, y e as vel. de camada

x=mod1(:,1)';
y=mod1(:,2)';

% conversao para unidades de grau/seg
u1=0.9009e-7*mod1(:,3)';
v1=0.9009e-7*mod1(:,4)';
u2=0.9009e-7*mod1(:,5)';
v2=0.9009e-7*mod1(:,6)';
n=length(x);

% parametros de AO

corrlen = input('Enter correlation length: ')
err = input('Enter error: ')


% eliminando pontos sem observaoes na camada 2

f0=find(u2);
ff=length(f0);
x2=zeros(ff,1);
y2=zeros(ff,1);
u2c=zeros(ff,1);
v2c=zeros(ff,1);
m=0;

x1=mod1(:,1)';
y1=mod1(:,2)';

for j=1:n
 if abs(u2(j)) > 0,
 m=m+1;
 u2c(m)=u2(j);
 v2c(m)=v2(j);
 x2(m)=x1(j);
 y2(m)=y1(j);
 end
end
pack

% aplica metodos das imagens

x1m=zeros(2*length(x1),1); x1m(1:length(x1))=x1;
y1m=zeros(2*length(x1),1); y1m(1:length(x1))=y1;
u1m=zeros(2*length(x1),1); u1m(1:length(x1))=u1';
v1m=zeros(2*length(x1),1); v1m(1:length(x1))=v1';

xol1=length(x1);

x2m=zeros(2*length(x2),1); x2m(1:length(x2))=x2;
y2m=zeros(2*length(x2),1); y2m(1:length(x2))=y2;
u2m=zeros(2*length(x2),1); u2m(1:length(x2))=u2c;
v2m=zeros(2*length(x2),1); v2m(1:length(x2))=v2c;

xol2=length(x2);

for i=1:xol1
[x1m(xol1+i),y1m(xol1+i),u1m(xol1+i),v1m(xol1+i)]=mirror(x1(i),y1(i),u1(i),v1(i));
end
for i=1:xol2
[x2m(xol2+i),y2m(xol2+i),u2m(xol2+i),v2m(xol2+i)]=mirror(x2(i),y2(i),u2c(i),v2c(i));
end


%----------------------------------------
% Mapeamento para funcao de corrente
%       e vorticidade relativa 
%             Camada 1
%---------------------------------------

xx=[-53:0.25:-44];
yy=[0:0.25:9];

% aproximacao linear da isobata de 200

iy=-45.5-xx(1:31);iy=[0 iy 0];
ix=xx(1:31);ix=[-53 ix -53];

[xg,yg]=meshgrid(xx,yy);
xc=reshape(xg,37*37,1);
yc=reshape(yg,37*37,1);

[xg,yg]=meshgrid(xx,yy);
xc=reshape(xg,37*37,1);
yc=reshape(yg,37*37,1);

[zeta1,V1]=vortoa(xc,yc,x1m,y1m,u1m,v1m,corrlen,err);
     psi1= vectoa(xc,yc,x1m,y1m,u1m,v1m,corrlen,err,0);

 ZETA1 =(reshape(zeta1,length(yy),length(xx)));
 PSI1 =1.2321e10*(reshape(psi1,length(yy),length(xx)));

%----------------------------------------
% Mapeamento para funcao de corrente 
%       e vorticidade relativa
%             Camada 2
%---------------------------------------

[zeta2,V2]=vortoa(xc,yc,x2m,y2m,u2m,v2m,corrlen,err);
     psi2= vectoa(xc,yc,x2m,y2m,u2m,v2m,corrlen,err,0);

ZETA2 = (reshape(zeta2,length(yy),length(xx)));
PSI2 = 1.2321e10*(reshape(psi2,length(yy),length(xx)));

 %------------------------------------------------------
 %   Mapeamento de vorticidade potencial
 %------------------------------------------------------

 load param_naolinear

 lat0=5;
 gprime=epsin*9.81;
 H1=hn(1);H2=hn(2);
 omg=7.292e-05;
 f0=2*omg*sin(5*pi/180);

 Fr1=(f0*f0/(gprime*H1));         
 Fr2=(f0*f0/(gprime*H2));  

 STR1= Fr1*( PSI2 - PSI1); % vorticidade de estiramento (s^-1)
 STR2= Fr2*( PSI1 - PSI2);

 beta=2*omg*cos(5*pi/180)/57.6577;   %beta is in 1/degree/sec

  
 betay=beta*(yg-lat0*ones(37));     %  termo beta (s^-1)
 
% soma as tres partes para obter vp 

 qg1= (ZETA1+betay+STR1);        
 qg2= (ZETA2+betay+STR2);


%------------------------------------------------
% plota linhas de corrente
%------------------------------------------------

load /home/ilson/cenpes_2004/dados/wx4/wx_coast


figure(1)
cs=contourf(xx,yy,1e-5*PSI1,20);
%clabel(cs,'manual')
     caxis([-0.5 0.5]);
axis([-53 -44 0 9])
axis('square')
%set(gca, 'fontsize',16)
hold on
fill(ix,iy,[.9 .9 .9])
f=find(isnan(lonc));lf=length(f);f=[1;f];
for i=1:lf
   fill(lonc(f(i)+1:f(i+1)-1),latc(f(i)+1:f(i+1)-1),[.6 .6 .6]);
end
titt=['\bf Funcao de Corrente da Camada Superior - Westrax 4 (1e-05 m2/s)'];
tt=title(titt)
xlabel('\bf Longitude ')
ylabel('\bf Latitude ')
hold on
%fill(cx,cy,[0.6 0.6 0.6])
colorbar
%text(-52.5,1.7,'unidades: 1e05 m^{2} s^{-1})
%print -depsc2 psi_superior.eps 

figure(2)
cs=contourf(xx,yy,1e-5*PSI2,20);
%clabel(cs,'manual')
     caxis([-0.5 0.5]);
axis([-53 -44 0 9])
axis('square')
titt1=['\bf Funcao de Corrente da Camada Inferior - Westrax 4 (1e-05 m2/s)'];
tt=title(titt1);
xlabel('\bf Longitude')
ylabel('\bf Latitude')
hold on
fill(ix,iy,[.9 .9 .9])
f=find(isnan(lonc));lf=length(f);f=[1;f];
for i=1:lf
   fill(lonc(f(i)+1:f(i+1)-1),latc(f(i)+1:f(i+1)-1),[.6 .6 .6]);
end
hold on
%fill(cx,cy,[.6 .6 .6])
caxis([-.4 .4])
colorbar
%text(-52.5,1.7,'unidades: 1e05 m^{2} s^{-1})

%print -depsc2 psi_inferior.eps 

%------------------------------------------------
% Plota termo beta
%------------------------------------------------
figure(3)
cs=contourf(xx,yy,1e5*betay,20);
caxis([-1 1])
axis([-53 -44 0 9])
axis('square')
tit1=['\bf Variacao de Vorticidade Planetaria - Westrax 4'];
tt=title(tit1)
xlabel('\bf Longitude')
ylabel('\bf Latitude')
hold on
fill(ix,iy,[.9 .9 .9])
f=find(isnan(lonc));lf=length(f);f=[1;f];
for i=1:lf
   fill(lonc(f(i)+1:f(i+1)-1),latc(f(i)+1:f(i+1)-1),[.6 .6 .6]);
end
colorbar

%------------------------------------------------
% plota vorticidade relativa
%------------------------------------------------

figure(4)
cs=contourf(xx,yy,1e5*ZETA1,20);
%clabel(cs,'manual')
caxis([-1 1])
axis([-53 -44 0 9])
axis('square')
tit1=['\bf Vorticidade Relativa da Camada Superior - Westrax 4'];
tt=title(tit1)
xlabel('\bf Longitude')
ylabel('\bf Latitude')
hold on
fill(ix,iy,[.9 .9 .9])
f=find(isnan(lonc));lf=length(f);f=[1;f];
for i=1:lf
   fill(lonc(f(i)+1:f(i+1)-1),latc(f(i)+1:f(i+1)-1),[.6 .6 .6]);
end
colorbar
%text(-52.5,1.7,'unidades: 1e-05 s^{-1}')
%print -depsc2 vort_rel_sup.eps 
 
figure(5)
cs=contourf(xx,yy,1e5*ZETA2,16);
caxis([-1 1])
axis([-53 -44 0 9])
axis('square')
tit1=['\bf Vorticidade Relativa da Camada Inferior - Westrax 4'];
tt=title( tit1)
xlabel('\bf Longitude')
ylabel('\bf Latitude')
hold on
fill(ix,iy,[.9 .9 .9])
f=find(isnan(lonc));lf=length(f);f=[1;f];
for i=1:lf
   fill(lonc(f(i)+1:f(i+1)-1),latc(f(i)+1:f(i+1)-1),[.6 .6 .6]);
end
hold on
%fill(cx,cy,[.6 .6 .6])
caxis([-.6 .5])
colorbar
%text(-52.5,1.7,'unidades: 1e-05 s^{-1}')

%print -depsc2 vort_rel_inf.eps 


%------------------------------------------------
% plota vorticidade de estiramento
%------------------------------------------------

figure(6)
cs=contourf(xx,yy,1e5*STR1,16);
%clabel(cs,'manual')
caxis([-1 1])
axis([-53 -44 0 9])
axis('square')
tit1=['\bf Vorticidade de Estiramento da Camada Superior - Westrax 4'];
tt=title( tit1)
xlabel('\bf Longitude')
ylabel('\bf Latitude')
hold on
fill(ix,iy,[.9 .9 .9])
f=find(isnan(lonc));lf=length(f);f=[1;f];
for i=1:lf
   fill(lonc(f(i)+1:f(i+1)-1),latc(f(i)+1:f(i+1)-1),[.6 .6 .6]);
end
hold on
%fill(cx,cy,[0.6 0.6 0.6])
colorbar
%text(-52.5,1.7,'unidades: 1e-05 s^{-1}')

%print -depsc2 vort_est_sup.eps 

figure(7)
cs=contourf(xx,yy,1e5*STR2,16);
%clabel(cs,'manual')
caxis([-1 1])
axis([-53 -44 0 9])
axis('square')
tit1=['\bf Vorticidade de Estiramento da Camada Inferior - Westrax 4'];
tt=title( tit1)
xlabel('\bf Longitude')
ylabel('\bf Latitude')
hold on
fill(ix,iy,[.9 .9 .9])
f=find(isnan(lonc));lf=length(f);f=[1;f];
for i=1:lf
   fill(lonc(f(i)+1:f(i+1)-1),latc(f(i)+1:f(i+1)-1),[.6 .6 .6]);
end
%fill(cx,cy,[0.6 0.6 0.6])
colorbar
%text(-52.5,1.7,'unidades: 1e-05 s^{-1}')

%print -depsc2 vort_est_inf.eps 

%------------------------------------------------
% plota vorticidade potencial
%------------------------------------------------

figure(8)
[cs,h]=contourf(xx,yy,1e5*qg1,20);
caxis([-1 1])
set(h,'linestyle','--')
%clabel(cs,'manual')
axis([-53 -44 0 9])
axis('square')
tit1=['\bf Vorticidade Potencial da Camada Superior - Westrax 4'];
tt=title( tit1)
xlabel('\bf Longitude')
ylabel('\bf Latitude')
hold on
fill(ix,iy,[.9 .9 .9])
f=find(isnan(lonc));lf=length(f);f=[1;f];
for i=1:lf
   fill(lonc(f(i)+1:f(i+1)-1),latc(f(i)+1:f(i+1)-1),[.6 .6 .6]);
end
%fill(cx,cy,[0.6 0.6 0.6])
colorbar
%text(-52.5,1.7,'unidades: 1e-05 s^{-1}')
%print -depsc2 vort_pot_sup.eps 


figure(9)
[cs,h]=contourf(xx,yy,1e5*qg2,20);
set(h,'linestyle','--')
caxis([-1 1])
axis([-53 -44 0 9])
axis('square')
tit1=['\bf Vorticidade Potencial da Camada Inferior - Westrax 4'];
tt=title( tit1)
xlabel('\bf Longitude')
ylabel('\bf Latitude')
hold on
fill(ix,iy,[.9 .9 .9])
f=find(isnan(lonc));lf=length(f);f=[1;f];
for i=1:lf
   fill(lonc(f(i)+1:f(i+1)-1),latc(f(i)+1:f(i+1)-1),[.6 .6 .6]);
end
%fill(cx,cy,[0.6 0.6 0.6])
colorbar
%text(-52.5,1.7,'unidades: 1e-05 s^{-1}')
%print -depsc2 vort_pot_inf.eps 

%------------------------------------------------
% plota superposicao dos campos
% de  vorticidade potencial e
% funcao de corrente
%------------------------------------------------

figure(10)
[cs,H]=contour(xx,yy,1e-5*PSI1,20,'--b');
set(H(:),'linewidth',2)
hold on
[cs1,H]=contour(xx,yy,1e5*qg1,20,'k')
set(H(:),'linewidth',2)
axis([-53 -44 0 9])
axis('square')
hold on
titt=['\bf Funcao de Corrente (Azul) e Vorticidade Potencial da Camada Superior (Preto) - Westrax 4'];
tt=title(titt)
xlabel('\bf Longitude')
ylabel('\bf Latitude')
hold on
fill(ix,iy,[.9 .9 .9])
f=find(isnan(lonc));lf=length(f);f=[1;f];
for i=1:lf
   fill(lonc(f(i)+1:f(i+1)-1),latc(f(i)+1:f(i+1)-1),[.6 .6 .6]);
end
%fill(cx,cy,[0.6 0.6 0.6])
%gtext('unidades (azul): \times 10^{5} m^{2} s^{-1}','fontsize',12)
%gtext('unidades (preto): \times 10^{-5} s^{-1}','fontsize',12)
%print -depsc2 corr_vort_sup.eps 


figure(11)
[cs,H]=contour(xx,yy,1e-5*PSI2,20,'b--');set(H(:),'linewidth',2);
%clabel(cs,'manual','color','b')
hold on
[cs1,H]=contour(xx,yy,1e5*qg2,20,'k');set(H(:),'linewidth',2);
%clabel(cs1,'manual','color','k')
axis([-53 -44 0 9])
axis('square')
titt1=['\bf Função de Corrente (azul) e Vorticidade Potencial da Camada Inferior (Preto) - Westrax 4'];
tt=title(titt1)
xlabel('\bf Longitude ')
ylabel('\bf Latitude')
hold on
fill(ix,iy,[.9 .9 .9])
f=find(isnan(lonc));lf=length(f);f=[1;f];
for i=1:lf
   fill(lonc(f(i)+1:f(i+1)-1),latc(f(i)+1:f(i+1)-1),[.6 .6 .6]);
end
%fill(cx,cy,[0.6 0.6 0.6])
%gtext('unidades (azul): \times 10^{5} m^{2} s^{-1}','fontsize',12)
%gtext('unidades (preto): \times 10^{-5} s^{-1}','fontsize',12)
%print -depsc2 corr_vort_inf.eps 





