% Programa para calculo de derivadas por diferencas finitas
clear all
close all

dx=pi/25;
x=0:dx:pi;
f=sin(x);

nudad=size(x,2);
deriv01=zeros(nudad,1);
deriv02=zeros(nudad,1);
deriv03=zeros(nudad,1);
deriv_seg=zeros(nudad,1);
deriv_analit=zeros(nudad,1);
deriv_seg_analit=zeros(nudad,1);

for j=1:nudad-1
deriv01(j)=(f(j+1)-f(j))/(x(j+1)-x(j));
end
for j=2:nudad
deriv02(j)=(f(j)-f(j-1))/(x(j)-x(j-1));
end
for j=2:nudad-1
deriv03(j)=(f(j+1)-f(j-1))/(x(j+1)-x(j-1));
end
for j=2:nudad-1
deriv_seg(j)=(f(j+1)-2*f(j)+f(j-1))/dx^2;
end

figure (1)
subplot(3,1,1)
plot(x,f)
axis([0 pi -1.5 1.5])
grid on
title('Funcao original')
ylabel('f')
subplot(3,1,2)
plot(x,deriv01,'r')
hold
plot(x,deriv02,'k')
plot(x,deriv03,'b')
axis([0 pi -1.5 1.5])
grid on
title('Derivadas: avan (verm), retard (preto), centr (azul)')
ylabel('df/dx')
subplot(3,1,3)
plot(x,deriv_seg)
axis([0 pi -1.5 1.5])
grid on
title('Derivada segunda (centrada)')
ylabel('d2f/dx2')
xlabel('x')

deriv_analit=cos(x);
deriv_analit=deriv_analit';
deriv_seg_analit=-sin(x);
deriv_seg_analit=deriv_seg_analit';

difer01=deriv01-deriv_analit;
difer02=deriv02-deriv_analit;
difer03=deriv03-deriv_analit;
difer_seg=deriv_seg-deriv_seg_analit;

max_min_dp_difer01=[max(difer01(2:nudad-1)) min(difer01(2:nudad-1)) std(difer01(2:nudad-1))]
max_min_dp_difer02=[max(difer02(2:nudad-1)) min(difer02(2:nudad-1)) std(difer02(2:nudad-1))]
max_min_dp_difer03=[max(difer03(2:nudad-1)) min(difer03(2:nudad-1)) std(difer03(2:nudad-1))]

figure (2)
subplot(3,1,1)
plot(x(2:nudad-1),difer01(2:nudad-1))
grid on
title('Erro da diferenca finita avancada')
ylabel('df/dx')
subplot(3,1,2)
plot(x(2:nudad-1),difer02(2:nudad-1))
grid on
title('Erro da diferenca finita retardada')
ylabel('df/dx')
subplot(3,1,3)
plot(x(2:nudad-1),difer03(2:nudad-1))
grid on
title('Erro da diferenca finita centrada')
ylabel('df/dx')
xlabel('x')

