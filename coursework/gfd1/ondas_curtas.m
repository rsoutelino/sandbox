%Declarar variaveis e ctes.
L=200; % Usar 100
N0=0.5;
K=2*pi/L;
g=9.8;
w=sqrt(g*K);
T=2*pi/w;

%%%%%%%%% Trajetorias de ondas curtas %%%%%%%%%

x=1;
t=1:0.1:max(T)+2;
z0=0;
z2=-2;
z4=-4;
z6=-6;
z8=-8;
z10=-10;

%% Equacoes Curtas%%



for i=1:length(t)
N(i)=N0*cos((K*x)-(2*pi*w*t(i)));    
X0(i)=N0*sin(K*x-w*t(i))*(cosh(K*z0)+sinh(K*z0));
X2(i)=N0*sin(K*x-w*t(i))*(cosh(K*z2)+sinh(K*z2));
X4(i)=N0*sin(K*x-w*t(i))*(cosh(K*z4)+sinh(K*z4));
X6(i)=N0*sin(K*x-w*t(i))*(cosh(K*z6)+sinh(K*z6));
X8(i)=N0*sin(K*x-w*t(i))*(cosh(K*z8)+sinh(K*z8));
X10(i)=N0*sin(K*x-w*t(i))*(cosh(K*z10)+sinh(K*z10));

Z0(i)=N0*cos(K*x-w*t(i))*(sinh(K*z0)+cosh(K*z0));
Z2(i)=N0*cos(K*x-w*t(i))*(sinh(K*z2)+cosh(K*z2))-2;
Z4(i)=N0*cos(K*x-w*t(i))*(sinh(K*z4)+cosh(K*z4))-4;
Z6(i)=N0*cos(K*x-w*t(i))*(sinh(K*z6)+cosh(K*z6))-6;
Z8(i)=N0*cos(K*x-w*t(i))*(sinh(K*z8)+cosh(K*z8))-8;
Z10(i)=N0*cos(K*x-w*t(i))*(sinh(K*z10)+cosh(K*z10))-10;
end

plot(X0,Z0,'k',X2,Z2,'k',X4,Z4,'k',X6,Z6,'k',X8,Z8,'k',X10,Z10,'k')
hold
x=-8:length(N)-9; % Eixo para Eta
plot(x,N,'k:') % Plota elevacao
legend('Eta')
title('Trajetoria da particula em agua profunda (Onda curta)')
axis([-8 8 -11 3])
print -dtiff ondas_curtas.tiff