%Declarar variaveis e ctes.
L=200; 
N0=10; % 0.5
K=2*pi/L;
H=500;
H(1)=0.01;
g=9.8;

w=sqrt(g*K^2*H);
T=2*pi/w;

%%%%%%%%% Trajetorias de ondas longas %%%%%%%%%

x=0.0;
t=1:T+5;
H=500;
z0=0.01;
z100=-100;
z200=-200;
z300=-300;
z400=-400;
z500=-500;

%% Equacoes Longas%%
for i=1:length(t)
N(i)=N0*cos((K*x)-(2*pi*w*t(i)));
X0(i)=(-(N0/K*H)*sin((K*x)-(w*t(i))));
X100(i)=(-(N0/K*H)*sin((K*x)-(w*t(i))));
X200(i)=(-(N0/K*H)*sin((K*x)-(w*t(i))));
X300(i)=(-(N0/K*H)*sin((K*x)-(w*t(i))));
X400(i)=(-(N0/K*H)*sin((K*x)-(w*t(i))));
X500(i)=(-(N0/K*H)*sin((K*x)-(w*t(i))));

Z0(i)=((N0/H)*(z0+H))*cos((K*x)-(w*t(i)));
Z100(i)=((N0/H)*(z100+H))*cos((K*x)-(w*t(i)))-100;
Z200(i)=((N0/H)*(z300+H))*cos((K*x)-(w*t(i)))-200;
Z300(i)=((N0/H)*(z500+H))*cos((K*x)-(w*t(i)))-300;
Z400(i)=((N0/H)*(z400+H))*cos((K*x)-(w*t(i)))-400;
Z500(i)=((N0/H)*(z500+H))*cos((K*x)-(w*t(i)))-500;
end


plot(X0,Z0,'k',X100,Z100,'k',X200,Z200,'k',X300,Z300,'k',X400,Z400,'k',X500,Z500,'k')
title('Trajetoria da particula em agua rasa (Onda longa)')
axis([-600000 600000 -600 100])
hold
x=-600000:440:600000;
%plot(x,N,'k:') % Plota elevacao
legend('Eta')
%print -dtiff ondas_longas.tiff