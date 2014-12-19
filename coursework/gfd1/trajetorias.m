%%%%%%%%% Trajetorias de ondas longas e curtas %%%%%%%%%

%%%% Ondas Longas %%%%
%%% Declaracao de Variaveis %%%


L=1000;
K=1/L;
H=8;
N0=2;
x=1;
t=1:1:2*pi/0.02;
w=0.02;

z0=0;
z4=-4;
z6=-6;
z8=-8;
z10=-10;


%% Equacoes %%
for i=1:length(t)
X(i)=-(N0/K*H)*sin((K*x)-(w*t(i)));

Z0(i)=((N0/H)*(z0+H))*cos((K*x)-(w*t(i)));
Z4(i)=((N0/H)*(z4+H))*cos((K*x)-(w*t(i)));
Z6(i)=((N0/H)*(z6+H))*cos((K*x)-(w*t(i)));
Z8(i)=((N0/H)*(z8+H))*cos((K*x)-(w*t(i)));
Z10(i)=((N0/H)*(z10+H))*cos((K*x)-(w*t(i)));
end

plot(X,Z0,'k',X,Z4,'b',X,Z6,'r',X,Z8,'g',X,Z10,'c')
subplot
title('Trajetoria da particula em agua rasa (Onda longa)')
xlabel('Eixo X')
ylabel('Eixo Z')
legend('Z0','Z4','Z8','Z10')