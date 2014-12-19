clear all
clf

% Isto nós já sabemos:
T=15; % Periodo da onda (ja sabemos de antemao)
t=0:.1:100; % tempo
t=t'; % importante
N=size(t);

% Isto nós não sabemos (faz-de-conta)
A0=5;
A1=8.75e-3;
A2=0.24;
A3=pi/6;
E=0.15;

% Simulamos uns dados de modo que saberemos quão perto do sinal
% verdadeiro estamos 
a=A0+A1*t;
s=A2*sin((2*pi/T)*t+A3);
n=E*randn(N);
Y=a+s+n;
%------------------------------------------------------------------
% Faz de conta que não sabemos como foi feito Y

% Tentemos algo simples

% X é um modelo linear: y = [ b1 ; b2 ] * [ 1 t]

X=[ones(N) t];

B=(X'*X)\(X'*Y);

figure(1)
plot(t,Y,'.b',t,B(1)+B(2)*t,'r',t,a,'--m')
title('Falta algo na nossa receita')
legend('Dados','Modelo Linear','Tendencia Verdadeira')

disp('Acertamos os coeficientes?')

[B';[A0 A1]]

disp('Quase em cheio')

% Tentemos algo oscilatório

% X é um modelo seno + linear: y = [ b1 ; b2 ; b3 ; b4 ] * [ 1 t sin(wt) cos(wt)]

wt=(2*pi/T)*t;

X=[ones(N) t sin(wt) cos(wt)];

B=(X'*X)\(X'*Y);

yh=B(1)+B(2)*t+B(3)*sin(wt)+B(4)*cos(wt);

figure(2)
plot(t,Y,'.b',t,yh,'r',t,a+s,'--m')
title('Fascinante...')
legend('Dados','Modelo Oscilatorio','Oscilacao Verdadeira')

disp('Acertamos os coeficientes?')

[B';[A0 A1 A2 A3]]

disp('Aham...bem...aí depende...cof cof...')

% tente isto:
C=[B(1) B(2) sqrt(B(3)^2 + B(4)^2) atan2(B(4),B(3))];

[C;[A0 A1 A2 A3]]

disp('Agora melhorou. Explique o que aconteceu.')

disp('Faltou checar o erro do ajuste... ')
% Consulte a formula 3.12.8 do Emery & Thomson e veras a luz
Se=sqrt( (1/(N(1)-2)) * sum((Y-yh).^2) )



vY=std(Y)^2;
dv=std(Y-yh)^2;
pv=(vY-dv)/vY;

disp(['O modelo explicou ' num2str(round(pv*100)) '% da variancia dos dados'])

% Olhe agora como o E&T faz, na seção 3.13
Cxy=(1/(N(1)-1))*sum((Y-mean(Y)).*(yh-mean(yh)))
r=Cxy/(std(Y)*std(yh))
r2=r^2