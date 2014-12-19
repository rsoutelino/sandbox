clear all;close all

x=[0.8133 2.7258 2.4117 6.1832 4.8636 6.1139 8.0668 8.0593 8.9044 9.1677]; 
t=[1 2 3 4 5 6 7 8 9 10]; 

n=length(x);

a=( n*sum(x.*t) - sum(x)*sum(t) ) / ( n*sum(x.^2) - sum(x)*sum(x) );

b=mean(t)-a*mean(x);

%%% montando a reta 

xx=-15:15;
tt=b*ones(size(xx)) + a*ones(size(xx)).*xx;

%%% plotando
figure
plot(x,t,'*',xx,tt,'r')
axis([-5 15 -5 15])
legend('Pontos Amostrais','Reta Ajustada')
title('Ajuste de Minimos Quadrados','fontsize',12,'fontweight','bold')