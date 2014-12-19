%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adveccao de uma concentracao na forma retangular                 
% com esquema QUICK                                                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%ler parametros do modelo
jmax=250;      %numero de pontos da grade em x
mmax=2400;     %numero de passos de tempo
T=12;         %passo de tempo, em segundos
h=10000;      %espacamento da grade, em m
u=50;          %velocidade da corrente em x, em m/s
ci=1;          %valor da concentracao inicial
cin=25;        %indice do ponto inicial da grade com ci
cfi=50;        %indice do ponto final da grade com ci
freqplot=100;  %frequencia de plotagem
conc=2;        %limite maximo do eixo da concentracao no grafico

%calculo das constantes do modelo
q=u*T/h;
catu=zeros(jmax,1);
cren=zeros(jmax,1);
a1=1/8;a2=-7/8;a3=3/8;a4=3/8;

%condicao inicial
catu(cin:cfi)=ci;
ccin=catu;

%loop no tempo
kplot=1;
%kfig=0;
for m=2:mmax
    tempo=m*T;
    kplot=kplot+1;
    % esquema avancado no tempo e retardado no espaco
    cren(3:jmax-1)=catu(3:jmax-1)-q*(catu(1:jmax-3)*a1+catu(2:jmax-2)*a2+catu(3:jmax-1)*a3+catu(4:jmax)*a4); 
        if(kplot==freqplot)
       	kplot=0;
         %kfig=kfig+1;
         %figure(kfig)
         plot(ccin,'r')
         hold
    		plot(cren)
    		axis([1 jmax -conc conc]);
         title(['Adveccao de sinal retangular (quick) - tempo ',num2str(tempo/3600),'horas'])
         xlabel('Pontos de grade')
         ylabel('Concentracao')
         pause
         hold
    end
    catu=cren;
end
    