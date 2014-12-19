%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adveccao de uma concentracao na forma retangular                 
% com esquema QUICKEST                                                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%ler parametros do modelo
jmax=250;        %numero de pontos da grade em x
mmax=2400;       %numero de passos de tempo
T=12;           %passo de tempo, em segundos
h=10000;        %espacamento da grade, em m
u=50;            %velocidade da corrente em x, em m/s
ci=1;            %valor da concentracao inicial
cin=25;          %indice do ponto inicial da grade com ci
cfi=50;          %indice do ponto final da grade com ci
freqplot=100;    %frequencia de plotagem
conc=2;          %limite maximo do eixo da concentracao no grafico

%calculo das constantes do modelo
q=u*T/h;
catu=zeros(jmax,1);
cren=zeros(jmax,1);
parc1=zeros(jmax,1);
parc2=zeros(jmax,1);
parc3=zeros(jmax,1);

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
    parc1(3:jmax-1)=(catu(4:jmax)-catu(2:jmax-2))/2;
    parc2(3:jmax-1)=(-2*catu(3:jmax-1)+catu(2:jmax-2)+catu(4:jmax))*q/2;
    parc3(3:jmax-1)=(catu(1:jmax-3)-3*catu(2:jmax-2)+3*catu(3:jmax-1)-catu(4:jmax))*(1-q*q)/6;
    cren(3:jmax-1)=catu(3:jmax-1)-q*(parc1(3:jmax-1)-parc2(3:jmax-1)+parc3(3:jmax-1));
    if(kplot==freqplot)
       	kplot=0;
         %kfig=kfig+1;
         %figure(kfig)
         plot(ccin,'r')
         hold
    		plot(cren)
    		axis([1 jmax -conc conc]);
         title(['Adveccao de sinal retangular (quickest) - tempo ',num2str(tempo/3600),'horas'])
         xlabel('Pontos de grade')
         ylabel('Concentracao')
         pause
         hold
    end
    catu=cren;
end