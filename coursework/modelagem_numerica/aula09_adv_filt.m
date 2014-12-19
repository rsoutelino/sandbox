%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adveccao de uma concentracao na forma retangular                 
% com esquema de segunda ordem                                                        
% (o modo computacional do esquema de segunda ordem
% eh filtrado)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%ler parametros do modelo
jmax=250;      %numero de pontos da grade em x
mmax=240;      %numero de passos de tempo
T=120;         %passo de tempo, em segundos
h=10000;       %espacamento da grade, em m
u=50;          %velocidade da corrente em x, em m/s
xin=1;         %valor da concentracao inicial
pin=25;        %indice do ponto inicial da grade com ci
pfi=50;        %indice do ponto final da grade com ci
freqplot=10;   %frequencia de plotagem
conc=2;        %limite maximo do eixo da concentracao no grafico
gama=0.5;      %coeficiente do filtro de 3 pontos

%calculo das constantes do modelo
q=u*T/h;
gama1=1-gama;
gama2=gama/2;
cant=zeros(jmax,1);
catu=zeros(jmax,1);
cren=zeros(jmax,1);
cfant=zeros(jmax,1);
cfatu=zeros(jmax,1);
cfren=zeros(jmax,1);
cfilt=zeros(jmax,1);

%condicao inicial
cant(pin:pfi)=xin;
cfant(pin:pfi)=xin;
ccin=cant;

%esquema avancado no tempo e centrado no espaco
%somente para o primeiro avanco no tempo
catu(2:jmax-1)=cant(2:jmax-1)-(q/2)*(cant(3:jmax)-cant(1:jmax-2)); 
cfatu(2:jmax-1)=cfant(2:jmax-1)-(q/2)*(cfant(3:jmax)-cfant(1:jmax-2)); 

%loop no tempo
kplot=2;
%kfig=0;
for m=3:mmax
    tempo=m*T;
    kplot=kplot+1;
    %esquema centrado no tempo e no espaco
	 cren(2:jmax-1)=cant(2:jmax-1)-q*(catu(3:jmax)-catu(1:jmax-2)); 
	 cfren(2:jmax-1)=cfant(2:jmax-1)-q*(cfatu(3:jmax)-cfatu(1:jmax-2)); 
     for j=2:jmax-1
         cfilt(j)=gama1*cfren(j)+gama2*(cfren(j-1)+cfren(j+1));
     end
     cfren=cfilt;
        if(kplot==freqplot)
       	kplot=0;
         %kfig=kfig+1;
         %figure(kfig)
         plot(ccin,'r')
         hold
         plot(cren,'g')
         plot(cfren)
         axis([1 jmax -conc conc]);
         title(['Adv de sinal retangular (esq. 2a. ordem) c/ filt - tempo ',num2str(tempo/3600),'horas'])
         xlabel('Pontos de grade')
         ylabel('Concentracao')
         pause
         hold
    end
    cant=catu;
    catu=cren;
    cfant=cfatu;
    cfatu=cfren;
end
    