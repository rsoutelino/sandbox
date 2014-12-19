%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Equacao de difusao unidimensional p/ sinal retangular%
%c(x)=0 para x < 225 e x > 275                        %
%c(x)=100 para 225 <=  x <= 275                       %
%solucao atraves de esquema explicito                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%parametros iniciais
jmax=500;  %numero de pontos de grade
mmax=3800; %numero de passos de tempo

%condicao inicial
j=1:jmax;
catu(j)=0;
catu(225:275)=100; 

%condicao de estabilidade
%q<=1/2  portanto    q=0.3125
Dh=10;   %em m2/s
T=2;     %em s
h=8;     %em m 
q=Dh*T/(h*h);

%esquema explicito avan t, centr em x
cren=catu;
freqplot=360;
kplot=359;
kfig=0;
for m=2:mmax
    tempo=(m-2)*T;
    kplot=kplot+1;
    cren(2:jmax-1)=catu(2:jmax-1)+ q*(catu(3:jmax)-2*catu(2:jmax-1)+catu(1:jmax-2));
    cren((cren<0))=0;   
    catu=cren;
        if(kplot==freqplot)
        kplot=0;
        XX=num2str(tempo);
        eval(['save -ascii fren_exp',XX,'.dat cren']);
        figure(1)
        plot(cren)
        grid on
        axis([1 jmax 0 100]);
        title(['tempo ',num2str(tempo/3600),' horas'])
        xlabel('PONTOS DE GRADE')
        grafico=['print -djpeg tempo_exp',XX];
        eval(grafico);
        end
end
