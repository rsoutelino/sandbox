%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Equacao de difusao unidimensional p/ sinal retangular%
%c(x)=0 para x < 225 e x > 275                        %
%c(x)=100 para 225 <=  x <= 275                       %
%solucao atraves de esquema pseudo-implicito          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%parametros iniciais
jmax=500;  %numero de pontos de grade
mmax=1900; %numero de passos de tempo

%condicao inicial
j=1:jmax;
catu(j)=0;
catu(225:275)=100; 

Dh=10;   %em m2/s
T=4;     %em s
h=8;     %em m 
q=Dh*T/(h*h);

%esquema pseudo-implicito cent t, centr em x
cant=catu;
cren=catu;
freqplot=180;
kplot=179;
kcig=0;
for m=3:mmax
    tempo=(m-3)*T;
    kplot=kplot+1;
    cren(2:jmax-1)=(cant(2:jmax-1)+ 2*q*(catu(3:jmax)-cant(2:jmax-1)+catu(1:jmax-2)))/(1+2*q);
    cren((cren<0))=0;
    cant=catu;
    catu=cren;
    if(kplot==freqplot)
    kplot=0;
    XX=num2str(tempo);
    eval(['save -ascii cren_pseu',XX,'.dat cren']);
    figure(1)
    plot(cren)
    grid on
    axis([1 jmax 0 100]);
    title(['tempo ',num2str(tempo/3600),' horas'])
    xlabel('PONTOS DE GRADE')
    grafico=['print -djpeg tempo_pseu',XX];
    eval(grafico);
    end
end
