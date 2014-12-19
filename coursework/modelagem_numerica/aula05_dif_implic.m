%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Equacao de difusao unidimensional p/ sinal retangular%
%c(x)=0 para x < 225 e x > 275                        %
%c(x)=100 para 225 <=  x <= 275                       %
%solucao atraves de esquema implicito                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%parametros iniciais
jmax=500;  %numero de pontos de grade
mmax=792;  %numero de passos de tempo

%condicao inicial
j=1:jmax;
catu(j)=0;
catu(225:275)=100; 
cant=catu;
cren=catu;

Dh=10;   %em m2/s
T=10;    %em s
h=8;     %em m 
q=Dh*T/(h*h);

%esquema implicito cent t, centr em x
freqplot=72;
kplot=71;
kfig=0;
sj=zeros(1,jmax);
ej=zeros(1,jmax);
dj=zeros(1,jmax);

aj=Dh/(h*h);
bj=2*Dh/(h*h)+1/(2*T);
cj=Dh/(h*h);
    
for m=2:mmax
   tempo=(m-2)*T;
   kplot=kplot+1;
   
   dj=cant*1/(2*T);
   %varredura ascendente
   for j=2:jmax-1   
   sj(j)=cj/(bj-aj*sj(j-1));
   ej(j)=(dj(j)+aj*ej(j-1))/(bj-aj*sj(j-1));
   end
   %varredura descendente
   for j=jmax-1:-1:2    
   cren(j)=sj(j)*cren(j+1)+ej(j);
   end
   cren((cren<0))=0;
   
   cant=catu;
   catu=cren;
        
   if(kplot==freqplot)
 	kplot=0;
   XX=num2str(tempo);
   eval(['save -ascii cren_imp',XX,'.dat cren']);
   figure(1)
   plot(cren)
   grid on
   axis([1 jmax 0 100]);
   title(['tempo ',num2str(tempo/3600),' horas'])
   xlabel('PONTOS DE GRADE')
   grafico=['print -djpeg tempo_imp',XX];
   eval(grafico);
   end
 end
