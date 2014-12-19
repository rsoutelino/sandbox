% SOLUCAO DO SISTEMA DE EQUACOES HIDRODINAMICAS 1D LINEAR 
% COM ESQUEMA DE PRIMEIRA ORDEM NO TEMPO E NO ESPACO. 
% CONDICOES DE CONTORNO: ELEVACAO SENOIDAL NO PONTO INICIAL DA GRADE   

% CONFIGURACAO DA GRADE
% E U E U E U E U E U E U
% + - + - + - + - + - + -
% j-1 (j) j+1 j+2
clear all
close all
% FORNECIMENTO DE PARAMETROS DA GRADE E DO PROCESSAMENTO (SI)
jmax=50;     % numero de pontos da grade
mmax=600;    % numero de passos de tempo
h=100;       % espacamento de grade
T=3;         % passo de tempo
prof=3;      % profundidade media
amp=0.5;     % amplitude da elevacao no inicio da grade
per=120;     % periodo da elevacao no inicio da grade
freqplot=10; % frequencia de plotagem
% DIMENSIONAMENTO DAS MATRIZES DE ELEVACAO E CORRENTES
eatu=zeros(jmax,1);
eren=zeros(jmax,1);
uatu=zeros(jmax,1);
uren=zeros(jmax,1);
% CONSTANTES DO MODELO	
g=9.8;
omega=2*pi/per;
ragh=sqrt(g/prof);
qu=g*T/h;
qe=T*prof/h;
amp2=2*amp;
cor2=2;
% CONDICOES INICIAIS (NO INICIO DA GRADE) 
eatu(1)=amp*sin(omega*T);
% LOOP NO TEMPO
kplot=1;
%kfig=0;
for m=2:mmax
   tempo=m*T;
   kplot=kplot+1;
   arg=omega*tempo;
% LOOPS NO ESPACO 
    uren(1:jmax-1)=uatu(1:jmax-1)-qu*(eatu(2:jmax)-eatu(1:jmax-1));
    eren(2:jmax-1)=eatu(2:jmax-1)-qe*(uren(2:jmax-1)-uren(1:jmax-2));
% CONDICOES DE CONTORNO DE CORRENTE
    uren(jmax)=2*uren(jmax-1)-uren(jmax-2);    
% CONDICOES DE CONTORNO DE ELEVACAO
    eren(1)=amp*sin(arg); 
    eren(jmax)=uren(jmax)/ragh;
% PLOTAGEM DE RESULTADOS DO MODELO 
         if(kplot==freqplot)
       	 kplot=0;
         %kfig=kfig+1;
         %figure(kfig)
         subplot(2,1,1)
         plot(eren)
         axis([1 jmax -amp2 amp2]);
         grid on
         title(['Elevação senoidal na borda (1a. ordem) - tempo = ',num2str(tempo/60),' minutos'])
         ylabel('Elevacao (m)')
         subplot(2,1,2)
         plot(uren)
         axis([1 jmax -cor2 cor2]);
         grid on
         xlabel('PONTOS DE GRADE')
         ylabel('Corrente (m/s)')
         pause
    	 end
% EVOLUCAO NO TEMPO DAS VARIAVEIS 
	eatu=eren;
	uatu=uren;
end