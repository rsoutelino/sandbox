% SOLUCAO DA EQUACAO DA ADVECCAO UNI-DIMENSIONAL LINEAR 
% POR METODO AVANCADO NO TEMPO E RETARDADO NO ESPACO. 
% SOLUCAO DA ADVECCAO DE UM SINAL SENOIDAL NO PONTO INICIAL DA GRADE   

% FORNECER NUMERO DE PONTOS DA GRADE, NUMERO DE PASSOS DE TEMPO, VELOCIDADE (m/s),
% ESPACAMENTO DE GRADE (m), PASSO DE TEMPO (s),
% AMPLITUDE DA OSCILACAO NA BORDA (m) E SEU PERIODO (s)
% E FREQUENCIA DE PLOTAGEM

% CONSTANTES DO MODELO
jmax=200;
mmax=360;
u=5;
h=10;
T=1;
amp=0.5;
per=100;
freqplo=5;

q=u*T/h;
omega=2*pi/per;
catu=zeros(jmax,1);
cren=zeros(mmax,1);

% CONDICOES INICIAIS (NA BORDA)
catu(1)=amp*sin(omega*2*T);
contplo=1;
amp2=amp*2;

% LOOP NO TEMPO
% CONDICOES DE CONTORNO
% FORMULA DE RECORRENCIA
% PLOTAGEM (PRESSIONE ENTER PARA EVOLUIR NO TEMPO)
% EVOLUCAO NO TEMPO DAS VARIAVEIS
for m=2:mmax
   tempo=m*T;
   cren(1)=amp*sin(omega*tempo);
   cren(2:jmax-1)=catu(2:jmax-1)-q*(catu(2:jmax-1)-catu(1:jmax-2));
   contplo=contplo+1;
   if(contplo==freqplo)
      contplo=0;
   	plot(cren)
   	axis([1 jmax -amp2 amp2]);
   	title(['Adveccao de sinal senoidal na borda (1a ordem) - tempo ',num2str(tempo),' segundos'])
    xlabel('PONTOS DE GRADE')
    ylabel('m')
   	pause
   end
   catu=cren;
end
