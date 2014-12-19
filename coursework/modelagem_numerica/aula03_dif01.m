% SOLUCAO DA EQUACAO DA DIFUSAO UNI-DIMENSIONAL LINEAR 
% PELO METODO AVANCADO NO TEMPO E CENTRADO NO ESPACO. 
% SOLUCAO DA DIFUSAO DE UMA CONCENTRACAO INICIAL    

% FORNECER NUMERO DE PONTOS DA GRADE, NUMERO DE PASSOS DE TEMPO,
% COEFICIENTE DE DIFUSAO (m^2/s),
% ESPACAMENTO DE GRADE (m), PASSO DE TEMPO (s)
% E FREQUENCIA DE PLOTAGEM

% CONSTANTES DO MODELO
jmax=25;
mmax=200;
d=10;
h=100;
T=10;
freqplo=5;

q=d*T/h/h;

catu=zeros(jmax,1);
catu(12)=99;
ccin=catu;
cren=zeros(mmax,1);
contplo=1;

% LOOP NO TEMPO
% FORMULA DE RECORRENCIA
% PLOTAGEM (PRESSIONE ENTER PARA EVOLUIR NO TEMPO)
% EVOLUCAO NO TEMPO DAS VARIAVEIS
for m=2:mmax
   tempo=m*T;
   cren(2:jmax-1)=catu(2:jmax-1)+q*(catu(3:jmax)-2*catu(2:jmax-1)+catu(1:jmax-2));
   for j=1:jmax
      if(cren(j)<0)
         cren(j)=0;
      end
   end
   contplo=contplo+1;
   if(contplo==freqplo)
      contplo=0;
      plot(ccin,'r')
      hold
   	  plot(cren)
   	  axis([1 jmax -1 100]);
   	  title(['Difusão (avan tempo, centr esp) - tempo ',num2str(tempo),' segundos'])
      pause
      hold
   end
   catu=cren;
end
