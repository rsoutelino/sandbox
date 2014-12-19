%%%%%%%%%%%%%%%%%%%%%%%%%
% Filtragem de Funcao F %                             
%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% Definindo a funcao
F=[3 5 2 7 3 9 5 10 8 11 9 14 11 16 13 18 14 20 18 25];

% Aplicacao de filtros de 3, 5 e 7 pontos
nudad=size(F,2);
Fil_3=F;
Fil_5=F;
Fil_7=F;
for j=2:nudad-1
Fil_3(j)=mean(F(j-1:j+1));
end
for j=3:nudad-2
Fil_5(j)=mean(F(j-2:j+2));
end
for j=4:nudad-3
Fil_7(j)=mean(F(j-3:j+3));
end

%Plotagem das funcoes original e filtradas
figure(1)
plot(F)
title(['Funcoes original (azul) e filtradas, com 3, 5 e 7 pontos'])
xlabel('Pontos de Grade j')
hold
plot(Fil_3,'r')
plot(Fil_5,'g')
plot(Fil_7,'k')


