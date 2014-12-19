%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filtragem de Funcao F 2D %                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% Definindo a funcao
jmax=8;
kmax=8;
j=1:jmax;
k=1:kmax;
v1=0.5;
v2=0.5;
F=[3 4 4 6 7 4 2 1; 6 7 9 9 6 5 5 4; 7 9 11 12 14 11 9 8;...
      8 12 13 18 14 8 6 6 ; 9 10 11 12 18 16 11 9; 6 9 11 13 14 12 9 6;...
      5 7 8 8 4 3 3 2; 2 4 3 6 5 4 2 1];

% Plotagem da funcao original
figure(1)
contourf(F);colorbar
axis([1 jmax 1 kmax]);
title(['Funcao Sem Filtragem'])
xlabel('Pontos de Grade j')
ylabel('Pontos de Grade k')

% Aplicacao do filtro de 9 pontos
Fil_9(2:jmax-1,2:kmax-1)=F(2:jmax-1,2:kmax-1)+0.5*v1*(1-v1)*(F(1:jmax-2,2:kmax-1)...
   +F(3:jmax,2:kmax-1)+F(2:jmax-1,1:kmax-2)+F(2:jmax-1,3:kmax)-4*F(2:jmax-1,2:kmax-1))...
   +0.25*v1*v1*(F(3:jmax,1:kmax-2)+F(1:jmax-2,1:kmax-2)...
   +F(1:jmax-2,3:kmax)+F(3:jmax,3:kmax)-4*F(2:jmax-1,2:kmax-1));

% Plotando os resultados do filtro de 9 pontos
figure(2)
contourf(Fil_9);colorbar
axis([1 jmax-1 1 kmax-1]);
title(['Funcao Com Filtragem 2D de 9 pontos'])
xlabel('Pontos de Grade j')
ylabel('Pontos de Grade k')

% Aplicacao do filtro de 5 pontos
Fil_5(2:jmax-1,2:kmax-1)=v2*F(2:jmax-1,2:kmax-1)+0.25*(1-v2)*(F(1:jmax-2,2:kmax-1)...
   +F(3:jmax,2:kmax-1)+F(2:jmax-1,1:kmax-2)+F(2:jmax-1,3:kmax));

% Plotando os resultados do filtro de 5 pontos
figure(3)
contourf(Fil_5);colorbar
axis([1 jmax-1 1 kmax-1]);
title(['Funcao Com Filtragem 2D de 5 pontos'])
xlabel('Pontos de Grade j')
ylabel('Pontos de Grade k')