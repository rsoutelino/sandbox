%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%              Cálculo de EOFs            %%%
%%%               WESTRAX 2                 %%%
%%%     Rafael Soutelino - agosto / 2006    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc

%%% LOOPING PARA LEITURA DOS DADOS

lon=[]; lat=[]; nest=[]; k=0; 

for i=[22:37 39:50];

    eval(['load ../ctd/w2cpz0',num2str(i),'.pro;']);
    eval(['data = w2cpz0',num2str(i),';'])
    p=data(:,1); 
    if data(10,8)==-9999 % tirando as estacoes que nao tem perfilagem pegasus 
        disp(['Estacao ',num2str(i),' nao tem pegasus!'])    
    else 
        disp('.')
        f=find(p==3000);
        if isempty(f)==0;
          k=k+1;
          lon=[lon data(1,6)*-1];
          lat=[lat data(1,5)];
          Usec(5:3002,k) = data(5:3002,8);
          Vsec(5:3002,k) = data(5:3002,9);
          psec(5:3002,k) = data(5:3002,1);
          nest=[nest i];
        else
          disp(['Estacao ',num2str(i),' nao atinge profundidade desejada'])
        end
    end
    eval(['clear w2cpz0',num2str(i)])    
end

clear data i f k p
 
Usec=Usec'/100;
Vsec=Vsec'/100;

uvsec = [Usec;Vsec];

%  dist=sw_dist(lat,lon,'km');
%  dist=[0 cumsum(dist)];

% calculando as EOFs

disp(' ')
disp('Calculando EOFs ...')
[perc,Fi,amp]=eoft(uvsec);

perc=flipud(perc); % o vetor sai de cabeca pra baixo

% separando as amplitudes U e V
ampU = amp(1:end/2,:); ampU=ampU';
ampV = amp((end/2+1):end,:); ampV=ampV';
















