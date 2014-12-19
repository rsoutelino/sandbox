%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%        PROGRAMA PARA CALCULO DE         %%%
%%%     FUNCAO DE CORRENTE GEOSTROFICA -    %%%
%%%               WESTRAX 2                 %%%
%%%     Rafael Soutelino - junho / 2006     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc

%%% LOOPING PARA LEITURA DOS DADOS

lon=[]; lat=[]; nest=[];

for i=[22:37 39:50]

    eval(['load ../ctd/w2cpz0',num2str(i),'.pro;']);
    eval(['data = w2cpz0',num2str(i),';'])
    if data(2,8)==9999 % tirando as estacoes que nao tem perfilagem pegasus
        x=0; % só pra constar no if
    else
        p=data(:,1);
        f=find(p==1500);
        if isempty(f)==0;
           lon=[lon data(1,6)*-1];
           lat=[lat data(1,5)];
           eval(['p',num2str(i),' = data(15:end,1);'])
           eval(['T',num2str(i),' = data(15:end,3);'])
           eval(['S',num2str(i),' = data(15:end,5);'])
           nest=[nest i];
        end
    end
    eval(['clear w2cpz0',num2str(i)])    
end

clear data

% macetes para plotar mapas

lonlim=[-53.5 -43]; latlim=[-0.5 10]; % limites de silveira_etal2000

m_proj('mercator','long',[lonlim(1) lonlim(2)],'lat',[latlim(1) latlim(2)],'on');

% comandos para plotar as variaveis
m_plot, m_contour, m_contour, m_quiver

% plotar linha de costa
m_usercoast('costa.mat','patch',[0 0 0])

