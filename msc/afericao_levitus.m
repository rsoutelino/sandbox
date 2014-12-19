%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  CONTROLE DE QUALIDADE DOS DADOS HIDROGRAFICOS
%                OCEANO LESTE 2
%      COMPARACAO COM OS CAMPOS DO LEVITUS
%      Rafael Soutelino - Mestrado - IOUSP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc;

%%% LIMITES OESTE-LESTE e SUL-NORTE DA AREA DE INTERESSE %%%

lonlim = [-41 -33];
latlim = [-20 -10];

%%% PERIODO DESEJADO %%%

disp(' ');
disp(['  ### Escolha do Periodo ###']);
disp(' ');
disp(['  - Para o campo anual: [anual]']);
disp(['  - Para o campo mensal: [jan] ou [fev] ou ...']);
disp(['  - Para media dos campos mensais: [jan,fev,mar,...]']);
disp(' ');
time = '[jan,fev,mar]'; %input(['  Periodo: '],'s');

time([1 end]) = [];

ft = find(time==',');

if isempty(ft) == 1
  Time = cellstr(time);
else
  for i = 1:length(ft)
    Time{i} = time(ft(i)-3:ft(i)-1);
  end
  Time{i+1} = time(ft(end)+1:ft(end)+3);
end

%%% LEITURA DA CLIMATOLOGIA DECLARADA %%%

for i = 1:length(Time)

  eval(['load ',Time{i},'/Twoa_',Time{i},'_brazilcoast.dat;']);
  eval(['load ',Time{i},'/Swoa_',Time{i},'_brazilcoast.dat;']);

  eval(['Twoa = Twoa_',Time{i},'_brazilcoast;']);
  eval(['Swoa = Swoa_',Time{i},'_brazilcoast;']);
  
  eval(['clear Twoa_',Time{i},'_brazilcoast;']);
  eval(['clear Swoa_',Time{i},'_brazilcoast;']);
  
  lat = Twoa(:,1)';
  lon = Twoa(:,2)';
  
  flat = find(lat >= min(latlim) & lat <= max(latlim));
  flon = find(lon >= min(lonlim) & lon <= max(lonlim));


  ii = 1;
  for jj = 1:length(flon)
    if isempty(find(flat == flon(jj))) == 0
      q(ii) = flon(jj);
      ii = ii+1;
    end
  end

  eval(['T',Time{i},' = Twoa(q,3:end);']);
  eval(['S',Time{i},' = Swoa(q,3:end);']);

end

Twoa = [Tjan; Tfev; Tmar];
Swoa = [Sjan; Sfev; Smar];

%%% Carregando os dados da Oceano Leste 2 %%%
Soe2 = []; Toe2 = [];

for i=[1:4 7:102 105:112];
     eval(['load ../leste2/ctd/filtrados/lesteII_ctd',num2str(i),'.dat;'])
     eval(['ctd = lesteII_ctd',num2str(i),';'])
     eval(['clear lesteII_ctd',num2str(i)])

     t = ctd(:,2);
     s = ctd(:,3);
     
     Soe2 = [Soe2; s]; Toe2 = [Toe2; t]; 
end


% plotando os diagramas TS espalhados

[Sg,Tg] = meshgrid(33.5:.25:38,-2:30);
Dg = sw_dens0(Sg,Tg)-1000;

figure(1)
set(1,'Color','w')
[c,h] = contour(Sg,Tg,Dg,20:1:40,'k');
hold on;
set(h,'Color',[0.7 0.7 0.7]);
cl = clabel(c,h,'VerticalAlignment','middle','fontsize',9,'labelspacing',400);
set(cl,'Color',[0.7 0.7 0.7]);
plot(Soe2,Toe2,'r.','markersize',.5)
plot(Swoa,Twoa,'b.','markersize',.5)
legend('Oceano Leste II','Levitus',0)
title('')
xlabel('Salinidade','fontsize',12)
ylabel('Temperatura Potencial ( \circ C)','fontsize',12)











































